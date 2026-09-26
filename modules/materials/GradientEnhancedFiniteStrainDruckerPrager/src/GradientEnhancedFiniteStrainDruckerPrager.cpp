/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
 *
 * Thomas Mader thomas.mader@boku.ac.at
 *
 * This file is part of the MAteRialMOdellingToolbox (marmot).
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of marmot.
 * ---------------------------------------------------------------------
 */
#include "Marmot/GradientEnhancedFiniteStrainDruckerPrager.h"
#include "Marmot/MarmotEigenSystems.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotNumericalDifferentiation.h"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <functional>
#include <stdexcept>

namespace Marmot::Materials {

  using namespace Fastor;
  using namespace FastorStandardTensors;
  namespace Differentiation = NumericalAlgorithms::Differentiation;

  namespace {

    /// absolute tolerance of the return-map residual (the yield function is scaled by the cohesion)
    constexpr double innerNewtonTol        = 1e-12;
    constexpr int    nMaxInnerNewtonCycles = 50;
    constexpr int    nMaxHalvings          = 10;
    /// the deviatoric Mandel stress counts as vanished (apex) below this fraction of the cohesive strength
    constexpr double apexTol = 1e-10;
    /// relative tolerance of the coaxiality of a cone solution with the trial state
    constexpr double coaxialityTol = 1e-6;

    /// the i-th material property, after checking that it exists: the references below are bound in the
    /// member initializer list, i.e. before the constructor body could validate the property count
    const double& property( const double* properties, int nProperties, int i )
    {
      constexpr int nRequired = 10;
      if ( nProperties < nRequired )
        throw std::invalid_argument( MakeString() << "GradientEnhancedFiniteStrainDruckerPrager: expected at least "
                                                  << nRequired << " material properties, got " << nProperties );
      return properties[i];
    }

    using Matrix9d = Eigen::Matrix< double, 9, 9 >;

    /// d(F Fp^-1) / dF in the row-major flattening of both
    Matrix9d dFeTrial_dF( const Tensor33d& FpInv )
    {
      Matrix9d D = Matrix9d::Zero();
      for ( int i = 0; i < 3; i++ )
        for ( int J = 0; J < 3; J++ )
          for ( int L = 0; L < 3; L++ )
            D( 3 * i + J, 3 * i + L ) = FpInv( L, J );
      return D;
    }

    /// d ln det F / dF = F^-T in the row-major flattening
    Eigen::RowVectorXd dLogDet_dF( const Tensor33d& F )
    {
      const Tensor33d    Finv = Fastor::inverse( F );
      Eigen::RowVectorXd d( 9 );
      for ( int k = 0; k < 3; k++ )
        for ( int L = 0; L < 3; L++ )
          d( 3 * k + L ) = Finv( L, k );
      return d;
    }

    /**
     * Newton's method for R( X ) = 0 with the Jacobian by the complex step, and a backtracking line search: a step
     * is halved while the residual cannot be evaluated, does not decrease, or leads to an inadmissible iterate. On
     * success, R and dR_dX belong to the converged X.
     */
    bool newtonWithBacktracking( const Differentiation::Complex::vector_to_vector_function_type& residual,
                                 Eigen::VectorXd&                                                X,
                                 Eigen::VectorXd&                                                R,
                                 Eigen::MatrixXd&                                                dR_dX,
                                 const std::function< bool( const Eigen::VectorXd& ) >&          admissible )
    {
      try {
        std::tie( R, dR_dX ) = Differentiation::Complex::forwardDifference( residual, X );
      }
      catch ( const std::exception& ) {
        return false;
      }
      for ( int counter = 0; counter <= nMaxInnerNewtonCycles; counter++ ) {
        if ( !R.allFinite() || !dR_dX.allFinite() )
          return false;
        if ( R.norm() < innerNewtonTol )
          return true;

        const Eigen::VectorXd dX   = -dR_dX.colPivHouseholderQr().solve( R );
        double                step = 1.0;
        for ( int halving = 0; halving <= nMaxHalvings; halving++, step *= 0.5 ) {
          const Eigen::VectorXd Xtrial = X + step * dX;
          if ( !admissible( Xtrial ) ) {
            if ( halving == nMaxHalvings )
              return false;
            continue;
          }
          try {
            const auto [Rtrial, Jtrial] = Differentiation::Complex::forwardDifference( residual, Xtrial );
            if ( Rtrial.allFinite() && ( Rtrial.norm() < R.norm() || halving == nMaxHalvings ) ) {
              X     = Xtrial;
              R     = Rtrial;
              dR_dX = Jtrial;
              break;
            }
          }
          catch ( const std::exception& ) {
            if ( halving == nMaxHalvings )
              return false;
          }
        }
      }
      return false;
    }

    Fastor::Tensor< double, 3 > principalValues( const Tensor33d& symmetric )
    {
      return Math::computeEigenSystemJacobi( symmetric ).first;
    }

  } // namespace

  GradientEnhancedFiniteStrainDruckerPrager::GradientEnhancedFiniteStrainDruckerPrager(
    const double* materialProperties,
    int           nMaterialProperties,
    int           materialNumber )
    : MarmotMaterialGradientEnhancedFiniteStrain( materialProperties, nMaterialProperties, materialNumber ),
      K( property( materialProperties, nMaterialProperties, 0 ) ),
      G( property( materialProperties, nMaterialProperties, 1 ) ),
      c0( property( materialProperties, nMaterialProperties, 2 ) ),
      frictionAngle( property( materialProperties, nMaterialProperties, 3 ) ),
      dilatancyAngle( property( materialProperties, nMaterialProperties, 4 ) ),
      H( property( materialProperties, nMaterialProperties, 5 ) ),
      softeningModulus( property( materialProperties, nMaterialProperties, 6 ) ),
      maxDamage( property( materialProperties, nMaterialProperties, 7 ) ),
      nonLocalRadius( property( materialProperties, nMaterialProperties, 8 ) ),
      weightingParameter( property( materialProperties, nMaterialProperties, 9 ) ),
      eta( outerConeParameters( frictionAngle ).first ),
      xi( outerConeParameters( frictionAngle ).second ),
      etaBar( outerConeParameters( dilatancyAngle ).first )
  {
    if ( K <= 0.0 || G <= 0.0 )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": bulk and shear modulus must be positive" );
    if ( maxDamage < 0.0 || maxDamage >= 1.0 )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": expected 0 <= maximum damage < 1" );
    if ( c0 <= 0.0 || softeningModulus <= 0.0 )
      throw std::invalid_argument( MakeString()
                                   << __PRETTY_FUNCTION__ << ": cohesion and softening modulus must be positive" );
    if ( dilatancyAngle > frictionAngle || frictionAngle < 0.0 || dilatancyAngle < 0.0 || frictionAngle >= 90.0 )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__
                                                << ": expected 0 <= dilatancy angle <= friction angle < 90 deg" );

    stateLayout.add( "Fp", 9 );     // plastic deformation gradient
    stateLayout.add( "alphaP", 1 ); // plastic hardening variable
    stateLayout.add( "alphaD", 1 ); // local damage variable, the source L of the nonlocal balance
    stateLayout.add( "kappa", 1 );  // damage history
    stateLayout.add( "omega", 1 );  // scalar damage
    stateLayout.finalize();
  }

  std::pair< double, double > GradientEnhancedFiniteStrainDruckerPrager::outerConeParameters( double angle )
  {
    const double s = std::sin( Math::degToRad( angle ) );
    const double c = std::cos( Math::degToRad( angle ) );
    const double d = std::sqrt( 3.0 ) * ( 3.0 - s );
    return { 6.0 * s / d, 6.0 * c / d };
  }

  double GradientEnhancedFiniteStrainDruckerPrager::getDensity( const double* stateVars ) const
  {
    if ( nMaterialProperties <= 10 )
      throw std::runtime_error( MakeString() << __PRETTY_FUNCTION__ << ": density not provided (property 10)" );
    return materialProperties[10];
  }

  void GradientEnhancedFiniteStrainDruckerPrager::initializeYourself( double* stateVars, int nStateVars )
  {
    for ( int i = 0; i < nStateVars; ++i )
      stateVars[i] = 0.0;

    const Tensor33d I = identity< double >();
    std::memcpy( stateLayout.getPtr( stateVars, "Fp" ), I.data(), 9 * sizeof( double ) );
  }

  GradientEnhancedFiniteStrainDruckerPrager::ReturnMapping GradientEnhancedFiniteStrainDruckerPrager::returnMapping(
    const Tensor33d& F,
    const Tensor33d& FpOld,
    double           alphaPOld ) const
  {
    const Tensor33d FpOldInv = Fastor::inverse( FpOld );
    const Tensor33d FeTrial  = F % FpOldInv;

    const double JeTrial = Fastor::determinant( FeTrial );
    if ( !( JeTrial > 0.0 ) || !std::isfinite( JeTrial ) )
      throw StressUpdateFailed( MakeString() << __PRETTY_FUNCTION__ << ": non-positive elastic volume ratio" );

    if ( yieldFunction( mandelStress( FeTrial ), alphaPOld ) <= 0.0 ) {
      ReturnMapping elastic;
      elastic.Fe              = FeTrial;
      elastic.FpNew           = FpOld;
      elastic.alphaP          = alphaPOld;
      elastic.dEpVol          = 0.0;
      elastic.plasticWork     = 0.0;
      elastic.dFe_dF          = dFeTrial_dF( FpOldInv );
      elastic.dDeltaAlphaD_dF = Eigen::RowVectorXd::Zero( 9 );
      return elastic;
    }

    // beyond the apex pressure, try the apex first: its admissibility decides between the two (the volumetric
    // plastic flow only lowers the pressure and the apex pressure only grows with the hardening, so
    // p_trial >= p_apex( alpha_n ) is necessary for the apex)
    const double pTrial       = Fastor::trace( mandelStress( FeTrial ) ) / 3.;
    const bool   apexPossible = eta > 0.0 && etaBar > 0.0 && pTrial >= xi * ( c0 + H * alphaPOld ) / eta;
    if ( apexPossible ) {
      try {
        return returnToApex( F, FpOld, alphaPOld );
      }
      catch ( const StressUpdateFailed& ) {
        // not admissible: the state belongs to the cone
      }
    }

    bool       converged = false;
    const auto cone      = returnToCone( F, FpOld, alphaPOld, converged );
    if ( converged )
      return cone;

    if ( !apexPossible && eta > 0.0 && etaBar <= 0.0 && pTrial >= xi * ( c0 + H * alphaPOld ) / eta )
      return returnToApex( F, FpOld, alphaPOld ); // reports that there is no apex without dilatancy

    throw StressUpdateFailed( MakeString() << __PRETTY_FUNCTION__ << ": return to the cone not successful" );
  }

  GradientEnhancedFiniteStrainDruckerPrager::ReturnMapping GradientEnhancedFiniteStrainDruckerPrager::returnToCone(
    const Tensor33d& F,
    const Tensor33d& FpOld,
    double           alphaPOld,
    bool&            converged ) const
  {
    using namespace Eigen;
    using complexDouble = std::complex< double >;

    const Tensor33d FpOldInv = Fastor::inverse( FpOld );
    const Tensor33d FeTrial  = F % FpOldInv;

    auto residual = [&]( const VectorXcd& X_ ) -> VectorXcd {
      return coneResidual< complexDouble >( X_, FeTrial, alphaPOld );
    };

    // initial guesses: the return of the linearized (Hencky) problem, exact for small strains, and fractions of it,
    // which keep the iterates away from the singular vertex of the cone
    const Tensor33d MTrial        = mandelStress( FeTrial );
    const Tensor33d devTrial      = deviatoric( MTrial );
    const double    dLambdaLinear = std::max( 0.0,
                                           yieldFunction( MTrial, alphaPOld ) /
                                             ( G + K * eta * etaBar + xi * xi * H ) );

    // an iterate whose deviator has turned against the trial deviator has crossed the vertex of the cone
    const auto admissible = [&]( const VectorXd& X_ ) {
      const Tensor33d dev_ = deviatoric( mandelStress( Tensor33d( X_.head( 9 ).eval().data() ) ) );
      return Fastor::inner( dev_, devTrial ) > 0.0;
    };

    VectorXd X( 11 ), R;
    MatrixXd dR_dX;
    converged = false;
    for ( const double fraction : { 1.0, 0.5, 0.25, 0.9, 0.0 } ) {
      const double    dLambda0 = fraction * dLambdaLinear;
      const Tensor33d dFp0     = exponentialMap( Tensor33d( dLambda0 * flowDirection( MTrial ) ) );
      const Tensor33d Fe0      = FeTrial % Fastor::inverse( dFp0 );
      X.head( 9 )              = Map< const Matrix< double, 9, 1 > >( Fe0.data() );
      X( 9 )                   = alphaPOld + xi * dLambda0;
      X( 10 )                  = dLambda0;
      if ( admissible( X ) && newtonWithBacktracking( residual, X, R, dR_dX, admissible ) ) {
        converged = true;
        break;
      }
    }
    if ( !converged )
      return {};
    converged = true;

    ReturnMapping r;
    r.plastic               = true;
    r.Fe                    = Tensor33d( X.head( 9 ).eval().data() );
    r.alphaP                = X( 9 );
    const double    dLambda = X( 10 );
    const Tensor33d M       = mandelStress( r.Fe );

    // a solution on the cone needs a positive multiplier and a deviatoric stress. For an isotropic model it is
    // moreover coaxial with the trial state (M and M_trial commute), with a deviator that points the same way as the
    // trial deviator (the finite-strain analogue of sqrt(J2_trial) - G dLambda >= 0); near the apex, Newton may
    // also converge to a spurious, non-coaxial solution. In all these cases the apex takes over.
    const Tensor33d dev        = deviatoric( M );
    const double    normDev    = std::sqrt( Fastor::inner( dev, dev ) );
    const double    normTr     = std::sqrt( Fastor::inner( devTrial, devTrial ) );
    const Tensor33d commutator = dev % devTrial - devTrial % dev;
    if ( dLambda < 0.0 || std::sqrt( 0.5 ) * normDev <= apexTol * xi * c0 || Fastor::inner( dev, devTrial ) <= 0.0 ||
         std::sqrt( Fastor::inner( commutator, commutator ) ) > coaxialityTol * normDev * normTr ) {
      converged = false;
      return {};
    }

    r.FpNew = Fastor::inverse( r.Fe ) % FeTrial % FpOld;

    const Tensor33d dEp = dLambda * flowDirection( M );
    r.dEpVol            = etaBar * dLambda; // tr dg/dM = etaBar
    r.plasticWork       = Fastor::inner( M, dEp );

    // implicit function theorem: R( X, FeTrial ) = 0 with dR/dFeTrial = [-I; 0; 0]
    MatrixXd dR_dFeTrial     = MatrixXd::Zero( 11, 9 );
    dR_dFeTrial.topRows( 9 ) = -Matrix9d::Identity();
    const MatrixXd dX_dF     = -dR_dX.colPivHouseholderQr().solve( dR_dFeTrial * dFeTrial_dF( FpOldInv ) );
    r.dFe_dF                 = dX_dF.topRows( 9 );

    // the local damage increment etaBar dLambda (dLambda >= 0)
    r.dDeltaAlphaD_dF = etaBar * dX_dF.row( 10 );

    return r;
  }

  GradientEnhancedFiniteStrainDruckerPrager::ReturnMapping GradientEnhancedFiniteStrainDruckerPrager::returnToApex(
    const Tensor33d& F,
    const Tensor33d& FpOld,
    double           alphaPOld ) const
  {
    using namespace Eigen;
    using complexDouble = std::complex< double >;

    if ( eta <= 0.0 || etaBar <= 0.0 )
      throw StressUpdateFailed( MakeString()
                                << __PRETTY_FUNCTION__ << ": no apex to return to without friction and dilatancy" );

    const Tensor33d FpOldInv   = Fastor::inverse( FpOld );
    const Tensor33d FeTrial    = F % FpOldInv;
    const double    thetaTrial = std::log( Fastor::determinant( FeTrial ) );

    auto residual = [&]( const VectorXcd& X_ ) -> VectorXcd {
      return apexResidual< complexDouble >( X_, thetaTrial, alphaPOld );
    };

    VectorXd X( 2 );
    X << thetaTrial, alphaPOld;
    VectorXd   R;
    MatrixXd   dR_dX;
    const bool converged = newtonWithBacktracking( residual, X, R, dR_dX, []( const VectorXd& ) { return true; } );
    if ( !converged || X( 1 ) < alphaPOld )
      throw StressUpdateFailed( MakeString() << __PRETTY_FUNCTION__ << ": return to the apex not successful" );

    const double theta = X( 0 );

    // by isotropy, Fp is determined up to a rotation: take Fe = Je^(1/3) I
    ReturnMapping r;
    r.plastic = true;
    r.Fe      = std::exp( theta / 3. ) * identity< double >();
    r.FpNew   = std::exp( -theta / 3. ) * F;
    r.alphaP  = X( 1 );

    // the plastic log strain increment: the trial elastic log strain minus the new, spherical one
    const auto dEpPrincipal = [&]( const Tensor33d& F_, double theta_ ) {
      const Tensor33d                   FeTrial_ = F_ % FpOldInv;
      const Fastor::Tensor< double, 3 > b = principalValues( Tensor33d( Fastor::transpose( FeTrial_ ) % FeTrial_ ) );
      Fastor::Tensor< double, 3 >       dEp;
      for ( int a = 0; a < 3; a++ )
        dEp( a ) = 0.5 * std::log( b( a ) ) - theta_ / 3.;
      return dEp;
    };
    r.dEpVol = thetaTrial - theta;

    // the apex is admissible only if the plastic increment lies in the subdifferential of g at the vertex:
    // sqrt(2) |dev dEp| <= dEp_v / etaBar (otherwise the state belongs to the cone)
    {
      const Fastor::Tensor< double, 3 > dEp    = dEpPrincipal( F, theta );
      const double                      dEpVol = dEp( 0 ) + dEp( 1 ) + dEp( 2 );
      double                            dev2   = 0.0;
      for ( int a = 0; a < 3; a++ )
        dev2 += std::pow( dEp( a ) - dEpVol / 3., 2 );
      if ( std::sqrt( 2.0 * dev2 ) > ( 1.0 + 1e-8 ) * dEpVol / etaBar + 1e-14 )
        throw StressUpdateFailed( MakeString() << __PRETTY_FUNCTION__ << ": the apex state is not admissible" );
    }
    const double p = Fastor::trace( mandelStress( r.Fe ) ) / 3.;
    r.plasticWork  = p * ( thetaTrial - theta );

    // implicit function theorem: R( X, thetaTrial ) = 0 with dR/dthetaTrial = [0; -xi/etaBar]
    const Vector2d           dR_dThetaTrial( 0.0, -xi / etaBar );
    const Vector2d           dX_dThetaTrial = -dR_dX.colPivHouseholderQr().solve( dR_dThetaTrial );
    const Eigen::RowVectorXd dTheta_dF      = dX_dThetaTrial( 0 ) * dLogDet_dF( F );

    r.dFe_dF = MatrixXd::Zero( 9, 9 );
    for ( int i = 0; i < 3; i++ )
      r.dFe_dF.row( 3 * i + i ) = std::exp( theta / 3. ) / 3. * dTheta_dF;

    // the local damage increment thetaTrial - theta (>= 0 at an admissible apex)
    r.dDeltaAlphaD_dF = dLogDet_dF( F ) - dTheta_dF;

    return r;
  }

  void GradientEnhancedFiniteStrainDruckerPrager::computeStress( ConstitutiveResponse< 3 >& response,
                                                                 AlgorithmicModuli< 3 >&    tangents,
                                                                 const Deformation< 3 >&    deformation,
                                                                 const TimeIncrement&       timeIncrement ) const
  {
    using namespace Eigen;
    using complexDouble = std::complex< double >;

    double* sv = response.stateVars;
    // read the plastic deformation gradient through the same Tensor33d( ptr ) construction that writes it
    const Tensor33d FpOld( stateLayout.getPtr( sv, "Fp" ) );
    double&         alphaP = stateLayout.getAs< double& >( sv, "alphaP" );
    double&         alphaD = stateLayout.getAs< double& >( sv, "alphaD" );
    double&         kappa  = stateLayout.getAs< double& >( sv, "kappa" );
    double&         omega  = stateLayout.getAs< double& >( sv, "omega" );

    const Tensor33d F( deformation.F );
    const double    N = deformation.N;

    ReturnMapping r;
    try {
      r = returnMapping( F, FpOld, alphaP );
    }
    catch ( const StressUpdateFailed& ) {
      throw;
    }
    catch ( const std::exception& e ) {
      // e.g. a failed tensor exponential: to the host, this is a stress update that needs a smaller increment
      throw StressUpdateFailed( MakeString() << __PRETTY_FUNCTION__ << ": " << e.what() );
    }

    std::memcpy( stateLayout.getPtr( sv, "Fp" ), r.FpNew.data(), 9 * sizeof( double ) );
    alphaP = r.alphaP;
    alphaD += deltaAlphaLocal( r.dEpVol );

    // implicit-gradient damage, irreversible through the history maximum kappa
    const double m             = weightingParameter;
    const double omegaOld      = omega;
    const double kappaOld      = kappa;
    const double alphaWeighted = m * N + ( 1.0 - m ) * alphaD;
    kappa                      = std::max( kappaOld, alphaWeighted );
    const double omegaFree     = kappa > 0.0 ? 1.0 - std::exp( -kappa / softeningModulus ) : 0.0;
    omega                      = std::min( omegaFree, maxDamage );
    const bool   loading       = alphaWeighted > kappaOld && kappa > 0.0 && omegaFree < maxDamage;
    const double dOmega_dKappa = loading ? std::exp( -kappa / softeningModulus ) / softeningModulus : 0.0;

    // effective Kirchhoff stress and its sensitivity to Fe, by the complex step
    const auto tauOfFe = [&]( const VectorXcd& Fe_ ) -> VectorXcd {
      const Tensor33t< complexDouble > tau = kirchhoffStress( Tensor33t< complexDouble >( Fe_.eval().data() ) );
      return Map< const Matrix< complexDouble, 9, 1 > >( tau.data() );
    };
    const auto [tauEffFlat,
                dTauEff_dFe] = Differentiation::Complex::forwardDifference( tauOfFe,
                                                                            Map< const Matrix< double, 9, 1 > >(
                                                                              r.Fe.data() ) );
    const Tensor33d tauEff( tauEffFlat.data() );
    const double    psiEff = ContinuumMechanics::EnergyDensityFunctions::
      PenceGouPotentialB( Tensor33d( ContinuumMechanics::DeformationMeasures::rightCauchyGreen( r.Fe ) ), K, G );

    response.tau                  = ( 1.0 - omega ) * tauEff;
    response.L                    = alphaD;
    response.nonLocalRadius       = nonLocalRadius;
    response.elasticEnergyDensity = ( 1.0 - omega ) * psiEff;
    // cumulative: the host carries the dissipation of the previous increments in
    response.dissipation += ( 1.0 - omega ) * r.plasticWork + psiEff * ( omega - omegaOld );

    // tangents: tau = (1 - omega) tauEff( Fe( F ) ), omega( kappa ), kappa = m N + (1 - m) alphaD( F )
    const Eigen::RowVectorXd dL_dF      = r.dDeltaAlphaD_dF;
    const Eigen::RowVectorXd dOmega_dF  = dOmega_dKappa * ( 1.0 - m ) * dL_dF;
    const double             dOmega_dN  = dOmega_dKappa * m;
    const MatrixXd           dTauEff_dF = dTauEff_dFe * r.dFe_dF;

    for ( int i = 0; i < 3; i++ )
      for ( int j = 0; j < 3; j++ ) {
        for ( int k = 0; k < 3; k++ )
          for ( int l = 0; l < 3; l++ )
            tangents.dTau_dF( i, j, k, l ) = ( 1.0 - omega ) * dTauEff_dF( 3 * i + j, 3 * k + l ) -
                                             tauEff( i, j ) * dOmega_dF( 3 * k + l );
        tangents.dTau_dN( i, j ) = -tauEff( i, j ) * dOmega_dN;
      }
    for ( int k = 0; k < 3; k++ )
      for ( int l = 0; l < 3; l++ )
        tangents.dL_dF( k, l ) = dL_dF( 3 * k + l );
    tangents.dL_dN = 0.0;
  }

} // namespace Marmot::Materials
