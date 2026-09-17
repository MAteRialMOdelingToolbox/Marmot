/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck,
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
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

#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialHughesWinget.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotVoigt.h"
#include <Eigen/Geometry>
#include <functional>
#include <vector>

using namespace Marmot;
using namespace Marmot::Materials;
using namespace Marmot::Testing;
using namespace Marmot::FastorStandardTensors;

namespace {

  /// Isotropic elastic stiffness in Voigt notation.
  Marmot::Matrix6d isotropicStiffness( double E, double nu )
  {
    Marmot::Matrix6d C   = Marmot::Matrix6d::Zero();
    const double     fac = E / ( ( 1.0 + nu ) * ( 1.0 - 2.0 * nu ) );
    for ( int i = 0; i < 3; i++ )
      for ( int j = 0; j < 3; j++ )
        C( i, j ) = fac * ( i == j ? ( 1.0 - nu ) : nu );
    for ( int i = 3; i < 6; i++ )
      C( i, i ) = fac * ( 1.0 - 2.0 * nu ) / 2.0;
    return C;
  }

  /**
   * @brief Minimal hypoelastic stub: linear elastic, plus one state variable accumulating the
   *        magnitude of the applied strain increments.
   *
   * The state variable makes it observable whether the wrapped material was driven at all, which the
   * rigid-rotation test relies on.
   */
  class StubLinearElastic : public MarmotMaterialHypoElastic {
  public:
    StubLinearElastic( const double* props, int nProps, int matNumber )
      : MarmotMaterialHypoElastic( props, nProps, matNumber )
    {
      stateLayout.add( "accumulatedStrain", 1 );
      stateLayout.finalize();
    }

    void computeStress( state3D&                state,
                        Marmot::Matrix6d&       dStress_dStrain,
                        const Marmot::Vector6d& dStrain,
                        const timeInfo& ) const override
    {
      const Marmot::Matrix6d C = isotropicStiffness( materialProperties[0], materialProperties[1] );
      state.stress += C * dStrain;
      dStress_dStrain = C;
      if ( state.stateVars )
        state.stateVars[0] += dStrain.norm();
    }

    double getDensity( const double* ) const override { return 1.0; }
  };

  /**
   * @brief Hypoelastic stub whose updated stress depends *nonlinearly* on the incoming stress.
   *
   * @f$ \sigma^{(n+1)} = \sigma_{\text{rot}} / (1 + k\,\|\sigma_{\text{rot}}\|) + C:\Delta\varepsilon @f$,
   * so @f$ \partial\sigma^{(n+1)}/\partial\sigma_{\text{rot}} \neq \boldsymbol{I} @f$. This is what
   * separates the Analytic tangent from the Exact one, in the same way a plastic return map does,
   * without pulling a real plasticity model into a core test.
   */
  class StubStressDependent : public MarmotMaterialHypoElastic {
  public:
    StubStressDependent( const double* props, int nProps, int matNumber )
      : MarmotMaterialHypoElastic( props, nProps, matNumber )
    {
      stateLayout.add( "accumulatedStrain", 1 );
      stateLayout.finalize();
    }

    void computeStress( state3D&                state,
                        Marmot::Matrix6d&       dStress_dStrain,
                        const Marmot::Vector6d& dStrain,
                        const timeInfo& ) const override
    {
      const Marmot::Matrix6d C = isotropicStiffness( materialProperties[0], materialProperties[1] );
      const double           k = materialProperties[2];

      state.stress    = state.stress / ( 1.0 + k * state.stress.norm() ) + C * dStrain;
      dStress_dStrain = C;
      if ( state.stateVars )
        state.stateVars[0] += dStrain.norm();
    }

    double getDensity( const double* ) const override { return 1.0; }
  };

  /// Rotation tensor from an axis-angle pair.
  Tensor33d rotationTensor( const Eigen::Vector3d& axis, double angle )
  {
    const Eigen::Matrix3d Q = Eigen::AngleAxisd( angle, axis.normalized() ).toRotationMatrix();
    Tensor33d             t;
    Marmot::mapEigenToFastor( t ) = Q;
    return t;
  }

  /// Matrix product of two Fastor 3x3 tensors.
  Tensor33d mul( const Tensor33d& a, const Tensor33d& b )
  {
    return a % b;
  }

  /// Largest absolute component-wise difference of two fourth-order tensors.
  double maxAbsDiff( const Tensor3333d& a, const Tensor3333d& b )
  {
    double m = 0.0;
    for ( size_t i = 0; i < a.size(); i++ )
      m = std::max( m, std::abs( a.data()[i] - b.data()[i] ) );
    return m;
  }

  /**
   * @brief Reference tangent: central differences of tau(F), computed *in the test* so that the
   *        production Numerical mode is never used to validate itself.
   */
  template < typename MaterialType >
  Tensor3333d referenceTangent( MaterialType&                material,
                                const Tensor33d&             F,
                                const std::vector< double >& stateIn,
                                double                       time,
                                double                       dT )
  {
    const int   n = material.getNumberOfRequiredStateVars();
    Tensor3333d dTau_dF( 0.0 );

    auto tauAt = [&]( const Tensor33d& Fp ) {
      std::vector< double >                                 scratch = stateIn;
      MarmotMaterialFiniteStrain::ConstitutiveResponse< 3 > r;
      r.stateVars = scratch.data();
      MarmotMaterialFiniteStrain::AlgorithmicModuli< 3 > t;
      material.computeStress( r, t, { Fp }, { time, dT } );
      return r.tau;
    };

    for ( int k = 0; k < 3; k++ )
      for ( int l = 0; l < 3; l++ ) {
        const double h  = 1e-7 * std::max( 1.0, std::abs( F( k, l ) ) );
        Tensor33d    Fp = F, Fm = F;
        Fp( k, l ) += h;
        Fm( k, l ) -= h;
        const Tensor33d taup = tauAt( Fp );
        const Tensor33d taum = tauAt( Fm );
        for ( int i = 0; i < 3; i++ )
          for ( int j = 0; j < 3; j++ )
            dTau_dF( i, j, k, l ) = ( taup( i, j ) - taum( i, j ) ) / ( 2.0 * h );
      }
    (void)n;
    return dTau_dF;
  }

} // namespace

using Wrapper          = HughesWingetWrapper< StubLinearElastic >;
using WrapperNonlinear = HughesWingetWrapper< StubStressDependent >;
using WrapperExact     = HughesWingetWrapper< StubStressDependent, HughesWingetTangent::Exact >;
using WrapperNumerical = HughesWingetWrapper< StubStressDependent, HughesWingetTangent::Numerical >;

namespace {
  const std::vector< double > elasticProps   = { 20000.0, 0.3 };
  const std::vector< double > nonlinearProps = { 20000.0, 0.3, 2e-3 };

  /// Allocate and initialise a state vector for a wrapper instance.
  template < typename W >
  std::vector< double > freshState( W& w )
  {
    std::vector< double > s( w.getNumberOfRequiredStateVars(), 0.0 );
    w.initializeYourself( s.data(), int( s.size() ) );
    return s;
  }

  /// Drive one increment, returning the Kirchhoff stress and tangent.
  template < typename W >
  std::pair< Tensor33d, Tensor3333d > step( W&                     w,
                                            std::vector< double >& state,
                                            const Tensor33d&       F,
                                            double                 time = 0.0,
                                            double                 dT   = 1.0 )
  {
    MarmotMaterialFiniteStrain::ConstitutiveResponse< 3 > r;
    r.stateVars = state.data();
    MarmotMaterialFiniteStrain::AlgorithmicModuli< 3 > t;
    w.computeStress( r, t, { F }, { time, dT } );
    return { r.tau, t.dTau_dF };
  }
} // namespace

/// The state layout must expose exactly the three documented slots.
void testStateLayout()
{
  Wrapper w( elasticProps.data(), int( elasticProps.size() ), 1 );

  throwExceptionOnFailure( w.getNumberOfRequiredStateVars() == 9 + 6 + 1,
                           "unexpected state layout size in " + std::string( __PRETTY_FUNCTION__ ) );

  auto state = freshState( w );
  throwExceptionOnFailure( w.getStateView( "HughesWinget_F_n", state.data() ).stateSize == 9,
                           "F_n slot has the wrong size in " + std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( w.getStateView( "HughesWinget_sigma_n", state.data() ).stateSize == 6,
                           "sigma_n slot has the wrong size in " + std::string( __PRETTY_FUNCTION__ ) );

  // initializeYourself must leave F_n = I, otherwise F_mid is singular on the first increment.
  const Tensor33d Fn( state.data() );
  throwExceptionOnFailure( checkIfEqual( Fn, Spatial3D::I, 1e-15 ),
                           "F_n was not initialised to the identity in " + std::string( __PRETTY_FUNCTION__ ) );
}

/// An undeformed increment must produce no stress at all.
void testUndeformed()
{
  Wrapper w( elasticProps.data(), int( elasticProps.size() ), 1 );
  auto    state = freshState( w );

  const auto [tau, dTau_dF] = step( w, state, Spatial3D::I );

  throwExceptionOnFailure( checkIfEqual( tau, Tensor33d( 0.0 ), 1e-15 ),
                           "undeformed state produced non-zero stress in " + std::string( __PRETTY_FUNCTION__ ) );
}

/**
 * @brief Superimposed rigid rotation must be reproduced exactly, in a single increment.
 *
 * For F^(n+1) = Q F^(n) the Hughes-Winget increment dl is exactly skew, hence dEps vanishes and the
 * Cayley transform returns dR = Q identically. The wrapped material must therefore not be strained at
 * all, and the stress must merely co-rotate. This holds for any rotation with Q + I invertible.
 */
void testRigidRotationIsExact()
{
  const auto axes = fibonacciLatticeHemisphere< 12 >();

  for ( int a = 0; a < axes.rows(); a++ )
    for ( double angleDeg : { 5.0, 30.0, 90.0, 150.0 } ) {

      Wrapper w( elasticProps.data(), int( elasticProps.size() ), 1 );
      auto    state = freshState( w );

      // Pre-load with a deviatoric stretch so that sigma_n is non-trivial.
      Tensor33d F1               = Spatial3D::I;
      F1( 0, 0 )                 = 1.02;
      F1( 1, 1 )                 = 0.99;
      F1( 0, 1 )                 = 0.015;
      const auto [tauOld, dummy] = step( w, state, F1 );

      const double accumulatedBefore = state.back();

      const double    theta = angleDeg * Constants::Pi / 180.0;
      const Tensor33d Q     = rotationTensor( Eigen::Vector3d( axes( a, 0 ),
                                                           axes( a, 1 ),
                                                           std::sqrt( std::max( 0.0,
                                                                                1.0 - axes( a, 0 ) * axes( a, 0 ) -
                                                                                  axes( a, 1 ) * axes( a, 1 ) ) ) ),
                                          theta );

      const auto [tauNew, dummy2] = step( w, state, mul( Q, F1 ) );

      // The wrapped material must not have seen any strain.
      throwExceptionOnFailure( checkIfEqual( state.back(), accumulatedBefore, 1e-14 ),
                               "rigid rotation strained the wrapped material in " +
                                 std::string( __PRETTY_FUNCTION__ ) );

      // The Kirchhoff stress must co-rotate: tau_new = Q tau_old Q^T.
      const Tensor33d expected = mul( mul( Q, tauOld ), Fastor::transpose( Q ) );
      throwExceptionOnFailure( checkIfEqual( tauNew, expected, 1e-11 ),
                               "rigid rotation did not co-rotate the stress in " + std::string( __PRETTY_FUNCTION__ ) );
    }
}

/// At small strain the wrapper must reproduce the wrapped material driven directly.
void testSmallStrainAgreement()
{
  Wrapper w( elasticProps.data(), int( elasticProps.size() ), 1 );
  auto    state = freshState( w );

  StubLinearElastic                  direct( elasticProps.data(), int( elasticProps.size() ), 1 );
  double                             directStateVar = 0.0;
  MarmotMaterialHypoElastic::state3D directState;
  directState.stateVars = &directStateVar;

  // 20 increments of a small symmetric strain.
  Eigen::Matrix3d H = Eigen::Matrix3d::Zero();
  H( 0, 0 )         = 1e-6;
  H( 1, 1 )         = -4e-7;
  H( 0, 1 ) = H( 1, 0 ) = 7e-7;

  Eigen::Matrix3d Ftotal = Eigen::Matrix3d::Identity();
  for ( int n = 0; n < 20; n++ ) {
    Ftotal += H;
    Tensor33d F;
    Marmot::mapEigenToFastor( F ) = Ftotal;
    step( w, state, F );

    Marmot::Matrix6d C = Marmot::Matrix6d::Zero();
    direct.computeStress( directState,
                          C,
                          ContinuumMechanics::VoigtNotation::voigtFromStrainMatrix< 3 >( H ),
                          { 0.0, 1.0 } );
  }

  Tensor33d Ffinal;
  Marmot::mapEigenToFastor( Ffinal ) = Ftotal;
  const auto [tau, dTau_dF]          = step( w, state, Ffinal );

  const Eigen::Matrix3d sigmaDirect = ContinuumMechanics::VoigtNotation::stressMatrixFromVoigt< 3 >(
    directState.stress );
  Tensor33d expected;
  Marmot::mapEigenToFastor( expected ) = sigmaDirect;

  // Compare the Cauchy stress: tau = J sigma, and the J factor alone is larger than the tolerance below.
  const Tensor33d sigma = Tensor33d( tau / Ftotal.determinant() );

  // The mid-step objective rate and the naive additive strain differ at O(total strain), so the
  // relative deviation is bounded by ~2e-5 here; that is a property of the formulation, not an error.
  throwExceptionOnFailure( checkIfEqual( sigma, expected, 1e-5 ),
                           "small-strain response deviates from the wrapped material in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

namespace {
  /// A few deformation gradients combining stretch, shear and rotation.
  std::vector< Tensor33d > probeDeformations()
  {
    std::vector< Tensor33d > out;

    Tensor33d a = Spatial3D::I;
    a( 0, 0 )   = 1.10;
    a( 1, 1 )   = 0.95;
    a( 2, 2 )   = 1.03;
    a( 0, 1 )   = 0.20;
    a( 1, 2 )   = 0.10;
    a( 2, 0 )   = 0.03;
    out.push_back( a );

    // The same, with a substantial rigid rotation superimposed.
    out.push_back( mul( rotationTensor( Eigen::Vector3d( 0.3, -0.7, 0.65 ), 0.9 ), a ) );

    Tensor33d c = Spatial3D::I;
    c( 0, 1 )   = 0.35; // simple shear
    out.push_back( c );

    return out;
  }

  /// Preload a wrapper so that sigma_n is non-zero, which is what activates the rotational tangent term.
  template < typename W >
  void preload( W& w, std::vector< double >& state )
  {
    Tensor33d F = Spatial3D::I;
    F( 0, 0 )   = 1.05;
    F( 1, 1 )   = 0.97;
    F( 0, 1 )   = 0.04;
    F( 2, 1 )   = 0.02;
    step( w, state, F );
  }
} // namespace

/// The analytic tangent must reproduce central differences of tau(F) for a path-independent material.
void testAnalyticTangentVsCentralDifference()
{
  for ( const auto& F : probeDeformations() ) {
    Wrapper w( elasticProps.data(), int( elasticProps.size() ), 1 );
    auto    state = freshState( w );
    preload( w, state );

    const std::vector< double > stateIn = state;
    const auto [tau, analytic]          = step( w, state, F );
    const Tensor3333d reference         = referenceTangent( w, F, stateIn, 0.0, 1.0 );

    const double scale = std::max( 1.0, Fastor::norm( reference ) );
    throwExceptionOnFailure( checkIfEqual( analytic, reference, 1e-6 * scale ),
                             "analytic tangent deviates from central differences in " +
                               std::string( __PRETTY_FUNCTION__ ) );
  }
}

/**
 * @brief With a stress-dependent wrapped material the Analytic mode is knowingly inexact, and the
 *        Exact and Numerical modes are not.
 */
void testExactAndNumericalTangentModes()
{
  bool analyticDeviatedSomewhere = false;

  for ( const auto& F : probeDeformations() ) {

    WrapperNonlinear wAnalytic( nonlinearProps.data(), int( nonlinearProps.size() ), 1 );
    WrapperExact     wExact( nonlinearProps.data(), int( nonlinearProps.size() ), 1 );
    WrapperNumerical wNumeric( nonlinearProps.data(), int( nonlinearProps.size() ), 1 );

    auto sAnalytic = freshState( wAnalytic );
    auto sExact    = freshState( wExact );
    auto sNumeric  = freshState( wNumeric );
    preload( wAnalytic, sAnalytic );
    preload( wExact, sExact );
    preload( wNumeric, sNumeric );

    const std::vector< double > stateIn = sAnalytic;

    const auto [tauA, tangentAnalytic] = step( wAnalytic, sAnalytic, F );
    const auto [tauE, tangentExact]    = step( wExact, sExact, F );
    const auto [tauN, tangentNumeric]  = step( wNumeric, sNumeric, F );

    // The stress itself must not depend on how the tangent is obtained.
    throwExceptionOnFailure( checkIfEqual( tauA, tauE, 1e-12 ) && checkIfEqual( tauA, tauN, 1e-12 ),
                             "tangent mode changed the stress in " + std::string( __PRETTY_FUNCTION__ ) );

    const Tensor3333d reference = referenceTangent( wAnalytic, F, stateIn, 0.0, 1.0 );
    const double      scale     = std::max( 1.0, Fastor::norm( reference ) );

    throwExceptionOnFailure( checkIfEqual( tangentExact, reference, 1e-5 * scale ),
                             "Exact tangent deviates from central differences in " +
                               std::string( __PRETTY_FUNCTION__ ) );
    throwExceptionOnFailure( checkIfEqual( tangentNumeric, reference, 1e-5 * scale ),
                             "Numerical tangent deviates from central differences in " +
                               std::string( __PRETTY_FUNCTION__ ) );

    // Guard against Exact silently falling back to Analytic: for at least one of these deformations
    // the shortcut must be measurably worse. It need not be worse for all of them, since the neglected
    // term scales with the rotational part of the increment.
    if ( maxAbsDiff( tangentAnalytic, reference ) > 1e-5 * scale )
      analyticDeviatedSomewhere = true;
  }

  throwExceptionOnFailure( analyticDeviatedSomewhere,
                           "the stress-dependent stub never exercised the d(sigma)/d(sigma_rot) term, so "
                           "this test cannot distinguish Exact from Analytic in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

/// The Kirchhoff stress must stay symmetric.
void testStressSymmetry()
{
  for ( const auto& F : probeDeformations() ) {
    Wrapper w( elasticProps.data(), int( elasticProps.size() ), 1 );
    auto    state = freshState( w );
    preload( w, state );
    const auto [tau, t] = step( w, state, F );

    throwExceptionOnFailure( checkIfEqual( tau, Tensor33d( Fastor::transpose( tau ) ), 1e-12 ),
                             "Kirchhoff stress is not symmetric in " + std::string( __PRETTY_FUNCTION__ ) );
  }
}

/**
 * @brief Purely volumetric deformation: tau = J sigma, with the Hughes-Winget normal strain increment
 *        dEps_ii = 2 (lambda - 1) / (lambda + 1) rather than log(lambda) or (lambda - 1).
 */
void testVolumetricScaling()
{
  const double lambda = 1.05;

  Wrapper w( elasticProps.data(), int( elasticProps.size() ), 1 );
  auto    state = freshState( w );

  const auto [tau, t] = step( w, state, Tensor33d( lambda * Spatial3D::I ) );

  const double     dEpsIi  = 2.0 * ( lambda - 1.0 ) / ( lambda + 1.0 );
  Marmot::Vector6d dStrain = Marmot::Vector6d::Zero();
  dStrain.head( 3 ).setConstant( dEpsIi );

  StubLinearElastic                  direct( elasticProps.data(), int( elasticProps.size() ), 1 );
  double                             sv = 0.0;
  MarmotMaterialHypoElastic::state3D ds;
  ds.stateVars       = &sv;
  Marmot::Matrix6d C = Marmot::Matrix6d::Zero();
  direct.computeStress( ds, C, dStrain, { 0.0, 1.0 } );

  const double J = lambda * lambda * lambda;
  Tensor33d    expected;
  Marmot::mapEigenToFastor( expected ) = J * ContinuumMechanics::VoigtNotation::stressMatrixFromVoigt< 3 >( ds.stress );

  throwExceptionOnFailure( checkIfEqual( tau, expected, 1e-10 ),
                           "volumetric response does not match J * sigma in " + std::string( __PRETTY_FUNCTION__ ) );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{ testStateLayout,
                                                       testUndeformed,
                                                       testRigidRotationIsExact,
                                                       testSmallStrainAgreement,
                                                       testAnalyticTangentVsCentralDifference,
                                                       testExactAndNumericalTangentModes,
                                                       testStressSymmetry,
                                                       testVolumetricScaling };

  executeTestsAndCollectExceptions( tests );
  return 0;
}
