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
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrainFactory.h"
#include "Marmot/MarmotMaterialGradientEnhancedHughesWinget.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include <Eigen/Geometry>
#include <functional>
#include <memory>
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

  /// Material properties for IncrementalGradientStub: E, nu, k (nonlocal degradation coefficient of
  /// the stress increment), R (nonlocal radius, c = R^2), c2 (linear coupling of KLocal on N).
  const std::vector< double > stubProps = { 20000.0, 0.25, 0.5, 0.02, 3.0 };

  /**
   * @brief A minimal, genuinely INCREMENTAL, hypoelastic, gradient-enhanced test-only material.
   *
   * @f$ \sigma^{(n+1)} = \sigma_{\text{in}} + (1 - k\,\bar N)\,C:\Delta\varepsilon @f$ -- it reads AND
   * updates the stress it is handed, which is the property GradientEnhancedHughesWingetWrapper
   * actually requires of a wrapped material (see the wrapper's doxygen): the degradation multiplies
   * only the newly-added strain-increment term, so @f$ \partial\sigma^{(n+1)}/\partial\sigma_{\text{in}}
   * = I @f$ exactly, regardless of the nonlocal field.
   *
   * A scalar "equivalent strain" state variable accumulates a fixed linear functional of
   * @f$ \Delta\varepsilon @f$ with distinct, non-zero weights on all six Voigt components, so
   * @c dKLocalddStrain -- and hence the wrapper's @c dL_dF -- is non-trivial and exercises every strain
   * component, including shear. @c KLocal also has an explicit linear dependence on the nonlocal field,
   * so @c dKLocalddK (the wrapper's @c dL_dN) is non-trivial too. The stress' dependence on the
   * nonlocal field (the same isotropic degradation) makes @c dStressddK -- the wrapper's @c dTau_dN --
   * non-trivial without disturbing the exact incremental structure above.
   */
  class IncrementalGradientStub : public MarmotMaterialGeneralGradientEnhancedHypoElastic< 1 > {
  public:
    IncrementalGradientStub( const double* props, int nProps, int matNumber )
      : MarmotMaterialGeneralGradientEnhancedHypoElastic< 1 >( props, nProps, matNumber )
    {
      stateLayout.add( "equivalentStrain", 1 );
      stateLayout.finalize();
    }

    void computeStress( response& res, tangents& tan, const increment& inc ) const override
    {
      const double E  = materialProperties[0];
      const double nu = materialProperties[1];
      const double k  = materialProperties[2];
      const double R  = materialProperties[3];
      const double c2 = materialProperties[4];

      const Marmot::Matrix6d C = isotropicStiffness( E, nu );
      const Marmot::Vector6d a = ( Marmot::Vector6d() << 1.0, 0.8, 0.6, 0.4, 0.3, 0.2 ).finished();

      const double phi         = inc.K( 0 );
      const double degradation = 1.0 - k * phi;

      double& eqStrain = stateLayout.getAs< double& >( res.stateVars, "equivalentStrain" );

      // Incremental: reads and updates the handed-in stress, as a genuine hypoelastic rate law must.
      res.stress = res.stress + degradation * ( C * inc.dStrain );

      eqStrain += a.dot( inc.dStrain );
      res.KLocal( 0 ) = eqStrain + c2 * phi;
      res.c( 0 )      = R * R;

      tan.dStressddStrain          = degradation * C;
      tan.dStressddK.col( 0 )      = -k * ( C * inc.dStrain );
      tan.dKLocalddStrain.row( 0 ) = a.transpose();
      tan.dKLocalddK( 0, 0 )       = c2;
    }

    double getDensity( const double* ) const override { return 1.0; }

    std::vector< double > getNonlocalViscosity( const double* ) const override { return { 1.0 }; }
  };

  using Wrapper          = GradientEnhancedHughesWingetWrapper< IncrementalGradientStub >;
  using WrapperExact     = GradientEnhancedHughesWingetWrapper< IncrementalGradientStub, HughesWingetTangent::Exact >;
  using WrapperNumerical = GradientEnhancedHughesWingetWrapper< IncrementalGradientStub,
                                                                HughesWingetTangent::Numerical >;

  using GEFiniteStrain  = MarmotMaterialGradientEnhancedFiniteStrain;
  using GEResponse      = GEFiniteStrain::ConstitutiveResponse< 3 >;
  using GETangents      = GEFiniteStrain::AlgorithmicModuli< 3 >;
  using GEDeformation   = GEFiniteStrain::Deformation< 3 >;
  using GETimeIncrement = GEFiniteStrain::TimeIncrement;

  /// Allocate and initialise a state vector for a gradient-enhanced finite-strain material instance.
  std::vector< double > freshState( GEFiniteStrain& w )
  {
    std::vector< double > s( w.getNumberOfRequiredStateVars(), 0.0 );
    w.initializeYourself( s.data(), int( s.size() ) );
    return s;
  }

  /// Result of a single computeStress call: the pieces every test case needs.
  struct StepResult {
    Tensor33d  tau;
    double     L;
    double     nonLocalRadius;
    GETangents tangents;
  };

  /// Drive one increment, returning the Kirchhoff stress, local driving force, nonlocal radius and tangents.
  StepResult step( GEFiniteStrain&        w,
                   std::vector< double >& state,
                   const Tensor33d&       F,
                   double                 N,
                   double                 time = 0.0,
                   double                 dT   = 1.0 )
  {
    GEResponse r;
    r.stateVars = state.data();
    GETangents t;
    w.computeStress( r, t, GEDeformation{ F, N }, GETimeIncrement{ time, dT } );
    return { r.tau, r.L, r.nonLocalRadius, t };
  }

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

  /**
   * @brief A minimal test-only gradient-enhanced hypoelastic material that records the nonlocal
   *        increment @c dK it was handed, and otherwise responds trivially.
   *
   * Used exclusively by testNonlocalIncrementBookkeeping() to verify that the wrapper hands the
   * wrapped material @f$ dK = \bar N - \bar N^{(n)} @f$, not the total field.
   */
  class RecordingIncrementMaterial : public MarmotMaterialGeneralGradientEnhancedHypoElastic< 1 > {
  public:
    RecordingIncrementMaterial( const double* props, int nProps, int matNumber )
      : MarmotMaterialGeneralGradientEnhancedHypoElastic< 1 >( props, nProps, matNumber )
    {
      stateLayout.add( "recordedIncrement", 1 );
      stateLayout.finalize();
    }

    void computeStress( response& res, tangents& tan, const increment& inc ) const override
    {
      double& recorded = stateLayout.getAs< double& >( res.stateVars, "recordedIncrement" );
      recorded         = inc.dK( 0 );

      // Trivial response: no stress, no coupling. Sufficient for the wrapper's kinematics, which
      // only needs a finite, well-defined nonlocal interaction to report a finite radius.
      res.stress      = Marmot::Vector6d::Zero();
      res.KLocal( 0 ) = 0.0;
      res.c( 0 )      = 1.0;
      (void)tan;
    }

    double getDensity( const double* ) const override { return 1.0; }

    std::vector< double > getNonlocalViscosity( const double* ) const override { return { 1.0 }; }
  };

  using RecordingWrapper = GradientEnhancedHughesWingetWrapper< RecordingIncrementMaterial >;

} // namespace

/// The wrapper's state layout must expose the four documented slots, sized 9 + 6 + 1 + (the wrapped
/// material's own state), in that order.
void testStateLayout()
{
  Wrapper w( stubProps.data(), int( stubProps.size() ), 1 );

  IncrementalGradientStub bare( stubProps.data(), int( stubProps.size() ), 1 );
  const int               nBase = bare.getNumberOfRequiredStateVars();

  throwExceptionOnFailure( w.getNumberOfRequiredStateVars() == 16 + nBase,
                           "unexpected total state layout size in " + std::string( __PRETTY_FUNCTION__ ) );

  auto state = freshState( w );
  throwExceptionOnFailure( w.getStateView( "HughesWinget_F_n", state.data() ).stateSize == 9,
                           "F_n slot has the wrong size in " + std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( w.getStateView( "HughesWinget_sigma_n", state.data() ).stateSize == 6,
                           "sigma_n slot has the wrong size in " + std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( w.getStateView( "HughesWinget_N_n", state.data() ).stateSize == 1,
                           "N_n slot has the wrong size in " + std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( w.getStateView( "materialstate", state.data() ).stateSize == nBase,
                           "materialstate slot has the wrong size in " + std::string( __PRETTY_FUNCTION__ ) );

  // initializeYourself must leave F_n = I, otherwise F_mid is singular on the first increment.
  const Tensor33d Fn( state.data() );
  throwExceptionOnFailure( checkIfEqual( Fn, Spatial3D::I, 1e-15 ),
                           "F_n was not initialised to the identity in " + std::string( __PRETTY_FUNCTION__ ) );
}

/// An undeformed increment must produce no stress, and a finite, positive nonlocal radius.
void testUndeformed()
{
  Wrapper w( stubProps.data(), int( stubProps.size() ), 1 );
  auto    state = freshState( w );

  const auto s = step( w, state, Spatial3D::I, 0.0 );

  throwExceptionOnFailure( checkIfEqual( s.tau, Tensor33d( 0.0 ), 1e-12 ),
                           "undeformed state produced non-zero stress in " + std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( std::isfinite( s.nonLocalRadius ) && s.nonLocalRadius > 0.0,
                           "nonlocal radius is not finite and positive in " + std::string( __PRETTY_FUNCTION__ ) );
}

/**
 * @brief Superimposed rigid rotation must be reproduced exactly, in a single increment.
 *
 * For F^(n+1) = Q F^(n) the Hughes-Winget increment dl is exactly skew, hence dEps vanishes and the
 * Cayley transform returns dR = Q identically -- see testRigidRotationIsExact in
 * TestMarmotMaterialHughesWinget.cpp for the derivation. With a genuinely incremental wrapped material
 * (sigma_new = sigma_in + ... with dStrain = 0) the wrapped material returns exactly the rotated
 * stress it was handed, so tau must co-rotate exactly, and the nonlocal field being fixed leaves L and
 * the nonlocal radius unchanged.
 */
void testRigidRotationIsExact()
{
  Wrapper w( stubProps.data(), int( stubProps.size() ), 1 );
  auto    state = freshState( w );

  // Preload with a deviatoric stretch so that the carried-over stress is non-trivial (and anisotropic).
  Tensor33d F1 = Spatial3D::I;
  F1( 0, 0 )   = 1.05;
  F1( 1, 1 )   = 0.97;
  F1( 0, 1 )   = 0.04;
  F1( 2, 1 )   = 0.02;

  const double N  = 0.2;
  const auto   s1 = step( w, state, F1, N );

  const double    theta = 30.0 * Constants::Pi / 180.0;
  const Tensor33d R     = rotationTensor( Eigen::Vector3d( 0.3, -0.7, 0.65 ), theta );
  const Tensor33d F2    = mul( R, F1 );

  const auto s2 = step( w, state, F2, N );

  const Tensor33d expectedTau = mul( mul( R, s1.tau ), Fastor::transpose( R ) );
  const double    scaleTau    = std::max( 1.0, Fastor::norm( s1.tau ) );
  throwExceptionOnFailure( checkIfEqual( s2.tau, expectedTau, 1e-10 * scaleTau ),
                           "rigid rotation did not co-rotate the Kirchhoff stress in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  throwExceptionOnFailure( checkIfEqual( s2.L, s1.L, 1e-10 * std::max( 1.0, std::abs( s1.L ) ) ),
                           "rigid rotation changed the local driving force in " + std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( checkIfEqual( s2.nonLocalRadius, s1.nonLocalRadius, 1e-12 ),
                           "rigid rotation changed the nonlocal radius in " + std::string( __PRETTY_FUNCTION__ ) );
}

/// At small strain the wrapper must reproduce the bare small-strain material driven along the same path.
void testSmallStrainAgreement()
{
  Wrapper w( stubProps.data(), int( stubProps.size() ), 1 );
  auto    state = freshState( w );

  IncrementalGradientStub direct( stubProps.data(), int( stubProps.size() ), 1 );
  std::vector< double >   directState( direct.getNumberOfRequiredStateVars(), 0.0 );

  IncrementalGradientStub::response directRes;
  IncrementalGradientStub::tangents directTan;
  directRes.stress      = Marmot::Vector6d::Zero();
  directRes.KLocal( 0 ) = 0.0;
  directRes.c( 0 )      = 0.0;
  directRes.stateVars   = directState.data();

  // 20 increments of a small symmetric strain, with the nonlocal field ramping up in step.
  Eigen::Matrix3d H = Eigen::Matrix3d::Zero();
  H( 0, 0 )         = 1e-6;
  H( 1, 1 )         = -4e-7;
  H( 0, 1 ) = H( 1, 0 ) = 7e-7;
  const double dN       = 5e-3;

  Eigen::Matrix3d Ftotal = Eigen::Matrix3d::Identity();
  double          Nrun   = 0.0;

  for ( int n = 0; n < 20; n++ ) {
    Ftotal += H;
    Nrun += dN;

    Tensor33d F;
    Marmot::mapEigenToFastor( F ) = Ftotal;
    step( w, state, F, Nrun );

    IncrementalGradientStub::increment inc;
    inc.dStrain = ContinuumMechanics::VoigtNotation::voigtFromStrainMatrix< 3 >( H );
    inc.K( 0 )  = Nrun;
    inc.dK( 0 ) = dN;
    inc.time    = double( n );
    inc.dT      = 1.0;

    direct.computeStress( directRes, directTan, inc );
  }

  Tensor33d Ffinal;
  Marmot::mapEigenToFastor( Ffinal ) = Ftotal;
  const auto s                       = step( w, state, Ffinal, Nrun );

  // Compare the Cauchy stress: tau = J sigma, and the J factor alone is larger than the tolerance below.
  const Tensor33d sigma = Tensor33d( s.tau / Ftotal.determinant() );

  Tensor33d expectedSigma;
  Marmot::mapEigenToFastor( expectedSigma ) = ContinuumMechanics::VoigtNotation::stressMatrixFromVoigt< 3 >(
    directRes.stress );

  // The mid-step objective rate and the naive additive strain differ at O(total strain), the same
  // O(1e-5) deviation documented in the small-strain HughesWingetWrapper test.
  throwExceptionOnFailure( checkIfEqual( sigma, expectedSigma, 1e-5 ),
                           "small-strain stress deviates from the wrapped material in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  const double scaleL = std::max( 1.0, std::abs( directRes.KLocal( 0 ) ) );
  throwExceptionOnFailure( checkIfEqual( s.L, directRes.KLocal( 0 ), 1e-8 * scaleL ),
                           "small-strain local driving force deviates from the wrapped material in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  throwExceptionOnFailure( checkIfEqual( s.nonLocalRadius, std::sqrt( directRes.c( 0 ) ), 1e-12 ),
                           "small-strain nonlocal radius deviates from the wrapped material in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

/**
 * @brief Analytic vs. numerical tangent, in a regime where the nonlocal field never leaves zero (no
 *        degradation), but the wrapper still carries a non-trivial stress forward across an increment
 *        with genuine shear and a rotational component.
 *
 * With IncrementalGradientStub's exact d(sigma^(n+1))/d(sigma_rot) = I, Analytic mode is now expected
 * to reproduce the numerical (forward-difference) tangent tightly on all four blocks.
 */
void testAnalyticTangentVsNumericalElastic()
{
  Wrapper          wAnalytic( stubProps.data(), int( stubProps.size() ), 1 );
  WrapperNumerical wNumerical( stubProps.data(), int( stubProps.size() ), 1 );

  auto sAnalytic  = freshState( wAnalytic );
  auto sNumerical = freshState( wNumerical );

  // Preload with a moderate stretch and shear, at N = 0 throughout.
  Tensor33d Fpre = Spatial3D::I;
  Fpre( 0, 0 )   = 1.03;
  Fpre( 1, 1 )   = 0.98;
  Fpre( 0, 1 )   = 0.02;
  Fpre( 1, 2 )   = 0.01;
  step( wAnalytic, sAnalytic, Fpre, 0.0 );
  step( wNumerical, sNumerical, Fpre, 0.0 );

  // A further increment with genuine shear and a rotational component, still at N = 0.
  Tensor33d Ffinal = Spatial3D::I;
  Ffinal( 0, 0 )   = 1.08;
  Ffinal( 1, 1 )   = 0.94;
  Ffinal( 0, 1 )   = 0.12;
  Ffinal( 1, 2 )   = 0.06;
  Ffinal( 2, 0 )   = 0.04;

  const auto sA = step( wAnalytic, sAnalytic, Ffinal, 0.0 );
  const auto sN = step( wNumerical, sNumerical, Ffinal, 0.0 );

  const double scaleTauF = std::max( 1.0, Fastor::norm( sN.tangents.dTau_dF ) );
  throwExceptionOnFailure( checkIfEqual( sA.tangents.dTau_dF, sN.tangents.dTau_dF, 1e-5 * scaleTauF ),
                           "dTau_dF deviates from the numerical tangent in the elastic regime in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  const double scaleTauN = std::max( 1.0, Fastor::norm( sN.tangents.dTau_dN ) );
  throwExceptionOnFailure( checkIfEqual( sA.tangents.dTau_dN, sN.tangents.dTau_dN, 1e-5 * scaleTauN ),
                           "dTau_dN deviates from the numerical tangent in the elastic regime in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  const double scaleLF = std::max( 1.0, Fastor::norm( sN.tangents.dL_dF ) );
  throwExceptionOnFailure( checkIfEqual( sA.tangents.dL_dF, sN.tangents.dL_dF, 1e-5 * scaleLF ),
                           "dL_dF deviates from the numerical tangent in the elastic regime in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  const double scaleLN = std::max( 1.0, std::abs( sN.tangents.dL_dN ) );
  throwExceptionOnFailure( checkIfEqual( sA.tangents.dL_dN, sN.tangents.dL_dN, 1e-5 * scaleLN ),
                           "dL_dN deviates from the numerical tangent in the elastic regime in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

/**
 * @brief All four tangent blocks of the Exact-tangent instance must agree with the Numerical-tangent
 *        instance, both driven along an identical history to a partially degraded, sheared state.
 */
void testExactTangentVsNumericalDamaged()
{
  WrapperExact     wExact( stubProps.data(), int( stubProps.size() ), 1 );
  WrapperNumerical wNumerical( stubProps.data(), int( stubProps.size() ), 1 );

  auto sExact     = freshState( wExact );
  auto sNumerical = freshState( wNumerical );

  // Drive both instances along an identical history, building up shear and a partially degraded
  // nonlocal field. The commit path is independent of the tangent mode, so both state-variable
  // arrays remain bit-identical throughout.
  struct HistoryPoint {
    double Fxx, Fyy, Fxy, Fyz, N;
  };
  const std::vector< HistoryPoint > history = { { 1.02, 0.99, 0.03, 0.00, 0.15 },
                                                { 1.05, 0.97, 0.06, 0.02, 0.30 },
                                                { 1.08, 0.95, 0.10, 0.05, 0.45 } };

  auto driveHistory = [&]( GEFiniteStrain& w, std::vector< double >& state ) {
    for ( const auto& h : history ) {
      Tensor33d F = Spatial3D::I;
      F( 0, 0 )   = h.Fxx;
      F( 1, 1 )   = h.Fyy;
      F( 0, 1 )   = h.Fxy;
      F( 1, 2 )   = h.Fyz;
      step( w, state, F, h.N );
    }
  };
  driveHistory( wExact, sExact );
  driveHistory( wNumerical, sNumerical );

  // A final increment with genuine shear -- dL_dF is reported against the engineering-shear Voigt
  // convention, so a pure diagonal stretch would never exercise the entries most likely to be wrong.
  Tensor33d Ffinal    = Spatial3D::I;
  Ffinal( 0, 0 )      = 1.10;
  Ffinal( 1, 1 )      = 0.94;
  Ffinal( 0, 1 )      = 0.15;
  Ffinal( 1, 2 )      = 0.08;
  Ffinal( 2, 0 )      = 0.05;
  const double Nfinal = 0.55;

  const auto sE = step( wExact, sExact, Ffinal, Nfinal );
  const auto sN = step( wNumerical, sNumerical, Ffinal, Nfinal );

  const double scaleTauF = std::max( 1.0, Fastor::norm( sN.tangents.dTau_dF ) );
  throwExceptionOnFailure( checkIfEqual( sE.tangents.dTau_dF, sN.tangents.dTau_dF, 1e-5 * scaleTauF ),
                           "dTau_dF deviates from the numerical tangent in " + std::string( __PRETTY_FUNCTION__ ) );

  const double scaleTauN = std::max( 1.0, Fastor::norm( sN.tangents.dTau_dN ) );
  throwExceptionOnFailure( checkIfEqual( sE.tangents.dTau_dN, sN.tangents.dTau_dN, 1e-5 * scaleTauN ),
                           "dTau_dN deviates from the numerical tangent in " + std::string( __PRETTY_FUNCTION__ ) );

  const double scaleLF = std::max( 1.0, Fastor::norm( sN.tangents.dL_dF ) );
  throwExceptionOnFailure( checkIfEqual( sE.tangents.dL_dF, sN.tangents.dL_dF, 1e-5 * scaleLF ),
                           "dL_dF deviates from the numerical tangent in " + std::string( __PRETTY_FUNCTION__ ) );

  const double scaleLN = std::max( 1.0, std::abs( sN.tangents.dL_dN ) );
  throwExceptionOnFailure( checkIfEqual( sE.tangents.dL_dN, sN.tangents.dL_dN, 1e-5 * scaleLN ),
                           "dL_dN deviates from the numerical tangent in " + std::string( __PRETTY_FUNCTION__ ) );
}

/**
 * @brief The wrapper must hand the wrapped material dK = N - N_n, the increment since the last
 *        accepted (committed) call -- not the total nonlocal field.
 */
void testNonlocalIncrementBookkeeping()
{
  const std::vector< double > dummyProps = { 0.0 };
  RecordingWrapper            w( dummyProps.data(), int( dummyProps.size() ), 1 );
  auto                        state = freshState( w );

  auto recordedIncrementAfter = [&]( double N ) {
    step( w, state, Spatial3D::I, N );
    return *w.getStateView( "recordedIncrement", state.data() ).stateLocation;
  };

  const double dK1 = recordedIncrementAfter( 0.3 );
  throwExceptionOnFailure( checkIfEqual( dK1, 0.3, 1e-14 ),
                           "first call did not see dK = N - 0 in " + std::string( __PRETTY_FUNCTION__ ) );

  const double dK2 = recordedIncrementAfter( 0.5 );
  throwExceptionOnFailure( checkIfEqual( dK2, 0.2, 1e-14 ),
                           "second call did not see dK = 0.5 - 0.3 in " + std::string( __PRETTY_FUNCTION__ ) );

  const double dK3 = recordedIncrementAfter( 0.5 );
  throwExceptionOnFailure( checkIfEqual( dK3, 0.0, 1e-14 ),
                           "repeating N did not see dK = 0 in " + std::string( __PRETTY_FUNCTION__ ) );
}

/**
 * @brief A zero eigen stress from an unstressed guess must converge immediately, to the identity.
 *
 * This is the commonest call there is -- geostatic initialisation of an unloaded region -- and it
 * used to THROW. The convergence test read
 * @c R.norm() / std::min( normalStress.norm(), 1.0 ) <= 1e-10, whose denominator is zero when the
 * stress is, so an exactly satisfied residual evaluated as 0/0 = NaN, NaN is not <= 1e-10, Newton
 * applied a zero correction five times over and gave up. Nothing caught it because this file is
 * the only caller in the tree and it did not cover the zero case.
 *
 * The second half is the other direction of the same mistake: min() capped the denominator at 1,
 * so a stress of order 100 was being held to 1e-12 relative. Asking for a non-zero eigen stress
 * exercises that branch.
 */
void testEigenDeformationAtZeroStressConverges()
{
  Wrapper w( stubProps.data(), int( stubProps.size() ), 1 );
  auto    state = freshState( w );

  const auto [F0, F1, F2] = w.findEigenDeformationForEigenStress( { 1.0, 1.0, 1.0 }, { 0.0, 0.0, 0.0 }, state.data() );

  throwExceptionOnFailure( checkIfEqual( F0, 1.0, 1e-12 ) && checkIfEqual( F1, 1.0, 1e-12 ) &&
                             checkIfEqual( F2, 1.0, 1e-12 ),
                           "a zero eigen stress did not return the identity deformation in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // And a non-zero target still converges, which is what the relative half of the scale is for.
  auto stateB             = freshState( w );
  const auto [G0, G1, G2] = w.findEigenDeformationForEigenStress( { 1.0, 1.0, 1.0 },
                                                                  { -10.0, -10.0, -10.0 },
                                                                  stateB.data() );

  throwExceptionOnFailure( G0 < 1.0 && checkIfEqual( G0, G1, 1e-12 ) && checkIfEqual( G1, G2, 1e-12 ),
                           "a hydrostatic compressive eigen stress did not give a uniform contraction in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

/**
 * @brief The wrapper must forward the wrapped material's micro-inertia, unchanged.
 *
 * This accessor is the whole reason the gradient-enhanced field can be made second order in time,
 * and a wrapped material reaches the solver only through it. Two ways to break it silently are the
 * base-state slice -- reading the wrapper's own state rather than the wrapped material's -- and the
 * first-vector entry, both of which would return a plausible-looking zero. Zero is also exactly
 * what "no micro-inertia" means, so the failure would look like a deck that simply did not ask for
 * the hyperbolic scheme.
 */
void testNonlocalMicroInertiaIsForwarded()
{
  // <= eta^2/4 for the base stub's eta of 1, so the value is an admissible one.
  constexpr double microInertia = 7.5e-9;

  /// A stub whose micro-inertia is a fixed, non-zero, recognisable number.
  class MicroInertiaStub : public IncrementalGradientStub {
  public:
    using IncrementalGradientStub::IncrementalGradientStub;

    std::vector< double > getNonlocalMicroInertia( const double* ) const override { return { 7.5e-9 }; }
  };

  GradientEnhancedHughesWingetWrapper< MicroInertiaStub > w( stubProps.data(), int( stubProps.size() ), 1 );
  auto                                                    state = freshState( w );

  throwExceptionOnFailure( checkIfEqual( w.getNonlocalMicroInertia( state.data() ), microInertia, 1e-20 ),
                           "the wrapper did not forward the wrapped material's micro-inertia in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // And a material that provides none still reports none -- the default this interface relies on.
  Wrapper plain( stubProps.data(), int( stubProps.size() ), 1 );
  auto    plainState = freshState( plain );

  throwExceptionOnFailure( plain.getNonlocalMicroInertia( plainState.data() ) == 0.0,
                           "a wrapped material providing no micro-inertia did not report zero in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

/**
 * @brief The energy/dissipation round trip must close at @f$ J \neq 1 @f$.
 *
 * This interface's densities are per unit REFERENCE volume; the wrapped small-strain material knows
 * only its own current volume, so the wrapper multiplies by @f$ J @f$ on the way out and must divide
 * by the previous @f$ J @f$ on the way in. Skipping the second half rescales the history once per
 * increment and compounds it -- invisible to every material in this tree, because none of them
 * touches these fields, and invisible to every other test here for the same reason.
 *
 * The stub below accumulates a FIXED amount per call, so after n increments at a steady volume
 * change the reported total must be exactly n times that amount times J, and nothing else.
 */
void testEnergyDensityRoundTripAtNonUnitJacobian()
{
  constexpr double perIncrement = 3.0;

  /// Accumulates a fixed dissipation per call, in its own (current-volume) density.
  class AccumulatingStub : public IncrementalGradientStub {
  public:
    using IncrementalGradientStub::IncrementalGradientStub;

    void computeStress( response& res, tangents& tan, const increment& inc ) const override
    {
      IncrementalGradientStub::computeStress( res, tan, inc );
      res.dissipation += perIncrement;
    }
  };

  GradientEnhancedHughesWingetWrapper< AccumulatingStub > w( stubProps.data(), int( stubProps.size() ), 1 );
  auto                                                    state = freshState( w );

  // A pure dilatation, held constant after the first increment, so J is the same on every entry.
  Tensor33d F = Spatial3D::I;
  F( 0, 0 ) = F( 1, 1 ) = F( 2, 2 ) = 1.1;
  const double J                    = 1.1 * 1.1 * 1.1;

  GEResponse r;
  r.stateVars = state.data();
  GETangents t;

  for ( int i = 1; i <= 3; ++i ) {
    w.computeStress( r, t, GEDeformation{ F, 0.0 }, GETimeIncrement{ double( i ), 1.0 } );

    // i increments of `perIncrement` in the wrapped material's own density, reported per reference
    // volume. Without the inverse conversion this grows by a further factor of J each time.
    throwExceptionOnFailure( checkIfEqual( r.dissipation, i * perIncrement * J, 1e-10 ),
                             MakeString() << __PRETTY_FUNCTION__ << ": after " << i << " increments the dissipation is "
                                          << r.dissipation << ", expected " << i * perIncrement * J );
  }
}

/**
 * @brief The eigen-deformation overloads must pre-stretch F, and the plane-strain ones must delegate.
 *
 * The base class' convenience surface -- computeStress and computePlaneStrain with an eigen
 * deformation, and their explicit counterparts -- is what a geostatic or thermally pre-strained
 * analysis reaches for, and none of it had a caller in the tree. An eigen deformation multiplies
 * the diagonal of F, so an element handed the identity with F0 must see exactly what an element
 * handed diag(F0) without one sees; that equality is what makes the overload worth having and is
 * what this checks.
 */
void testEigenDeformationAndPlaneStrainOverloads()
{
  constexpr double F0 = 1.05;

  Tensor33d stretched = Spatial3D::I;
  stretched( 0, 0 ) = stretched( 1, 1 ) = stretched( 2, 2 ) = F0;

  // Reference: the stretch applied directly, no eigen deformation.
  Wrapper    wRef( stubProps.data(), int( stubProps.size() ), 1 );
  auto       refState = freshState( wRef );
  const auto ref      = step( wRef, refState, stretched, 0.1 );

  // The same thing expressed as an eigen deformation on an undeformed element.
  Wrapper    wEig( stubProps.data(), int( stubProps.size() ), 1 );
  auto       eigState = freshState( wEig );
  GEResponse rEig;
  rEig.stateVars = eigState.data();
  GETangents tEig;
  /* Through a base reference, which is how the element holds it. The wrapper overrides the
   * four-argument computeStress and thereby HIDES the five-argument eigen-deformation overload for
   * anyone calling on the derived type -- harmless in practice, since every consumer works through
   * MarmotMaterialGradientEnhancedFiniteStrain*, but worth knowing.
   */
  GEFiniteStrain& baseEig = wEig;
  baseEig.computeStress( rEig,
                         tEig,
                         GEDeformation{ Spatial3D::I, 0.1 },
                         GETimeIncrement{ 0.0, 1.0 },
                         std::make_tuple( F0, F0, F0 ) );

  throwExceptionOnFailure( checkIfEqual( rEig.tau, ref.tau, 1e-12 ),
                           "an eigen deformation did not reproduce the directly stretched response in " +
                             std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( checkIfEqual( rEig.L, ref.L, 1e-12 ),
                           "the local driving force differs under an eigen deformation in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // dTau_dN is NOT chain-ruled: the nonlocal field is not deformed. dTau_dF is.
  throwExceptionOnFailure( checkIfEqual( tEig.dTau_dN, ref.tangents.dTau_dN, 1e-12 ),
                           "dTau_dN was scaled by the eigen deformation, which does not deform the "
                           "nonlocal field, in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // computePlaneStrain delegates to computeStress, in both overloads.
  Wrapper    wPs( stubProps.data(), int( stubProps.size() ), 1 );
  auto       psState = freshState( wPs );
  GEResponse rPs;
  rPs.stateVars = psState.data();
  GETangents tPs;
  wPs.computePlaneStrain( rPs, tPs, GEDeformation{ stretched, 0.1 }, GETimeIncrement{ 0.0, 1.0 } );

  throwExceptionOnFailure( checkIfEqual( rPs.tau, ref.tau, 1e-12 ),
                           "computePlaneStrain did not delegate to computeStress in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  Wrapper    wPsE( stubProps.data(), int( stubProps.size() ), 1 );
  auto       psEState = freshState( wPsE );
  GEResponse rPsE;
  rPsE.stateVars = psEState.data();
  GETangents tPsE;
  wPsE.computePlaneStrain( rPsE,
                           tPsE,
                           GEDeformation{ Spatial3D::I, 0.1 },
                           GETimeIncrement{ 0.0, 1.0 },
                           std::make_tuple( F0, F0, F0 ) );

  throwExceptionOnFailure( checkIfEqual( rPsE.tau, ref.tau, 1e-12 ),
                           "the eigen-deformation computePlaneStrain did not delegate in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // And the explicit counterparts, which discard the tangent.
  Wrapper    wEx( stubProps.data(), int( stubProps.size() ), 1 );
  auto       exState = freshState( wEx );
  GEResponse rEx;
  rEx.stateVars = exState.data();
  wEx.computeStressExplicit( rEx, GEDeformation{ stretched, 0.1 }, GETimeIncrement{ 0.0, 1.0 } );

  throwExceptionOnFailure( checkIfEqual( rEx.tau, ref.tau, 1e-12 ),
                           "computeStressExplicit disagrees with computeStress in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  Wrapper    wPsEx( stubProps.data(), int( stubProps.size() ), 1 );
  auto       psExState = freshState( wPsEx );
  GEResponse rPsEx;
  rPsEx.stateVars = psExState.data();
  wPsEx.computePlaneStrainExplicit( rPsEx, GEDeformation{ stretched, 0.1 }, GETimeIncrement{ 0.0, 1.0 } );

  throwExceptionOnFailure( checkIfEqual( rPsEx.tau, ref.tau, 1e-12 ),
                           "computePlaneStrainExplicit disagrees with computeStress in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // The value constructor, which nothing else builds.
  double           dummyState = 0.0;
  const GEResponse built( ref.tau, 1.5, 2.5, 3.5, 4.5, &dummyState );
  throwExceptionOnFailure( checkIfEqual( built.L, 1.5, 1e-15 ) && checkIfEqual( built.nonLocalRadius, 2.5, 1e-15 ) &&
                             checkIfEqual( built.elasticEnergyDensity, 3.5, 1e-15 ) &&
                             checkIfEqual( built.dissipation, 4.5, 1e-15 ) && built.stateVars == &dummyState,
                           "the ConstitutiveResponse value constructor did not store its arguments in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

/**
 * @brief The factory must round-trip a registered name, and name what it knows when it cannot.
 *
 * The factory is how GCDP reaches this wrapper -- it registers GCDP/HUGHES-WINGET from its own
 * repository and never constructs the wrapper directly, which is the opposite of what every other
 * test in this file does. A registration or lookup regression would therefore pass this suite while
 * breaking the integration the wrapper exists for.
 */
void testFactoryRoundTripAndUnknownName()
{
  using Factory = MarmotLibrary::MarmotMaterialGradientEnhancedFiniteStrainFactory;

  const std::string name = "TESTONLY/GRADIENT-ENHANCED-HUGHES-WINGET";

  const bool registered = Factory::registerMaterial< Wrapper >( name );
  throwExceptionOnFailure( registered,
                           "registerMaterial did not report success in " + std::string( __PRETTY_FUNCTION__ ) );

  std::unique_ptr< MarmotMaterialGradientEnhancedFiniteStrain > created(
    Factory::createMaterial( name, stubProps.data(), int( stubProps.size() ), 1 ) );

  throwExceptionOnFailure( created != nullptr,
                           "createMaterial returned nullptr for a registered name in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // Constructed through the factory it must behave as the directly constructed one does.
  throwExceptionOnFailure( created->getNumberOfRequiredStateVars() ==
                             Wrapper( stubProps.data(), int( stubProps.size() ), 1 ).getNumberOfRequiredStateVars(),
                           "the factory-created material has a different state layout in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  bool threw = false;
  try {
    Factory::createMaterial( "NO/SUCH-MATERIAL", stubProps.data(), int( stubProps.size() ), 1 );
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw,
                           "createMaterial accepted an unregistered name in " + std::string( __PRETTY_FUNCTION__ ) );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{ testStateLayout,
                                                       testUndeformed,
                                                       testEigenDeformationAtZeroStressConverges,
                                                       testNonlocalMicroInertiaIsForwarded,
                                                       testEnergyDensityRoundTripAtNonUnitJacobian,
                                                       testEigenDeformationAndPlaneStrainOverloads,
                                                       testFactoryRoundTripAndUnknownName,
                                                       testRigidRotationIsExact,
                                                       testSmallStrainAgreement,
                                                       testAnalyticTangentVsNumericalElastic,
                                                       testExactTangentVsNumericalDamaged,
                                                       testNonlocalIncrementBookkeeping };

  executeTestsAndCollectExceptions( tests );
  return 0;
}
