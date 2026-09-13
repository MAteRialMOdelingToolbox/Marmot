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
#include "Marmot/MarmotMaterialGradientEnhancedHughesWinget.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
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

int main()
{
  auto tests = std::vector< std::function< void() > >{ testStateLayout,
                                                       testUndeformed,
                                                       testRigidRotationIsExact,
                                                       testSmallStrainAgreement,
                                                       testAnalyticTangentVsNumericalElastic,
                                                       testExactTangentVsNumericalDamaged,
                                                       testNonlocalIncrementBookkeeping };

  executeTestsAndCollectExceptions( tests );
  return 0;
}
