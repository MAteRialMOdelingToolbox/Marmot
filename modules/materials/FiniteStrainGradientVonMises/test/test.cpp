#include "Marmot/FiniteStrainGradientVonMises.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>
#include <cmath>
#include <functional>
#include <vector>

using namespace Marmot;
using namespace Marmot::Materials;
using namespace Marmot::Testing;
using namespace Marmot::FastorStandardTensors;

using Mat = MarmotMaterialGradientPlasticityFiniteStrain< 1 >;
using Res = Mat::response;
using Tan = Mat::tangents;
using Inc = Mat::increment;

namespace {

  FiniteStrainGradientVonMises makeMaterial( const std::vector< double >& props )
  {
    return FiniteStrainGradientVonMises( props.data(), static_cast< int >( props.size() ), 1 );
  }

  // Reference material properties: K, G, fy0, H, g, density, nonlocalViscosity
  std::vector< double > referenceProps()
  {
    return { 3500.0, 1500.0, 250.0, 1500.0, 25.0, 2000.0, 0.0 };
  }

  // Returns a freshly-initialized state vector (Fp = I, kappa = laplaceKappa = 0)
  std::vector< double > freshState( const FiniteStrainGradientVonMises& mat )
  {
    std::vector< double > stateVars( mat.getNumberOfRequiredStateVars(), 0.0 );
    // initializeYourself is non-const by design (it writes into the state array),
    // but does not depend on/modify the material's own members.
    const_cast< FiniteStrainGradientVonMises& >( mat ).initializeYourself( stateVars.data(),
                                                                            static_cast< int >( stateVars.size() ) );
    return stateVars;
  }

  struct EvalResult {
    Tensor33d tau;
    double    f;
    Tan       tan;
  };

  // Evaluate the material starting from a fixed "old" state (stateVarsOld is copied,
  // never mutated), so repeated calls with perturbed inputs are fully independent --
  // this is what makes the central-difference tangent check below clean.
  EvalResult evaluate( const FiniteStrainGradientVonMises& mat,
                       const std::vector< double >&        stateVarsOld,
                       const Tensor33d&                     F,
                       const double                         dLambda,
                       const double                         laplaceDLambda )
  {
    std::vector< double > stateVars = stateVarsOld;

    Res res;
    res.tau                  = Tensor33d( 0.0 );
    res.f                    = Eigen::Vector< double, 1 >::Zero();
    res.stateVars            = stateVars.data();
    res.elasticEnergyDensity = 0.0;
    res.dissipation          = 0.0;

    Tan tan;
    Inc inc;
    inc.deformation.F      = F;
    inc.dLambda( 0 )        = dLambda;
    inc.laplaceDLambda( 0 ) = laplaceDLambda;
    inc.time                = 0.0;
    inc.dT                  = 1.0;

    mat.computeStress( res, tan, inc );

    return { res.tau, res.f( 0 ), tan };
  }

  Tensor33d simpleShearF( const double gamma )
  {
    Tensor33d F = Spatial3D::I;
    F( 0, 1 )   = gamma;
    return F;
  }

} // namespace

void testFiniteStrainGradientVonMisesElasticStep()
{
  const auto mat        = makeMaterial( referenceProps() );
  const auto stateOld    = freshState( mat );

  // Small shear -> trial stress well below fy0=250 -> elastic point.
  const Tensor33d F = simpleShearF( 1.0e-4 );

  const auto result = evaluate( mat, stateOld, F, 0.0, 0.0 );

  throwExceptionOnFailure( std::isfinite( result.f ), "Non-finite FB residual in elastic step test" );
  throwExceptionOnFailure( std::isfinite( result.tau( 0, 1 ) ) && std::isfinite( result.tau( 1, 1 ) ),
                           "Non-finite Kirchhoff stress in elastic step test" );

  // Fischer-Burmeister residual must vanish (up to the 1e-12 regularization) when
  // dLambda = 0 and the trial state is truly elastic (a > 0, b = 0 => FB(a,0) = 0).
  throwExceptionOnFailure( checkIfEqual( result.f, 0.0, 1.0e-6 ),
                           "FB residual not consistent with dLambda = 0 for elastic trial state" );

  // Kirchhoff stress must be symmetric.
  for ( int i = 0; i < 3; ++i )
    for ( int j = 0; j < 3; ++j )
      throwExceptionOnFailure( checkIfEqual( result.tau( i, j ), result.tau( j, i ), 1.0e-8 ),
                               "Kirchhoff stress tensor not symmetric in elastic step test" );
}

void testFiniteStrainGradientVonMisesPlasticConsistency()
{
  const auto mat     = makeMaterial( referenceProps() );
  const auto stateOld = freshState( mat );

  // Large shear -> elastic trial stress well above fy0=250 -> yielding point.
  const Tensor33d F = simpleShearF( 0.15 );

  // At dLambda = 0 the trial state should register as violating the yield surface:
  // a = -f_tr < 0, b = 0 => FB(a,0) = -2*a = 2*f_tr > 0.
  const auto atZeroLambda = evaluate( mat, stateOld, F, 0.0, 0.0 );
  throwExceptionOnFailure( std::isfinite( atZeroLambda.f ), "Non-finite FB residual at dLambda = 0" );
  throwExceptionOnFailure( atZeroLambda.f > 0.0,
                           "Expected positive FB residual for a yielding trial state with dLambda = 0" );

  // A small positive dLambda should move the response into the plastic branch and
  // reduce the (unsigned) FB residual compared to the dLambda = 0 evaluation.
  const double dLambda = 1.0e-3;
  const auto   result  = evaluate( mat, stateOld, F, dLambda, 0.0 );

  throwExceptionOnFailure( std::isfinite( result.f ), "Non-finite FB residual in plastic consistency test" );
  throwExceptionOnFailure( std::abs( result.f ) < atZeroLambda.f,
                           "FB residual did not decrease when a positive dLambda was supplied for a yielding point" );

  for ( int i = 0; i < 3; ++i )
    for ( int j = 0; j < 3; ++j ) {
      throwExceptionOnFailure( std::isfinite( result.tau( i, j ) ), "Non-finite Kirchhoff stress entry in plastic test" );
      throwExceptionOnFailure( checkIfEqual( result.tau( i, j ), result.tau( j, i ), 1.0e-8 ),
                               "Kirchhoff stress tensor not symmetric in plastic consistency test" );
    }
}

namespace {

  // Central-difference check of tan.dTau_dF against computeStress() re-evaluated at
  // perturbed F, always starting from the SAME fixed old state so the comparison is clean.
  void checkTangentByCentralDifference( const FiniteStrainGradientVonMises& mat,
                                        const std::vector< double >&        stateOld,
                                        const Tensor33d&                     F,
                                        const double                         dLambda,
                                        const std::string&                   label )
  {
    const auto reference = evaluate( mat, stateOld, F, dLambda, 0.0 );

    // h chosen to balance central-difference truncation error (~h^2) against
    // floating-point cancellation error (~eps/h); 1e-5 sits well inside the
    // sweet spot for the O(0.1-1) stress/deformation magnitudes used here.
    const double h = 1.0e-5;

    double maxAbsRef = 0.0;
    for ( int i = 0; i < 3; ++i )
      for ( int j = 0; j < 3; ++j )
        maxAbsRef = std::max( maxAbsRef, std::abs( reference.tau( i, j ) ) );

    for ( int p = 0; p < 3; ++p ) {
      for ( int q = 0; q < 3; ++q ) {

        Tensor33d Fp = F;
        Fp( p, q )   = F( p, q ) + h;
        Tensor33d Fm = F;
        Fm( p, q )   = F( p, q ) - h;

        const auto resPlus  = evaluate( mat, stateOld, Fp, dLambda, 0.0 );
        const auto resMinus = evaluate( mat, stateOld, Fm, dLambda, 0.0 );

        for ( int i = 0; i < 3; ++i ) {
          for ( int j = 0; j < 3; ++j ) {
            const double dTau_dF_numeric = ( resPlus.tau( i, j ) - resMinus.tau( i, j ) ) / ( 2.0 * h );
            const double dTau_dF_analytic = reference.tan.dTau_dF( i, j, p, q );

            const double scale = std::max( { std::abs( dTau_dF_numeric ), std::abs( dTau_dF_analytic ), 1.0 } );
            const double relErr = std::abs( dTau_dF_numeric - dTau_dF_analytic ) / scale;

            throwExceptionOnFailure( relErr < 1.0e-6,
                                     MakeString() << label << ": mismatch in dTau_dF(" << i << "," << j << "," << p
                                                  << "," << q << ") analytic=" << dTau_dF_analytic
                                                  << " numeric=" << dTau_dF_numeric << " relErr=" << relErr );
          }
        }
      }
    }
  }

} // namespace

// Regression test: the Fischer-Burmeister residual must be evaluated at the UPDATED
// (returned) stress, not the trial stress. With a trial-based residual, plastic flow
// cannot reduce the residual at all -- and with softening (H < 0) the residual even grows
// with dLambda -- so the global Newton stalls at the onset of plasticity. Verified here by
// requiring that increasing dLambda monotonically reduces the FB residual for a yielding
// state of a SOFTENING material.
void testFiniteStrainGradientVonMisesSofteningConsistency()
{
  // K, G, fy0, H (negative = softening), g, density, nonlocalViscosity
  const auto mat      = makeMaterial( { 3500.0, 1500.0, 250.0, -400.0, 25.0, 2000.0, 0.0 } );
  const auto stateOld = freshState( mat );

  const Tensor33d F = simpleShearF( 0.15 ); // trial state well beyond yield

  const auto atZero  = evaluate( mat, stateOld, F, 0.0, 0.0 );
  const auto atSmall = evaluate( mat, stateOld, F, 1.0e-3, 0.0 );
  const auto atLarge = evaluate( mat, stateOld, F, 5.0e-3, 0.0 );

  throwExceptionOnFailure( atZero.f > 0.0, "Expected positive FB residual at dLambda = 0 for a yielding state" );
  throwExceptionOnFailure( std::abs( atSmall.f ) < atZero.f,
                           "FB residual did not decrease with dLambda for a softening material "
                           "(residual evaluated at the trial instead of the updated stress?)" );
  throwExceptionOnFailure( std::abs( atLarge.f ) < atZero.f,
                           "FB residual did not decrease further with larger dLambda for a softening material" );
}

// Regression test: at the virgin state (F = I exactly, zero stress) the von Mises norm
// sqrt(3 J2) is non-differentiable and its unguarded dual gradient is 0/0 = NaN, which
// poisoned df_dF / dFddLambda / dFddLaplacian (and thus the global stiffness of every
// virgin material point) before the guard was added in FiniteStrainGradientVonMises.cpp.
void testFiniteStrainGradientVonMisesVirginStateTangentsFinite()
{
  const auto mat      = makeMaterial( referenceProps() );
  const auto stateOld = freshState( mat );

  const Tensor33d F      = Spatial3D::I;
  const auto      result = evaluate( mat, stateOld, F, 0.0, 0.0 );

  for ( int i = 0; i < 9; ++i ) {
    throwExceptionOnFailure( std::isfinite( result.tan.df_dF( 0, i ) ),
                             MakeString() << "Non-finite df_dF(0," << i << ") at virgin state (F = I)" );
    throwExceptionOnFailure( std::isfinite( result.tan.dTau_ddLambda( i, 0 ) ),
                             MakeString() << "Non-finite dTau_ddLambda(" << i << ",0) at virgin state (F = I)" );
  }
  throwExceptionOnFailure( std::isfinite( result.tan.dFddLambda( 0, 0 ) ),
                           "Non-finite dFddLambda at virgin state (F = I)" );
  throwExceptionOnFailure( std::isfinite( result.tan.dFddLaplacian( 0, 0 ) ),
                           "Non-finite dFddLaplacian at virgin state (F = I)" );
  for ( int i = 0; i < 3; ++i )
    for ( int j = 0; j < 3; ++j )
      for ( int k = 0; k < 3; ++k )
        for ( int l = 0; l < 3; ++l )
          throwExceptionOnFailure( std::isfinite( result.tan.dTau_dF( i, j, k, l ) ),
                                   "Non-finite dTau_dF at virgin state (F = I)" );
}

void testFiniteStrainGradientVonMisesTangentElastic()
{
  const auto mat     = makeMaterial( referenceProps() );
  const auto stateOld = freshState( mat );

  const Tensor33d F = simpleShearF( 1.0e-4 );
  checkTangentByCentralDifference( mat, stateOld, F, 0.0, "elastic" );
}

void testFiniteStrainGradientVonMisesTangentPlastic()
{
  const auto mat     = makeMaterial( referenceProps() );
  const auto stateOld = freshState( mat );

  const Tensor33d F = simpleShearF( 0.15 );
  checkTangentByCentralDifference( mat, stateOld, F, 1.0e-3, "plastic" );
}

// Central-difference check of the complementarity-residual couplings df_dF and dFddLambda
// (the K_LU / K_LL blocks of the global Newton matrix) in the plastic branch.
void testFiniteStrainGradientVonMisesYieldResidualTangents()
{
  const auto mat      = makeMaterial( referenceProps() );
  const auto stateOld = freshState( mat );

  const Tensor33d F       = simpleShearF( 0.15 );
  const double    dLambda = 1.0e-3;
  const double    h       = 1.0e-6;

  const auto reference = evaluate( mat, stateOld, F, dLambda, 0.0 );

  // df / dF_pq
  for ( int p = 0; p < 3; ++p )
    for ( int q = 0; q < 3; ++q ) {
      Tensor33d Fp = F;
      Fp( p, q ) += h;
      Tensor33d Fm = F;
      Fm( p, q ) -= h;

      const double numeric  = ( evaluate( mat, stateOld, Fp, dLambda, 0.0 ).f -
                                evaluate( mat, stateOld, Fm, dLambda, 0.0 ).f ) /
                              ( 2.0 * h );
      const double analytic = reference.tan.df_dF( 0, 3 * p + q );
      const double scale    = std::max( { std::abs( numeric ), std::abs( analytic ), 1.0 } );
      throwExceptionOnFailure( std::abs( numeric - analytic ) / scale < 1.0e-5,
                               MakeString() << "mismatch in df_dF(0," << 3 * p + q << ") analytic=" << analytic
                                            << " numeric=" << numeric );
    }

  // df / d dLambda
  {
    const double numeric  = ( evaluate( mat, stateOld, F, dLambda + h, 0.0 ).f -
                              evaluate( mat, stateOld, F, dLambda - h, 0.0 ).f ) /
                            ( 2.0 * h );
    const double analytic = reference.tan.dFddLambda( 0, 0 );
    const double scale    = std::max( { std::abs( numeric ), std::abs( analytic ), 1.0 } );
    throwExceptionOnFailure( std::abs( numeric - analytic ) / scale < 1.0e-5,
                             MakeString() << "mismatch in dFddLambda analytic=" << analytic
                                          << " numeric=" << numeric );
  }

  // df / d laplaceDLambda
  {
    const double numeric  = ( evaluate( mat, stateOld, F, dLambda, h ).f -
                              evaluate( mat, stateOld, F, dLambda, -h ).f ) /
                            ( 2.0 * h );
    const double analytic = reference.tan.dFddLaplacian( 0, 0 );
    const double scale    = std::max( { std::abs( numeric ), std::abs( analytic ), 1.0 } );
    throwExceptionOnFailure( std::abs( numeric - analytic ) / scale < 1.0e-5,
                             MakeString() << "mismatch in dFddLaplacian analytic=" << analytic
                                          << " numeric=" << numeric );
  }
}

int main()
{
  const std::vector< std::function< void() > > tests = {
    testFiniteStrainGradientVonMisesElasticStep,
    testFiniteStrainGradientVonMisesPlasticConsistency,
    testFiniteStrainGradientVonMisesSofteningConsistency,
    testFiniteStrainGradientVonMisesVirginStateTangentsFinite,
    testFiniteStrainGradientVonMisesTangentElastic,
    testFiniteStrainGradientVonMisesTangentPlastic,
    testFiniteStrainGradientVonMisesYieldResidualTangents,
  };

  executeTestsAndCollectExceptions( tests );
  return 0;
}
