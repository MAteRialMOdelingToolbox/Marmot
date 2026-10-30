#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotMaterialGeneralGradientEnhancedHypoElastic.h"
#include "Marmot/MarmotMaterialGeneralGradientEnhancedHypoElasticFactory.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotNumericalDifferentiation.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include "memory"
#include <algorithm>
#include <array>
#include <stdexcept>
#include <string>

using namespace Marmot::Testing;
using response  = typename MarmotMaterialGeneralGradientEnhancedHypoElastic< 6 >::response;
using tangents  = typename MarmotMaterialGeneralGradientEnhancedHypoElastic< 6 >::tangents;
using increment = typename MarmotMaterialGeneralGradientEnhancedHypoElastic< 6 >::increment;

namespace {

  // {E, nu, k1, k2, k3, k4, k5, c, m}; c3/fc/transitionWidth/transitionExponent fall back to their defaults
  const std::array< double, 9 > defaultMaterialProperties = { 2.5e4, 0.18, 1.5e-4, 110, 30, 100, 1e-4, 25.0, 1.0 };

  std::unique_ptr< MarmotMaterialGeneralGradientEnhancedHypoElastic< 6 > > createGM7( const double* materialProperties,
                                                                                      const int nMaterialProperties )
  {
    const int                                                                elLabel = 1;
    std::unique_ptr< MarmotMaterialGeneralGradientEnhancedHypoElastic< 6 > > mat(
      MarmotLibrary::MarmotMaterialGeneralGradientEnhancedHypoElasticFactory< 6 >::createMaterial( "GM7",
                                                                                                   materialProperties,
                                                                                                   nMaterialProperties,
                                                                                                   elLabel ) );
    return mat;
  }

  // Runs a strain path with a fixed local strain increment over nIncrements steps, using GM7's staggered
  // (Picard-style) update of the nonlocal driving variable to a single material point (nonlocal = local),
  // and returns the final macroscopic stress.
  Marmot::Vector6d runStaggeredStrainPath( const Marmot::Vector6d& dStrain, int nIncrements )
  {
    auto mat = createGM7( defaultMaterialProperties.data(), defaultMaterialProperties.size() );

    response  res;
    tangents  tan;
    increment inc;

    res.stress  = { 0., 0., 0., 0., 0., 0. };
    inc.dStrain = dStrain;
    inc.K       = Eigen::Vector< double, 6 >::Zero();
    inc.dK      = inc.dStrain; // nonlocal strain increment equals local increment (single material point)

    const int             nStateVars = mat->getNumberOfRequiredStateVars();
    std::vector< double > stateVars( nStateVars, 0.0 );
    res.stateVars = stateVars.data();

    for ( int i = 0; i < nIncrements; i++ ) {
      // small increments needed for reasonable results (staggered approach)
      mat->computeStress( res, tan, inc );
      inc.K = res.KLocal; // explicit (staggered) update of nonlocal variable (nonlocal = local if only one qp)
    }

    return res.stress;
  }

  // Verifies the analytical dStressddStrain against a forwardDifference numerical tangent, after settling
  // the material onto a strain path over nIncrements staggered steps.
  void checkTangentForStrainPath( const std::string&      testName,
                                  const Marmot::Vector6d& dStrain,
                                  int                     nIncrements,
                                  double                  tangentTol )
  {
    auto mat = createGM7( defaultMaterialProperties.data(), defaultMaterialProperties.size() );

    response  res;
    tangents  tan;
    increment inc;

    res.stress  = { 0., 0., 0., 0., 0., 0. };
    inc.dStrain = dStrain;
    inc.K       = Eigen::Vector< double, 6 >::Zero();
    inc.dK      = inc.dStrain;

    const int             nStateVars = mat->getNumberOfRequiredStateVars();
    std::vector< double > stateVars( nStateVars, 0.0 );
    res.stateVars = stateVars.data();

    for ( int i = 0; i < nIncrements; i++ ) {
      mat->computeStress( res, tan, inc );
      inc.K = res.KLocal;
    }

    std::vector< double > stateVarsOld( nStateVars );
    stateVarsOld    = stateVars;
    response resOld = res;

    Eigen::MatrixXd C;

    const Marmot::Vector6d dEps = inc.dStrain;
    C                           = Marmot::NumericalAlgorithms::Differentiation::forwardDifference(
      [&]( const Marmot::Vector6d& dE_ ) {
        // reset stateVars to old state
        stateVars  = stateVarsOld;
        response r = resOld;

        increment i = inc;
        i.dStrain   = dE_;

        // compute stress
        mat->computeStress( r, tan, i );

        return r.stress;
      },
      dEps );

    stateVars = stateVarsOld;
    res       = resOld;
    mat->computeStress( res, tan, inc );

    throwExceptionOnFailure( checkIfEqual< double >( C, tan.dStressddStrain, tangentTol ),
                             MakeString() << testName << " tangent failed" );
  }

  // Coordinate invariance checks. Marmot::Testing::spinTurbokreisel only accepts a
  // MarmotMaterialPointSolverHypoElastic, so GM7 gets its own version of the same check: the staggered
  // strain path of runStaggeredStrainPath is applied in rotated coordinate systems (the rotated
  // increment is used for both the local and the nonlocal strain). Stress, dStressddStrain,
  // dStressddK and KLocal are then transformed back to the reference frame and compared.
  struct GM7Result {
    Marmot::Vector6d stress;
    Marmot::Matrix6d dStressddStrain;
    Marmot::Matrix6d dStressddK;
    Marmot::Vector6d KLocal;
  };

  GM7Result runStaggeredStrainPathWithTangents( const Marmot::Vector6d& dStrain, int nIncrements )
  {
    auto mat = createGM7( defaultMaterialProperties.data(), defaultMaterialProperties.size() );

    response  res;
    tangents  tan;
    increment inc;

    res.stress  = { 0., 0., 0., 0., 0., 0. };
    inc.dStrain = dStrain;
    inc.K       = Eigen::Vector< double, 6 >::Zero();
    inc.dK      = inc.dStrain;

    const int             nStateVars = mat->getNumberOfRequiredStateVars();
    std::vector< double > stateVars( nStateVars, 0.0 );
    res.stateVars = stateVars.data();

    for ( int i = 0; i < nIncrements; i++ ) {
      mat->computeStress( res, tan, inc );
      inc.K = res.KLocal;
    }

    return { res.stress, tan.dStressddStrain, tan.dStressddK, res.KLocal };
  }

  // The same rotations as Marmot::Testing::spinTurbokreisel: 100 coordinate systems from a
  // Fibonacci lattice on the hemisphere.
  std::vector< Eigen::Matrix3d > turbokreiselRotations()
  {
    using namespace Eigen;
    constexpr int                N   = 100;
    const Matrix< double, N, 2 > pts = fibonacciLatticeHemisphere< N >();

    std::vector< Matrix3d > rotations;
    for ( int i = 0; i < N; i++ )
      rotations.push_back(
        ( AngleAxisd( pts( i, 1 ), Vector3d::UnitX() ) * AngleAxisd( pts( i, 0 ), Vector3d::UnitY() ) )
          .matrix()
          .transpose() );
    return rotations;
  }

  // The 24 proper rotations of the cube (signed permutation matrices with det = +1). The N37
  // quadrature has full cubic symmetry, so these rotations map the set of microplanes onto itself.
  std::vector< Eigen::Matrix3d > cubicSymmetryRotations()
  {
    std::vector< Eigen::Matrix3d > rotations;
    std::array< int, 3 >           perm = { 0, 1, 2 };
    do {
      for ( int signs = 0; signs < 8; signs++ ) {
        Eigen::Matrix3d Q = Eigen::Matrix3d::Zero();
        for ( int i = 0; i < 3; i++ )
          Q( i, perm[i] ) = ( signs >> i & 1 ) ? -1. : 1.;
        if ( Q.determinant() > 0 )
          rotations.push_back( Q );
      }
    } while ( std::next_permutation( perm.begin(), perm.end() ) );
    return rotations;
  }

  void checkCoordinateInvarianceForStrainPath( const std::string&                    testName,
                                               const Marmot::Vector6d&               dStrain,
                                               int                                   nIncrements,
                                               const std::vector< Eigen::Matrix3d >& rotations,
                                               double                                stressTol,
                                               double                                stiffnessTol )
  {
    using namespace Marmot::ContinuumMechanics::VoigtNotation::Transformations;

    const GM7Result ref = runStaggeredStrainPathWithTangents( dStrain, nIncrements );

    for ( const auto& Q : rotations ) {
      const GM7Result rot = runStaggeredStrainPathWithTangents( transformStrainToLocalSystem( dStrain, Q ),
                                                                nIncrements );

      throwExceptionOnFailure( checkIfEqual< double >( transformStressToGlobalSystem( rot.stress, Q ),
                                                       ref.stress,
                                                       stressTol ),
                               MakeString() << testName << " stress failed" );
      throwExceptionOnFailure( checkIfEqual< double >( transformStiffnessToGlobalSystem( rot.dStressddStrain, Q ),
                                                       ref.dStressddStrain,
                                                       stiffnessTol ),
                               MakeString() << testName << " dStressddStrain failed" );
      throwExceptionOnFailure( checkIfEqual< double >( transformStiffnessToGlobalSystem( rot.dStressddK, Q ),
                                                       ref.dStressddK,
                                                       stiffnessTol ),
                               MakeString() << testName << " dStressddK failed" );
      // KLocal is the accumulated local strain, a purely kinematic quantity, so it must match exactly
      throwExceptionOnFailure( checkIfEqual< double >( transformStrainToGlobalSystem( rot.KLocal, Q ),
                                                       ref.KLocal,
                                                       1e-12 ),
                               MakeString() << testName << " KLocal failed" );
    }
  }

} // namespace

void testOedometricCompression()
{
  const Marmot::Vector6d stress = runStaggeredStrainPath( { -1e-5, 0., 0., 0., 0., 0. }, 100 );

  std::cout.precision( 10 );
  std::cout << "oedometricCompression stress = " << stress.transpose() << std::endl;

  Marmot::Vector6d stressTarget = { -27.01211634, -6.219820639, -6.219820639, 0., 0., 0. };

  throwExceptionOnFailure( checkIfEqual< double >( stress, stressTarget, 1e-6 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testOedometricExtension()
{
  const Marmot::Vector6d stress = runStaggeredStrainPath( { 1e-5, 0., 0., 0., 0., 0. }, 100 );

  std::cout << "oedometricExtension stress = " << stress.transpose() << std::endl;

  Marmot::Vector6d stressTarget = { 2.105194143, 2.569889942, 2.569889942, 0., 0., 0. };

  throwExceptionOnFailure( checkIfEqual< double >( stress, stressTarget, 1e-6 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testTriaxialCompression()
{
  const Marmot::Vector6d stress = runStaggeredStrainPath( { -1e-5, -2e-6, -2e-6, 0., 0., 0. }, 100 );

  std::cout << "triaxialCompression stress = " << stress.transpose() << std::endl;

  Marmot::Vector6d stressTarget = { -29.90014834, -12.77615032, -12.77615032, 0., 0., 0. };

  throwExceptionOnFailure( checkIfEqual< double >( stress, stressTarget, 1e-6 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testGeneralTriaxialCompression()
{
  const Marmot::Vector6d stress = runStaggeredStrainPath( { -1e-5, -6e-6, -2e-6, 0., 0., 0. }, 100 );

  std::cout << "generalTriaxialCompression stress = " << stress.transpose() << std::endl;

  Marmot::Vector6d stressTarget = { -32.44738362, -23.86011106, -15.27283851, 0., 0., 0. };

  throwExceptionOnFailure( checkIfEqual< double >( stress, stressTarget, 1e-6 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testShear()
{
  const Marmot::Vector6d stress = runStaggeredStrainPath( { 0., 0., 0., 5e-6, -3e-6, 2e-6 }, 100 );

  std::cout << "shear stress = " << stress.transpose() << std::endl;

  Marmot::Vector6d stressTarget =
    { -1.024015932, -0.8022645829, -0.2246639328, 3.759950781, -2.437000706, 1.969265099 };

  throwExceptionOnFailure( checkIfEqual< double >( stress, stressTarget, 1e-6 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testCombinedNormalShear()
{
  const Marmot::Vector6d stress = runStaggeredStrainPath( { -5e-6, -2e-6, -2e-6, 3e-6, -2e-6, 1e-6 }, 100 );

  std::cout << "combinedNormalShear stress = " << stress.transpose() << std::endl;

  Marmot::Vector6d stressTarget = { -16.08907703, -9.691132544, -9.691132544, 3.198972243, -2.132648162, 1.066324081 };

  throwExceptionOnFailure( checkIfEqual< double >( stress, stressTarget, 1e-6 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

// With m = 0 the over-nonlocal strain equals the local strain, so GM7 must reproduce MicroplaneM7: same stress,
// and dStressddStrain + dStressddK equal to the MicroplaneM7 tangent. The strain path is given as a list of
// (increment, number of steps) segments, so that load reversals are covered too.
void checkEquivalenceWithMicroplaneM7( const std::string&                                       testName,
                                       const std::vector< std::pair< Marmot::Vector6d, int > >& segments )
{
  // {E, nu, k1, k2, k3, k4, k5, c, m = 0}, and the same properties without c and m for MicroplaneM7
  const std::array< double, 9 > propertiesGM7 = { 2.5e4, 0.18, 1.5e-4, 110, 30, 100, 1e-4, 25.0, 0.0 };
  const std::array< double, 7 > propertiesM7  = { 2.5e4, 0.18, 1.5e-4, 110, 30, 100, 1e-4 };

  auto                                         gm7 = createGM7( propertiesGM7.data(), propertiesGM7.size() );
  std::unique_ptr< MarmotMaterialHypoElastic > m7(
    MarmotLibrary::MarmotMaterialHypoElasticFactory::createMaterial( "MICROPLANEM7",
                                                                     propertiesM7.data(),
                                                                     propertiesM7.size(),
                                                                     1 ) );

  std::vector< double > stateVarsGM7( gm7->getNumberOfRequiredStateVars(), 0.0 );
  std::vector< double > stateVarsM7( m7->getNumberOfRequiredStateVars(), 0.0 );

  response  res;
  tangents  tan;
  increment inc;
  res.stress    = { 0., 0., 0., 0., 0., 0. };
  res.stateVars = stateVarsGM7.data();
  inc.K         = Eigen::Vector< double, 6 >::Zero();

  MarmotMaterialHypoElastic::state3D state;
  state.stress    = Marmot::Vector6d::Zero();
  state.stateVars = stateVarsM7.data();
  MarmotMaterialHypoElastic::timeInfo timeInfo{ 0.0, 1.0 };
  Marmot::Matrix6d                    tangentM7;

  const auto relativeDifference = []( const auto& a, const auto& b ) {
    return ( a - b ).norm() / std::max( b.norm(), 1.0 );
  };

  int step = 0;
  for ( const auto& [dStrain, nSteps] : segments ) {
    for ( int i = 0; i < nSteps; i++, step++ ) {
      inc.dStrain = dStrain;
      inc.dK      = dStrain;
      gm7->computeStress( res, tan, inc );
      inc.K = res.KLocal;

      m7->computeStress( state, tangentM7, dStrain, timeInfo );
      timeInfo.time += timeInfo.dT;

      const Marmot::Vector6d stressGM7  = res.stress;
      const Marmot::Matrix6d tangentGM7 = tan.dStressddStrain + tan.dStressddK;

      throwExceptionOnFailure( relativeDifference( stressGM7, state.stress ) < 1e-10 &&
                                 relativeDifference( tangentGM7, tangentM7 ) < 1e-10,
                               MakeString() << testName << ": GM7 (m = 0) differs from MicroplaneM7 in step " << step
                                            << "\nstress GM7: " << stressGM7.transpose()
                                            << "\nstress M7:  " << state.stress.transpose() << "\ntangent GM7:\n"
                                            << tangentGM7 << "\ntangent M7:\n"
                                            << tangentM7 );
    }
  }
}

void testEquivalenceWithMicroplaneM7()
{
  using V = Marmot::Vector6d;
  checkEquivalenceWithMicroplaneM7( "oedometricCompression", { { V{ -1e-5, 0., 0., 0., 0., 0. }, 100 } } );
  checkEquivalenceWithMicroplaneM7( "triaxialCompression", { { V{ -1e-5, -6e-6, -2e-6, 0., 0., 0. }, 100 } } );
  checkEquivalenceWithMicroplaneM7( "shear", { { V{ 0., 0., 0., 5e-6, -3e-6, 2e-6 }, 100 } } );
  checkEquivalenceWithMicroplaneM7( "combinedNormalShear", { { V{ -5e-6, -2e-6, -2e-6, 3e-6, -2e-6, 1e-6 }, 100 } } );
  // tension beyond the peak, then unloading into compression: covers the unloading transition and the sign
  // changes of the microplane normal stresses
  checkEquivalenceWithMicroplaneM7( "tensionUnloadingCompression",
                                    { { V{ 1e-5, 0., 0., 0., 0., 0. }, 100 },
                                      { V{ -1e-5, 0., 0., 0., 0., 0. }, 200 } } );
}

void testTangent()
{
  // initialize material parameters
  auto mat = createGM7( defaultMaterialProperties.data(), defaultMaterialProperties.size() );

  response  res;
  tangents  tan;
  increment inc;

  res.stress  = { 0., 0., 0., 0., 0., 0. };
  inc.dStrain = { 1e-6, 1e-6, 1e-6, 0., 0., 0. };
  inc.K       = Eigen::Vector< double, 6 >::Zero();
  inc.dK      = inc.dStrain; // nonlocal strain increment equals local increment (single material point)

  const int             nStateVars = mat->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVars, 0.0 );
  res.stateVars = stateVars.data();

  Marmot::Matrix6d Cel, Ctarget;
  Ctarget = Marmot::ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( defaultMaterialProperties[0],
                                                                                defaultMaterialProperties[1] );

  // compute material response
  for ( int i = 0; i < 5; i++ ) {
    mat->computeStress( res, tan, inc );
    inc.K = res.KLocal;
    if ( i == 0 )
      Cel = tan.dStressddStrain;
  }

  throwExceptionOnFailure( checkIfEqual< double >( Cel, Ctarget, 1e-7 ),
                           MakeString() << __PRETTY_FUNCTION__ << " elastic tangent failed" );

  checkTangentForStrainPath( __PRETTY_FUNCTION__, { 1e-6, 1e-6, 1e-6, 0., 0., 0. }, 5, 0.5 );
}

void testTangentOedometricExtension()
{
  checkTangentForStrainPath( __PRETTY_FUNCTION__, { 1e-6, 0., 0., 0., 0., 0. }, 5, 0.5 );
}

void testTangentCombinedNormalShear()
{
  checkTangentForStrainPath( __PRETTY_FUNCTION__, { -5e-7, -2e-7, -2e-7, 3e-7, -2e-7, 1e-7 }, 5, 0.5 );
}

// Coordinate invariance ("spinTurbokreisel"). GM7 is isotropic in theory, but the discrete N37
// microplane quadrature is only approximately isotropic: the response is exactly invariant only in the
// elastic range, and once the boundaries are active it differs by the quadrature error. So the
// tolerances are loose. Measured on this path: up to 5.3 % component-wise stress deviation (0.8 % in
// norm). The tangents deviate by about 9 % in norm, but small off-diagonal entries differ by up to a
// factor of ~100, and checkIfEqual compares every entry relative to itself. So the tangent tolerance
// only guards against NaN/inf and gross errors; the testCubicSymmetryInvariance*() tests below are
// the sharp check.
void testCoordinateInvarianceTurbokreisel()
{
  // same total strain as the Turbokreisel test of MicroplaneM7, beyond the elastic range
  checkCoordinateInvarianceForStrainPath( __PRETTY_FUNCTION__,
                                          { -5e-4, -2e-4, -2e-4, 3e-4, -2e-4, 1e-4 },
                                          10,
                                          turbokreiselRotations(),
                                          1e-1,
                                          1e3 );
}

// Coordinate invariance under the symmetry group of the quadrature. A rotation from the cube group
// only permutes the microplanes and turns their in-plane shear directions, so a correct implementation
// must give the same response up to round-off in every regime, including softening. This catches
// errors in the projection tensors, in the treatment of the two in-plane shear directions and in the
// homogenization, which the loose Turbokreisel tolerances above would hide.
void testCubicSymmetryInvarianceCombinedNormalShear()
{
  // direction of testCombinedNormalShear(), but 10x the strain to go beyond the elastic range
  checkCoordinateInvarianceForStrainPath( __PRETTY_FUNCTION__,
                                          { -5e-5, -2e-5, -2e-5, 3e-5, -2e-5, 1e-5 },
                                          100,
                                          cubicSymmetryRotations(),
                                          1e-10,
                                          1e-8 );
}

void testCubicSymmetryInvarianceTensionShear()
{
  // tension-dominated path, well beyond the peak (tensile softening)
  checkCoordinateInvarianceForStrainPath( __PRETTY_FUNCTION__,
                                          { 3e-5, 6e-6, -3e-6, 1.5e-5, -9e-6, 6e-6 },
                                          100,
                                          cubicSymmetryRotations(),
                                          1e-10,
                                          1e-8 );
}

// Fewer than the 9 required material properties (E, nu, k1 ... k5, c, m) must be rejected with
// std::invalid_argument from the property check, before c and m are read past the end of the property array.
void testTooFewMaterialProperties()
{
  const std::array< double, 8 > materialProperties = { 2.5e4, 0.18, 1.5e-4, 110, 30, 100, 1e-4, 25.0 };

  bool rejected = false;
  try {
    createGM7( materialProperties.data(), materialProperties.size() );
  }
  catch ( const std::invalid_argument& e ) {
    rejected = std::string( e.what() ).find( "material properties required" ) != std::string::npos;
  }

  throwExceptionOnFailure( rejected, MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

int main()
{
  testOedometricCompression();
  testOedometricExtension();
  testTriaxialCompression();
  testGeneralTriaxialCompression();
  testShear();
  testCombinedNormalShear();
  testEquivalenceWithMicroplaneM7();
  testTangent();
  testTangentOedometricExtension();
  testTangentCombinedNormalShear();
  testCoordinateInvarianceTurbokreisel();
  testCubicSymmetryInvarianceCombinedNormalShear();
  testCubicSymmetryInvarianceTensionShear();
  testTooFewMaterialProperties();
  return 0;
}
