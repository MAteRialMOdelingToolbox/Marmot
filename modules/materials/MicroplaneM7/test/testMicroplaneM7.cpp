#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotMaterialPointSolverHypoElastic.h"
#include "Marmot/MarmotNumericalDifferentiation.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include "memory"
#include <algorithm>
#include <array>

using namespace Marmot::Testing;
using namespace Marmot::Solvers;

std::unique_ptr< MarmotMaterialHypoElastic > createMarmotMaterialHypoElastic( std::string   materialName,
                                                                              const double* materialProperties,
                                                                              const int     nMaterialProperties )
{

  const int                                    elLabel = 1;
  std::unique_ptr< MarmotMaterialHypoElastic > mat(
    MarmotLibrary::MarmotMaterialHypoElasticFactory::createMaterial( materialName,
                                                                     materialProperties,
                                                                     nMaterialProperties,
                                                                     elLabel ) );
  return mat;
}

void testOedometricCompression()
{
  // initialize material parameters
  const std::array< double, 7 > materialProperties  = { 2.5e4, 0.18, 1.5e-4, 110, 30, 100, 1e-4 };
  const int                     nMaterialProperties = 7;

  // create the material
  auto mat = createMarmotMaterialHypoElastic( "MICROPLANEM7", materialProperties.data(), nMaterialProperties );

  // set increment
  MarmotMaterialHypoElastic::state3D state;
  Marmot::Vector6d                   dStrain = { -1e-3, 0., 0., 0., 0., 0. };
  Marmot::Matrix6d                   dStressDDStrain;
  const int                          nStateVars = mat->getNumberOfRequiredStateVars();
  std::vector< double >              stateVars( nStateVars, 0.0 );
  state.stateVars      = stateVars.data();
  double       timeOld = 0.0;
  const double dT      = 0.1;

  // compute material response
  for ( int i = 0; i < 10; i++ ) {
    MarmotMaterialHypoElastic::timeInfo timeInfo{ timeOld + dT, dT };
    mat->computeStress( state, dStressDDStrain, dStrain, timeInfo );
    timeOld += dT;
  }

  std::cout.precision( 10 );
  std::cout << "oedometricCompression stress = " << state.stress.transpose() << std::endl;

  // set the correct result
  Marmot::Vector6d stressTarget = { -198.4872406, -59.19080234, -59.19080234, 0., 0., 0. };

  throwExceptionOnFailure( checkIfEqual< double >( state.stress, stressTarget, 1e-6 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testTriaxialCompression()
{
  // initialize material parameters
  const std::array< double, 7 > materialProperties  = { 2.5e4, 0.18, 1.5e-4, 110, 30, 100, 1e-4 };
  const int                     nMaterialProperties = 7;

  // create the material
  auto mat = createMarmotMaterialHypoElastic( "MICROPLANEM7", materialProperties.data(), nMaterialProperties );

  // set increment
  MarmotMaterialHypoElastic::state3D state;
  Marmot::Vector6d                   dStrain = { -1e-3, -2e-4, -2e-4, 0., 0., 0. };
  Marmot::Matrix6d                   dStressDDStrain;
  const int                          nStateVars = mat->getNumberOfRequiredStateVars();
  std::vector< double >              stateVars( nStateVars, 0.0 );
  state.stateVars      = stateVars.data();
  double       timeOld = 0.0;
  const double dT      = 0.05;

  // compute material response
  for ( int i = 0; i < 20; i++ ) {
    MarmotMaterialHypoElastic::timeInfo timeInfo{ timeOld + dT, dT };
    mat->computeStress( state, dStressDDStrain, dStrain, timeInfo );
    timeOld += dT;
  }

  std::cout << "triaxialCompression stress = " << state.stress.transpose() << std::endl;

  // set the correct result
  Marmot::Vector6d stressTarget = { -337.6493166, -231.1318392, -231.1318392, 0., 0., 0. };

  throwExceptionOnFailure( checkIfEqual< double >( state.stress, stressTarget, 1e-6 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testShear()
{
  // initialize material parameters
  const std::array< double, 7 > materialProperties  = { 2.5e4, 0.18, 1.5e-4, 110, 30, 100, 1e-4 };
  const int                     nMaterialProperties = 7;

  // create the material
  auto mat = createMarmotMaterialHypoElastic( "MICROPLANEM7", materialProperties.data(), nMaterialProperties );

  // set increment
  MarmotMaterialHypoElastic::state3D state;
  Marmot::Vector6d                   dStrain = { 0., 0., 0., 5e-4, -3e-4, 2e-4 };
  Marmot::Matrix6d                   dStressDDStrain;
  const int                          nStateVars = mat->getNumberOfRequiredStateVars();
  std::vector< double >              stateVars( nStateVars, 0.0 );
  state.stateVars      = stateVars.data();
  double       timeOld = 0.0;
  const double dT      = 0.1;

  // compute material response
  for ( int i = 0; i < 10; i++ ) {
    MarmotMaterialHypoElastic::timeInfo timeInfo{ timeOld + dT, dT };
    mat->computeStress( state, dStressDDStrain, dStrain, timeInfo );
    timeOld += dT;
  }

  std::cout << "shear stress = " << state.stress.transpose() << std::endl;

  // set the correct result
  Marmot::Vector6d stressTarget = { -27.1304588, -25.03769381, -17.82531926, 25.04060546, -18.02493041, 16.15387587 };

  throwExceptionOnFailure( checkIfEqual< double >( state.stress, stressTarget, 1e-6 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testTangent()
{
  // initialize material parameters
  const std::array< double, 7 > materialProperties  = { 2.5e4, 0.18, 1.5e-4, 110, 30, 100, 1e-4 };
  const int                     nMaterialProperties = 7;

  // create the material
  auto mat = createMarmotMaterialHypoElastic( "MICROPLANEM7", materialProperties.data(), nMaterialProperties );

  // set increment
  MarmotMaterialHypoElastic::state3D state;
  Marmot::Vector6d                   dStrain = { 1e-5, 1e-5, 1e-5, 0., 0., 0. };
  Marmot::Matrix6d                   dStressDDStrain;
  const int                          nStateVars = mat->getNumberOfRequiredStateVars();
  std::vector< double >              stateVars( nStateVars, 0.0 );
  state.stateVars      = stateVars.data();
  double       timeOld = 0.0;
  const double dT      = 0.1;

  Marmot::Matrix6d Cel, Ctarget;
  Ctarget = Marmot::ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( materialProperties[0],
                                                                                materialProperties[1] );

  // compute material response
  for ( int i = 0; i < 5; i++ ) {
    MarmotMaterialHypoElastic::timeInfo timeInfo{ timeOld + dT, dT };
    mat->computeStress( state, dStressDDStrain, dStrain, timeInfo );
    timeOld += dT;
    if ( i == 0 )
      Cel = dStressDDStrain;
  }

  std::vector< double > stateVarsOld( nStateVars );
  stateVarsOld               = stateVars;
  Marmot::Vector6d stressOld = state.stress;

  Eigen::MatrixXd C;

  const Eigen::VectorXd dEps = dStrain;
  // ----------------------------------------
  C = Marmot::NumericalAlgorithms::Differentiation::forwardDifference(
    [&]( const Eigen::VectorXd& dE_ ) {
      // reset stateVars to old state
      stateVars = stateVarsOld;

      MarmotMaterialHypoElastic::state3D localState;
      localState.stress    = stressOld;
      localState.stateVars = stateVars.data();

      MarmotMaterialHypoElastic::timeInfo timeInfo{ timeOld + dT, dT };
      // compute stress
      mat->computeStress( localState, dStressDDStrain, dE_, timeInfo );

      return localState.stress;
    },
    dEps );
  // ----------------------------------------

  stateVars    = stateVarsOld;
  state.stress = stressOld;
  MarmotMaterialHypoElastic::timeInfo timeInfo{ timeOld + dT, dT };
  mat->computeStress( state, dStressDDStrain, dStrain, timeInfo );

  throwExceptionOnFailure( checkIfEqual< double >( Cel, Ctarget, 1e-7 ),
                           MakeString() << __PRETTY_FUNCTION__ << " elastic tangent failed" );

  throwExceptionOnFailure( checkIfEqual< double >( C, dStressDDStrain, 0.5 ),
                           MakeString() << __PRETTY_FUNCTION__ << " plastic tangent failed" );
}

void testOedometricExtension()
{
  // initialize material parameters
  const std::array< double, 7 > materialProperties  = { 2.5e4, 0.18, 1.5e-4, 110, 30, 100, 1e-4 };
  const int                     nMaterialProperties = 7;

  // create the material
  auto mat = createMarmotMaterialHypoElastic( "MICROPLANEM7", materialProperties.data(), nMaterialProperties );

  // set increment
  MarmotMaterialHypoElastic::state3D state;
  Marmot::Vector6d                   dStrain = { 1e-4, 0., 0., 0., 0., 0. };
  Marmot::Matrix6d                   dStressDDStrain;
  const int                          nStateVars = mat->getNumberOfRequiredStateVars();
  std::vector< double >              stateVars( nStateVars, 0.0 );
  state.stateVars      = stateVars.data();
  double       timeOld = 0.0;
  const double dT      = 0.1;

  // compute material response
  for ( int i = 0; i < 10; i++ ) {
    MarmotMaterialHypoElastic::timeInfo timeInfo{ timeOld + dT, dT };
    mat->computeStress( state, dStressDDStrain, dStrain, timeInfo );
    timeOld += dT;
  }

  std::cout << "oedometricExtension stress = " << state.stress.transpose() << std::endl;

  // set the correct result
  Marmot::Vector6d stressTarget = { 2.105194144, 2.569889942, 2.569889942, 0., 0., 0. };

  throwExceptionOnFailure( checkIfEqual< double >( state.stress, stressTarget, 1e-6 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testGeneralTriaxialCompression()
{
  // initialize material parameters
  const std::array< double, 7 > materialProperties  = { 2.5e4, 0.18, 1.5e-4, 110, 30, 100, 1e-4 };
  const int                     nMaterialProperties = 7;

  // create the material
  auto mat = createMarmotMaterialHypoElastic( "MICROPLANEM7", materialProperties.data(), nMaterialProperties );

  // set increment: three distinct principal normal strains (non-axisymmetric)
  MarmotMaterialHypoElastic::state3D state;
  Marmot::Vector6d                   dStrain = { -1e-3, -6e-4, -2e-4, 0., 0., 0. };
  Marmot::Matrix6d                   dStressDDStrain;
  const int                          nStateVars = mat->getNumberOfRequiredStateVars();
  std::vector< double >              stateVars( nStateVars, 0.0 );
  state.stateVars      = stateVars.data();
  double       timeOld = 0.0;
  const double dT      = 0.05;

  // compute material response
  for ( int i = 0; i < 20; i++ ) {
    MarmotMaterialHypoElastic::timeInfo timeInfo{ timeOld + dT, dT };
    mat->computeStress( state, dStressDDStrain, dStrain, timeInfo );
    timeOld += dT;
  }

  std::cout << "generalTriaxialCompression stress = " << state.stress.transpose() << std::endl;

  // set the correct result
  Marmot::Vector6d stressTarget = { -378.1410099, -334.6771595, -274.3581159, 0., 0., 0. };

  throwExceptionOnFailure( checkIfEqual< double >( state.stress, stressTarget, 1e-6 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testCombinedNormalShear()
{
  // initialize material parameters
  const std::array< double, 7 > materialProperties  = { 2.5e4, 0.18, 1.5e-4, 110, 30, 100, 1e-4 };
  const int                     nMaterialProperties = 7;

  // create the material
  auto mat = createMarmotMaterialHypoElastic( "MICROPLANEM7", materialProperties.data(), nMaterialProperties );

  // set increment: combined normal compression and shear in the same increment
  MarmotMaterialHypoElastic::state3D state;
  Marmot::Vector6d                   dStrain = { -5e-4, -2e-4, -2e-4, 3e-4, -2e-4, 1e-4 };
  Marmot::Matrix6d                   dStressDDStrain;
  const int                          nStateVars = mat->getNumberOfRequiredStateVars();
  std::vector< double >              stateVars( nStateVars, 0.0 );
  state.stateVars      = stateVars.data();
  double       timeOld = 0.0;
  const double dT      = 0.1;

  // compute material response
  for ( int i = 0; i < 10; i++ ) {
    MarmotMaterialHypoElastic::timeInfo timeInfo{ timeOld + dT, dT };
    mat->computeStress( state, dStressDDStrain, dStrain, timeInfo );
    timeOld += dT;
  }

  std::cout << "combinedNormalShear stress = " << state.stress.transpose() << std::endl;

  // set the correct result
  Marmot::Vector6d stressTarget = { -159.0734483, -100.9032682, -102.2851853, 28.43759883, -19.00010476, 9.570438297 };

  throwExceptionOnFailure( checkIfEqual< double >( state.stress, stressTarget, 1e-6 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

// Shared tangent-verification routine: drives the material through `nIncrements` steps of
// `dStrain` (the same path used by the corresponding golden-value stress test, so the
// boundary/transition code actually engaged there is what gets checked here), then compares
// the analytical tangent returned by computeStress against a forward-difference numerical
// tangent for one more increment of the same size. Used to guard the NormalResponse
// struct-unification refactor (mp_response2 -> NormalResponse, projection tensor now applied
// inside normal_elasticLaw/tensile_boundary): a mistake there would show up as the analytical
// tangent diverging from the numerical one.
void checkTangentForStrainPath( const std::string&      testName,
                                const Marmot::Vector6d& dStrain,
                                const int               nIncrements,
                                const double            tangentTol )
{
  const std::array< double, 7 > materialProperties  = { 2.5e4, 0.18, 1.5e-4, 110, 30, 100, 1e-4 };
  const int                     nMaterialProperties = 7;

  auto mat = createMarmotMaterialHypoElastic( "MICROPLANEM7", materialProperties.data(), nMaterialProperties );

  MarmotMaterialHypoElastic::state3D state;
  Marmot::Matrix6d                   dStressDDStrain;
  const int                          nStateVars = mat->getNumberOfRequiredStateVars();
  std::vector< double >              stateVars( nStateVars, 0.0 );
  state.stateVars      = stateVars.data();
  double       timeOld = 0.0;
  const double dT      = 0.1;

  for ( int i = 0; i < nIncrements; i++ ) {
    MarmotMaterialHypoElastic::timeInfo timeInfo{ timeOld + dT, dT };
    mat->computeStress( state, dStressDDStrain, dStrain, timeInfo );
    timeOld += dT;
  }

  std::vector< double > stateVarsOld( nStateVars );
  stateVarsOld               = stateVars;
  Marmot::Vector6d stressOld = state.stress;

  Eigen::MatrixXd       C;
  const Eigen::VectorXd dEps = dStrain;
  C                          = Marmot::NumericalAlgorithms::Differentiation::forwardDifference(
    [&]( const Eigen::VectorXd& dE_ ) {
      stateVars = stateVarsOld;

      MarmotMaterialHypoElastic::state3D localState;
      localState.stress    = stressOld;
      localState.stateVars = stateVars.data();

      MarmotMaterialHypoElastic::timeInfo timeInfo{ timeOld + dT, dT };
      mat->computeStress( localState, dStressDDStrain, dE_, timeInfo );

      return localState.stress;
    },
    dEps );

  stateVars    = stateVarsOld;
  state.stress = stressOld;
  MarmotMaterialHypoElastic::timeInfo timeInfo{ timeOld + dT, dT };
  mat->computeStress( state, dStressDDStrain, dStrain, timeInfo );

  throwExceptionOnFailure( checkIfEqual< double >( C, dStressDDStrain, tangentTol ),
                           MakeString() << testName << " tangent failed" );
}

void testTangentOedometricExtension()
{
  // Small strain increment (same order of magnitude as testTangent()) so the tensile boundary
  // is only just engaged on some microplanes, avoiding the strongly non-smooth shear-boundary
  // branch switching (shear_lawCaner's R27/R28 selection on sign(sig_N)) seen at larger strains.
  Marmot::Vector6d dStrain = { 1e-5, 0., 0., 0., 0., 0. };
  checkTangentForStrainPath( __PRETTY_FUNCTION__, dStrain, 5, 0.5 );
}

void testTangentCombinedNormalShear()
{
  Marmot::Vector6d dStrain = { -2e-5, -1e-5, -1e-5, 1.5e-5, -1e-5, 0.5e-5 };
  checkTangentForStrainPath( __PRETTY_FUNCTION__, dStrain, 5, 0.5 );
}

// Coordinate invariance ("spinTurbokreisel"): the same strain path is applied in 100 rotated
// coordinate systems, and stress and tangent are transformed back and compared with the reference
// frame. M7 is isotropic in theory. The discrete N37 microplane quadrature is only approximately
// isotropic, though: stress and tangent are exactly invariant only in the elastic range, and once the
// boundaries are active they differ by the quadrature error (about 1 % on this compression-shear path,
// much more in tension softening). So the tolerances are loose. The testCubicSymmetryInvariance*()
// tests below are the sharp check.
void testCoordinateInvarianceTurbokreisel()
{
  std::vector< double > materialProperties = { 2.5e4, 0.18, 1.5e-4, 110, 30, 100, 1e-4 };
  auto                  solveropts         = MarmotMaterialPointSolverHypoElastic::SolverOptions();
  std::string           matName            = "MICROPLANEM7";

  auto solver = MarmotMaterialPointSolverHypoElastic( matName,
                                                      materialProperties.data(),
                                                      materialProperties.size(),
                                                      solveropts );

  MarmotMaterialPointSolverHypoElastic::Step step;
  step.isStrainComponentControlled = { true, true, true, true, true, true };
  step.isStressComponentControlled = { false, false, false, false, false, false };
  step.stressIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  // same total strain as testCombinedNormalShear(), beyond the elastic range
  step.strainIncrementTarget = { -5e-3, -2e-3, -2e-3, 3e-3, -2e-3, 1e-3 };
  // fixed increment size, so that every rotated frame follows the same incremental path
  step.dTStart = 0.1;
  step.dTMax   = 0.1;

  solver.addStep( step );

  // Measured on this path: up to 5.3 % component-wise stress deviation (0.8 % in norm). The tangent
  // deviates by about 9 % in norm, but its small off-diagonal entries differ by up to a factor of ~100,
  // and spinTurbokreisel checks every entry relative to itself. So the tangent tolerance only guards
  // against NaN/inf and gross errors; the tangent is checked sharply in testCubicSymmetryInvariance*().
  throwExceptionOnFailure( spinTurbokreisel( solver, 1e-1, 1e3 ),
                           MakeString() << __PRETTY_FUNCTION__ << " Turbokreisel failed!" );
}

// The 24 proper rotations of the cube (signed permutation matrices with det = +1). The N37 quadrature
// has full cubic symmetry, so these rotations map the set of microplanes onto itself.
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

// Drives a fresh material point through nIncrements constant strain increments and returns the final
// stress and tangent.
std::pair< Marmot::Vector6d, Marmot::Matrix6d > runStrainPath( const Marmot::Vector6d& dStrain, const int nIncrements )
{
  const std::array< double, 7 > materialProperties  = { 2.5e4, 0.18, 1.5e-4, 110, 30, 100, 1e-4 };
  const int                     nMaterialProperties = 7;

  auto mat = createMarmotMaterialHypoElastic( "MICROPLANEM7", materialProperties.data(), nMaterialProperties );

  MarmotMaterialHypoElastic::state3D state;
  Marmot::Matrix6d                   dStressDDStrain = Marmot::Matrix6d::Zero();
  std::vector< double >              stateVars( mat->getNumberOfRequiredStateVars(), 0.0 );
  state.stateVars      = stateVars.data();
  double       timeOld = 0.0;
  const double dT      = 0.1;

  for ( int i = 0; i < nIncrements; i++ ) {
    MarmotMaterialHypoElastic::timeInfo timeInfo{ timeOld + dT, dT };
    mat->computeStress( state, dStressDDStrain, dStrain, timeInfo );
    timeOld += dT;
  }

  return { state.stress, dStressDDStrain };
}

// Coordinate invariance under the symmetry group of the quadrature. A rotation from the cube group
// only permutes the microplanes and turns their in-plane shear directions, so a correct implementation
// must give the same response up to round-off in every regime, including softening. This catches
// errors in the projection tensors, in the treatment of the two in-plane shear directions and in the
// homogenization, which the loose Turbokreisel tolerances above would hide.
void checkCubicSymmetryInvarianceForStrainPath( const std::string&      testName,
                                                const Marmot::Vector6d& dStrain,
                                                const int               nIncrements )
{
  using namespace Marmot::ContinuumMechanics::VoigtNotation::Transformations;

  const auto [refStress, refStiffness] = runStrainPath( dStrain, nIncrements );

  for ( const auto& Q : cubicSymmetryRotations() ) {
    const auto [stress, stiffness] = runStrainPath( transformStrainToLocalSystem( dStrain, Q ), nIncrements );

    throwExceptionOnFailure( checkIfEqual< double >( transformStressToGlobalSystem( stress, Q ), refStress, 1e-10 ),
                             MakeString() << testName << " stress failed" );
    throwExceptionOnFailure( checkIfEqual< double >( transformStiffnessToGlobalSystem( stiffness, Q ),
                                                     refStiffness,
                                                     1e-8 ),
                             MakeString() << testName << " tangent failed" );
  }
}

void testCubicSymmetryInvarianceCombinedNormalShear()
{
  // same strain path as testCombinedNormalShear()
  checkCubicSymmetryInvarianceForStrainPath( __PRETTY_FUNCTION__, { -5e-4, -2e-4, -2e-4, 3e-4, -2e-4, 1e-4 }, 10 );
}

void testCubicSymmetryInvarianceTensionShear()
{
  // tension-dominated path, well beyond the peak (tensile softening)
  checkCubicSymmetryInvarianceForStrainPath( __PRETTY_FUNCTION__, { 3e-5, 6e-6, -3e-6, 1.5e-5, -9e-6, 6e-6 }, 100 );
}

int main()
{
  testOedometricCompression();
  testTriaxialCompression();
  testShear();
  testTangent();
  testOedometricExtension();
  testGeneralTriaxialCompression();
  testCombinedNormalShear();
  testTangentOedometricExtension();
  testTangentCombinedNormalShear();
  testCoordinateInvarianceTurbokreisel();
  testCubicSymmetryInvarianceCombinedNormalShear();
  testCubicSymmetryInvarianceTensionShear();
  return 0;
}
