#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotMaterialPointSolverHypoElastic.h"
#include "Marmot/MarmotNumericalDifferentiation.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include "Marmot/RankineSofteningPlasticity.h"
#include <Eigen/Eigenvalues>
#include <stdexcept>

using namespace Marmot;
using namespace Testing;
using namespace Materials;
using namespace Materials;

namespace {

  // E, nu, fy0, HLin, deltaYieldStress, delta (linear softening)
  std::vector< double > linearSofteningProperties = { 30000., 0.2, 3., -3000., 0., 0. };
  // combined linear and exponential softening
  std::vector< double > expSofteningProperties = { 30000., 0.2, 3., -500., -1., 5000. };
  // perfect plasticity
  std::vector< double > perfectlyPlasticProperties = { 30000., 0.2, 3., 0., 0., 0. };

  const double E      = 30000.;
  const double nu     = 0.2;
  const double fy0    = 3.;
  const double H      = -3000.;
  const double G      = E / ( 2. * ( 1. + nu ) );
  const double lambda = E * nu / ( ( 1. + nu ) * ( 1. - 2. * nu ) );

  struct MaterialPointState {
    Math::Vector6d  stress    = Math::Vector6d::Zero();
    Eigen::Vector2d stateVars = Eigen::Vector2d::Zero(); // kappa, acYield

    double kappa() const { return stateVars( 0 ); }
    double acYield() const { return stateVars( 1 ); }
  };

  // single stress update from a given state, returning the updated state and the algorithmic tangent
  std::pair< MaterialPointState, Math::Matrix6d > update( const RankineSofteningPlasticity& material,
                                                          const MaterialPointState&         stateOld,
                                                          const Math::Vector6d&             dStrain )
  {
    MaterialPointState           stateNew = stateOld;
    MaterialHypoElastic::state3D state;
    state.stress    = stateNew.stress;
    state.stateVars = stateNew.stateVars.data();
    Math::Matrix6d                      tangent;
    const MaterialHypoElastic::timeInfo tInfo{ 0.0, 1.0 };
    material.computeStress( state, tangent, dStrain, tInfo );
    stateNew.stress = state.stress;
    return { stateNew, tangent };
  }

  // true if f throws an exception of type Ex
  template < typename Ex, typename F >
  bool throws( F&& f )
  {
    try {
      f();
    }
    catch ( const Ex& ) {
      return true;
    }
    return false;
  }

  // compare the algorithmic tangent against a central difference approximation
  void checkTangent( const RankineSofteningPlasticity& material,
                     const MaterialPointState&         stateOld,
                     const Math::Vector6d&             dStrain )
  {
    using namespace Math::Differentiation::Numerical;
    const Math::Matrix6d tangent = update( material, stateOld, dStrain ).second;

    // differentiate wrt the scaled strain increment, since the default step is relative to max(1, |x|)
    const double          strainScale = 1e-4;
    const Eigen::MatrixXd numericalTangent =
      centralDifference( [&]( const Eigen::VectorXd& y )
                           -> Eigen::VectorXd { return update( material, stateOld, strainScale * y ).first.stress; },
                         dStrain / strainScale ) /
      strainScale;

    // relative to the tangent, floored by the elastic scale since the tangent vanishes at the apex for perfect
    // plasticity
    const double relError = ( tangent - numericalTangent ).cwiseAbs().maxCoeff() /
                            std::max( tangent.cwiseAbs().maxCoeff(), 1e-3 * E );
    throwExceptionOnFailure( relError < 1e-6,
                             Utilities::MakeString() << "tangent check failed, relative error " << relError );
  }

  // check that the returned stress satisfies all three Rankine yield conditions
  void checkAdmissibility( const MaterialPointState& state, const std::vector< double >& p )
  {
    const double fy = p[2] + p[3] * state.kappa() + p[4] * ( 1. - std::exp( -p[5] * state.kappa() ) );
    const Eigen::SelfAdjointEigenSolver< Eigen::Matrix3d > eigenSolver(
      Marmot::ContinuumMechanics::Voigt::voigtToStress( state.stress ) );
    throwExceptionOnFailure( eigenSolver.eigenvalues().maxCoeff() <= fy + 1e-10, "returned stress is not admissible" );
    throwExceptionOnFailure( state.kappa() >= 0., "negative hardening variable" );
  }

} // namespace

void testUniaxialStrain()
{
  const RankineSofteningPlasticity material( linearSofteningProperties.data(), 6, 1 );

  const double         eps = 2e-4;
  const Math::Vector6d dE  = ( Math::Vector6d() << eps, 0, 0, 0, 0, 0 ).finished();
  const auto [state, _]    = update( material, MaterialPointState(), dE );

  // single active surface in direction 1
  const double         dLambda      = ( ( lambda + 2. * G ) * eps - fy0 ) / ( lambda + 2. * G + H );
  const Math::Vector6d stressTarget = ( Math::Vector6d() << fy0 + H * dLambda,
                                        lambda * ( eps - dLambda ),
                                        lambda * ( eps - dLambda ),
                                        0,
                                        0,
                                        0 )
                                        .finished();

  throwExceptionOnFailure( checkIfEqual( state.kappa(), dLambda, 1e-15 ), "uniaxial strain: kappa mismatch" );
  throwExceptionOnFailure( checkIfEqual< double >( state.stress, stressTarget, 1e-11 ),
                           "uniaxial strain: comparison with analytical solution failed" );
  checkTangent( material, MaterialPointState(), dE );
}

void testEqualBiaxialStrain()
{
  const RankineSofteningPlasticity material( linearSofteningProperties.data(), 6, 1 );

  const double         eps = 1e-4;
  const Math::Vector6d dE  = ( Math::Vector6d() << eps, eps, 0, 0, 0, 0 ).finished();
  const auto [state, _]    = update( material, MaterialPointState(), dE );

  // return to the corner of surfaces 1 and 2 with equal consistency parameters
  const double         sigTrial     = ( 2. * lambda + 2. * G ) * eps;
  const double         dLambda      = ( sigTrial - fy0 ) / ( 2. * lambda + 2. * G + 2. * H );
  const double         sig          = fy0 + 2. * H * dLambda;
  const Math::Vector6d stressTarget = ( Math::Vector6d() << sig, sig, lambda * ( 2. * eps - 2. * dLambda ), 0, 0, 0 )
                                        .finished();

  throwExceptionOnFailure( checkIfEqual( state.kappa(), 2. * dLambda, 1e-15 ), "biaxial strain: kappa mismatch" );
  throwExceptionOnFailure( checkIfEqual< double >( state.stress, stressTarget, 1e-11 ),
                           "biaxial strain: comparison with analytical solution failed" );
  checkTangent( material, MaterialPointState(), dE );
}

void testIsotropicTension()
{
  const RankineSofteningPlasticity material( linearSofteningProperties.data(), 6, 1 );

  // isotropic strain increments, applied twice to continue softening at the apex
  const double         eps = 1e-4;
  const Math::Vector6d dE  = ( Math::Vector6d() << eps, eps, eps, 0, 0, 0 ).finished();

  MaterialPointState state;
  for ( int i = 1; i <= 2; i++ ) {
    checkTangent( material, state, dE );
    state = update( material, state, dE ).first;

    // return to the apex, i.e., the intersection of all three surfaces with equal consistency parameters:
    // sigma = 3K ( eps - dLambda ) = fy0 + 3 H dLambda
    const double         threeK       = 3. * lambda + 2. * G;
    const double         dLambda      = ( threeK * i * eps - fy0 ) / ( threeK + 3. * H );
    const double         sig          = fy0 + 3. * H * dLambda;
    const Math::Vector6d stressTarget = ( Math::Vector6d() << sig, sig, sig, 0, 0, 0 ).finished();

    throwExceptionOnFailure( checkIfEqual( state.kappa(), 3. * dLambda, 1e-15 ), "isotropic tension: kappa mismatch" );
    throwExceptionOnFailure( checkIfEqual< double >( state.stress, stressTarget, 1e-11 ),
                             "isotropic tension: comparison with analytical solution failed" );
  }
}

void testBiaxialTension()
{
  std::string matName = "RANKINESOFTENINGPLASTICITY";
  auto        solver  = Materials::MaterialPointSolverHypoElastic( matName,
                                                           linearSofteningProperties.data(),
                                                           static_cast< int >( linearSofteningProperties.size() ),
                                                           MaterialPointSolverHypoElastic::SolverOptions() );

  const double eps = 2e-4;

  // equal biaxial tension under plane stress
  MaterialPointSolverHypoElastic::Step step;
  step.isStrainComponentControlled = { true, true, false, false, false, false };
  step.isStressComponentControlled = { false, false, true, true, true, true };
  step.stressIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step.strainIncrementTarget       = { eps, eps, 0.0, 0.0, 0.0, 0.0 };
  step.dTStart                     = 0.25;
  solver.addStep( step );
  solver.solve();

  // return to the intersection of surfaces 1 and 2:
  // eps - dLambda = (1 - nu) sigma / E, sigma = fy0 + 2 H dLambda
  const double         dLambda      = ( eps - ( 1. - nu ) * fy0 / E ) / ( 1. + 2. * H * ( 1. - nu ) / E );
  const double         sig          = fy0 + 2. * H * dLambda;
  const Math::Vector6d stressTarget = ( Math::Vector6d() << sig, sig, 0, 0, 0, 0 ).finished();

  const auto result = solver.getHistory().back();
  throwExceptionOnFailure( checkIfEqual( result.stateVars( 0 ), 2. * dLambda, 1e-14 ),
                           "biaxial tension: kappa mismatch" );
  throwExceptionOnFailure( checkIfEqual< double >( result.stress, stressTarget, 1e-9 ),
                           "biaxial tension: comparison with analytical solution failed" );

  // tangent for continued biaxial loading from the corner
  const MaterialPointState state{ .stress = result.stress, .stateVars = result.stateVars };
  checkTangent( RankineSofteningPlasticity( linearSofteningProperties.data(), 6, 1 ),
                state,
                ( Math::Vector6d() << 1e-5, 1e-5, -( 2. * nu / ( 1. - nu ) ) * 1e-5, 0, 0, 0 ).finished() );
}

void testActiveSetCorrection()
{
  const RankineSofteningPlasticity material( linearSofteningProperties.data(), 6, 1 );

  // both trial principal stresses 1 and 2 exceed fy0, but only surface 1 is active at the solution
  const Math::Vector6d dE = ( Math::Vector6d() << 1.036e-4, 0.656e-4, 0, 0, 0, 0 ).finished();
  const auto [state, _]   = update( material, MaterialPointState(), dE );

  const double dLambda = ( state.stress( 0 ) - fy0 ) / H;
  throwExceptionOnFailure( checkIfEqual( state.kappa(), dLambda, 1e-15 ), "active set: kappa mismatch" );
  throwExceptionOnFailure( state.stress( 1 ) < state.stress( 0 ), "active set: surface 2 must be inactive" );
  checkAdmissibility( state, linearSofteningProperties );
  checkTangent( material, MaterialPointState(), dE );
}

void testGeneralLoadingPath()
{
  const RankineSofteningPlasticity material( expSofteningProperties.data(), 6, 1 );

  // rotated, multiaxial loading with shear, applied in several increments
  const Math::Vector6d dE = ( Math::Vector6d() << 4e-5, 2e-5, -1e-5, 3e-5, -2e-5, 1e-5 ).finished();

  MaterialPointState state;
  for ( int i = 0; i < 10; i++ ) {
    checkTangent( material, state, dE );
    state = update( material, state, dE ).first;
    checkAdmissibility( state, expSofteningProperties );
  }
  throwExceptionOnFailure( state.kappa() > 0., "general loading path: no plastic flow occurred" );
}

void testPerfectPlasticity()
{
  const RankineSofteningPlasticity material( perfectlyPlasticProperties.data(), 6, 1 );

  const double eps = 1e-4;

  // uniaxial strain, return to surface 1
  {
    const Math::Vector6d dE           = ( Math::Vector6d() << 2. * eps, 0, 0, 0, 0, 0 ).finished();
    const auto [state, _]             = update( material, MaterialPointState(), dE );
    const double         dLambda      = 2. * eps - fy0 / ( lambda + 2. * G );
    const Math::Vector6d stressTarget = ( Math::Vector6d() << fy0,
                                          lambda * ( 2. * eps - dLambda ),
                                          lambda * ( 2. * eps - dLambda ),
                                          0,
                                          0,
                                          0 )
                                          .finished();

    throwExceptionOnFailure( checkIfEqual( state.kappa(), dLambda, 1e-15 ), "perfect plasticity uniaxial: kappa" );
    throwExceptionOnFailure( checkIfEqual< double >( state.stress, stressTarget, 1e-11 ),
                             "perfect plasticity uniaxial: comparison with analytical solution failed" );
    checkTangent( material, MaterialPointState(), dE );
    checkTangent( material, state, dE );
  }

  // equal biaxial strain, return to the intersection of surfaces 1 and 2
  {
    const Math::Vector6d dE           = ( Math::Vector6d() << eps, eps, 0, 0, 0, 0 ).finished();
    const auto [state, _]             = update( material, MaterialPointState(), dE );
    const double         dLambda      = eps - fy0 / ( 2. * lambda + 2. * G );
    const Math::Vector6d stressTarget = ( Math::Vector6d() << fy0, fy0, lambda * ( 2. * eps - 2. * dLambda ), 0, 0, 0 )
                                          .finished();

    throwExceptionOnFailure( checkIfEqual( state.kappa(), 2. * dLambda, 1e-15 ), "perfect plasticity biaxial: kappa" );
    throwExceptionOnFailure( checkIfEqual< double >( state.stress, stressTarget, 1e-11 ),
                             "perfect plasticity biaxial: comparison with analytical solution failed" );
    checkTangent( material, MaterialPointState(), dE );
    checkTangent( material, state, dE );
  }

  // isotropic strain, return to the apex; the tangent vanishes
  {
    const Math::Vector6d dE           = ( Math::Vector6d() << eps, eps, eps, 0, 0, 0 ).finished();
    const auto [state, tangent]       = update( material, MaterialPointState(), dE );
    const double         dLambda      = eps - fy0 / ( 3. * lambda + 2. * G );
    const Math::Vector6d stressTarget = ( Math::Vector6d() << fy0, fy0, fy0, 0, 0, 0 ).finished();

    throwExceptionOnFailure( checkIfEqual( state.kappa(), 3. * dLambda, 1e-15 ), "perfect plasticity apex: kappa" );
    throwExceptionOnFailure( checkIfEqual< double >( state.stress, stressTarget, 1e-11 ),
                             "perfect plasticity apex: comparison with analytical solution failed" );
    throwExceptionOnFailure( tangent.isZero( 1e-8 * E ), "perfect plasticity apex: tangent must vanish" );
    checkTangent( material, MaterialPointState(), dE );
    checkTangent( material, state, dE );
  }

  // rotated, multiaxial loading with shear; the maximum principal stress stays on the yield surface
  {
    const Math::Vector6d dE = ( Math::Vector6d() << 4e-5, 2e-5, -1e-5, 3e-5, -2e-5, 1e-5 ).finished();

    MaterialPointState state;
    for ( int i = 0; i < 10; i++ ) {
      checkTangent( material, state, dE );
      state = update( material, state, dE ).first;
      checkAdmissibility( state, perfectlyPlasticProperties );
    }
    const Eigen::SelfAdjointEigenSolver< Eigen::Matrix3d > eigenSolver(
      Marmot::ContinuumMechanics::Voigt::voigtToStress( state.stress ) );
    throwExceptionOnFailure( state.kappa() > 0. && checkIfEqual( eigenSolver.eigenvalues().maxCoeff(), fy0, 1e-11 ),
                             "perfect plasticity general loading path: stress not on the yield surface" );
  }
}

void testYieldFlag()
{
  const RankineSofteningPlasticity material( linearSofteningProperties.data(), 6, 1 );

  const Math::Vector6d dE = ( Math::Vector6d() << 5e-5, 0, 0, 0, 0, 0 ).finished();

  // elastic loading, plastic loading, elastic unloading
  MaterialPointState state = update( material, MaterialPointState(), dE ).first;
  throwExceptionOnFailure( state.acYield() == 0., "yield flag: elastic step flagged as plastic" );
  state = update( material, state, 2. * dE ).first;
  throwExceptionOnFailure( state.acYield() == 1., "yield flag: plastic step not flagged" );
  state = update( material, state, -dE ).first;
  throwExceptionOnFailure( state.acYield() == 0., "yield flag: unloading flagged as plastic" );
}

void testUniaxialStress()
{
  std::string matName = "RANKINESOFTENINGPLASTICITY";
  auto        solver  = Materials::MaterialPointSolverHypoElastic( matName,
                                                           linearSofteningProperties.data(),
                                                           static_cast< int >( linearSofteningProperties.size() ),
                                                           MaterialPointSolverHypoElastic::SolverOptions() );

  const double eps = 2e-4;

  MaterialPointSolverHypoElastic::Step step;
  step.isStrainComponentControlled = { true, false, false, false, false, false };
  step.isStressComponentControlled = { false, true, true, true, true, true };
  step.stressIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step.strainIncrementTarget       = { eps, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step.dTStart                     = 0.25;
  solver.addStep( step );
  solver.solve();

  // sigma = E ( eps - kappa ) = fy0 + H kappa
  const double         kappa        = ( E * eps - fy0 ) / ( E + H );
  const Math::Vector6d stressTarget = ( Math::Vector6d() << fy0 + H * kappa, 0, 0, 0, 0, 0 ).finished();

  throwExceptionOnFailure( checkIfEqual< double >( solver.getHistory().back().stress, stressTarget, 1e-9 ),
                           "uniaxial stress: comparison with analytical solution failed" );
}

void testCoordinateInvariance()
{
  std::string matName = "RANKINESOFTENINGPLASTICITY";
  auto        solver  = Materials::MaterialPointSolverHypoElastic( matName,
                                                           expSofteningProperties.data(),
                                                           static_cast< int >( expSofteningProperties.size() ),
                                                           MaterialPointSolverHypoElastic::SolverOptions() );

  MaterialPointSolverHypoElastic::Step step;
  step.isStrainComponentControlled = { true, true, true, true, true, true };
  step.isStressComponentControlled = { false, false, false, false, false, false };
  step.stressIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step.strainIncrementTarget       = { 2e-4, 0., 0., 0., 3e-4, 0.0 };
  step.dTStart                     = 1.0;
  solver.addStep( step );

  throwExceptionOnFailure( spinTurbokreisel( solver, 1e-10, 1e-8 ), "Turbokreisel failed!" );
}

/// Strain increments of a diverged iterate (about 1e288, e.g. of an arc-length iteration) overflow the elastic
/// energy: the material must request a cutback instead of returning non-finite results.
void testDivergedIterateRequestsCutback()
{
  const RankineSofteningPlasticity material( perfectlyPlasticProperties.data(),
                                             static_cast< int >( perfectlyPlasticProperties.size() ),
                                             1 );
  const Math::Vector6d             dStrainDiverged = ( Math::Vector6d() << -2.1304181664190551e+287,
                                           -1.0379772584363431e+287,
                                           1.1922728518060270e+045,
                                           -2.9804698180504283e+287,
                                           -3.2636090085991650e+287,
                                           -2.4320172494955940e+287 )
                                           .finished();
  bool cutback = false;
  try {
    update( material, MaterialPointState(), dStrainDiverged );
  }
  catch ( const Utilities::StressUpdateFailed& ) {
    cutback = true;
  }
  throwExceptionOnFailure( cutback, "diverged iterate did not request a cutback" );
  throwExceptionOnFailure( update( material, MaterialPointState(), Math::Vector6d::Constant( 1e-5 ) )
                             .first.stress.allFinite(),
                           "regular increment gave non-finite results" );
}

void testPropertiesAndZeroIncrement()
{
  throwExceptionOnFailure( throws< std::runtime_error >(
                             [] { RankineSofteningPlasticity( linearSofteningProperties.data(), 5, 1 ); } ),
                           "properties: missing properties not detected" );

  const RankineSofteningPlasticity material( linearSofteningProperties.data(), 6, 1 );
  throwExceptionOnFailure( throws< std::runtime_error >( [&] { material.getDensity( nullptr ); } ),
                           "properties: missing density not detected" );
  std::vector< double > withDensity = linearSofteningProperties;
  withDensity.push_back( 2.4e-9 );
  const RankineSofteningPlasticity materialWithDensity( withDensity.data(),
                                                        static_cast< int >( withDensity.size() ),
                                                        1 );
  throwExceptionOnFailure( materialWithDensity.getDensity( nullptr ) == 2.4e-9, "properties: wrong density" );

  // a zero increment returns the elastic stiffness and leaves the state unchanged
  MaterialPointState stateOld;
  stateOld.stress << 1., 0.5, 0., 0.2, 0., 0.;
  const auto [state, tangent] = update( material, stateOld, Math::Vector6d::Zero() );
  throwExceptionOnFailure( state.stress == stateOld.stress && state.stateVars == stateOld.stateVars,
                           "zero increment: state changed" );
  throwExceptionOnFailure( tangent == ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu ),
                           "zero increment: tangent is not the elastic stiffness" );
}

int main()
{
  std::vector< std::function< void( void ) > > tests = { testUniaxialStrain,
                                                         testEqualBiaxialStrain,
                                                         testIsotropicTension,
                                                         testBiaxialTension,
                                                         testActiveSetCorrection,
                                                         testGeneralLoadingPath,
                                                         testPerfectPlasticity,
                                                         testYieldFlag,
                                                         testUniaxialStress,
                                                         testCoordinateInvariance,
                                                         testDivergedIterateRequestsCutback,
                                                         testPropertiesAndZeroIncrement };

  executeTestsAndCollectExceptions( tests );
  return 0;
}
