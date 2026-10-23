#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotMaterialPointSolverHypoElastic.h"
#include "Marmot/MarmotNumericalDifferentiation.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include "Marmot/RankineDamagePlasticity.h"
#include "Marmot/RankineSofteningPlasticity.h"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <limits>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::Solvers;
using namespace Marmot::Materials;

namespace {

  const double E        = 30000.;
  const double nu       = 0.2;
  const double ft       = 3.;
  const double epsilonF = 5e-4;
  const double G        = E / ( 2. * ( 1. + nu ) );
  const double lambda   = E * nu / ( ( 1. + nu ) * ( 1. - 2. * nu ) );

  // E, nu, ft, epsilonF
  std::vector< double > properties = { E, nu, ft, epsilonF };
  // without damage
  std::vector< double > noDamageProperties = { E, nu, ft, std::numeric_limits< double >::infinity() };
  // perfectly plastic RankineSofteningPlasticity
  std::vector< double > perfectlyPlasticProperties = { E, nu, ft, 0., 0., 0. };

  struct MaterialPointState {
    Vector6d        stress    = Vector6d::Zero();
    Eigen::Vector3d stateVars = Eigen::Vector3d::Zero(); // kappaD, acYield, omega

    double kappaD() const { return stateVars( 0 ); }
    double acYield() const { return stateVars( 1 ); }
    double omega() const { return stateVars( 2 ); }
  };

  struct UpdateResult {
    MaterialPointState state;
    Matrix6d           tangent;
    // increments of elastic energy density and dissipation, and their derivatives wrt the strain increment
    double   dPsi, dDissipation;
    Vector6d dPsi_dE, dDissipation_dE;
  };

  // single stress update from a given state
  UpdateResult update( const MarmotMaterialHypoElastic& material,
                       const MaterialPointState&        stateOld,
                       const Vector6d&                  dStrain )
  {
    UpdateResult                              result{ .state           = stateOld,
                                                      .tangent         = Matrix6d::Zero(),
                                                      .dPsi            = 0.,
                                                      .dDissipation    = 0.,
                                                      .dPsi_dE         = Vector6d::Zero(),
                                                      .dDissipation_dE = Vector6d::Zero() };
    MarmotMaterialHypoElastic::state3D        state( stateOld.stress, 0., 0., result.state.stateVars.data() );
    const MarmotMaterialHypoElastic::timeInfo tInfo{ 0.0, 1.0 };
    material.computeStress( state, result.tangent, dStrain, tInfo );
    result.state.stress    = state.stress;
    result.dPsi            = state.elasticEnergyDensity;
    result.dDissipation    = state.dissipation;
    result.dPsi_dE         = state.dElasticEnergyDensity_dStrain;
    result.dDissipation_dE = state.dDissipation_dStrain;
    return result;
  }

  // compare the algorithmic tangent and the energy tangents against central difference approximations
  void checkTangent( const RankineDamagePlasticity& material,
                     const MaterialPointState&      stateOld,
                     const Vector6d&                dStrain )
  {
    using namespace NumericalAlgorithms::Differentiation;
    const auto result = update( material, stateOld, dStrain );

    Eigen::Matrix< double, 8, 6 > tangent;
    tangent << result.tangent, result.dPsi_dE.transpose(), result.dDissipation_dE.transpose();

    // differentiate wrt the scaled strain increment, since the default step is relative to max(1, |x|)
    const double strainScale = 1e-4;
    const Eigen::MatrixXd
      numericalTangent = centralDifference(
                           [&]( const Eigen::VectorXd& y ) -> Eigen::VectorXd {
                             const auto r = update( material, stateOld, strainScale * y );
                             return ( Eigen::VectorXd( 8 ) << r.state.stress, r.dPsi, r.dDissipation ).finished();
                           },
                           dStrain / strainScale ) /
                         strainScale;

    // moduli and each energy gradient relative to their own scale, floored since the tangent may vanish
    const Eigen::MatrixXd error    = tangent - numericalTangent;
    const double          relError = std::max(
      { error.topRows( 6 ).cwiseAbs().maxCoeff() / std::max( tangent.topRows( 6 ).cwiseAbs().maxCoeff(), 1e-3 * E ),
        error.row( 6 ).cwiseAbs().maxCoeff() / std::max( tangent.row( 6 ).cwiseAbs().maxCoeff(), 1e-3 * ft ),
        error.row( 7 ).cwiseAbs().maxCoeff() / std::max( tangent.row( 7 ).cwiseAbs().maxCoeff(), 1e-3 * ft ) } );
    throwExceptionOnFailure( relError < 1e-6, MakeString() << "tangent check failed, relative error " << relError );
  }

  // check that the effective stress satisfies all three Rankine yield conditions
  void checkAdmissibility( const MaterialPointState& state )
  {
    const Eigen::SelfAdjointEigenSolver< Eigen::Matrix3d > eigenSolver(
      ContinuumMechanics::VoigtNotation::voigtToStress( state.stress ) );
    throwExceptionOnFailure( eigenSolver.eigenvalues().maxCoeff() <= ( 1. - state.omega() ) * ft * ( 1. + 1e-10 ),
                             "returned stress is not admissible" );
    throwExceptionOnFailure( checkIfEqual( state.omega(), 1. - std::exp( -state.kappaD() / epsilonF ), 1e-15 ),
                             "damage inconsistent with kappaD" );
  }

  MarmotMaterialPointSolverHypoElastic::Step uniaxialStressStep( double dEps, double timeStart )
  {
    MarmotMaterialPointSolverHypoElastic::Step step;
    step.timeStart                   = timeStart;
    step.timeEnd                     = timeStart + 1.;
    step.isStrainComponentControlled = { true, false, false, false, false, false };
    step.isStressComponentControlled = { false, true, true, true, true, true };
    step.stressIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    step.strainIncrementTarget       = { dEps, 0.0, 0.0, 0.0, 0.0, 0.0 };
    step.dTStart                     = 0.02;
    step.dTMax                       = 0.02;
    return step;
  }

} // namespace

void testUniaxialStress()
{
  std::string matName = "RANKINEDAMAGEPLASTICITY";
  auto        solver  = MarmotMaterialPointSolverHypoElastic( matName,
                                                              properties.data(),
                                                              properties.size(),
                                                              MarmotMaterialPointSolverHypoElastic::SolverOptions() );

  // loading far into the softening regime, unloading, reloading
  const double epsLoad = 2e-3, epsUnload = -1e-3;
  solver.addStep( uniaxialStressStep( epsLoad, 0. ) );
  solver.addStep( uniaxialStressStep( epsUnload, 1. ) );
  solver.addStep( uniaxialStressStep( -epsUnload, 2. ) );
  solver.solve();

  double kappaDMax = 0;
  for ( const auto& entry : solver.getHistory() ) {
    const double eps    = entry.strain( 0 );
    const double kappaD = entry.stateVars( 0 );
    const double omega  = entry.stateVars( 2 );
    throwExceptionOnFailure( kappaD >= kappaDMax - 1e-15, "uniaxial stress: kappaD decreased" );
    kappaDMax = std::max( kappaD, kappaDMax );

    // softening branch: constant effective stress ft, sigma = ft exp(-eps^p / epsilonF) with eps^p = kappaD
    if ( entry.time <= 1.0 + 1e-12 && eps > ft / E ) {
      throwExceptionOnFailure( checkIfEqual( kappaD, eps - ft / E, 1e-13 ), "uniaxial stress: kappaD mismatch" );
      throwExceptionOnFailure( checkIfEqual( entry.stress( 0 ), ft * std::exp( -kappaD / epsilonF ), 1e-10 ),
                               "uniaxial stress: softening law violated" );
    }
    // elastic unloading and reloading with the damaged stiffness (1 - omega) E
    if ( entry.time > 1.0 + 1e-12 ) {
      throwExceptionOnFailure( checkIfEqual( entry.stress( 0 ), ( 1. - omega ) * E * ( eps - kappaD ), 1e-10 ),
                               "uniaxial stress: damaged unloading stiffness mismatch" );
      throwExceptionOnFailure( checkIfEqual( kappaD, epsLoad - ft / E, 1e-13 ),
                               "uniaxial stress: damage evolved during unloading" );
    }
    throwExceptionOnFailure( entry.stress.tail( 5 ).cwiseAbs().maxCoeff() < 1e-9,
                             "uniaxial stress: lateral stress not zero" );
  }
  throwExceptionOnFailure( kappaDMax > 3. * epsilonF, "uniaxial stress: softening regime not reached" );

  // the unloading tangent is the damaged elastic stiffness
  const MaterialPointState      state{ .stress    = solver.getHistory().back().stress,
                                       .stateVars = solver.getHistory().back().stateVars };
  const RankineDamagePlasticity material( properties.data(), properties.size(), 1 );
  const auto [_, tangent, dPsi, dD, dPsi_dE, dD_dE] = update( material,
                                                              state,
                                                              ( Vector6d() << -1e-5, 0, 0, 0, 0, 0 ).finished() );
  const Matrix6d Cel = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );
  throwExceptionOnFailure( checkIfEqual< double >( tangent, ( 1. - state.omega() ) * Cel, 1e-12 ),
                           "uniaxial stress: unloading tangent mismatch" );
}

void testUniaxialStrain()
{
  const RankineDamagePlasticity material( properties.data(), properties.size(), 1 );

  // single plastic increment from the virgin state: the effective stress is the perfectly plastic solution
  const double   eps          = 4e-4;
  const Vector6d dE           = ( Vector6d() << eps, 0, 0, 0, 0, 0 ).finished();
  const auto     result       = update( material, MaterialPointState(), dE );
  const double   dLambda      = eps - ft / ( lambda + 2. * G );
  const double   omega        = 1. - std::exp( -dLambda / epsilonF );
  const Vector6d stressTarget = ( 1. - omega ) *
                                ( Vector6d() << ft, lambda * ( eps - dLambda ), lambda * ( eps - dLambda ), 0, 0, 0 )
                                  .finished();

  throwExceptionOnFailure( checkIfEqual( result.state.kappaD(), dLambda, 1e-15 ), "uniaxial strain: kappaD mismatch" );
  throwExceptionOnFailure( checkIfEqual( result.state.omega(), omega, 1e-15 ), "uniaxial strain: omega mismatch" );
  throwExceptionOnFailure( checkIfEqual< double >( result.state.stress, stressTarget, 1e-11 ),
                           "uniaxial strain: comparison with analytical solution failed" );
  checkTangent( material, MaterialPointState(), dE );
  checkTangent( material, result.state, dE );
}

void testNoDamageLimit()
{
  const RankineDamagePlasticity    material( noDamageProperties.data(), noDamageProperties.size(), 1 );
  const RankineSofteningPlasticity reference( perfectlyPlasticProperties.data(), 6, 1 );

  const std::vector< Vector6d > increments = { ( Vector6d() << 1e-4, 0, 0, 0, 0, 0 ).finished(),
                                               ( Vector6d() << 5e-5, 5e-5, 0, 0, 0, 0 ).finished(),
                                               ( Vector6d() << 4e-5, 4e-5, 4e-5, 0, 0, 0 ).finished(),
                                               ( Vector6d() << 4e-5, 2e-5, -1e-5, 3e-5, -2e-5, 1e-5 ).finished() };

  for ( const auto& dE : increments ) {
    MaterialPointState state;
    Vector6d           stressReference = Vector6d::Zero();
    Eigen::Vector2d    stateVarsReference( 0., 0. );

    for ( int i = 0; i < 5; i++ ) {
      const auto result = update( material, state, dE );
      state             = result.state;

      MarmotMaterialHypoElastic::state3D stateReference( stressReference, 0., 0., stateVarsReference.data() );
      Matrix6d                           tangentReference;
      reference.computeStress( stateReference, tangentReference, dE, { 0.0, 1.0 } );
      stressReference = stateReference.stress;

      throwExceptionOnFailure( state.stress == stressReference && result.tangent == tangentReference,
                               "no damage: stress or tangent differs from RankineSofteningPlasticity" );
      throwExceptionOnFailure( state.kappaD() == stateVarsReference( 0 ) && state.omega() == 0. &&
                                 state.acYield() == stateVarsReference( 1 ),
                               "no damage: state differs from RankineSofteningPlasticity" );
    }
    throwExceptionOnFailure( state.acYield() == 1., "no damage: no plastic flow occurred" );
  }
}

void testTangentSingleAndMultiSurface()
{
  const RankineDamagePlasticity material( properties.data(), properties.size(), 1 );

  // one (uniaxial), two (biaxial) and three (isotropic, apex) active surfaces, continued over several increments
  for ( const auto& dE : { ( Vector6d() << 1e-4, 0, 0, 0, 0, 0 ).finished(),
                           ( Vector6d() << 6e-5, 6e-5, 0, 0, 0, 0 ).finished(),
                           ( Vector6d() << 5e-5, 5e-5, 5e-5, 0, 0, 0 ).finished(),
                           ( Vector6d() << 6e-5, 5e-5, -1e-5, 0, 0, 0 ).finished() } ) {
    MaterialPointState state;
    for ( int i = 0; i < 4; i++ ) {
      checkTangent( material, state, dE );
      state = update( material, state, dE ).state;
      checkAdmissibility( state );
    }
    throwExceptionOnFailure( state.omega() > 0.1, "multi surface: no damage occurred" );
  }

  // active set correction: both trial stresses exceed ft, only surface 1 is active at the solution
  const Vector6d dE    = ( Vector6d() << 1.036e-4, 0.656e-4, 0, 0, 0, 0 ).finished();
  const auto     state = update( material, MaterialPointState(), dE ).state;
  checkAdmissibility( state );
  checkTangent( material, MaterialPointState(), dE );
}

void testTangentGeneralLoadingPath()
{
  const RankineDamagePlasticity material( properties.data(), properties.size(), 1 );

  // rotated, multiaxial loading with shear, followed by unloading and reloading in a different direction
  const Vector6d dE         = ( Vector6d() << 4e-5, 2e-5, -1e-5, 3e-5, -2e-5, 1e-5 ).finished();
  const Vector6d dEReloaded = ( Vector6d() << 1e-5, 5e-5, 1e-5, -2e-5, 0, 2e-5 ).finished();

  MaterialPointState state;
  for ( const Vector6d& increment : std::vector< Vector6d >( 8, dE ) ) {
    checkTangent( material, state, increment );
    state = update( material, state, increment ).state;
    checkAdmissibility( state );
  }

  // elastic unloading: damaged elastic tangent, no evolution of damage
  const double omega = state.omega();
  checkTangent( material, state, -0.1 * dE );
  state = update( material, state, -0.1 * dE ).state;
  throwExceptionOnFailure( state.acYield() == 0. && state.omega() == omega, "general path: unloading not elastic" );

  for ( int i = 0; i < 6; i++ ) {
    checkTangent( material, state, dEReloaded );
    state = update( material, state, dEReloaded ).state;
    checkAdmissibility( state );
  }
  throwExceptionOnFailure( state.omega() > omega, "general path: no further damage on reloading" );
}

void testEnergyBalance()
{
  const RankineDamagePlasticity material( properties.data(), properties.size(), 1 );

  // uniaxial strain from first yield: the effective stress is constant, so the backward Euler update satisfies the
  // energy balance sigma_n+1 : dE = dPsi + dDissipation exactly
  const double   epsYield = ft / ( lambda + 2. * G );
  const Vector6d dE       = ( Vector6d() << 1e-4, 0, 0, 0, 0, 0 ).finished();

  auto result = update( material, MaterialPointState(), ( Vector6d() << epsYield, 0, 0, 0, 0, 0 ).finished() );
  MaterialPointState state = result.state;
  double             work = 0, psi = result.dPsi, dissipation = 0;
  for ( int i = 0; i < 30; i++ ) {
    result = update( material, state, dE );
    state  = result.state;
    work += state.stress.dot( dE );
    psi += result.dPsi;
    dissipation += result.dDissipation;
  }

  const double psiTotal = 0.5 * ( 1. - state.omega() ) * ( lambda + 2. * G ) * epsYield * epsYield;
  throwExceptionOnFailure( checkIfEqual( psi, psiTotal, 1e-13 ), "energy balance: elastic energy mismatch" );
  throwExceptionOnFailure( checkIfEqual( work + 0.5 * ft * epsYield, psi + dissipation, 1e-12 ),
                           "energy balance violated" );
  throwExceptionOnFailure( dissipation > 0.9 * ( ft * epsilonF ), "energy balance: dissipation too small" );
}

void testCoordinateInvariance()
{
  std::string matName = "RANKINEDAMAGEPLASTICITY";
  auto        solver  = MarmotMaterialPointSolverHypoElastic( matName,
                                                              properties.data(),
                                                              properties.size(),
                                                              MarmotMaterialPointSolverHypoElastic::SolverOptions() );

  MarmotMaterialPointSolverHypoElastic::Step step;
  step.isStrainComponentControlled = { true, true, true, true, true, true };
  step.isStressComponentControlled = { false, false, false, false, false, false };
  step.stressIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step.strainIncrementTarget       = { 2e-4, 0., 0., 0., 3e-4, 0.0 };
  step.dTStart                     = 1.0;
  solver.addStep( step );

  throwExceptionOnFailure( spinTurbokreisel( solver, 1e-10, 1e-8 ), "Turbokreisel failed!" );
}

int main()
{
  std::vector< std::function< void( void ) > > tests = { testUniaxialStress,
                                                         testUniaxialStrain,
                                                         testNoDamageLimit,
                                                         testTangentSingleAndMultiSurface,
                                                         testTangentGeneralLoadingPath,
                                                         testEnergyBalance,
                                                         testCoordinateInvariance };

  executeTestsAndCollectExceptions( tests );
  return 0;
}
