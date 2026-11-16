#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotExceptions.h"
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
#include <stdexcept>

using namespace Marmot;
using namespace Testing;
using namespace Materials;
using namespace Materials;

namespace {

  const double E        = 30000.;
  const double nu       = 0.2;
  const double ft       = 3.;
  const double epsilonF = 5e-4;
  const double G        = E / ( 2. * ( 1. + nu ) );
  const double lambda   = E * nu / ( ( 1. + nu ) * ( 1. - 2. * nu ) );

  // E, nu, ft, epsilonF, damage threshold omegaMax
  std::vector< double > propertiesWith( double omegaMax )
  {
    return { E, nu, ft, epsilonF, omegaMax };
  }
  std::vector< double > properties = propertiesWith( 1. );
  // without damage
  std::vector< double > noDamageProperties = { E, nu, ft, std::numeric_limits< double >::infinity(), 1. };

  // damage law with threshold omegaMax
  double damage( double kappaD, double omegaMax = 1. )
  {
    return omegaMax * ( 1. - std::exp( -kappaD / ( omegaMax * epsilonF ) ) );
  }

  // perfectly plastic RankineSofteningPlasticity
  std::vector< double > perfectlyPlasticProperties = { E, nu, ft, 0., 0., 0. };

  struct MaterialPointState {
    Math::Vector6d  stress    = Math::Vector6d::Zero();
    Eigen::Vector3d stateVars = Eigen::Vector3d::Zero(); // kappaD, acYield, omega

    double kappaD() const { return stateVars( 0 ); }
    double acYield() const { return stateVars( 1 ); }
    double omega() const { return stateVars( 2 ); }
  };

  struct UpdateResult {
    MaterialPointState state;
    Math::Matrix6d     tangent;
    // increments of elastic energy density and dissipation
    double dPsi, dDissipation;
  };

  // single stress update from a given state
  UpdateResult update( const MaterialHypoElastic& material,
                       const MaterialPointState&  stateOld,
                       const Math::Vector6d&      dStrain )
  {
    UpdateResult result{ .state = stateOld, .tangent = Math::Matrix6d::Zero(), .dPsi = 0., .dDissipation = 0. };
    MaterialHypoElastic::state3D        state( stateOld.stress, 0., 0., result.state.stateVars.data() );
    const MaterialHypoElastic::timeInfo tInfo{ 0.0, 1.0 };
    material.computeStress( state, result.tangent, dStrain, tInfo );
    result.state.stress = state.stress;
    result.dPsi         = state.elasticEnergyDensity;
    result.dDissipation = state.dissipation;
    return result;
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
  void checkTangent( const RankineDamagePlasticity& material,
                     const MaterialPointState&      stateOld,
                     const Math::Vector6d&          dStrain )
  {
    using namespace Math::Differentiation::Numerical;
    const auto result = update( material, stateOld, dStrain );

    const Math::Matrix6d& tangent = result.tangent;

    // differentiate wrt the scaled strain increment, since the default step is relative to max(1, |x|)
    const double          strainScale      = 1e-4;
    const Eigen::MatrixXd numericalTangent = centralDifference(
                                               [&]( const Eigen::VectorXd& y ) -> Eigen::VectorXd {
                                                 const auto r = update( material, stateOld, strainScale * y );
                                                 return r.state.stress;
                                               },
                                               dStrain / strainScale ) /
                                             strainScale;

    // relative to the scale of the moduli, floored since the tangent may vanish
    const double relError = ( tangent - numericalTangent ).cwiseAbs().maxCoeff() /
                            std::max( tangent.cwiseAbs().maxCoeff(), 1e-3 * E );
    throwExceptionOnFailure( relError < 1e-6,
                             Utilities::MakeString() << "tangent check failed, relative error " << relError );
  }

  // check that the effective stress satisfies all three Rankine yield conditions
  void checkAdmissibility( const MaterialPointState& state, double omegaMax = 1. )
  {
    const Eigen::SelfAdjointEigenSolver< Eigen::Matrix3d > eigenSolver(
      Marmot::ContinuumMechanics::Voigt::voigtToStress( state.stress ) );
    throwExceptionOnFailure( eigenSolver.eigenvalues().maxCoeff() <= ( 1. - state.omega() ) * ft * ( 1. + 1e-10 ),
                             "returned stress is not admissible" );
    throwExceptionOnFailure( checkIfEqual( state.omega(), damage( state.kappaD(), omegaMax ), 1e-14 ),
                             "damage inconsistent with kappaD" );
  }

  MaterialPointSolverHypoElastic::Step uniaxialStressStep( double dEps, double timeStart )
  {
    MaterialPointSolverHypoElastic::Step step;
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
  auto        solver  = Materials::MaterialPointSolverHypoElastic( matName,
                                                           properties.data(),
                                                           static_cast< int >( properties.size() ),
                                                           MaterialPointSolverHypoElastic::SolverOptions() );

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
  const RankineDamagePlasticity material( properties.data(), static_cast< int >( properties.size() ), 1 );
  const auto [_, tangent, dPsi, dD] = update( material,
                                              state,
                                              ( Math::Vector6d() << -1e-5, 0, 0, 0, 0, 0 ).finished() );
  const Math::Matrix6d Cel          = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );
  throwExceptionOnFailure( checkIfEqual< double >( tangent, ( 1. - state.omega() ) * Cel, 1e-12 ),
                           "uniaxial stress: unloading tangent mismatch" );
}

void testUniaxialStrain()
{
  const RankineDamagePlasticity material( properties.data(), static_cast< int >( properties.size() ), 1 );

  // single plastic increment from the virgin state: the effective stress is the perfectly plastic solution
  const double         eps          = 4e-4;
  const Math::Vector6d dE           = ( Math::Vector6d() << eps, 0, 0, 0, 0, 0 ).finished();
  const auto           result       = update( material, MaterialPointState(), dE );
  const double         dLambda      = eps - ft / ( lambda + 2. * G );
  const double         omega        = 1. - std::exp( -dLambda / epsilonF );
  const Math::Vector6d stressTarget = ( 1. - omega ) * ( Math::Vector6d() << ft,
                                                         lambda * ( eps - dLambda ),
                                                         lambda * ( eps - dLambda ),
                                                         0,
                                                         0,
                                                         0 )
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
  const RankineDamagePlasticity    material( noDamageProperties.data(),
                                          static_cast< int >( noDamageProperties.size() ),
                                          1 );
  const RankineSofteningPlasticity reference( perfectlyPlasticProperties.data(), 6, 1 );

  const std::vector< Math::Vector6d > increments =
    { ( Math::Vector6d() << 1e-4, 0, 0, 0, 0, 0 ).finished(),
      ( Math::Vector6d() << 5e-5, 5e-5, 0, 0, 0, 0 ).finished(),
      ( Math::Vector6d() << 4e-5, 4e-5, 4e-5, 0, 0, 0 ).finished(),
      ( Math::Vector6d() << 4e-5, 2e-5, -1e-5, 3e-5, -2e-5, 1e-5 ).finished() };

  for ( const auto& dE : increments ) {
    MaterialPointState state;
    Math::Vector6d     stressReference = Math::Vector6d::Zero();
    Eigen::Vector2d    stateVarsReference( 0., 0. );

    for ( int i = 0; i < 5; i++ ) {
      const auto result = update( material, state, dE );
      state             = result.state;

      MaterialHypoElastic::state3D stateReference( stressReference, 0., 0., stateVarsReference.data() );
      Math::Matrix6d               tangentReference;
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
  const RankineDamagePlasticity material( properties.data(), static_cast< int >( properties.size() ), 1 );

  // one (uniaxial), two (biaxial) and three (isotropic, apex) active surfaces, continued over several increments
  for ( const auto& dE : { ( Math::Vector6d() << 1e-4, 0, 0, 0, 0, 0 ).finished(),
                           ( Math::Vector6d() << 6e-5, 6e-5, 0, 0, 0, 0 ).finished(),
                           ( Math::Vector6d() << 5e-5, 5e-5, 5e-5, 0, 0, 0 ).finished(),
                           ( Math::Vector6d() << 6e-5, 5e-5, -1e-5, 0, 0, 0 ).finished() } ) {
    MaterialPointState state;
    for ( int i = 0; i < 4; i++ ) {
      checkTangent( material, state, dE );
      state = update( material, state, dE ).state;
      checkAdmissibility( state );
    }
    throwExceptionOnFailure( state.omega() > 0.1, "multi surface: no damage occurred" );
  }

  // active set correction: both trial stresses exceed ft, only surface 1 is active at the solution
  const Math::Vector6d dE    = ( Math::Vector6d() << 1.036e-4, 0.656e-4, 0, 0, 0, 0 ).finished();
  const auto           state = update( material, MaterialPointState(), dE ).state;
  checkAdmissibility( state );
  checkTangent( material, MaterialPointState(), dE );
}

void testTangentGeneralLoadingPath()
{
  const RankineDamagePlasticity material( properties.data(), static_cast< int >( properties.size() ), 1 );

  // rotated, multiaxial loading with shear, followed by unloading and reloading in a different direction
  const Math::Vector6d dE         = ( Math::Vector6d() << 4e-5, 2e-5, -1e-5, 3e-5, -2e-5, 1e-5 ).finished();
  const Math::Vector6d dEReloaded = ( Math::Vector6d() << 1e-5, 5e-5, 1e-5, -2e-5, 0, 2e-5 ).finished();

  MaterialPointState state;
  for ( const Math::Vector6d& increment : std::vector< Math::Vector6d >( 8, dE ) ) {
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

void checkEnergyBalance( double omegaMax )
{
  const auto                    props = propertiesWith( omegaMax );
  const RankineDamagePlasticity material( props.data(), static_cast< int >( props.size() ), 1 );

  // uniaxial strain from first yield: the effective stress is constant, so the backward Euler update satisfies the
  // energy balance sigma_n+1 : dE = dPsi + dDissipation exactly
  const double         epsYield = ft / ( lambda + 2. * G );
  const Math::Vector6d dE       = ( Math::Vector6d() << 1e-4, 0, 0, 0, 0, 0 ).finished();

  auto result = update( material, MaterialPointState(), ( Math::Vector6d() << epsYield, 0, 0, 0, 0, 0 ).finished() );
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
  throwExceptionOnFailure( dissipation > 0.9 * omegaMax * ( ft * epsilonF ), "energy balance: dissipation too small" );
}

void testEnergyBalance()
{
  for ( const double omegaMax : { 1., 0.9 } )
    checkEnergyBalance( omegaMax );
}

/// Damage threshold omegaMax < 1: initial slope 1 / epsilonF, damage bounded by omegaMax with the residual stress
/// (1 - omegaMax) ft under uniaxial stress, consistent tangents, and no complete damage for large kappaD.
void testDamageThreshold()
{
  const double omegaMax = 0.9;
  auto         props    = propertiesWith( omegaMax );

  // uniaxial stress far into the softening regime: sigma = (1 - omega(kappaD)) ft with kappaD = eps - ft / E
  std::string matName = "RANKINEDAMAGEPLASTICITY";
  auto        solver  = Materials::MaterialPointSolverHypoElastic( matName,
                                                           props.data(),
                                                           static_cast< int >( props.size() ),
                                                           MaterialPointSolverHypoElastic::SolverOptions() );
  solver.addStep( uniaxialStressStep( 5e-3, 0. ) );
  solver.solve();
  for ( const auto& entry : solver.getHistory() ) {
    const double eps    = entry.strain( 0 );
    const double kappaD = entry.stateVars( 0 );
    if ( eps > ft / E ) {
      throwExceptionOnFailure( checkIfEqual( kappaD, eps - ft / E, 1e-13 ), "threshold: kappaD mismatch" );
      throwExceptionOnFailure( checkIfEqual( entry.stress( 0 ), ( 1. - damage( kappaD, omegaMax ) ) * ft, 1e-10 ),
                               "threshold: softening law violated" );
    }
    throwExceptionOnFailure( entry.stateVars( 2 ) < omegaMax, "threshold: damage exceeds omegaMax" );
  }
  const double residualStress = solver.getHistory().back().stress( 0 );
  throwExceptionOnFailure( residualStress > ( 1. - omegaMax ) * ft && residualStress < 1.001 * ( 1. - omegaMax ) * ft,
                           "threshold: no residual stress (1 - omegaMax) ft" );

  // initial slope d omega / d kappaD = 1 / epsilonF, independent of omegaMax
  const double epsYield = ft / ( lambda + 2. * G );
  for ( const double w : { 1., omegaMax } ) {
    const auto                    p = propertiesWith( w );
    const RankineDamagePlasticity material( p.data(), static_cast< int >( p.size() ), 1 );
    const auto                    state = update( material,
                               MaterialPointState(),
                               ( Math::Vector6d() << epsYield + 1e-9, 0, 0, 0, 0, 0 ).finished() )
                         .state;
    throwExceptionOnFailure( checkIfEqual( state.omega() / state.kappaD(), 1. / epsilonF, 1e-5 ),
                             "threshold: initial slope is not 1 / epsilonF" );
  }

  // consistent tangents with one, two and three active surfaces
  const RankineDamagePlasticity material( props.data(), static_cast< int >( props.size() ), 1 );
  for ( const auto& dE : { ( Math::Vector6d() << 1e-4, 0, 0, 0, 0, 0 ).finished(),
                           ( Math::Vector6d() << 6e-5, 6e-5, 0, 0, 0, 0 ).finished(),
                           ( Math::Vector6d() << 5e-5, 5e-5, 5e-5, 0, 0, 0 ).finished() } ) {
    MaterialPointState state;
    for ( int i = 0; i < 6; i++ ) {
      checkTangent( material, state, dE );
      state = update( material, state, dE ).state;
      checkAdmissibility( state, omegaMax );
    }
  }

  // large kappaD: the damage saturates at omegaMax, whereas omegaMax = 1 means complete damage (cutback)
  const Math::Vector6d dELarge = ( Math::Vector6d() << 1., 0, 0, 0, 0, 0 ).finished();
  const auto           large   = update( material, MaterialPointState(), dELarge );
  throwExceptionOnFailure( checkIfEqual( large.state.omega(), omegaMax, 1e-15 ) &&
                             checkIfEqual( large.state.stress( 0 ), ( 1. - omegaMax ) * ft, 1e-12 ),
                           "threshold: damage does not saturate at omegaMax" );
  const RankineDamagePlasticity unbounded( properties.data(), static_cast< int >( properties.size() ), 1 );
  throwExceptionOnFailure( throws< Utilities::StressUpdateFailed >(
                             [&] { update( unbounded, MaterialPointState(), dELarge ); } ),
                           "threshold: complete damage without threshold did not request a cutback" );
}

void testCoordinateInvariance()
{
  std::string matName = "RANKINEDAMAGEPLASTICITY";
  auto        solver  = Materials::MaterialPointSolverHypoElastic( matName,
                                                           properties.data(),
                                                           static_cast< int >( properties.size() ),
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

/// Quadrature point states of a diverged iterate (dissipation control, tensile bar at the onset of dissipation, strain
/// increments of about 1e288): the update overflows (the elastic energy to inf, complete
/// damage at the plastic point) and must request a cutback instead of returning non-finite results.
void testDivergedIterateRequestsCutback()
{
  const std::vector< double >   props = { 210000., 0.3, 693., 3.3, 1. };
  const RankineDamagePlasticity material( props.data(), static_cast< int >( props.size() ), 1 );

  struct Case {
    MaterialPointState state;
    Math::Vector6d     dStrain;
  };
  std::vector< Case > cases( 3 );
  // elastic quadrature point
  cases[0].state.stress << 6.5087488797518449e+02, -1.1460284989614635e+00, -9.2737855777728662e-01,
    1.2095787026234459e+00, 5.0457762116698635e-01, -5.0845805470899919e-01;
  cases[0].dStrain << -2.1304181664190551e+287, -1.0379772584363431e+287, 1.1922728518060270e+045,
    -2.9804698180504283e+287, -3.2636090085991650e+287, -2.4320172494955940e+287;
  // plastic quadrature point
  cases[1].state.stress << 6.6713904544172669e+02, 3.5094066506532569e+00, 5.5710345496692897e+00,
    5.0537225453201051e-01, -2.2062866091754861e-02, -3.2671228863947905e+00;
  cases[1].state.stateVars << 1.0908213453838e-04, 1., 3.3054645967478535e-05;
  cases[1].dStrain << -7.9390800699293370e+287, 2.4129632418183885e+287, 2.9152508981520328e+044,
    -2.6781038693586638e+287, -1.0266542866855759e+288, -7.6505516681351224e+287;
  // elastic quadrature point, larger increment
  cases[2].state.stress << 6.3594567924872888e+02, 4.6210448040177383e+01, 4.5978543623505857e+01,
    -7.7865816008572999e+00, -6.7477652659058602e+00, -3.9672606952376455e-01;
  cases[2].dStrain << -7.0908769355913750e+288, 8.0118534558425438e+287, -3.0696396257343197e+009,
    4.6872184635987446e+288, 2.7759733987994883e+288, -1.8229633857534692e+287;

  for ( const auto& c : cases ) {
    bool cutback = false;
    try {
      update( material, c.state, c.dStrain );
    }
    catch ( const Utilities::StressUpdateFailed& ) {
      cutback = true;
    }
    throwExceptionOnFailure( cutback, "diverged iterate did not request a cutback" );
  }

  // the converged states themselves with a regular increment remain valid
  for ( const auto& c : cases ) {
    const auto result = update( material, c.state, Math::Vector6d::Constant( 1e-4 ) );
    throwExceptionOnFailure( result.state.stress.allFinite() && std::isfinite( result.dPsi ) &&
                               std::isfinite( result.dDissipation ),
                             "regular increment gave non-finite results" );
  }
}

void testPropertiesAndZeroIncrement()
{
  throwExceptionOnFailure( throws< std::runtime_error >( [] { RankineDamagePlasticity( properties.data(), 4, 1 ); } ),
                           "properties: missing properties not detected" );
  const std::vector< double > noSoftening = { E, nu, ft, 0., 1. };
  throwExceptionOnFailure( throws< std::invalid_argument >( [&] {
                             RankineDamagePlasticity( noSoftening.data(), static_cast< int >( noSoftening.size() ), 1 );
                           } ),
                           "properties: nonpositive epsilonF not detected" );
  for ( const double omegaMax : { 0., -0.1, 1.1 } ) {
    const auto p = propertiesWith( omegaMax );
    throwExceptionOnFailure( throws< std::invalid_argument >(
                               [&] { RankineDamagePlasticity( p.data(), static_cast< int >( p.size() ), 1 ); } ),
                             "properties: damage threshold outside (0, 1] not detected" );
  }

  const RankineDamagePlasticity material( properties.data(), 5, 1 );
  throwExceptionOnFailure( throws< std::runtime_error >( [&] { material.getDensity( nullptr ); } ),
                           "properties: missing density not detected" );
  const std::vector< double >   withDensity = { E, nu, ft, epsilonF, 1., 2.4e-9 };
  const RankineDamagePlasticity materialWithDensity( withDensity.data(), static_cast< int >( withDensity.size() ), 1 );
  throwExceptionOnFailure( materialWithDensity.getDensity( nullptr ) == 2.4e-9, "properties: wrong density" );

  // a zero increment from a damaged plastic state returns the damaged elastic stiffness, leaves the state unchanged
  // and is not flagged as plastic
  const auto damaged = update( material, MaterialPointState(), ( Math::Vector6d() << 4e-4, 0, 0, 0, 0, 0 ).finished() )
                         .state;
  const auto result = update( material, damaged, Math::Vector6d::Zero() );
  throwExceptionOnFailure( damaged.omega() > 0. && damaged.acYield() == 1. && result.state.stress == damaged.stress &&
                             result.state.kappaD() == damaged.kappaD() && result.state.omega() == damaged.omega(),
                           "zero increment: state changed" );
  throwExceptionOnFailure( result.state.acYield() == 0., "zero increment: flagged as plastic" );
  throwExceptionOnFailure( checkIfEqual< double >( result.tangent,
                                                   ( 1. - damaged.omega() ) *
                                                     ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E,
                                                                                                                 nu ),
                                                   1e-12 ),
                           "zero increment: tangent is not the damaged elastic stiffness" );

  // a completely damaged old state has no effective stress
  MaterialPointState completelyDamaged;
  completelyDamaged.stateVars << 1e6, 1., 1.;
  throwExceptionOnFailure( throws< Utilities::StressUpdateFailed >(
                             [&] { update( material, completelyDamaged, Math::Vector6d::Zero() ); } ),
                           "complete damage of the old state did not request a cutback" );

  // complete damage within the increment
  throwExceptionOnFailure( throws< Utilities::StressUpdateFailed >( [&] {
                             update( material,
                                     MaterialPointState(),
                                     ( Math::Vector6d() << 10., 0, 0, 0, 0, 0 ).finished() );
                           } ),
                           "complete damage within the increment did not request a cutback" );
}

int main()
{
  std::vector< std::function< void( void ) > > tests = { testUniaxialStress,
                                                         testUniaxialStrain,
                                                         testNoDamageLimit,
                                                         testTangentSingleAndMultiSurface,
                                                         testTangentGeneralLoadingPath,
                                                         testEnergyBalance,
                                                         testDamageThreshold,
                                                         testCoordinateInvariance,
                                                         testDivergedIterateRequestsCutback,
                                                         testPropertiesAndZeroIncrement };

  executeTestsAndCollectExceptions( tests );
  return 0;
}
