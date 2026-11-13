#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotMaterialPointSolverHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/VonMises.h"

using namespace Marmot::Testing;
using namespace Marmot::Materials;

void testVonMisesCoordinateInvariance()
{
  // material properties
  std::vector< double > materialProperties = { 210000., 0.3, 200., 2100., 20., 20 };
  auto                  solveropts         = Marmot::Materials::MaterialPointSolverHypoElastic::SolverOptions();
  std::string           matName            = "VONMISES";

  // create material point solver instance
  auto solver = Marmot::Materials::MaterialPointSolverHypoElastic( matName,
                                                                   &materialProperties[0],
                                                                   static_cast< int >( materialProperties.size() ),
                                                                   solveropts );

  // define a step
  Marmot::Materials::MaterialPointSolverHypoElastic::Step step;

  // define step parameters
  step.isStrainComponentControlled = { true, true, true, true, true, true };
  step.isStressComponentControlled = { false, false, false, false, false, false };
  step.stressIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step.strainIncrementTarget       = { 0.01, 0., 0., 0., 0.03, 0.0 };
  step.dTStart                     = 1.0;

  // add step to solver
  solver.addStep( step );
  // check coordinate invariance with Turbokreisel
  throwExceptionOnFailure( spinTurbokreisel( solver, 1e-10, 1e-8 ), "Turbokreisel failed!" );
}

void testVonMises()
{
  // material properties
  std::vector< double > materialProperties = { 210000., 0.3, 200., 2100., 20., 20 };
  auto                  solveropts         = Marmot::Materials::MaterialPointSolverHypoElastic::SolverOptions();
  std::string           matName            = "VONMISES";

  // create material point solver instance
  auto solver = Marmot::Materials::MaterialPointSolverHypoElastic( matName,
                                                                   &materialProperties[0],
                                                                   static_cast< int >( materialProperties.size() ),
                                                                   solveropts );

  // define a step
  Marmot::Materials::MaterialPointSolverHypoElastic::Step step;

  // define step parameters
  step.isStrainComponentControlled = { true, true, true, true, true, true };
  step.isStressComponentControlled = { false, false, false, false, false, false };
  step.stressIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step.strainIncrementTarget       = { 0.00839244, 0.00089344, -0.00703916, 0.00013635, 0.00160548, 0.00572825 };
  step.dTStart                     = 1.0;

  // add step to solver
  solver.addStep( step );

  solver.solve();

  // read history
  auto history = solver.getHistory();
  // reference solution
  const Marmot::Math::Vector6d stressTarget =
    { 511.262747695, 395.408929527, 272.856322779, 1.0532516407, 12.4017194288, 44.2485420671 };

  // compare solutions
  throwExceptionOnFailure( checkIfEqual< double >( history.back().stress, stressTarget, 1e-9 ),
                           "comparison with reference solution failed" );
}

void checkVonMisesEnergies( const std::vector< double >& materialProperties, double kappaMin )
{
  using namespace Marmot;
  using namespace Marmot::Math;
  using namespace Marmot::Utilities;
  using namespace Marmot::Materials;

  const double        E = materialProperties[0], nu = materialProperties[1];
  const VonMisesModel material( materialProperties.data(), static_cast< int >( materialProperties.size() ), 0 );
  const Matrix6d      CelInv = ContinuumMechanics::Elasticity::Isotropic::complianceTensor( E, nu );

  struct MaterialPointState {
    MaterialHypoElastic::state3D state;
    double                       kappa   = 0.;
    Vector6d                     dPsi_dE = Vector6d::Zero();
    Vector6d                     dD_dE   = Vector6d::Zero();
  };

  // single stress update from a given state
  auto update = [&]( const MaterialPointState& old, const Vector6d& dStrain ) {
    MaterialPointState                  updated = old;
    Matrix6d                            tangent;
    const MaterialHypoElastic::timeInfo tInfo{ 0.0, 1.0 };
    updated.state.stateVars = &updated.kappa;
    material.computeStressWithEnergyTangents( updated.state, tangent, dStrain, tInfo, updated.dPsi_dE, updated.dD_dE );

    // the standard update gives the same state
    MaterialPointState plain = old;
    Matrix6d           tangentPlain;
    plain.state.stateVars = &plain.kappa;
    material.computeStress( plain.state, tangentPlain, dStrain, tInfo );
    throwExceptionOnFailure( plain.state.stress == updated.state.stress && tangentPlain == tangent &&
                               plain.kappa == updated.kappa &&
                               plain.state.elasticEnergyDensity == updated.state.elasticEnergyDensity &&
                               plain.state.dissipation == updated.state.dissipation,
                             "computeStress() and computeStressWithEnergyTangents() differ" );
    return updated;
  };

  // loading into softening, then elastic unloading
  const Vector6d     direction = ( Vector6d() << 1., -0.3, -0.2, 0.1, 0.4, -0.2 ).finished();
  MaterialPointState current;
  double             dissipationReference = 0.;
  for ( int i = 0; i < 60; i++ ) {
    const Vector6d           dStrain = ( i < 50 ? 1e-4 : -2e-5 ) * direction;
    const MaterialPointState updated = update( current, dStrain );
    const auto&              state   = updated.state;

    // plastic work sigma_{n+1} : dEpsP, with the plastic strain increment from the elastic law
    dissipationReference += state.stress.dot( dStrain - CelInv * ( state.stress - current.state.stress ) );

    throwExceptionOnFailure( checkIfEqual( state.elasticEnergyDensity,
                                           0.5 * state.stress.dot( CelInv * state.stress ),
                                           1e-10 ),
                             MakeString() << "elastic energy density wrong in increment " << i );
    throwExceptionOnFailure( checkIfEqual( state.dissipation, dissipationReference, 1e-9 ),
                             MakeString() << "dissipation wrong in increment " << i );

    // energy tangents against central differences wrt the strain increment
    Vector6d     dPsiNumerical, dDNumerical;
    const double h = 1e-9;
    for ( int j = 0; j < 6; j++ ) {
      Vector6d perturbation = Vector6d::Zero();
      perturbation( j )     = h;
      const auto plus       = update( current, dStrain + perturbation ).state;
      const auto minus      = update( current, dStrain - perturbation ).state;
      dPsiNumerical( j )    = ( plus.elasticEnergyDensity - minus.elasticEnergyDensity ) / ( 2. * h );
      dDNumerical( j )      = ( plus.dissipation - minus.dissipation ) / ( 2. * h );
    }
    const double scale = state.stress.cwiseAbs().maxCoeff();
    throwExceptionOnFailure( ( updated.dPsi_dE - dPsiNumerical ).cwiseAbs().maxCoeff() < 1e-5 * scale,
                             MakeString() << "elastic energy tangent wrong in increment " << i );
    throwExceptionOnFailure( ( updated.dD_dE - dDNumerical ).cwiseAbs().maxCoeff() < 1e-5 * scale,
                             MakeString() << "dissipation tangent wrong in increment " << i );

    current = updated;
  }
  throwExceptionOnFailure( current.kappa > kappaMin && current.state.dissipation > 0.,
                           "insufficient plastic flow in the energy test" );
}

void testVonMisesEnergies()
{
  // exponential softening fy = 200 - 150 (1 - exp(-500 kappa))
  checkVonMisesEnergies( { 210000., 0.3, 200., 0., -150., 500. }, 0. );
  // linear softening fy = 200 - 70000 kappa, beyond zero strength at kappa = 200 / 70000
  checkVonMisesEnergies( { 210000., 0.3, 200., -70000., 0., 1. }, 200. / 70000. );
}

// the explicit update gives the stress, energies and hardening of the implicit one; a zero increment has the stress as
// elastic energy tangent and no dissipation
void testVonMisesExplicitUpdateAndZeroIncrement()
{
  using namespace Marmot;
  using namespace Marmot::Math;
  using namespace Marmot::Materials;

  const std::vector< double >         props = { 210000., 0.3, 200., 0., -150., 500. };
  const VonMisesModel                 material( props.data(), static_cast< int >( props.size() ), 0 );
  const MaterialHypoElastic::timeInfo tInfo{ 0.0, 1.0 };
  const Vector6d                      direction = ( Vector6d() << 1., -0.3, -0.2, 0.1, 0.4, -0.2 ).finished();

  MaterialHypoElastic::state3D implicitState, explicitState;
  double                       kappaImplicit = 0., kappaExplicit = 0.;
  implicitState.stateVars = &kappaImplicit;
  explicitState.stateVars = &kappaExplicit;
  // elastic, then plastic increments
  for ( const double scale : { 1e-4, 1e-3, 1e-3 } ) {
    Matrix6d tangent;
    material.computeStress( implicitState, tangent, scale * direction, tInfo );
    material.computeStressExplicit( explicitState, scale * direction, tInfo );
    // the explicit return mapping iterates a fixed number of times: equal up to round-off
    const double stressScale = implicitState.stress.cwiseAbs().maxCoeff();
    throwExceptionOnFailure( checkIfEqual< double >( implicitState.stress,
                                                     explicitState.stress,
                                                     1e-12 * stressScale ) &&
                               checkIfEqual( kappaImplicit, kappaExplicit, 1e-15 ) &&
                               checkIfEqual( implicitState.elasticEnergyDensity,
                                             explicitState.elasticEnergyDensity,
                                             1e-12 * implicitState.elasticEnergyDensity ) &&
                               checkIfEqual( implicitState.dissipation,
                                             explicitState.dissipation,
                                             1e-12 * std::max( implicitState.dissipation, 1e-12 ) ),
                             "computeStressExplicit() differs from computeStress()" );
  }
  throwExceptionOnFailure( kappaExplicit > 0., "no plastic flow in the explicit update" );

  Matrix6d                           tangent;
  Vector6d                           dPsi_dE, dD_dE;
  const MaterialHypoElastic::state3D stateOld = implicitState;
  material.computeStressWithEnergyTangents( implicitState, tangent, Vector6d::Zero(), tInfo, dPsi_dE, dD_dE );
  throwExceptionOnFailure( implicitState.stress == stateOld.stress && dPsi_dE == stateOld.stress && dD_dE.isZero( 0. ),
                           "energy tangents of a zero increment wrong" );
}

int main()
{
  std::vector< std::function< void( void ) > > tests = { testVonMises,
                                                         testVonMisesCoordinateInvariance,
                                                         testVonMisesEnergies,
                                                         testVonMisesExplicitUpdateAndZeroIncrement };

  executeTestsAndCollectExceptions( tests );
  return 0;
}
