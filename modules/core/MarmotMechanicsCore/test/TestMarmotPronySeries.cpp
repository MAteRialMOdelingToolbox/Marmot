#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotPronySeries.h"
#include "Marmot/MarmotTesting.h"

using namespace Marmot::Testing;

void testPronySeriesWithZeroMaxwellElements()
{
  /*
   * This test checks if the Elastic case is recovered
   * when the number of Maxwell elements is set to zero.
   */

  // define the number of Maxwell elements
  const int nMaxwell = 0;

  using namespace Marmot::ContinuumMechanics;
  const auto C0 = Elasticity::Isotropic::stiffnessTensor( 100, 0.3 );

  using namespace Marmot::ContinuumMechanics::Viscoelasticity::PronySeries;
  Properties     props = { nMaxwell, C0, Marmot::Vector6d::Zero(), Marmot::Matrix6d::Zero() };
  StateVarMatrix stateVars;

  // define the strain increment, stress, and stiffness
  Marmot::Vector6d dStrain   = Vector6d::Zero();
  Marmot::Vector6d stress    = Vector6d::Zero();
  Marmot::Matrix6d stiffness = Matrix6d::Zero();

  // set strain increment
  dStrain( 0 ) = 1.0e-3;
  dStrain( 1 ) = 2.0e-3;
  dStrain( 2 ) = 4.0e-3;
  dStrain( 3 ) = 3.0e-3;
  dStrain( 4 ) = 5.0e-3;
  dStrain( 5 ) = 6.0e-3;

  evaluatePronySeries( props, stress, stiffness, stateVars, dStrain, 1.0, false );

  throwExceptionOnFailure( checkIfEqual< double >( stress, C0 * dStrain ),
                           MakeString() << __PRETTY_FUNCTION__ << "stress computation failed" );

  throwExceptionOnFailure( checkIfEqual< double >( stiffness, C0 ),
                           MakeString() << __PRETTY_FUNCTION__ << "stiffness computation failed" );
}

void testPronySeriesWithOneMaxwellElement()
{
  /*
   * This test checks if the Prony series implementation
   * is correct for the case of one Maxwell element.
   */

  // define the number of Maxwell elements
  const int nMaxwell = 1;

  using namespace Marmot::ContinuumMechanics;
  const auto C0 = Elasticity::Isotropic::stiffnessTensor( 100, 0.0 ); // zero Poisson's ratio for simplicity

  using namespace Marmot::ContinuumMechanics::Viscoelasticity::PronySeries;
  Properties        props        = { nMaxwell, C0, Marmot::Vector6d::Zero(), Marmot::Matrix6d::Zero() };
  StateVarMatrix    stateVars    = Marmot::Matrix6d::Zero();
  mapStateVarMatrix mapStateVars = Eigen::Map< StateVarMatrix >( stateVars.data(), stateVars.rows(), stateVars.cols() );

  // define the strain increment, stress, and stiffness
  Marmot::Vector6d dStrain   = Vector6d::Zero();
  Marmot::Vector6d stress    = Vector6d::Zero();
  Marmot::Matrix6d stiffness = Matrix6d::Zero();

  // set strain increment
  dStrain( 0 ) = 1.0e-3;
  dStrain( 1 ) = 2.0e-3;
  dStrain( 2 ) = 3.0e-3;
  dStrain( 3 ) = 4.0e-3;
  dStrain( 4 ) = 5.0e-3;
  dStrain( 5 ) = 6.0e-3;

  // define the Prony series parameters
  Marmot::Matrix6d relaxationTimes = Matrix6d::Zero();
  relaxationTimes.setConstant( 1.0 );
  Marmot::Matrix6d C1 = C0 * 0.1;

  props.pronyRelaxationTimes = relaxationTimes;
  props.pronyStiffnesses     = C1;

  // test for a huge time step, where the Prony series should converge to the Elastic case
  evaluatePronySeries( props, stress, stiffness, stateVars, dStrain, 1e10, false );
  throwExceptionOnFailure( checkIfEqual< double >( stress, C0 * dStrain, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__ << "stress computation failed (t -> inf)" );

  // testfor zero time step, where the Prony series should converge to the Elastic case
  stress.setZero();
  evaluatePronySeries( props, stress, stiffness, stateVars, dStrain, 0, true );
  throwExceptionOnFailure( checkIfEqual< double >( stress, C0 * dStrain, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__ << "stress computation failed (dT = 0)" );

  // test for a small time step
  const double     dt           = 0.2;
  Marmot::Vector6d stressTarget = 2 * C0 * dStrain + ( stress - C0 * dStrain ) * std::exp( -dt ) +
                                  ( 1. - std::exp( -dt ) ) / dt * C1 * dStrain;
  evaluatePronySeries( props, stress, stiffness, stateVars, dStrain, dt, true );

  throwExceptionOnFailure( checkIfEqual< double >( stress, stressTarget, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__ << "stress computation failed (loading dT != 0)" );

  // new stress target for zero strain increment (relaxation)
  stressTarget = 2 * C0 * dStrain + ( stress - 2 * C0 * dStrain ) * std::exp( -dt );
  dStrain.setZero();
  evaluatePronySeries( props, stress, stiffness, stateVars, Vector6d::Zero(), dt, false );

  throwExceptionOnFailure( checkIfEqual< double >( stress, stressTarget, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__ << "stress computation failed (relaxation)" );
}

int main()
{
  testPronySeriesWithZeroMaxwellElements();

  testPronySeriesWithOneMaxwellElement();

  return 0;
}
