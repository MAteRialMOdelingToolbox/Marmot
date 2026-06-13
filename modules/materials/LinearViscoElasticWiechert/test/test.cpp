#include "Marmot/LinearViscoElasticWiechert.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotWiechert.h"

#include <Eigen/Dense>

using namespace Marmot::Testing;

namespace {

  void testGeneralizedMaxwellPowerLawRelaxation()
  {
    // [equilibrium E, nu, power-law m, power-law n, nMaxwell, minTau, timeToDays]
    const double                                  properties[7] = { 1e8, 0.3, 2e7, 0.25, 6., 1e-4, 1. };
    Marmot::Materials::LinearViscoElasticWiechert material( properties, 7, 1 );

    const int nStateVars = material.getNumberOfRequiredStateVars();
    throwExceptionOnFailure( nStateVars == 36, "Incorrect generalized-Maxwell state-variable count." );

    Eigen::VectorXd stateVars( nStateVars );
    material.initializeYourself( stateVars.data(), nStateVars );

    MarmotMaterialHypoElastic::state3D state{ Marmot::Vector6d::Zero(), 0.0, 0.0, stateVars.data() };
    Marmot::Matrix6d                   tangent = Marmot::Matrix6d::Zero();

    Marmot::Vector6d strainIncrement = Marmot::Vector6d::Zero();
    strainIncrement[0]               = 1e-3;
    material.computeStress( state, tangent, strainIncrement, { 0.0, 1e-10 } );
    const double stressAfterLoading = state.stress[0];

    const Eigen::Map< const Marmot::Materials::Wiechert::StateVarMatrix > branchStates( stateVars.data(), 6, 6 );
    throwExceptionOnFailure( !branchStates.col( 0 ).isApprox( branchStates.col( 5 ) ),
                             "Power-law approximation produced identical Maxwell branches." );

    material.computeStress( state, tangent, Marmot::Vector6d::Zero(), { 1e-10, 1e-2 } );
    const double stressAfterShortHold = state.stress[0];
    material.computeStress( state, tangent, Marmot::Vector6d::Zero(), { 1e-2, 10.0 } );
    const double stressAfterLongHold = state.stress[0];

    throwExceptionOnFailure( stressAfterLoading > stressAfterShortHold && stressAfterShortHold > stressAfterLongHold &&
                               stressAfterLongHold > 0.0,
                             "Generalized Maxwell material does not relax monotonically." );
  }

  void testDensity()
  {
    const double                                  properties[8] = { 1e8, 0.3, 2e7, 0.25, 4., 1e-4, 1., 2400. };
    Marmot::Materials::LinearViscoElasticWiechert material( properties, 8, 1 );
    throwExceptionOnFailure( checkIfEqual( material.getDensity( nullptr ), properties[7] ),
                             "LinearViscoElasticWiechert density is incorrect." );
  }

} // namespace

int main()
{
  executeTestsAndCollectExceptions( { testGeneralizedMaxwellPowerLawRelaxation, testDensity } );
  return 0;
}
