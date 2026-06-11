#include "Marmot/LinearViscoElasticWiechert.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotTesting.h"

#include <Eigen/Dense>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

// Use namespaces for brevity
using namespace Marmot::Testing;

std::unique_ptr< MarmotMaterialHypoElastic > createMarmotMaterialHypoElastic( const std::string& materialName,
                                                                              const double*      materialProperties,
                                                                              int                nMaterialProperties )
{
  (void)materialName;
  return std::make_unique< Marmot::Materials::LinearViscoElasticWiechert >( materialProperties,
                                                                            nMaterialProperties,
                                                                            1 );
}

// Function to test the viscoelastic material response
void testStressMaterialResponse()
{
  // Define material parameters matching LinearViscoElasticWiechert constructor
  // Indices:
  // [0]: E_0
  // [1]: nu_0
  // [2]: m
  // [3]: n
  // [4]: nMaxwell
  // [5]: minTau
  // [6]: timeToDays
  const double materialProperties[7] = { 1e8, 0.3, 1e-2, 1e-8, 1, 1e-2, 1e0 };
  const int    nMaterialProperties   = 7;

  // Create the material object
  auto mat = createMarmotMaterialHypoElastic( "LINEARVISCOELASTICWIECHERT", materialProperties, nMaterialProperties );

  // Assign state variables
  const int nStateVars = mat->getNumberOfRequiredStateVars();
  if ( nStateVars != 6 ) {
    throw std::runtime_error( "Unexpected number of required state variables in " +
                              std::string( __PRETTY_FUNCTION__ ) );
  }

  Eigen::VectorXd stateVar( nStateVars );
  stateVar.setZero();

  MarmotMaterialHypoElastic::state3D state{ Marmot::Vector6d::Zero(), 0.0, stateVar.data() };

  MarmotMaterialHypoElastic::timeInfo timeInfo;

  // Tangent matrix required by the new MarmotMaterialHypoElastic API
  Marmot::Matrix6d D_ijkl = Marmot::Matrix6d::Zero();

  // Stress vector used for checking the final response
  Marmot::Vector6d stress = Marmot::Vector6d::Zero();

  // Time bookkeeping
  const double    timeOld = 0.0;
  Eigen::VectorXd time( 2 );
  time.setZero();

  // ---------------------------------------------------------------------------
  // First increment: load-free initialization
  // ---------------------------------------------------------------------------
  double dT = 28.0;
  time[1] += dT;

  Marmot::Vector6d dstrain1 = Marmot::Vector6d::Zero();

  timeInfo.time = timeOld;
  timeInfo.dT   = dT;

  mat->computeStress( state, D_ijkl, dstrain1, timeInfo );
  stress = state.stress;

  // ---------------------------------------------------------------------------
  // Second increment: load application
  // ---------------------------------------------------------------------------
  dT = 1e-6;
  time[1] += dT;

  Marmot::Vector6d dstrain2 = Marmot::Vector6d::Zero();
  dstrain2[5]               = 1e-1;

  timeInfo.time = timeOld;
  timeInfo.dT   = dT;

  mat->computeStress( state, D_ijkl, dstrain2, timeInfo );
  stress = state.stress;

  // ---------------------------------------------------------------------------
  // Third increment: constant strain, relaxation
  // ---------------------------------------------------------------------------
  dT = 100.0;
  time[1] += dT;

  Marmot::Vector6d dstrain3 = Marmot::Vector6d::Zero();

  timeInfo.time = timeOld;
  timeInfo.dT   = dT;

  mat->computeStress( state, D_ijkl, dstrain3, timeInfo );
  stress = state.stress;

  // Expected stress value
  Marmot::Vector6d stressTarget = Marmot::Vector6d::Zero();
  stressTarget[5]               = 3846153.831362;

  throwExceptionOnFailure( checkIfEqual< double >( stress, stressTarget, 1e-6 ),
                           "stress computation failed for displacement jump in " + std::string( __PRETTY_FUNCTION__ ) );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{
    testStressMaterialResponse,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
