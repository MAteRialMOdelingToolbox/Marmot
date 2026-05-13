#include "Marmot/LinearViscoElasticWiechert.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotWiechert.h"
#include <Eigen/Dense>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <memory>

// Use namespaces for brevity
using namespace Marmot::Testing;
using namespace Marmot::ContinuumMechanics::Elasticity::Isotropic;
using namespace Marmot::ContinuumMechanics::Elasticity::TransverseIsotropic;

// Function to create a MarmotMaterialHypoElastic object
// Inputs:
// - materialName: The name of the material (e.g., "LINEARELASTIC")
// - materialProperties: Array of material parameters (e.g., Young's modulus, Poisson's ratio)
// - nMaterialProperties: Number of parameters in the materialProperties array
std::unique_ptr< MarmotMaterialHypoElastic > createMarmotMaterialHypoElastic( const std::string& materialName,
                                                                              const double*      materialProperties,
                                                                              int                nMaterialProperties )
{
  (void)materialName;
  return std::make_unique< Marmot::Materials::LinearViscoElasticWiechert >( materialProperties,
                                                                            nMaterialProperties,
                                                                            1 );
}

// Function to test the viscoelastic interface material response for given surface strain
void testStressMaterialResponse()
{
  // Define material parameters matching LinearViscoElasticInterface constructor
  // Indices: [0]: E_0, [1]: nu_0, [2]: m, [3]: n, [4]: nMaxwell, [5]: minTau, [6]: timeToDays
  // Padded to 16 elements for compatibility
  const double materialProperties[7] = { 1e8, 0.3, 1e-2, 1e-8, 1, 1e-2, 1e0 };
  const int    nMaterialProperties   = 7;

  // Create the material object
  auto mat = createMarmotMaterialHypoElastic( "LINEARVISCOELASTICWIECHERT", materialProperties, nMaterialProperties );

  // Assign state variables
  // number of required state vars
  int nStateVars = mat->getNumberOfRequiredStateVars();
  if ( nStateVars != 6 ) {
    throw std::runtime_error( "Unexpected number of required state variables in " +
                              std::string( __PRETTY_FUNCTION__ ) );
  }

  // initialize state vars
  Eigen::VectorXd stateVar( nStateVars );
  stateVar.setZero();
  MarmotMaterialHypoElastic::state3D state{ Marmot::Vector6d::Zero(), 0.0, stateVar.data() };

  MarmotMaterialHypoElastic::timeInfo timeInfo;

  // first increment ( load free )
  const double    timeOld = 0.0; // Previous time step
  Eigen::VectorXd time( 2 );
  time.setZero();
  double dT = 28.0; // time increment
  time[1] += dT;

  // Define initial force/stress state (set to zero) and strain increment
  double stress[6] = { 0 };
  // Define matrices to store the tangent components
  double D_ijkl[36] = { 0 };
  // Define zero displacement and zero surface strain initial increments
  const double dstrain1[6] = { 0 };

  // compute material response
  timeInfo.time = timeOld;
  timeInfo.dT   = dT;
  mat->computeStress( state, D_ijkl, dstrain1, timeInfo );
  for ( int i = 0; i < 6; ++i )
    stress[i] = state.stress[i];

  // second increment ( load application )
  dT = 1e-6;
  time[1] += dT;
  const double dstrain2[6] = { 0, 0, 0, 0, 0, 1e-1 };

  // compute material response
  timeInfo.time = timeOld;
  timeInfo.dT   = dT;
  mat->computeStress( state, D_ijkl, dstrain2, timeInfo );
  for ( int i = 0; i < 6; ++i )
    stress[i] = state.stress[i];

  // third increment ( constant strain, relaxation )
  dT = 100.;
  time[1] += dT;
  const double dstrain3[6] = { 0 };

  // compute material response
  timeInfo.time = timeOld;
  timeInfo.dT   = dT;
  mat->computeStress( state, D_ijkl, dstrain3, timeInfo );
  for ( int i = 0; i < 6; ++i )
    stress[i] = state.stress[i];
  // Set the expected force and surface stress explicitly
  // Use the actual value previously printed by the test
  double stressTarget[6] = { 0., 0., 0., 0., 0., 3807692.305921 };

  // Convert to Eigen maps for easier comparison
  Eigen::Map< Eigen::VectorXd > stressVec( stress, 6 );
  Eigen::Map< Eigen::VectorXd > stressTargetVec( stressTarget, 6 );
  // Compare the computed stress to the expected stress and throw an exception if they differ
  throwExceptionOnFailure( checkIfEqual< double >( stressVec, stressTargetVec, 1e-6 ),
                           "stress computation failed for displacement jump in " + std::string( __PRETTY_FUNCTION__ ) );
}

int main()
{

  auto tests = std::vector< std::function< void() > >{
    testStressMaterialResponse, // test for surface stress response
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
