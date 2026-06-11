#include "Marmot/LinearViscoElasticInterface.h"
#include "Marmot/MarmotInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotTesting.h"

#include <Eigen/Dense>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

// Use namespaces for brevity
using namespace Marmot::Testing;

// Function to create a MarmotMaterialHypoElastic object
// Inputs:
// - materialName: The name of the material (e.g., "LINEARELASTIC")
// - materialProperties: Array of material parameters (e.g., Young's modulus, Poisson's ratio)
// - nMaterialProperties: Number of parameters in the materialProperties array
std::unique_ptr< MarmotInterfaceMaterialHypoElastic > createMarmotInterfaceMaterialHypoElastic(
  const std::string& materialName,
  const double*      materialProperties,
  int                nMaterialProperties )
{
  // Element label (arbitrary value)
  const int elLabel = 1;

  // Create the material object using Marmot's factory method
  auto mat = std::unique_ptr< MarmotInterfaceMaterialHypoElastic >( dynamic_cast< MarmotInterfaceMaterialHypoElastic* >(
    MarmotLibrary::MarmotInterfaceMaterialHypoElasticFactory::createMaterial( materialName,
                                                                              materialProperties,
                                                                              nMaterialProperties,
                                                                              elLabel ) ) );

  return mat; // Return the created material object
}

void computeStress( MarmotInterfaceMaterialHypoElastic& mat,
                    double*                             stateVars,
                    double*                             force,
                    double*                             surfaceStress,
                    double*                             Q_ij,
                    double*                             Z_ijkl,
                    double*                             H_ijk,
                    double*                             Y_ijkl,
                    const double*                       dU,
                    const double*                       dSurfaceStrain,
                    const double*                       normal,
                    const double                        timeOld,
                    const double                        dT )
{
  MarmotInterfaceMaterialHypoElastic::State         state{ force, surfaceStress, stateVars };
  MarmotInterfaceMaterialHypoElastic::Tangents      tangents{ Q_ij, Z_ijkl, H_ijk, Y_ijkl };
  MarmotInterfaceMaterialHypoElastic::Deformation   deformation{ dU, dSurfaceStrain, normal };
  MarmotInterfaceMaterialHypoElastic::TimeIncrement timeIncrement{ timeOld, dT };
  mat.computeStress( state, tangents, deformation, timeIncrement );
}

void testRejectsMultipleMaxwellElements()
{
  const double materialProperties[8] = { 1e4, 0.3, 1e-7, 1e-2, 1e-8, 2, 1e-2, 1e0 };

  bool rejected = false;
  try {
    createMarmotInterfaceMaterialHypoElastic( "LINEARVISCOELASTICINTERFACE", materialProperties, 8 );
  }
  catch ( const std::invalid_argument& ) {
    rejected = true;
  }

  throwExceptionOnFailure( rejected, "LinearViscoElasticInterface accepted more than one Maxwell element." );
}

// Function to test the viscoelastic interface material response for a displacement jump
void testForceMaterialResponse()
{
  // Define material parameters matching LinearViscoElasticInterface constructor
  // Indices: [0]: E_0, [1]: nu_0, [2]: h, [3]: m, [4]: n, [5]: nMaxwell, [6]: minTau, [7]: timeToDays
  // Padded to 16 elements for compatibility
  const double materialProperties[8] = { 1e4, 0.3, 1e-7, 1e-2, 1e-8, 1, 1e-2, 1e0 };
  const int    nMaterialProperties   = 8;

  // Create the material object
  auto mat = createMarmotInterfaceMaterialHypoElastic( "LINEARVISCOELASTICINTERFACE",
                                                       materialProperties,
                                                       nMaterialProperties );

  if ( !mat ) {
    throw std::runtime_error( "Material creation failed" );
  }

  // Assign state variables
  // number of required state vars
  int nStateVars = mat->getNumberOfRequiredStateVars();

  // initialize state vars
  Eigen::VectorXd stateVar( nStateVars );
  stateVar.setZero();
  mat->initializeYourself( stateVar.data(), nStateVars );

  // first increment ( load free )
  const double    timeOld = 0.0; // Previous time step
  Eigen::VectorXd time( 2 );
  time.setZero();
  double dT = 28.0; // time increment
  time[1] += dT;

  // Define initial force/stress state (set to zero) and strain increment
  double force[3]          = { 0, 0, 0 };
  double surface_stress[9] = { 0 };
  // Define matrices to store the tangent components
  double H_inv_ij[9]          = { 0 };
  double Z_ijkl[81]           = { 0 };
  double H_inv_nF_ijk[27]     = { 0 };
  double Yn_H_inv_Fn_ijkl[81] = { 0 };
  // Define zero displacement and zero surface strain initial increments
  const double dU1[6]               = { 0 };
  const double dSurface_strain1[18] = { 0 };
  // Define normal vector
  const double normal[3] = { 0, 0, 1 };

  // compute material response
  computeStress( *mat,
                 stateVar.data(),
                 force,
                 surface_stress,
                 H_inv_ij,
                 Z_ijkl,
                 H_inv_nF_ijk,
                 Yn_H_inv_Fn_ijkl,
                 dU1,
                 dSurface_strain1,
                 normal,
                 timeOld,
                 dT );

  // second increment ( load application )
  dT = 1e-6;
  time[1] += dT;
  const double dU2[6]               = { 0, 1e-3, 0, 0, 0, 0 };
  const double dSurface_strain2[18] = { 0 };

  // compute material response
  computeStress( *mat,
                 stateVar.data(),
                 force,
                 surface_stress,
                 H_inv_ij,
                 Z_ijkl,
                 H_inv_nF_ijk,
                 Yn_H_inv_Fn_ijkl,
                 dU2,
                 dSurface_strain2,
                 normal,
                 timeOld,
                 dT );

  // third increment ( constant strain, relaxation )
  dT = 100.;
  time[1] += dT;
  const double dU3[6]               = { 0, 0, 0, 0, 0, 0 };
  const double dSurface_strain3[18] = { 0 };
  // compute material response
  computeStress( *mat,
                 stateVar.data(),
                 force,
                 surface_stress,
                 H_inv_ij,
                 Z_ijkl,
                 H_inv_nF_ijk,
                 Yn_H_inv_Fn_ijkl,
                 dU3,
                 dSurface_strain3,
                 normal,
                 timeOld,
                 dT );

  // expected force and surface stress
  double forceTarget[3]          = { 0, 38461538.461538, 0 };
  double surface_stressTarget[9] = { 0, 0, 0, 0, 0, -3.84615384615385, 0, -3.84615384615385, 0 };

  // Convert to Eigen maps for easier comparison
  Eigen::Map< Eigen::Vector3d > forceVec( force );
  Eigen::Map< Eigen::Vector3d > forceTargetVec( forceTarget );
  Eigen::Map< Eigen::VectorXd > surface_stressVec( surface_stress, 9 );
  Eigen::Map< Eigen::VectorXd > surface_stressTargetVec( surface_stressTarget, 9 );

  // Compare the computed stress to the expected stress and throw an exception if they differ
  throwExceptionOnFailure( checkIfEqual< double >( forceVec, forceTargetVec, 1e-6 ),
                           "force computation failed for displacement jump in " + std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( checkIfEqual< double >( surface_stressVec, surface_stressTargetVec, 1e-6 ),
                           "surface stress computation failed for surface shear strain in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Function to test the viscoelastic interface material response for given surface strain
void testSurfaceStressMaterialResponse()
{
  // Define material parameters matching LinearViscoElasticInterface constructor
  // Indices: [0]: E_0, [1]: nu_0, [2]: h, [3]: m, [4]: n, [5]: nMaxwell, [6]: minTau, [7]: timeToDays
  // Padded to 16 elements for compatibility
  const double materialProperties[8] = { 1e8, 0.3, 1e-7, 1e-2, 1e-8, 1, 1e-2, 1e0 };
  const int    nMaterialProperties   = 8;

  // Create the material object
  auto mat = createMarmotInterfaceMaterialHypoElastic( "LINEARVISCOELASTICINTERFACE",
                                                       materialProperties,
                                                       nMaterialProperties );

  // Assign state variables
  // number of required state vars
  int nStateVars = mat->getNumberOfRequiredStateVars();

  // initialize state vars
  Eigen::VectorXd stateVar( nStateVars );
  stateVar.setZero();
  mat->initializeYourself( stateVar.data(), nStateVars );

  // first increment ( load free )
  const double    timeOld = 0.0; // Previous time step
  Eigen::VectorXd time( 2 );
  time.setZero();
  double dT = 28.0; // time increment
  time[1] += dT;

  // Define initial force/stress state (set to zero) and strain increment
  double force[3]          = { 0, 0, 0 };
  double surface_stress[9] = { 0 };
  // Define matrices to store the tangent components
  double H_inv_ij[9]          = { 0 };
  double Z_ijkl[81]           = { 0 };
  double H_inv_nF_ijk[27]     = { 0 };
  double Yn_H_inv_Fn_ijkl[81] = { 0 };
  // Define zero displacement and zero surface strain initial increments
  const double dU1[6]               = { 0 };
  const double dSurface_strain1[18] = { 0 };
  // Define normal vector
  const double normal[3] = { 0, 0, 1 };

  // compute material response
  computeStress( *mat,
                 stateVar.data(),
                 force,
                 surface_stress,
                 H_inv_ij,
                 Z_ijkl,
                 H_inv_nF_ijk,
                 Yn_H_inv_Fn_ijkl,
                 dU1,
                 dSurface_strain1,
                 normal,
                 timeOld,
                 dT );

  // second increment ( load application )
  dT = 1e-6;
  time[1] += dT;
  const double dU2[6]               = { 0 };
  const double dSurface_strain2[18] = { 0, 1e-1, 0, 1e-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

  // compute material response
  computeStress( *mat,
                 stateVar.data(),
                 force,
                 surface_stress,
                 H_inv_ij,
                 Z_ijkl,
                 H_inv_nF_ijk,
                 Yn_H_inv_Fn_ijkl,
                 dU2,
                 dSurface_strain2,
                 normal,
                 timeOld,
                 dT );

  // third increment ( constant strain, relaxation )
  dT = 100.;
  time[1] += dT;
  const double dU3[6]               = { 0 };
  const double dSurface_strain3[18] = { 0 };
  // compute material response
  computeStress( *mat,
                 stateVar.data(),
                 force,
                 surface_stress,
                 H_inv_ij,
                 Z_ijkl,
                 H_inv_nF_ijk,
                 Yn_H_inv_Fn_ijkl,
                 dU3,
                 dSurface_strain3,
                 normal,
                 timeOld,
                 dT );

  // expected force and surface stress
  double forceTarget[3]          = { 0 };
  double surface_stressTarget[9] = { 0., 0.384615, 0., 0.384615, 0., 0., 0., 0., 0. };

  // Convert to Eigen maps for easier comparison
  Eigen::Map< Eigen::Vector3d > forceVec( force );
  Eigen::Map< Eigen::Vector3d > forceTargetVec( forceTarget );
  Eigen::Map< Eigen::VectorXd > surface_stressVec( surface_stress, 9 );
  Eigen::Map< Eigen::VectorXd > surface_stressTargetVec( surface_stressTarget, 9 );

  // Compare the computed stress to the expected stress and throw an exception if they differ
  throwExceptionOnFailure( checkIfEqual< double >( forceVec, forceTargetVec, 1e-6 ),
                           "force computation failed for displacement jump in " + std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( checkIfEqual< double >( surface_stressVec, surface_stressTargetVec, 1e-6 ),
                           "surface stress computation failed for surface shear strain in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

int main()
{

  auto tests = std::vector< std::function< void() > >{
    testRejectsMultipleMaxwellElements,
    testForceMaterialResponse,         // test for force response
    testSurfaceStressMaterialResponse, // test for surface stress response
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
