#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotInterfaceMaterialHelperFunctions.h"
#include "Marmot/MarmotInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <iterator>

// Use namespaces for brevity
using namespace Marmot::Testing;
using namespace Marmot::ContinuumMechanics::Elasticity::Isotropic;
using namespace Marmot::ContinuumMechanics::Elasticity::TransverseIsotropic;
using namespace Marmot::Materials::InterfaceMaterialHelperFunctions;

// Function to create a MarmotInterfaceMaterialHypoElastic object via the
// interface-material factory.
// Inputs:
// - materialName: Registered interface material name (e.g. "LINEARELASTICINTERFACEBMGU")
// - materialProperties: Array of interface material parameters
// - nMaterialProperties: Number of entries in materialProperties
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

// Function to test the viscoelastic interface material response for a displacement jump
void testForceMaterialResponse()
{
  // Define material parameters (Young's modulus and Poisson's ratio)
  // E_0: Youngs modulus of interphase
  // nu_0: Poisson's ratio of interphase
  // h : thickness of the interphase
  // dummy : placeholder for 8th parameter
  //                                     E_0, nu_0,     h, dummy
  const double materialProperties[8] = { 2e4, 0.3, 2e4, 0.3, 1e4, 0.3, 1e-7, 0.0 };
  const int    nMaterialProperties   = 8;

  // Create the material object
  auto mat = createMarmotInterfaceMaterialHypoElastic( "LINEARELASTICINTERFACEBMGU",
                                                       materialProperties,
                                                       nMaterialProperties );
  // Assign state variables
  // number of required state vars
  int nStateVars = mat->getNumberOfRequiredStateVars();

  // initialize state vars
  Eigen::VectorXd stateVar( nStateVars );
  stateVar.setZero();
  mat->assignStateVars( stateVar.data(), nStateVars );

  // Define initial force/stress state (set to zero) and strain increment
  double force[3]          = { 0, 0, 0 };
  double surface_stress[9] = { 0 };
  // Define a matrix to store the tangent stiffness (stress-strain relation)
  double H_inv_ij[3 * 3]                 = { { 0 } };
  double Z_ijkl[3 * 3 * 3 * 3]           = { { 0 } };
  double H_inv_nF_ijk[3 * 3 * 3]         = { { 0 } };
  double Yn_H_inv_Fn_ijkl[3 * 3 * 3 * 3] = { { 0 } };
  // Define displacement and surface strain increments
  //  Apply a small displacement increment on the top surface
  const double dU[6]               = { 0, 1e-3, 0, 0, 0, 0 };
  const double dSurface_strain[18] = { 0 };
  // Define normal vector
  const double normal[3] = { 0, 0, 1 };

  // Define time parameters for the material response calculation
  const double timeOld = 0.0; // Previous time step
  const double dT      = 1.0; // Time increment
  double       pNewDT;        // Placeholder for the new time increment

  // Compute the stress response of the material
  mat->computeStress( force,
                      surface_stress,
                      H_inv_ij,
                      Z_ijkl,
                      H_inv_nF_ijk,
                      Yn_H_inv_Fn_ijkl,
                      dU,
                      dSurface_strain,
                      normal,
                      &timeOld,
                      dT,
                      pNewDT );

  // Define the expected stress values for the applied strain increment
  double forceTarget[3]          = { 0, 76923076.9230769, 0 };
  double surface_stressTarget[9] = { 0., 0., 0., 0., 0., 0., 0., 0., 0. };
  // Convert to Eigen maps for easier comparison
  Eigen::Map< Eigen::Vector3d > forceVec( force );
  Eigen::Map< Eigen::Vector3d > forceTargetVec( forceTarget );
  Eigen::Map< Eigen::VectorXd > surface_stressVec( surface_stress, 9 );
  Eigen::Map< Eigen::VectorXd > surface_stressTargetVec( surface_stressTarget, 9 );

  // Print computed values with full precision for verification
  std::cout << std::setprecision( 15 ) << "Computed force: [" << force[0] << ", " << force[1] << ", " << force[2] << "]"
            << std::endl;
  std::cout << std::setprecision( 15 ) << "Computed surface_stress: [";
  for ( int i = 0; i < 9; i++ ) {
    std::cout << surface_stress[i];
    if ( i < 8 )
      std::cout << ", ";
  }
  std::cout << "]" << std::endl;

  // Compare the computed stress to the expected stress and throw an exception if they differ
  throwExceptionOnFailure( checkIfEqual< double >( forceVec, forceTargetVec, 1e-10 ),
                           "force computation failed for displacement jump in " + std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( checkIfEqual< double >( surface_stressVec, surface_stressTargetVec, 1e-10 ),
                           "surface stress computation failed for surface shear strain in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Function to test the viscoelastic interface material response for given surface strain
void testSurfaceStressMaterialResponse()
{
  // Define material parameters (Young's modulus and Poisson's ratio)
  // E_0: Youngs modulus of interphase
  // nu_0: Poisson's ratio of interphase
  // h : thickness of the interphase
  // dummy : placeholder for 8th parameter
  //                                     E_0, nu_0,     h, dummy
  const double materialProperties[8] = { 2e4, 0.3, 2e4, 0.3, 1e4, 0.3, 1e-7, 0.0 };
  const int    nMaterialProperties   = 8;

  // Create the material object
  auto mat = createMarmotInterfaceMaterialHypoElastic( "LINEARELASTICINTERFACEBMGU",
                                                       materialProperties,
                                                       nMaterialProperties );
  // Assign state variables
  // number of required state vars
  int nStateVars = mat->getNumberOfRequiredStateVars();

  // initialize state vars
  Eigen::VectorXd stateVar( nStateVars );
  stateVar.setZero();
  mat->assignStateVars( stateVar.data(), nStateVars );
  // Define initial force/stress state (set to zero) and strain increment
  double force[3]          = { 0, 0, 0 };
  double surface_stress[9] = { 0 };
  // Define matrices to store the material response outputs
  double H_inv_ij[21]              = { 0 };
  double Z_ijkl[21 * 21]           = { 0 };
  double H_inv_nF_ijk[21 * 3]      = { 0 };
  double Yn_H_inv_Fn_ijkl[21 * 21] = { 0 };
  // Define displacement and surface strain increments
  //  Apply a small displacement increment on the top surface
  const double dU[6]               = { 0, 0, 0, 0, 0, 0 };
  const double dSurface_strain[18] = { 0, 1e-3, 0, 1e-3, 0, 0, 0, 0, 0, 0, 1e-3, 0, 1e-3, 0, 0, 0, 0, 0 };
  // Define normal vector
  const double normal[3] = { 0, 0, 1 };

  // Define time parameters for the material response calculation
  const double timeOld = 0.0; // Previous time step
  const double dT      = 1.0; // Time increment
  double       pNewDT;        // Placeholder for the new time increment

  // Compute the stress response of the material
  mat->computeStress( force,
                      surface_stress,
                      H_inv_ij,
                      Z_ijkl,
                      H_inv_nF_ijk,
                      Yn_H_inv_Fn_ijkl,
                      dU,
                      dSurface_strain,
                      normal,
                      &timeOld,
                      dT,
                      pNewDT );

  double forceTarget[3]          = { 0, 0, 0 };
  double surface_stressTarget[9] = { 0, 7.69230769230769e-07, 0, 7.69230769230769e-07, 0, 0, 0, 0, 0 };
  // Convert to Eigen maps for easier comparison
  Eigen::Map< Eigen::Vector3d > forceVec( force );
  Eigen::Map< Eigen::Vector3d > forceTargetVec( forceTarget );
  Eigen::Map< Eigen::VectorXd > surface_stressVec( surface_stress, 9 );
  Eigen::Map< Eigen::VectorXd > surface_stressTargetVec( surface_stressTarget, 9 );

  // Print computed values with full precision for verification
  std::cout << std::setprecision( 15 ) << "Computed force (test2): [" << force[0] << ", " << force[1] << ", "
            << force[2] << "]" << std::endl;
  std::cout << std::setprecision( 15 ) << "Computed surface_stress (test2): [";
  for ( int i = 0; i < 9; i++ ) {
    std::cout << surface_stress[i];
    if ( i < 8 )
      std::cout << ", ";
  }
  std::cout << "]" << std::endl;

  // Compare the computed stress to the expected stress and throw an exception if they differ
  throwExceptionOnFailure( checkIfEqual< double >( forceVec, forceTargetVec, 1e-10 ),
                           "force computation failed for displacement jump in " + std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( checkIfEqual< double >( surface_stressVec, surface_stressTargetVec, 1e-10 ),
                           "surface stress computation failed for surface shear strain in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

int main()
{

  auto tests = std::vector< std::function< void() > >{
    testForceMaterialResponse,         // test for force response
    testSurfaceStressMaterialResponse, // test for surface stress response
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
