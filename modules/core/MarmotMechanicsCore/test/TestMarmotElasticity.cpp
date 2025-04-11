#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>

using namespace Marmot::Testing;
using namespace Marmot::ContinuumMechanics::Elasticity::Isotropic;
using namespace Marmot::ContinuumMechanics::Elasticity::TransverseIsotropic;
using namespace Marmot::ContinuumMechanics::Elasticity::Orthotropic;

void testStiffnessTensor_Isotropic()
{
  // Test parameters
  double E  = 1000; // Young's modulus in Pa
  double nu = 0.25; // Poisson's ratio

  // Compute the computed stiffness matrix using the function
  Marmot::Matrix6d computedC = stiffnessTensor( E, nu );

  // The expected stiffness matrix (expectedC)
  Marmot::Matrix6d expectedC;
  // clang-format off
  expectedC << 
    1200,   400,    400,    0,      0,      0,
    400,    1200,   400,    0,      0,      0,
    400,    400,    1200,   0,      0,      0,
    0,      0,      0,      400,    0,      0,
    0,      0,      0,      0,      400,    0,
    0,      0,      0,      0,      0,      400;
  // clang-format on

  // Check if the computed stiffness matrix is equal to the expected stiffness matrix
  throwExceptionOnFailure( checkIfEqual< double >( computedC, expectedC, 1e-15 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testStiffnessTensorKG_Isotropic()
{
  // Test parameters for K (bulk modulus) and G (shear modulus)
  double K = 666.66666666666667; // Bulk modulus in Pa
  double G = 400;                // Shear modulus in Pa
  // E = 1000 and nu = 0.25, the same as before, if these value are used for K and G

  // Compute the stiffness matrix using stiffnessTensorKG
  Marmot::Matrix6d computedC_KG = stiffnessTensorKG( K, G );

  // Now compute the expected stiffness matrix manually
  // The same as before
  Marmot::Matrix6d expectedC_KG;
  // clang-format off
  expectedC_KG << 
    1200,   400,    400,    0,      0,      0,
    400,    1200,   400,    0,      0,      0,
    400,    400,    1200,   0,      0,      0,
    0,      0,      0,      400,    0,      0,
    0,      0,      0,      0,      400,    0,
    0,      0,      0,      0,      0,      400;
  // clang-format on

  // Check if the computed stiffness matrix using stiffnessTensorKG is equal to the expected stiffness matrix
  throwExceptionOnFailure( checkIfEqual< double >( computedC_KG, expectedC_KG, 1e-15 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testComplianceTensor_Isotropic()
{
  // Test parameters
  double E  = 1000; // Young's modulus in Pa
  double nu = 0.25; // Poisson's ratio

  // Compute the compliance matrix (inverse of the stiffness matrix)
  Marmot::Matrix6d computedCInv = complianceTensor( E, nu );

  // The expected compliance matrix (calcualted manually)
  Marmot::Matrix6d expectedCInv;
  // clang-format off
  expectedCInv << 
    0.001,      -0.00025,   -0.00025,   0,          0,          0,
    -0.00025,   0.001,      -0.00025,   0,          0,          0,
    -0.00025,   -0.00025,   0.001,      0,          0,          0,
    0,          0,          0,          0.0025,     0,          0,
    0,          0,          0,          0,          0.0025,     0,
    0,          0,          0,          0,          0,          0.0025;
  // clang-format on

  // Check if the computed compliance matrix is equal to the expected compliance matrix
  throwExceptionOnFailure( checkIfEqual< double >( computedCInv, expectedCInv, 1e-15 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testComplianceTensor_TransverseIsotropic()
{
  // Test parameters for transverse isotropic material
  double E1   = 1000; // Young's modulus in the 1-direction (Pa)
  double E2   = 500;  // Young's modulus in the 2-direction (Pa)
  double nu12 = 0.5;  // Poisson's ratio (12-direction)
  double nu23 = 0.25; // Poisson's ratio (23-direction)
  double G12  = 400;  // Shear modulus in the 12-plane (Pa)

  // Compute the compliance matrix (inverse of the stiffness matrix) using the function
  Marmot::Matrix6d computedCInv = complianceTensor( E1, E2, nu12, nu23, G12 );

  // The expected compliance matrix (calculated manually)
  // double G23 = 200;
  Marmot::Matrix6d expectedCInv;
  // clang-format off
  expectedCInv << 
    0.001,      -0.001,     -0.001,     0,          0,          0,
    -0.001,     0.002,      -0.0005,    0,          0,          0,
    -0.001,     -0.0005,    0.002,      0,          0,          0,
    0,          0,          0,          0.0025,     0,          0,
    0,          0,          0,          0,          0.0025,     0,
    0,          0,          0,          0,          0,          0.005;
  // clang-format on

  // Check if the computed compliance matrix is equal to the expected compliance matrix
  throwExceptionOnFailure( checkIfEqual< double >( computedCInv, expectedCInv, 1e-15 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testStiffnessTensor_TransverseIsotropic()
{
  // Test parameters for transverse isotropic material
  double E1   = 1000; // Young's modulus in the 1-direction (Pa)
  double E2   = 500;  // Young's modulus in the 2-direction (Pa)
  double nu12 = 0.5;  // Poisson's ratio (12-direction)
  double nu23 = 0.25; // Poisson's ratio (23-direction)
  double G12  = 400;  // Shear modulus in the 12-plane (Pa)

  // Compute the stiffness matrix using the function
  Marmot::Matrix6d computedC = stiffnessTensor( E1, E2, nu12, nu23, G12 );

  // The expected stiffness matrix (calculated manually)
  Marmot::Matrix6d expectedC;
  // clang-format off
  expectedC << 
    -3000,      -2000,      -2000,      0,      0,      0,
    -2000,      -800,       -1200,      0,      0,      0,
    -2000,      -1200,      -800,       0,      0,      0,
    0,          0,          0,          400,    0,      0,
    0,          0,          0,          0,      400,    0,
    0,          0,          0,          0,      0,      200;
  // clang-format on

  // Check if the computed stiffness matrix is equal to the expected stiffness matrix
  throwExceptionOnFailure( checkIfEqual< double >( computedC, expectedC, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testComplianceTensor_Orthotropic()
{
  // Test parameters for orthotropic material
  double E1   = 1000; // Young's modulus in the 1-direction (Pa)
  double E2   = 500;  // Young's modulus in the 2-direction (Pa)
  double E3   = 200;  // Young's modulus in the 3-direction (Pa)
  double nu12 = 0.5;  // Poisson's ratio (12-direction)
  double nu23 = 0.25; // Poisson's ratio (23-direction)
  double nu13 = 0.1;  // Poisson's ratio (13-direction)
  double G12  = 400;  // Shear modulus in the 12-plane (Pa)
  double G23  = 250;  // Shear modulus in the 23-plane (Pa)
  double G31  = 200;  // Shear modulus in the 31-plane (Pa)

  // Compute the compliance matrix (inverse of the stiffness matrix) using the function
  Marmot::Matrix6d computedCInv = complianceTensor( E1, E2, E3, nu12, nu23, nu13, G12, G23, G31 );

  // The expected compliance matrix (calculated manually)
  Marmot::Matrix6d expectedCInv;
  // clang-format off
  expectedCInv << 
    0.001,      -0.001,     -0.0005,    0,          0,          0,
    -0.001,     0.002,      -0.00125,   0,          0,          0,
    -0.0005,    -0.00125,   0.005,      0,          0,          0,
    0,          0,          0,          0.0025,     0,          0,
    0,          0,          0,          0,          0.005,      0,
    0,          0,          0,          0,          0,          0.004;
  // clang-format on

  // Check if the computed compliance matrix is equal to the expected compliance matrix
  throwExceptionOnFailure( checkIfEqual< double >( computedCInv, expectedCInv, 1e-15 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testStiffnessTensor_Orthotropic()
{
  // Test parameters for orthotropic material
  double E1   = 1000; // Young's modulus in the 1-direction (Pa)
  double E2   = 500;  // Young's modulus in the 2-direction (Pa)
  double E3   = 200;  // Young's modulus in the 3-direction (Pa)
  double nu12 = 0.5;  // Poisson's ratio (12-direction)
  double nu23 = 0.25; // Poisson's ratio (23-direction)
  double nu13 = 0.1;  // Poisson's ratio (13-direction)
  double G12  = 400;  // Shear modulus in the 12-plane (Pa)
  double G23  = 250;  // Shear modulus in the 23-plane (Pa)
  double G31  = 200;  // Shear modulus in the 31-plane (Pa)

  // Compute the stiffness matrix using the function
  Marmot::Matrix6d computedC = stiffnessTensor( E1, E2, E3, nu12, nu23, nu13, G12, G23, G31 );

  // The expected stiffness matrix (calculated manually)
  Marmot::Matrix6d expectedC;
  // clang-format off
  expectedC << 
    5000,                       3333.3333333333333332,      1333.3333333333333333,      0,      0,      0,
    3333.3333333333333332,      2814.8148148148148147,      1037.037037037037037,       0,      0,      0,
    1333.3333333333333333,      1037.037037037037037,       592.59259259259259259,      0,      0,      0,
    0,                          0,                          0,                          400,    0,      0,
    0,                          0,                          0,                          0,      200,    0,
    0,                          0,                          0,                          0,      0,      250;
  // clang-format on

  // Check if the computed stiffness matrix is equal to the expected stiffness matrix
  throwExceptionOnFailure( checkIfEqual< double >( computedC, expectedC, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

int main()
{
  testStiffnessTensor_Isotropic();
  testStiffnessTensorKG_Isotropic();
  testComplianceTensor_Isotropic();
  testComplianceTensor_TransverseIsotropic();
  testStiffnessTensor_TransverseIsotropic();
  testComplianceTensor_Orthotropic();
  testStiffnessTensor_Orthotropic();
  return 0;
}
