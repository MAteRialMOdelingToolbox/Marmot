#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotLowerDimensionalStress.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotVoigt.h"

/*
 * ASSIGNEE: Johannes Thiel
 * * TODO: Test LowerDimensionalStress implementation for multiple cases and edge cases
 */

using namespace Eigen;
using namespace Marmot::ContinuumMechanics::Elasticity::Isotropic;
using namespace Marmot::ContinuumMechanics::PlaneStrain;
using namespace Marmot::ContinuumMechanics::PlaneStress;
using namespace Marmot::ContinuumMechanics::UniaxialStress;
using namespace Marmot::Testing;

Marmot::Matrix6d create_C_Matrix()
{
  // Definition of C (isotropic material)
  double           E  = 1000;
  double           nu = 0.25;
  Marmot::Matrix6d C  = stiffnessTensor( E, nu );
  return C;
}

void test_getUniaxialStressTangent()
{
  Marmot::Matrix6d C = create_C_Matrix();

  // Compute the result using the getUniaxialStressTangent function
  double computedResult = getUniaxialStressTangent( C );

  // Expected uniaxial stress tangent for isotropic material = E
  double expectedResult = 1000.0;

  throwExceptionOnFailure( checkIfEqual( computedResult, expectedResult, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void test_reduce3D_dStress_dDeformationGradient()
{
  // Create Tensor
  Marmot::EigenTensors::Tensor633d inputTensor;
  // Set all values to zero initially
  inputTensor.setZero();
  // Populate the tensor with some non-zero values to represent derivatives
  inputTensor( 0, 0, 0 ) = 1.0;  // σ_xx
  inputTensor( 0, 1, 0 ) = 2.0;  // σ_xx
  inputTensor( 0, 0, 1 ) = 3.0;  // σ_xx
  inputTensor( 0, 1, 1 ) = 4.0;  // σ_xx
  inputTensor( 1, 0, 0 ) = 5.0;  // σ_yy
  inputTensor( 1, 1, 0 ) = 6.0;  // σ_yy
  inputTensor( 1, 0, 1 ) = 7.0;  // σ_yy
  inputTensor( 1, 1, 1 ) = 8.0;  // σ_yy
  inputTensor( 3, 0, 0 ) = 9.0;  // σ_xy
  inputTensor( 3, 1, 0 ) = 10.0; // σ_xy
  inputTensor( 3, 0, 1 ) = 11.0; // σ_xy
  inputTensor( 3, 1, 1 ) = 12.0; // σ_xy
  inputTensor( 2, 0, 0 ) = 13.0; // σ_zz
  inputTensor( 2, 1, 0 ) = 14.0; // σ_zz
  inputTensor( 2, 0, 1 ) = 15.0; // σ_zz
  inputTensor( 2, 1, 1 ) = 16.0; // σ_zz
  inputTensor( 4, 0, 0 ) = 17.0; // σ_xz
  inputTensor( 4, 1, 2 ) = 18.0; // σ_xz
  inputTensor( 5, 0, 1 ) = 19.0; // σ_yz
  inputTensor( 5, 1, 2 ) = 20.0; // σ_yz

  /*
  inputTensor:
   1  2  0  3  4  0  0  0  0
   5  6  0  7  8  0  0  0  0
  13 14  0 15 16  0  0  0  0
   9 10  0 11 12  0  0  0  0
  17  0  0  0  0  0  0 18  0
   0  0  0 19  0  0  0 20  0
  */

  // Compute the result using the reduce3D_dStress_dDeformationGradient function
  Marmot::EigenTensors::Tensor322d computedResult = reduce3D_dStress_dDeformationGradient( inputTensor );

  // Expected results
  Marmot::EigenTensors::Tensor322d expectedResult;
  expectedResult( 0, 0, 0 ) = 1.0;  // σ_xx
  expectedResult( 0, 1, 0 ) = 2.0;  // σ_xx
  expectedResult( 0, 0, 1 ) = 3.0;  // σ_xx
  expectedResult( 0, 1, 1 ) = 4.0;  // σ_xx
  expectedResult( 1, 0, 0 ) = 5.0;  // σ_yy
  expectedResult( 1, 1, 0 ) = 6.0;  // σ_yy
  expectedResult( 1, 0, 1 ) = 7.0;  // σ_yy
  expectedResult( 1, 1, 1 ) = 8.0;  // σ_yy
  expectedResult( 2, 0, 0 ) = 9.0;  // σ_xy
  expectedResult( 2, 1, 0 ) = 10.0; // σ_xy
  expectedResult( 2, 0, 1 ) = 11.0; // σ_xy
  expectedResult( 2, 1, 1 ) = 12.0; // σ_xy

  /*
  expectedResult:
  1  2  3  4
  5  6  7  8
  9 10 11 12
  */

  throwExceptionOnFailure( checkIfEqual< double, 3, 2, 2 >( computedResult, expectedResult ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void test_getPlaneStrainTangent()
{
  Marmot::Matrix6d C = create_C_Matrix();

  // Compute the result using the getPlaneStrainTangent function
  Matrix3d computedResult = getPlaneStrainTangent( C );

  // Expected results
  Matrix3d expectedResult;
  // clang-format off
  expectedResult << 
        1200, 400,  0,
        400,  1200, 0,
        0,    0,    400;
  // clang-format on

  throwExceptionOnFailure( checkIfEqual< double >( computedResult, expectedResult ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void test_dStrainDStrainPlaneStrain()
{
  // Compute the result using the dStrainDStrainPlaneStrain function
  Matrix< double, 6, 3 > computedResult = dStrainDStrainPlaneStrain();

  // Expected results
  Matrix< double, 6, 3 > expectedResult;
  // clang-format off
  expectedResult << 1, 0, 0,
                    0, 1, 0,
                    0, 0, 0,
                    0, 0, 1,
                    0, 0, 0,
                    0, 0, 0;
  // clang-format on

  throwExceptionOnFailure( checkIfEqual< double >( computedResult, expectedResult ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void test_compute_dStress_dDeformationGradient()
{
  // Create Tensor
  Marmot::EigenTensors::Tensor633d inputTensor;
  // Set all values to 1 initially
  inputTensor.setConstant( 1.0 );
  // Set some different values
  inputTensor( 0, 0, 0 ) = 1.0;
  inputTensor( 0, 1, 0 ) = 2.0;
  inputTensor( 0, 0, 1 ) = 3.0;
  inputTensor( 0, 1, 1 ) = 4.0;
  inputTensor( 1, 0, 0 ) = 5.0;
  inputTensor( 1, 1, 0 ) = 6.0;
  inputTensor( 1, 0, 1 ) = 7.0;
  inputTensor( 1, 1, 1 ) = 8.0;
  inputTensor( 3, 0, 0 ) = 9.0;
  inputTensor( 3, 1, 0 ) = 10.0;
  inputTensor( 3, 0, 1 ) = 11.0;
  inputTensor( 3, 1, 1 ) = 12.0;
  inputTensor( 2, 0, 0 ) = 13.0;
  inputTensor( 2, 1, 0 ) = 14.0;
  inputTensor( 2, 0, 1 ) = 15.0;
  inputTensor( 2, 1, 1 ) = 16.0;
  inputTensor( 4, 0, 0 ) = 17.0;
  inputTensor( 4, 1, 2 ) = 18.0;
  inputTensor( 5, 0, 1 ) = 19.0;
  inputTensor( 5, 1, 2 ) = 20.0;

  /*
  inputTensor:
   1  2  1  3  4  1  1  1  1
   5  6  1  7  8  1  1  1  1
  13 14  1 15 16  1  1  1  1
   9 10  1 11 12  1  1  1  1
  17  1  1  1  1  1  1 18  1
   1  1  1 19  1  1  1 20  1
  */

  // Compute the result using the compute_dStress_dDeformationGradient function
  Marmot::EigenTensors::Tensor322d computedResult = compute_dStress_dDeformationGradient( inputTensor );

  // Expected results
  using namespace Marmot::ContinuumMechanics::TensorUtility::IndexNotation;

  Marmot::EigenTensors::Tensor322d expectedResult;
  // clang-format off
  for ( int m = 0; m < 2; m ++ )
    for ( int n = 0; n < 2; n ++ )
      for ( int k = 0; k < 2; k ++ )
        for ( int l = 0; l < 2; l ++ )
          expectedResult(    toVoigt<2> (m,n), k, l ) 
          =  -  inputTensor( toVoigt<3> (m,n), 2, 2 ) 
          * 1./ inputTensor( toVoigt<3> (2,2), 2, 2 )
          *     inputTensor( toVoigt<3> (2,2), k, l );
  // clang-format on

  /*
  expectedResult:
  -13 -14 -15 -16
  -13 -14 -15 -16
  -13 -14 -15 -16
  */

  throwExceptionOnFailure( checkIfEqual< double, 3, 2, 2 >( computedResult, expectedResult ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void test_dStressPlaneStressDStress()
{
  // Compute the result using the dStressPlaneStressDStress function
  Matrix< double, 3, 6 > computedResult = dStressPlaneStressDStress();

  // Expected results
  Matrix< double, 3, 6 > expectedResult;
  // clang-format off
  expectedResult << 1, 0, 0, 0, 0, 0,
                    0, 1, 0, 0, 0, 0,
                    0, 0, 0, 1, 0, 0;
  // clang-format on

  throwExceptionOnFailure( checkIfEqual< double >( computedResult, expectedResult ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void test_getPlaneStressTangent()
{
  Marmot::Matrix6d C = create_C_Matrix();

  // Compute the result using the getPlaneStressTangent function
  Matrix3d computedResult = getPlaneStressTangent( C );

  // Expected results
  Matrix< double, 3, 6 > dStressPlaneStressDStress;
  // clang-format off
  dStressPlaneStressDStress <<
    1, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0;
  // clang-format on

  Matrix< double, 6, 3 > dStrainDStrainPlaneStress;
  // clang-format off
  dStrainDStrainPlaneStress <<
    1,        0,        0,
    0,        1,        0,
    -1.0/3.0, -1.0/3.0, 0,
    0,        0,        1,
    0,        0,        0,
    0,        0,        0;
  // clang-format on

  Matrix3d expectedResult = dStressPlaneStressDStress * C * dStrainDStrainPlaneStress;

  throwExceptionOnFailure( checkIfEqual< double >( computedResult, expectedResult ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void test_planeStressCompensationStrain()
{
  Marmot::Vector6d strain;
  strain << 0.04, 0.06, 0.0, 0.08, 0.0, 0.0;
  double nu = 0.2;

  // Compute the result using the planeStressCompensationStrain function
  Marmot::Vector6d computedResult = planeStressCompensationStrain( strain, nu );

  // Expected results
  Marmot::Vector6d expectedResult;
  expectedResult << 0, 0, -0.025, 0, 0, 0;

  throwExceptionOnFailure( checkIfEqual< double >( computedResult, expectedResult ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void test_planeStressTangentTransformationMatrix()
{
  Marmot::Matrix6d tangent = create_C_Matrix();

  // Compute the result using the planeStressTangentTransformationMatrix function
  Marmot::Matrix6d computedResult = planeStressTangentTransformationMatrix( tangent );

  // Expected results
  Marmot::Matrix6d expectedResult;
  // clang-format off
  expectedResult <<
    1,        0,        0,  0,  0,  0,
    0,        1,        0,  0,  0,  0,
    -1.0/3.0, -1.0/3.0, 0,  0,  0,  0,
    0,        0,        0,  1,  0,  0,
    0,        0,        0,  0,  1,  0,
    0,        0,        0,  0,  0,  1;
  // clang-format on

  throwExceptionOnFailure( checkIfEqual< double >( computedResult, expectedResult ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void test_dStrainDStrainPlaneStress()
{
  Marmot::Matrix6d tangent = create_C_Matrix();

  // Compute the result using the dStrainDStrainPlaneStress function
  Matrix< double, 6, 3 > computedResult = dStrainDStrainPlaneStress( tangent );

  // Expected results
  Matrix< double, 6, 3 > expectedResult;
  // clang-format off
  expectedResult <<
    1,        0,        0,
    0,        1,        0,
    -1.0/3.0, -1.0/3.0, 0,
    0,        0,        1,
    0,        0,        0,
    0,        0,        0;
  // clang-format on

  throwExceptionOnFailure( checkIfEqual< double >( computedResult, expectedResult ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

int main()
{

  auto tests = std::vector< std::function< void() > >{ test_getUniaxialStressTangent,
                                                       test_reduce3D_dStress_dDeformationGradient,
                                                       test_getPlaneStrainTangent,
                                                       test_dStrainDStrainPlaneStrain,
                                                       test_compute_dStress_dDeformationGradient,
                                                       test_dStressPlaneStressDStress,
                                                       test_getPlaneStressTangent,
                                                       test_planeStressCompensationStrain,
                                                       test_planeStressTangentTransformationMatrix,
                                                       test_dStrainDStrainPlaneStress };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
