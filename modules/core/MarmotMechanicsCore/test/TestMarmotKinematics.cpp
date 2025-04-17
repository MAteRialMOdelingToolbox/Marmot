#include "Marmot/MarmotKinematics.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotTesting.h"

using namespace Eigen;
using namespace Marmot::Testing;
using namespace Marmot::ContinuumMechanics::Kinematics::VelocityGradient;
using namespace Marmot::ContinuumMechanics::Kinematics::Strain;
using namespace Marmot::ContinuumMechanics::Kinematics::DeformationGradient;

int transformVoigtIndices( int i, int j )
{
  int ij;
  if ( i == j ) {
    ij = i;
  }
  else if ( i * j == 0 ) {
    ij = i + j + 2;
  }
  else {
    ij = 5;
  }
  return ij;
}

void testDifferentiationsWithRegardsToVelocityGradient()
{
  // check dOmega_dVelocityGradient and dStretchingRate_dVelocityGradient
  double dwdlexpect;
  double dddlexpect;
  // loop through all indices/dimensions
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ )
      for ( int k = 0; k < 3; k++ )
        for ( int l = 0; l < 3; l++ ) {
          // calculate the expected differentiation of w wrt the velocity gradient
          dwdlexpect = 0.5 * ( ( i == k ? 1 : 0 ) * ( j == l ? 1 : 0 ) - ( k == j ? 1 : 0 ) * ( i == l ? 1 : 0 ) );
          // calculate the expected differentiation of d wrt the velocity gradient
          dddlexpect = 0.5 * ( ( i == k ? 1 : 0 ) * ( j == l ? 1 : 0 ) + ( j == k ? 1 : 0 ) * ( i == l ? 1 : 0 ) ) *
                       ( i == j ? 1 : 2 );
          // get the transformed voigt indices
          int ij = transformVoigtIndices( i, j );
          // check if expected and calculated values are the same
          throwExceptionOnFailure( checkIfEqual( dwdlexpect, dOmega_dVelocityGradient( i, j, k, l ) ),
                                   MakeString() << __PRETTY_FUNCTION__ << " Error in dOmega_dVelocityGradient." );
          throwExceptionOnFailure( checkIfEqual( dddlexpect, dStretchingRate_dVelocityGradient( ij, k, l ) ),
                                   MakeString()
                                     << __PRETTY_FUNCTION__ << " Error in dStretchingRate_dVelocityGradient." );
        }
}

void testDeformationGradientAnddEdFCalculation()
{
  // Test for 0 deformation
  Eigen::Matrix3d  F       = Eigen::Matrix3d::Identity(); // assign identity test value
  Marmot::Vector6d calcE   = GreenLagrange( F );          // calculate the value for 0 strain by function
  Marmot::Vector6d expectE = Marmot::Vector6d::Zero();    // assign the expected value for 0 strain
  // check if expected and calculated values are the same
  throwExceptionOnFailure( checkIfEqual< double >( calcE, expectE ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " Error for Green-Lagrange strain with 0 deformation." );
  // calculate the value for the differentiation of strain wrt deformation gradient by function
  Marmot::EigenTensors::Tensor633d dEdFcalc = dGreenLagrangedDeformationGradient( F );
  double                           dEdFexpect;
  // loop through all indices/dimensions
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ )
      for ( int k = 0; k < 3; k++ )
        for ( int l = 0; l < 3; l++ ) {
          // get the transformed voigt indices
          int ij = transformVoigtIndices( i, j );
          // calculate the expected value of strain wrt deformation gradient
          dEdFexpect = 0.5 * ( ( i == l ? 1 : 0 ) * F( k, j ) + ( j == l ? 1 : 0 ) * F( k, i ) ) * ( i == j ? 1 : 2 );
          // check if expected and calculated values are the same
          throwExceptionOnFailure( checkIfEqual( dEdFexpect, dEdFcalc( ij, k, l ) ),
                                   MakeString() << __PRETTY_FUNCTION__ << " Error in dEdF for 0 deformation." );
        }

  // Test for high deformation
  Eigen::Matrix3d Fh;
  Fh << 1.2, 0.1, 0.2, 0.4, 1.4, 0.3, 0.5, 0.6, 1.6;    // assign high test value
  Marmot::Vector6d calcEh = GreenLagrange( Fh );        // calculate the value for high strain by function
  Marmot::Vector6d expectEh;
  expectEh << 0.425, 0.665, 0.845, 0.980, 1.160, 1.400; // assign the expected value for high strain
  // check if expected and calculated values are the same
  throwExceptionOnFailure( checkIfEqual< double >( calcEh, expectEh, 1e-15 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " Error for Green-Lagrange strain with high deformation." );
  // calculate the value for the differentiation of strain wrt deformation gradient by function
  Marmot::EigenTensors::Tensor633d dEdFcalch = dGreenLagrangedDeformationGradient( Fh );
  // loop through all indices/dimensions
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ )
      for ( int k = 0; k < 3; k++ )
        for ( int l = 0; l < 3; l++ ) {
          // get the transformed voigt indices
          int ij = transformVoigtIndices( i, j );
          // calculate the expected value of strain wrt deformation gradient
          dEdFexpect = 0.5 * ( ( i == l ? 1 : 0 ) * Fh( k, j ) + ( j == l ? 1 : 0 ) * Fh( k, i ) ) * ( i == j ? 1 : 2 );
          // check if expected and calculated values are the same
          throwExceptionOnFailure( checkIfEqual( dEdFexpect, dEdFcalch( ij, k, l ) ),
                                   MakeString() << __PRETTY_FUNCTION__ << " Error in dEdF with high deformation." );
        }
}

void testDimensionalTransformationTo3DOfTheDeformationGradient()
{
  // Test 3D making 1x1
  const Eigen::Matrix< double, 1, 1 >                     t1D{ 5 };      // assign 1x1 deformation gradient
  const Eigen::Ref< const Eigen::Matrix< double, 1, 1 > > t1Dref( t1D ); // make it a Ref object
  Eigen::Matrix3d calcF1D = make3D( t1Dref ); // get value for 3D deformation gradient by function
  Eigen::Matrix3d expectF1D;
  expectF1D << 5, 0, 0, 0, 1, 0, 0, 0, 1;     // assign expected value of 3D deformation gradient
  // check if expected and calculated values are the same
  throwExceptionOnFailure( checkIfEqual< double >( calcF1D, expectF1D, 1e-15 ),
                           MakeString() << __PRETTY_FUNCTION__ << " Making 1D F 3D failed." );

  // Test 3D making 2x2
  const Eigen::Matrix< double, 2, 2 >                     t2D{ { 5, 2 }, { 3, 4 } }; // assign 2x2 deformation gradient
  const Eigen::Ref< const Eigen::Matrix< double, 2, 2 > > t2Dref( t2D );             // make it a Ref object
  Eigen::Matrix3d calcF2D = make3D( t2Dref ); // get value for 3D deformation gradient by function
  Eigen::Matrix3d expectF2D;
  expectF2D << 5, 2, 0, 3, 4, 0, 0, 0, 1;     // assign expected value of 3D deformation gradient
  // check if expected and calculated values are the same
  throwExceptionOnFailure( checkIfEqual< double >( calcF2D, expectF2D, 1e-15 ),
                           MakeString() << __PRETTY_FUNCTION__ << " Making 2D F 3D failed." );

  // Test 3D making 3x3
  const Eigen::Matrix< double, 3, 3 > t3D{ { 5, 2, 3 }, { 4, 10, 6 }, { 7, 8, 9 } }; // assign 3x3 deformation gradient
  const Eigen::Ref< const Eigen::Matrix< double, 3, 3 > > t3Dref( t3D );             // make it a Ref object
  Eigen::Matrix3d calcF3D = make3D( t3Dref ); // get value for 3D deformation gradient by function
  Eigen::Matrix3d expectF3D;
  expectF3D << 5, 2, 3, 4, 10, 6, 7, 8, 9;    // assign expected value of 3D deformation gradient
  // check if expected and calculated values are the same
  throwExceptionOnFailure( checkIfEqual< double >( calcF3D, expectF3D, 1e-15 ),
                           MakeString() << __PRETTY_FUNCTION__ << " Making 3D F 3D failed." );
}

int main()
{

  auto tests = std::vector< std::function< void() > >{ testDifferentiationsWithRegardsToVelocityGradient,
                                                       testDeformationGradientAnddEdFCalculation,
                                                       testDimensionalTransformationTo3DOfTheDeformationGradient };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
