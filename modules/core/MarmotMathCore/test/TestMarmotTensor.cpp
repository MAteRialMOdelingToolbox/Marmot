#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotTesting.h"

using namespace Marmot::ContinuumMechanics::CommonTensors;
using namespace Marmot::ContinuumMechanics::TensorUtility;
using namespace Marmot::Testing;
using namespace Eigen;

auto testInitialize_IFourthOrder()
{
  // variable to use for check
  double ikjl;
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ )
      for ( int k = 0; k < 3; k++ )
        for ( int l = 0; l < 3; l++ ) {
          // assign value to variable
          ikjl = ( i == k ? 1 : 0 ) * ( j == l ? 1 : 0 );
          // check if they are the same
          throwExceptionOnFailure( checkIfEqual( ikjl, IFourthOrder( i, j, k, l ) ),
                                   MakeString() << __PRETTY_FUNCTION__ << " Error in initializing IFourthOrder." );
        }
}

auto testInitialize_IFourthOrderTranspose()
{
  // variable to use for check
  double iljk;
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ )
      for ( int k = 0; k < 3; k++ )
        for ( int l = 0; l < 3; l++ ) {
          // assign value to variable
          iljk = ( i == l ? 1 : 0 ) * ( j == k ? 1 : 0 );
          // check if they are the same
          throwExceptionOnFailure( checkIfEqual( iljk, IFourthOrderTranspose( i, j, k, l ) ),
                                   MakeString()
                                     << __PRETTY_FUNCTION__ << " Error in initializing IFourthOrderTranspose." );
        }
}

auto testInitialize_I2xI2()
{
  // variable to use for check
  double ijkl;
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ )
      for ( int k = 0; k < 3; k++ )
        for ( int l = 0; l < 3; l++ ) {
          // assign value to variable
          ijkl = ( i == j ? 1 : 0 ) * ( k == l ? 1 : 0 );
          // check if they are the same
          throwExceptionOnFailure( checkIfEqual( ijkl, I2xI2( i, j, k, l ) ),
                                   MakeString() << __PRETTY_FUNCTION__ << " Error in initializing I2xI2." );
        }
}

auto testInitialize_Isym()
{
  // variable to use for check
  double ijklSym;
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ )
      for ( int k = 0; k < 3; k++ )
        for ( int l = 0; l < 3; l++ ) {
          // assign value to variable
          ijklSym = 0.5 * ( ( i == k ? 1 : 0 ) * ( j == l ? 1 : 0 ) + ( i == l ? 1 : 0 ) * ( j == k ? 1 : 0 ) );
          // check if they are the same
          throwExceptionOnFailure( checkIfEqual( ijklSym, Isym( i, j, k, l ) ),
                                   MakeString() << __PRETTY_FUNCTION__ << " Error in initializing Isym." );
        }
}

auto testInitialize_Iskew()
{
  // variable to use for check
  double ijklSkew;
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ )
      for ( int k = 0; k < 3; k++ )
        for ( int l = 0; l < 3; l++ ) {
          // assign value to variable
          ijklSkew = 0.5 * ( ( i == k ? 1 : 0 ) * ( j == l ? 1 : 0 ) - ( i == l ? 1 : 0 ) * ( j == k ? 1 : 0 ) );
          // check if they are the same
          throwExceptionOnFailure( checkIfEqual( ijklSkew, Iskew( i, j, k, l ) ),
                                   MakeString() << __PRETTY_FUNCTION__ << " Error in initializing Iskew." );
        }
}

auto testInitialize_dDeviatoricStress_dStress()
{
  // variable to use for check
  double Idev;
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ )
      for ( int k = 0; k < 3; k++ )
        for ( int l = 0; l < 3; l++ ) {
          // assign value to variable
          Idev = ( i == k ? 1 : 0 ) * ( j == l ? 1 : 0 ) - 1 / 3. * ( i == j ? 1 : 0 ) * ( k == l ? 1 : 0 );
          // check if they are the same
          throwExceptionOnFailure( checkIfEqual( Idev, dDeviatoricStress_dStress( i, j, k, l ) ),
                                   MakeString()
                                     << __PRETTY_FUNCTION__ << " Error in initializing dDeviatoricStress_dStress." );
        }
}

auto testInitialize_LeviCivita3D()
{
  // variable to use for check
  double e;
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ )
      for ( int k = 0; k < 3; k++ ) {
        // assign value to variable
        e = ( ( i != j && j != k && k != i )
                ? ( ( ( i > j && ( j > k || k == 2 ) ) || ( i == 0 && k == 1 ) ) ? -1.0 : 1.0 )
                : 0.0 );
        // check if they are the same
        throwExceptionOnFailure( checkIfEqual( e, LeviCivita3D( i, j, k ) ),
                                 MakeString() << __PRETTY_FUNCTION__ << " Error in initializing LeviCivita3D." );
      }
}

auto testInitialize_LeviCivita2D()
{
  // variable to use for check
  double e;
  for ( int j = 0; j < 2; j++ )
    for ( int k = 0; k < 2; k++ ) {
      // assign value to variable
      e = ( ( j != k ) ? ( ( j > k ) ? -1.0 : 1.0 ) : 0.0 );
      // check if they are the same
      throwExceptionOnFailure( checkIfEqual( e, LeviCivita2D( 0, j, k ) ),
                               MakeString() << __PRETTY_FUNCTION__ << " Error in initializing LeviCivita2D." );
    }
}

auto testDyadicProduct()
{
  Eigen::Matrix3d dyadicExpect;
  dyadicExpect << 4, 5, 6, 8, 10, 12, 12, 15, 18;
  Eigen::Vector3d a;
  a << 1, 2, 3;
  Eigen::Vector3d b;
  b << 4, 5, 6;
  throwExceptionOnFailure( checkIfEqual< double >( dyadicExpect, dyadicProduct( a, b ), 1e-15 ),
                           MakeString() << __PRETTY_FUNCTION__ << " Error in in dyadic product." );
}

auto testTensorMatrixConversions()
{
  Fastor::Tensor< double, 3, 3 >       tensor2D;
  Fastor::Tensor< double, 3, 3, 3 >    tensor3D;
  Fastor::Tensor< double, 3, 3, 3, 3 > tensor4D;

  for ( int i = 0; i < 3; ++i )
    for ( int j = 0; j < 3; ++j ) {
      tensor2D( i, j ) = 10 * i + j;

      for ( int k = 0; k < 3; ++k ) {
        tensor3D( i, j, k ) = 100 * i + 10 * j + k;

        for ( int l = 0; l < 3; ++l )
          tensor4D( i, j, k, l ) = 1000 * i + 100 * j + 10 * k + l;
      }
    }

  const auto matrix2D_3x3 = convert2ndOrderTensorToMatrix_3x3( tensor2D );
  const auto matrix3D_9x3 = convert3rdOrderTensorToMatrix_9x3( tensor3D );
  const auto matrix3D_3x9 = convert3rdOrderTensorToMatrix_3x9( tensor3D );
  const auto matrix4D_9x9 = convert4thOrderTensorToMatrix_9x9( tensor4D );

  for ( int i = 0; i < 3; ++i )
    for ( int j = 0; j < 3; ++j ) {
      throwExceptionOnFailure( checkIfEqual( tensor2D( i, j ), matrix2D_3x3( i, j ) ),
                               MakeString() << __PRETTY_FUNCTION__ << " Error in 2nd-order tensor conversion." );

      for ( int k = 0; k < 3; ++k ) {
        const int row = 3 * i + j;
        const int col = 3 * j + k;

        throwExceptionOnFailure( checkIfEqual( tensor3D( i, j, k ), matrix3D_9x3( row, k ) ),
                                 MakeString() << __PRETTY_FUNCTION__ << " Error in 3rd-order (9x3) conversion." );

        throwExceptionOnFailure( checkIfEqual( tensor3D( i, j, k ), matrix3D_3x9( i, col ) ),
                                 MakeString() << __PRETTY_FUNCTION__ << " Error in 3rd-order (3x9) conversion." );

        for ( int l = 0; l < 3; ++l ) {
          const int row4 = 3 * i + j;
          const int col4 = 3 * k + l;

          throwExceptionOnFailure( checkIfEqual( tensor4D( i, j, k, l ), matrix4D_9x9( row4, col4 ) ),
                                   MakeString() << __PRETTY_FUNCTION__ << " Error in 4th-order (9x9) conversion." );
        }
      }
    }
}

int main()
{

  auto tests = std::vector< std::function< void() > >{ testInitialize_IFourthOrder,
                                                       testInitialize_IFourthOrderTranspose,
                                                       testInitialize_I2xI2,
                                                       testInitialize_Isym,
                                                       testInitialize_Iskew,
                                                       testInitialize_dDeviatoricStress_dStress,
                                                       testInitialize_LeviCivita3D,
                                                       testInitialize_LeviCivita2D,
                                                       testDyadicProduct,
                                                       testTensorMatrixConversions };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
