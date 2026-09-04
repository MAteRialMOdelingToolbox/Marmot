#include "Marmot/MarmotTensor.h"

using namespace Marmot::TensorUtility;

namespace Marmot {
  namespace TensorUtility::EigenTensors {

    Tensor3333d Initialize_IFourthOrder()
    {
      Tensor3333d I;

      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          for ( int k = 0; k < 3; k++ )
            for ( int l = 0; l < 3; l++ ) {
              I( i, j, k, l ) = d( i, k ) * d( j, l );
            }
      return I;
    }

    Tensor3333d Initialize_IFourthOrderTranspose()
    {
      Tensor3333d IT;

      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          for ( int k = 0; k < 3; k++ )
            for ( int l = 0; l < 3; l++ ) {
              IT( i, j, k, l ) = d( i, l ) * d( j, k );
            }
      return IT;
    }

    Tensor3333d Initialize_I2xI2()
    {
      Tensor3333d I2xI2;

      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          for ( int k = 0; k < 3; k++ )
            for ( int l = 0; l < 3; l++ ) {
              I2xI2( i, j, k, l ) = d( i, j ) * d( k, l );
            }
      return I2xI2;
    }

    Tensor33d Initialize_I2()
    {
      Tensor33d I2;

      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          I2( i, j ) = d( i, j );

      return I2;
    }

    Tensor3333d Initialize_Isym()
    {
      Tensor3333d Isym;

      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          for ( int k = 0; k < 3; k++ )
            for ( int l = 0; l < 3; l++ ) {
              Isym( i, j, k, l ) = 0.5 * ( d( i, k ) * d( j, l ) + d( i, l ) * d( j, k ) );
            }
      return Isym;
    }

    Tensor3333d Initialize_Iskew()
    {
      Tensor3333d Iskew;

      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          for ( int k = 0; k < 3; k++ )
            for ( int l = 0; l < 3; l++ ) {
              Iskew( i, j, k, l ) = 0.5 * ( d( i, k ) * d( j, l ) - d( i, l ) * d( j, k ) );
            }
      return Iskew;
    }

    Tensor3333d Initialize_dDeviatoricStress_dStress()
    {
      Tensor3333d dsdsigma;
      dsdsigma.setZero();
      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          for ( int k = 0; k < 3; k++ )
            for ( int l = 0; l < 3; l++ ) {
              dsdsigma( i, j, k, l ) = d( i, k ) * d( j, l ) - 1. / 3 * d( i, j ) * d( k, l );
            }
      return dsdsigma;
    }

    Tensor333d Initialize_LeviCivita3D()
    {
      Tensor333d e;
      e.setConstant( 0.0 );

      e( 0, 1, 2 ) = 1.0;
      e( 1, 2, 0 ) = 1.0;
      e( 2, 0, 1 ) = 1.0;

      e( 2, 1, 0 ) = -1.0;
      e( 0, 2, 1 ) = -1.0;
      e( 1, 0, 2 ) = -1.0;

      return e;
    }

    Tensor122d Initialize_LeviCivita2D()
    {
      Tensor122d e;
      e.setConstant( 0.0 );

      e( 0, 0, 1 ) = 1.0;
      e( 0, 1, 0 ) = -1.0;

      return e;
    }

  } // namespace TensorUtility::EigenTensors

  namespace TensorUtility {
    Eigen::Matrix3d dyadicProduct( const Eigen::Vector3d& vector1, const Eigen::Vector3d& vector2 )
    {
      Eigen::Matrix3d dyade;

      for ( int i = 0; i < vector1.rows(); i++ )
        for ( int j = 0; j < vector1.rows(); j++ )
          dyade( i, j ) = vector1( i ) * vector2( j );

      return dyade;
    }
  } // namespace TensorUtility
} // namespace Marmot
