#include "Marmot/MarmotTensor.h"

using namespace Marmot::ContinuumMechanics::TensorUtility;

namespace Marmot {
  namespace ContinuumMechanics::CommonTensors {

    EigenTensors::Tensor3333d Initialize_IFourthOrder()
    {
      EigenTensors::Tensor3333d I;

      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          for ( int k = 0; k < 3; k++ )
            for ( int l = 0; l < 3; l++ ) {
              I( i, j, k, l ) = d( i, k ) * d( j, l );
            }
      return I;
    }

    EigenTensors::Tensor3333d Initialize_IFourthOrderTranspose()
    {
      EigenTensors::Tensor3333d IT;

      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          for ( int k = 0; k < 3; k++ )
            for ( int l = 0; l < 3; l++ ) {
              IT( i, j, k, l ) = d( i, l ) * d( j, k );
            }
      return IT;
    }

    EigenTensors::Tensor3333d Initialize_I2xI2()
    {
      EigenTensors::Tensor3333d I2xI2;

      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          for ( int k = 0; k < 3; k++ )
            for ( int l = 0; l < 3; l++ ) {
              I2xI2( i, j, k, l ) = d( i, j ) * d( k, l );
            }
      return I2xI2;
    }

    EigenTensors::Tensor33d Initialize_I2()
    {
      EigenTensors::Tensor33d I2;

      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          I2( i, j ) = d( i, j );

      return I2;
    }

    EigenTensors::Tensor3333d Initialize_Isym()
    {
      EigenTensors::Tensor3333d Isym;

      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          for ( int k = 0; k < 3; k++ )
            for ( int l = 0; l < 3; l++ ) {
              Isym( i, j, k, l ) = 0.5 * ( d( i, k ) * d( j, l ) + d( i, l ) * d( j, k ) );
            }
      return Isym;
    }

    EigenTensors::Tensor3333d Initialize_Iskew()
    {
      EigenTensors::Tensor3333d Iskew;

      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          for ( int k = 0; k < 3; k++ )
            for ( int l = 0; l < 3; l++ ) {
              Iskew( i, j, k, l ) = 0.5 * ( d( i, k ) * d( j, l ) - d( i, l ) * d( j, k ) );
            }
      return Iskew;
    }

    EigenTensors::Tensor3333d Initialize_dDeviatoricStress_dStress()
    {
      EigenTensors::Tensor3333d dsdsigma;
      dsdsigma.setZero();
      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          for ( int k = 0; k < 3; k++ )
            for ( int l = 0; l < 3; l++ ) {
              dsdsigma( i, j, k, l ) = d( i, k ) * d( j, l ) - 1. / 3 * d( i, j ) * d( k, l );
            }
      return dsdsigma;
    }

    EigenTensors::Tensor333d Initialize_LeviCivita3D()
    {
      EigenTensors::Tensor333d e;
      e.setConstant( 0.0 );

      e( 0, 1, 2 ) = 1.0;
      e( 1, 2, 0 ) = 1.0;
      e( 2, 0, 1 ) = 1.0;

      e( 2, 1, 0 ) = -1.0;
      e( 0, 2, 1 ) = -1.0;
      e( 1, 0, 2 ) = -1.0;

      return e;
    }

    EigenTensors::Tensor122d Initialize_LeviCivita2D()
    {
      EigenTensors::Tensor122d e;
      e.setConstant( 0.0 );

      e( 0, 0, 1 ) = 1.0;
      e( 0, 1, 0 ) = -1.0;

      return e;
    }

  } // namespace ContinuumMechanics::CommonTensors

  namespace ContinuumMechanics::TensorUtility {
    Eigen::Matrix3d dyadicProduct( const Eigen::Vector3d& vector1, const Eigen::Vector3d& vector2 )
    {
      Eigen::Matrix3d dyade;

      for ( int i = 0; i < vector1.rows(); i++ )
        for ( int j = 0; j < vector1.rows(); j++ )
          dyade( i, j ) = vector1( i ) * vector2( j );

      return dyade;
    }
  } // namespace ContinuumMechanics::TensorUtility
} // namespace Marmot
