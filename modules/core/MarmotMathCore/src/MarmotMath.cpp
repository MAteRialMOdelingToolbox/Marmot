#include "Marmot/MarmotMath.h"

namespace Marmot {
  namespace Math {
    // return linear interpolation of polynom y at given coordinates (x0, y0) and (x1, y1) at
    // point x
    double linearInterpolation( double x, double x0, double x1, double y0, double y1 )
    {
      return y0 + ( x - x0 ) * ( y1 - y0 ) / ( x1 - x0 );
    }

    // bounded version of std::exp
    double exp( double x )
    {
      if ( x <= -64 ) // underflow if arg < -708.4 (type double)
        return 0.0;
      if ( x >= 64 )  // overflow if arg > 709.8 (type double), leave ample margin (e.g. for
                      // squaring)
        return std::exp( 64 );
      return std::exp( x );
    }

    unsigned long factorial( unsigned int n )
    {
      if ( n == 0 || n == 1 )
        return 1;
      return n * factorial( n - 1 );
    }

    double makeReal( const complexDouble& value )
    {
      return value.real();
    }

    double makeReal( const autodiff::real& value )
    {
      return double( value );
    }
    double makeReal( const autodiff::dual& value )
    {
      return double( value );
    }

    double makeReal( const double& value )
    {
      return value;
    }

    // return the exponent to the power of ten of an expression like 5*10^5 --> return 5
    int getExponentPowerTen( const double x )
    {
      if ( x >= 1e-16 )      // positive number
        return floor( log10( x ) );
      else if ( x <= 1e-16 ) // negative number
        return floor( log10( abs( x ) ) );
      else                   // number close to 0
        return 0;
    }

    Matrix3d orthonormalCoordinateSystem( Vector3d& normalVector )
    {
      normalVector.normalize();
      Matrix3d coordinateSystem = Eigen::MatrixXd::Zero( 3, 3 );
      coordinateSystem.col( 0 ) = normalVector;

      if ( coordinateSystem( 0, 0 ) == 0 && coordinateSystem( 1, 0 ) == 0 ) {
        coordinateSystem( 1, 1 ) = 1.0;
      }
      else {
        coordinateSystem( 0, 1 ) = -coordinateSystem( 1, 0 );
        coordinateSystem( 1, 1 ) = coordinateSystem( 0, 0 );
      }
      coordinateSystem.col( 1 ).normalize();

      coordinateSystem.col( 2 ) = coordinateSystem.col( 0 ).cross( coordinateSystem.col( 1 ) );
      coordinateSystem.col( 2 ).normalize();

      return coordinateSystem;
    }

    Matrix3d orthonormalCoordinateSystem( const Vector3d& n1, const Vector3d& n2 )
    {

      // check if n1 and n2 are orthogonal

      if ( std::abs( n1.dot( n2 ) ) > 1e-15 )
        throw std::invalid_argument( "n1 and n2 not orthogonal" );

      Matrix3d coordinateSystem = Eigen::MatrixXd::Zero( 3, 3 );
      coordinateSystem.col( 0 ) = n1;
      coordinateSystem.col( 1 ) = n2;

      // make sure n1 and n2 are normalized
      coordinateSystem.col( 0 ).normalize();
      coordinateSystem.col( 1 ).normalize();

      coordinateSystem.col( 2 ) = coordinateSystem.col( 0 ).cross( coordinateSystem.col( 1 ) );
      coordinateSystem.col( 2 ).normalize();

      return coordinateSystem;
    }

    Matrix3d directionCosines( const Matrix3d& transformedCoordinateSystem )
    {
      return transformedCoordinateSystem.transpose();
    }

    Matrix3d transformToLocalSystem( const Matrix3d& T, const Matrix3d& transformedCoordinateSystem )
    {
      const Matrix3d N = transformedCoordinateSystem;
      return N.transpose() * T * N;
    }

    Matrix3d transformToGlobalSystem( const Matrix3d& T, const Matrix3d& transformedCoordinateSystem )
    {
      const Matrix3d N = transformedCoordinateSystem;
      return N * T * N.transpose();
    }

  } // namespace Math
} // namespace Marmot
