#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotConstants.h"
#include <cmath>
#include <math.h>

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
            if ( x >= 64 ) // overflow if arg > 709.8 (type double), leave ample margin (e.g. for
                           // squaring)
                return std::exp( 64 );
            return std::exp( x );
        }

        // return the exponent to the power of ten of an expression like 5*10^5 --> return 5
        int getExponentPowerTen( const double x )
        {
            if ( x >= 1e-16 ) // positive number
                return floor( log10( x ) );
            else if ( x <= 1e-16 ) // negative number
                return floor( log10( abs( x ) ) );
            else // number close to 0
                return 0;
        }

        double radToDeg( const double alpha ) { return alpha * 180 / Marmot::Constants::Pi; }

        double degToRad( const double alpha ) { return alpha / 180 * Marmot::Constants::Pi; }

        double macauly( double scalar ) { return scalar >= 0 ? scalar : 0.0; }

        int heaviside( double scalar ) { return scalar >= 0 ? 1 : 0; }

	Matrix3d OrthonomCoordSystem(Vector3d& NormVec)
	{
	   NormVec.normalize();
	   Matrix3d Coordsystem = Eigen::MatrixXd::Zero(3,3);
	   Coordsystem.col(0) = NormVec;
	  	   	
	   if (Coordsystem(0,0) == 0 && Coordsystem(1,0) == 0){
	   	Coordsystem(1,1) = 1.0;
	   } else{
	   	Coordsystem(0,1) = -Coordsystem(1,0);
	   	Coordsystem(1,1) =  Coordsystem(0,0);	
	   }
	   
	   Coordsystem.col(2) = Coordsystem.col(0).cross(Coordsystem.col(1));
	   Coordsystem.col(2).normalize();

	   return Coordsystem;
	}


        Matrix3d DirectionCosLocabalGlobal(const Matrix3d& LocalCoordinateSystem)
        {
	   Vector3d UnitVectorX1 = Eigen::MatrixXd::Zero(3,1);
	   UnitVectorX1(0) = 1;
	   
	   Matrix3d GlobalCoordinateSystem = OrthonormalCoordinateSystem(UnitVectorX1);
	   Matrix3d DirectionCos;
	   
	   for (int i = 0;i<=2;i++){
	   	for (int j = 0;j<=2;j++){
	   		DirectionCos(i,j) = LocalCoordinateSystem.col(i).dot(GlobalCoordinateSystem.col(j));
	   	}
	   }
	   
	   return DirectionCos;
	   
        }

	Eigen::MatrixXi DyadicProduct(const Eigen::VectorXd& Vector1, const Eigen::VectorXd& Vector2)
	{
		Eigen::MatrixXi Dyade;

		for (int i = 0;i<Vector1.rows();i++){
			for (int j = 0;j<Vector1.rows();j++){
				Dyade(i,j) = Vector1(i)*Vector2(j);
			}
		}

		return Dyade;

	}

        //int const kronecker (int i, int j)
        //{
            //if (i == j)
                //return 1;
            //else
                //return 0;
        //}

    } // namespace Math
} // namespace Marmot
