#pragma once
#include <string>
#include "bftTypedefs.h"

namespace bft {
    namespace Math {
        double linearInterpolation( double x, double x0, double x1, double y0, double y1 );
        double exp( double x );
        int    getExponentPowerTen( const double x );
        double radToDeg( const double alpha );
        double degToRad( const double alpha );
        double macauly( double scalar );
        int    heaviside( double scalar );
        //int const kronecker ( int i, int j);
        
        
        template <int nRows, int nCols>
        Eigen::Matrix<double, nRows, nCols> macaulyMatrix(const Eigen::Matrix<double, nRows, nCols>& mat )
        {
            Eigen::Matrix<double, nRows, nCols> positivePart = Eigen::Matrix<double, nRows, nCols>::Zero();
            for (int i=0; i<nRows; i++){
                for (int j=0; j<nCols; j++){
                    positivePart(i,j) = macauly( mat(i,j));
                }
            }
            return positivePart;
        }

    } // namespace Math
} // namespace bft
