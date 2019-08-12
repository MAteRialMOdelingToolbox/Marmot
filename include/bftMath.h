#pragma once
#include "bftConstants.h"
#include "bftTypedefs.h"
#include <string>

using namespace Eigen;

namespace bft {
    namespace Math {
        double linearInterpolation( double x, double x0, double x1, double y0, double y1 );
        double exp( double x );
        int    getExponentPowerTen( const double x );
        double radToDeg( const double alpha );
        double degToRad( const double alpha );
        double macauly( double scalar );
        int    heaviside( double scalar );

        template <int nRows, int nCols>
        Eigen::Matrix<double, nRows, nCols> macaulyMatrix( const Eigen::Matrix<double, nRows, nCols>& mat )
        {
            Eigen::Matrix<double, nRows, nCols> positivePart = Eigen::Matrix<double, nRows, nCols>::Zero();
            for ( int i = 0; i < nRows; i++ ) {
                for ( int j = 0; j < nCols; j++ ) {
                    positivePart( i, j ) = macauly( mat( i, j ) );
                }
            }
            return positivePart;
        }

        template <typename functionType, typename yType, typename... Args>
        yType explicitEuler( yType yN, const double dt, functionType fRate, Args&&... fRateArgs )
        {
            return yN + fRate( yN, fRateArgs... ) * dt;
        }

        template <int ySize, typename functionType, typename... Args>
        Eigen::Matrix<double, ySize, 1> semiImplicitEuler( Eigen::Matrix<double, ySize, 1> yN, const double dt, functionType fRate, Args&&... fRateArgs )
        {
            Eigen::Matrix<double, ySize, ySize> fS = Eigen::Matrix<double, ySize, ySize>::Zero();
            Eigen::Matrix<double, ySize, ySize> Iy = Eigen::Matrix<double, ySize, ySize>::Identity();

            VectorXd leftX( ySize );
            VectorXd rightX( ySize );

            for ( size_t i = 0; i < ySize; i++ ) {
                double volatile h = std::max( 1.0, std::abs( yN( i ) ) ) * Constants::cubicRootEps();
                leftX             = yN;
                leftX( i ) -= h;
                rightX = yN;
                rightX( i ) += h;
                fS.col( i ) = 1. / ( 2. * h ) * ( fRate( rightX, fRateArgs... ) - fRate( leftX, fRateArgs... ) );
            }

            return yN + ( Iy - dt * fS ).inverse() * dt * fRate( yN, fRateArgs... );
        }

        // returns central numerical differentiation of a vector-valued function f(x) with respect to vector x
        // --> use std::bind to create function f(x) from function with multiple arguments
        template <int nRows, int nCols, typename functionType>
        Eigen::Matrix<double, nRows, nCols> centralDiff( functionType f, const Eigen::Matrix<double, nCols, 1>& X )
        {
            Eigen::Matrix<double, nRows, nCols> dXdY = Eigen::Matrix<double, nRows, nCols>::Zero();

            VectorXd leftX( nCols );
            VectorXd rightX( nCols );

            for ( size_t i = 0; i < nCols; i++ ) {
                double volatile h = std::max( 1.0, std::abs( X( i ) ) ) * Constants::cubicRootEps();
                leftX             = X;
                leftX( i )       -= h;
                rightX            = X;
                rightX( i )      += h;
                dXdY.col( i ) = 1. / ( 2. * h ) * ( f( rightX ) - f( leftX ) );
            }

            return dXdY;
        }

        template <typename functionType, typename yType, typename... Args>
        yType explicitEulerRichardson( yType yN, const double dt, functionType fRate, Args&&... fRateArgs )
        {
            yType fN = fRate( yN, fRateArgs... );
            yType u  = yN + fN * dt;
            yType v  = yN + fN * dt / 2.;
            yType w  = v + fRate( v, fRateArgs... ) * dt / 2.;

            return 2. * w - u;
        }

    } // namespace Math
} // namespace bft
