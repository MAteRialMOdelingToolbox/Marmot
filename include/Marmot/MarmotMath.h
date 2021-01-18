/* ---------------------------------------------------------------------
 *                                       _   
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_ 
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_ 
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 * 
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck, 
 * 2020 - today
 * 
 * festigkeitslehre@uibk.ac.at
 * 
 * Matthias Neuner matthias.neuner@uibk.ac.at
 * Magdalena Schreter magdalena.schreter@uibk.ac.at
 * 
 * This file is part of the MAteRialMOdellingToolbox (marmot).
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of marmot.
 * ---------------------------------------------------------------------
 */

#pragma once
#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotTypedefs.h"
#include <algorithm> 

namespace Marmot {
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

            Eigen::VectorXd leftX( ySize );
            Eigen::VectorXd rightX( ySize );

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

            Eigen::VectorXd leftX( nCols );
            Eigen::VectorXd rightX( nCols );

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
        
        template <int ySize, typename functionType, typename... Args>
        std::tuple<Eigen::Matrix<double, ySize, 1>, double> explicitEulerRichardsonWithErrorEstimator( Eigen::Matrix<double, ySize, 1> yN, const double dt, const double TOL, functionType fRate, Args&&... fRateArgs )
        {

            typedef Eigen::Matrix<double, ySize, 1> ySized;
            ySized fN = fRate( yN, fRateArgs... );
            ySized u  = yN + fN * dt;
            ySized v  = yN + fN * dt / 2.;
            ySized w  = v + fRate( v, fRateArgs... ) * dt / 2.;
            ySized yNew = 2. * w - u;

            // error estimator
            const double AERR = 1.0;
            const double aI = AERR/TOL;
            const double rI = 1.0;
            double scaling = 0;
            ySized ESTVec = ySized::Zero();

            for (int i =0; i < ySize; i++){
                scaling =  aI + rI * std::max( abs(yNew(i)), abs(yN(i)) );
                ESTVec(i) = abs(w(i) - u(i) ) / abs( scaling );
            }

            const double EST = ESTVec.maxCoeff();
            const double tauNew = dt * std::min(2., std::max (0.2, 0.9 * std::sqrt(TOL/EST)));

            return std::make_tuple( yNew, tauNew);
        }

       	Matrix3d DirectionCosLocalToGlobal(const Matrix3d& LocalCoordinateSystem);
	Matrix3d OrthonormalCoordinateSystem(const Vector3d& normalVector);	
        double Polyfit(const Eigen::Matrix<double,4,1>& Xdata, const Eigen::Matrix<double,4,1>& Ydata, int& angle );
	Eigen::MatrixXd DyadicProduct(const Eigen::VectorXd& Vector1, const Eigen::VectorXd& Vector2);

    } // namespace Math
} // namespace Marmot
