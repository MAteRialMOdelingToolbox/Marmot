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
 * Alexander Dummer alexander.dummer@uibk.ac.at
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
#include <complex>

namespace Marmot {
  namespace Math {
    /**
     * Check if value \ref x is a valid floating point number */
    template < typename T >
    bool isNaN( T x )
    {
      return x != x;
    }

    /**
     * Linear interpolation at location \ref x between two points (\ref x0|\ref y0) and (\ref x1|\ref y1) */
    double linearInterpolation( double x, double x0, double x1, double y0, double y1 );

    /**
     * exponential of value \ref x with numerical limits check */
    double exp( double x );

    /**
     * compute the exponent to the power of ten of an expression, e.g., 5*10^5 --> return 5 */
    int    getExponentPowerTen( const double x );

    /**
     * convert angle \ref alpha in radiant to degree */
    constexpr double radToDeg( const double alpha ) { return alpha * 180 / Marmot::Constants::Pi; }

    /**
     * convert angle \ref alpha in degree to radiant */
    constexpr double degToRad( const double alpha ) { return alpha / 180 * Marmot::Constants::Pi; }

    /**
     * Macaulay function applied to \ref scalar */
    constexpr double macauly( double scalar ) { return scalar >= 0 ? scalar : 0.0; }

    /**
     * Heaviside function applied to \ref scalar */
    constexpr int heaviside( double scalar ) { return scalar >= 0 ? 1 : 0; }

    /**
     * Extract sign of value \ref val*/
    template < typename T >
    constexpr T sgn( T val )
    {
      return val / std::abs( val );
    }

    double makeReal( const double& value );
    double makeReal( const complexDouble& value );

    /**
     * apply Macaulay function to a matrix
     * @todo: Can be replaced easily with Eigen's array() functionality ??? */
    template < int nRows, int nCols >
    Eigen::Matrix< double, nRows, nCols > macaulyMatrix( const Eigen::Matrix< double, nRows, nCols >& mat )
    {
      Eigen::Matrix< double, nRows, nCols > positivePart = Eigen::Matrix< double, nRows, nCols >::Zero();
      for ( int i = 0; i < nRows; i++ ) {
        for ( int j = 0; j < nCols; j++ ) {
          positivePart( i, j ) = macauly( mat( i, j ) );
        }
      }
      return positivePart;
    }

    /**
     * Explicit Euler integration of function \ref fRate taking arguments \ref fRateArgs and initial value \ref a yN*/
    template < typename functionType, typename yType, typename... Args >
    yType explicitEuler( yType yN, const double dt, functionType fRate, Args&&... fRateArgs )
    {
      return yN + fRate( yN, fRateArgs... ) * dt;
    }

    /**
     * Semi-implicit Euler integration of function \ref fRate taking arguments \ref fRateArgs and initial value \ref yN
     * using central difference scheme for computing
     * @todo: Replace inverse bei solving equation system?
     * @todo: Use external central difference function? */
    template < int ySize, typename functionType, typename... Args >
    Eigen::Matrix< double, ySize, 1 > semiImplicitEuler( Eigen::Matrix< double, ySize, 1 > yN,
                                                         const double                      dt,
                                                         functionType                      fRate,
                                                         Args&&... fRateArgs )
    {
      Eigen::Matrix< double, ySize, ySize > fS = Eigen::Matrix< double, ySize, ySize >::Zero();
      Eigen::Matrix< double, ySize, ySize > Iy = Eigen::Matrix< double, ySize, ySize >::Identity();

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

    /**
     * returns central numerical differentiation of a vector-valued function \ref f with respect to vector \ref X
     * --> use std::bind to create function f(x) from function with multiple arguments */
    template < int nRows, int nCols, typename functionType >
    Eigen::Matrix< double, nRows, nCols > centralDiff( functionType f, const Eigen::Matrix< double, nCols, 1 >& X )
    {
      Eigen::Matrix< double, nRows, nCols > dXdY = Eigen::Matrix< double, nRows, nCols >::Zero();

      Eigen::VectorXd leftX( nCols );
      Eigen::VectorXd rightX( nCols );

      for ( size_t i = 0; i < nCols; i++ ) {
        double volatile h = std::max( 1.0, std::abs( X( i ) ) ) * Constants::cubicRootEps();
        leftX             = X;
        leftX( i ) -= h;
        rightX = X;
        rightX( i ) += h;
        dXdY.col( i ) = 1. / ( 2. * h ) * ( f( rightX ) - f( leftX ) );
      }

      return dXdY;
    }

    /**
     *  Explicit Euler integration with Richardson extrapolation of function \ref fRate taking arguments \ref fRateArgs
     * and initial value \ref yN
     * */
    template < typename functionType, typename yType, typename... Args >
    yType explicitEulerRichardson( yType yN, const double dt, functionType fRate, Args&&... fRateArgs )
    {
      yType fN = fRate( yN, fRateArgs... );
      yType u  = yN + fN * dt;
      yType v  = yN + fN * dt / 2.;
      yType w  = v + fRate( v, fRateArgs... ) * dt / 2.;

      return 2. * w - u;
    }

    /**
     *  Explicit Euler integration with error estimation based on Richardson extrapolation of function \ref fRate taking arguments \ref fRateArgs
     * and initial value \ref yN .
     * */
    template < int ySize, typename functionType, typename... Args >
    std::tuple< Eigen::Matrix< double, ySize, 1 >, double > explicitEulerRichardsonWithErrorEstimator(
      Eigen::Matrix< double, ySize, 1 > yN,
      const double                      dt,
      const double                      TOL,
      functionType                      fRate,
      Args&&... fRateArgs )
    {

      typedef Eigen::Matrix< double, ySize, 1 > ySized;
      ySized                                    fN   = fRate( yN, fRateArgs... );
      ySized                                    u    = yN + fN * dt;
      ySized                                    v    = yN + fN * dt / 2.;
      ySized                                    w    = v + fRate( v, fRateArgs... ) * dt / 2.;
      ySized                                    yNew = 2. * w - u;

      // error estimator
      const double AERR    = 1.0;
      const double aI      = AERR / TOL;
      const double rI      = 1.0;
      double       scaling = 0;
      ySized       ESTVec  = ySized::Zero();

      for ( int i = 0; i < ySize; i++ ) {
        scaling     = aI + rI * std::max( abs( yNew( i ) ), abs( yN( i ) ) );
        ESTVec( i ) = abs( w( i ) - u( i ) ) / abs( scaling );
      }

      const double EST    = ESTVec.maxCoeff();
      const double tauNew = dt * std::min( 2., std::max( 0.2, 0.9 * std::sqrt( TOL / EST ) ) );

      return std::make_tuple( yNew, tauNew );
    }

    /**
     * Computes the directional cosines between a transformed and the global cartesian coordinate system.
     */
    Matrix3d directionCosines( const Matrix3d& transformedCoordinateSystem );

    /**
     * Computes an orthonormal coordinate system from an unit normal vector as \f$ x_1 \f$ - axis.
     */
    Matrix3d orthonormalCoordinateSystem( Vector3d& normalVector );
  } // namespace Math
} // namespace Marmot
