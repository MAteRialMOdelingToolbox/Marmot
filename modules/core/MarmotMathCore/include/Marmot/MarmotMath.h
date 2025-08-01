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
#include "autodiff/forward/dual.hpp"
#include "autodiff/forward/real.hpp"
#include <algorithm>
#include <autodiff/forward/dual/dual.hpp>
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

    unsigned long factorial( unsigned int n );

    /**
     * compute the exponent to the power of ten of an expression, e.g., 5*10^5 --> return 5 */
    int getExponentPowerTen( const double x );

    /**
     * convert angle \ref alpha in radiant to degree */
    constexpr double radToDeg( const double alpha )
    {
      return alpha * 180 / Marmot::Constants::Pi;
    }

    /**
     * convert angle \ref alpha in degree to radiant */
    constexpr double degToRad( const double alpha )
    {
      return alpha / 180 * Marmot::Constants::Pi;
    }

    /**
     * Macaulay function applied to \ref scalar */
    constexpr double macauly( double scalar )
    {
      return scalar >= 0 ? scalar : 0.0;
    }

    /**
     * Heaviside function applied to \ref scalar */
    constexpr int heaviside( double scalar )
    {
      return scalar >= 0 ? 1 : 0;
    }

    constexpr int heavisideExclude0( double scalar )
    {
      return scalar > 0 ? 1 : 0;
    }

    /**
     * Extract sign of value \ref val*/
    template < typename T >
    constexpr int sgn( const T& val ) noexcept
    {
      return ( T( 0 ) < val ) - ( val < T( 0 ) );
    }

    double makeReal( const double& value );
    double makeReal( const std::complex< double >& value );
    double makeReal( const autodiff::real& value );

    template < typename T, typename G >
    double makeReal( const autodiff::detail::Dual< T, G >& number )
    {
      return double( number );
    }

    template < typename T, int... Rest >
    Eigen::Matrix< double, Rest... > makeReal( const Eigen::Matrix< T, Rest... > mat )
    {
      Eigen::Matrix< double, Rest... > out;
      const size_t                     m = static_cast< size_t >( mat.rows() );
      const size_t                     n = static_cast< size_t >( mat.cols() );

      for ( size_t i = 0; i < m; i++ ) {
        for ( size_t j = 0; j < n; j++ ) {
          out( i, j ) = makeReal( mat( i, j ) );
        }
      }
      return out;
    }
    template < typename T >
    Eigen::VectorXd makeReal( Eigen::Vector< T, Eigen::Dynamic > in )
    {

      int             inSize = in.size();
      Eigen::VectorXd out( inSize );
      for ( int i = 0; i < inSize; i++ ) {
        out( i ) = double( in( i ) );
      }
      return out;
    }

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
     *  Explicit Euler integration with error estimation based on Richardson extrapolation of function \ref fRate taking
     * arguments \ref fRateArgs and initial value \ref yN .
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

    /**
     * Computes an orthonormal coordinate system from two unit normal vectors as \f$ x_1 \f$ and \f$ x_2 \f$ - axis.
     */
    Matrix3d orthonormalCoordinateSystem( const Vector3d& n1, const Vector3d& n2 );

    Matrix3d transformToLocalSystem( const Matrix3d& T, const Matrix3d& transformedCoordinateSystem );

    Matrix3d transformToGlobalSystem( const Matrix3d& T, const Matrix3d& transformedCoordinateSystem );

  } // namespace Math
} // namespace Marmot
