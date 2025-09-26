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

    /** @brief Checks if a scalar is NaN
     *  @param x scalar to be checked
     *  @return true if x is NaN, false otherwise
     *
     *  Checks if value x is NaN (not a number).
     *  This is done by checking if x is unequal to itself, which is only true for NaN values.
     */
    template < typename T >
    bool isNaN( T x )
    {
      return x != x;
    }

    /** @brief Performs linear interpolation
     * @param x Point at which to interpolate
     * @param x0 First known point
     * @param x1 Second known point
     * @param y0 Value at first known point
     * @param y1 Value at second known point
     * @return Interpolated value at point x
     *
     * Performs linear interpolation to find the value at point x given two known points (x0, y0) and (x1, y1).
     * If x is outside the range [x0, x1], the function will perform linear extrapolation.
     */
    double linearInterpolation( double x, double x0, double x1, double y0, double y1 );

    /** @brief Computes the exponential of value \ref x with numerical limits check
     *  @param x Exponent to which e is raised
     *  @return exponential of x
     *
     * If x is larger than the maximum limit of double precision floating point numbers,
     * the maximum limit is returned. If \ref x is smaller than the minimum limit, the minimum limit is returned.
     */
    double exp( double x );

    /** @brief Computes the power of base \ref x with exponent \ref y with numerical limits check
     *  @param x Base
     *  @param y Exponent
     *  @return x raised to the power of y
     *
     *  If the result of \f$ x^y \f$ is larger than the maximum limit of double precision floating point numbers,
     *  the maximum limit is returned. If the result is smaller than the minimum limit, the minimum limit is returned.
     *  If \ref x is negative, the absolute value is used and a warning is printed to the console.
     */
    int getExponentPowerTen( const double x );

    /** @brief Convert angle from radiant to degree
     *  @param alpha Angle in radiant
     *  @return Angle in degree
     *
     *  Converts angle alpha given in radiant to degree.
     */
    constexpr double radToDeg( const double alpha )
    {
      return alpha * 180 / Marmot::Constants::Pi;
    }

    /** @brief Convert angle from degree to radiant
     *  @param alpha Angle in degree
     *  @return Angle in radiant
     *
     *  Converts angle alpha given in degree to radiant.
     */
    constexpr double degToRad( const double alpha )
    {
      return alpha / 180 * Marmot::Constants::Pi;
    }

    /** @brief Macaulay function applied to a scalar
     *  @param scalar Input value
     *  @return Positive part of scalar,
     *
     *  The Macaulay function, also as positive part operator, is defined as:
     *  \f[
     *  \langle x \rangle =
     *  \begin{cases}
     *  x, & x \geq 0 \\
     *  0, & x < 0
     *  \end{cases}
     *  \f]
     */
    constexpr double macauly( double scalar )
    {
      return scalar >= 0 ? scalar : 0.0;
    }

    /** @brief Heaviside step function applied to a scalar
     *  @param scalar Input value
     *  @return 1 if scalar >= 0, 0 otherwise
     *
     *  The Heaviside step function is defined as:
     *  \f[
     *  H(x) =
     *  \begin{cases}
     *  1, & x \geq 0 \\
     *  0, & x < 0
     *  \end{cases}
     *  \f]
     */
    constexpr int heaviside( double scalar )
    {
      return scalar >= 0 ? 1 : 0;
    }

    /** @brief Heaviside step function excluding zero applied to a scalar
     *  @param scalar Input value
     *  @return 1 if scalar > 0, 0 otherwise
     *
     *  The Heaviside step function excluding zero is defined as:
     *  \f[
     *  H(x) =
     *  \begin{cases}
     *  1, & x > 0 \\
     *  0, & x \leq 0
     *  \end{cases}
     *  \f]
     */
    constexpr int heavisideExclude0( double scalar )
    {
      return scalar > 0 ? 1 : 0;
    }

    /** @brief Sign function applied to a scalar
     *  @param val Input value
     *  @return 1 if val > 0, -1 if val < 0, 0 if val == 0
     *
     *  The sign function is defined as:
     *  \f[
     *  \sgn(x) =
     *  \begin{cases}
     *  1, & x > 0 \\
     *  0, & x = 0 \\
     * -1, & x < 0
     *  \end{cases}
     *  \f]
     */
    template < typename T >
    constexpr int sgn( const T& val ) noexcept
    {
      return ( T( 0 ) < val ) - ( val < T( 0 ) );
    }

    /** @brief Converts various scalar types to double precision floating point numbers
     *  @param value Input value of various numeric types
     *  @return Converted value as double
     *
     *  This function converts input values of different numeric types, including:
     *  - double
     *  - std::complex<double>
     *  - autodiff::real
     *  - autodiff::dual
     *  - Eigen::Matrix with elements of any type convertible to double
     *
     *  The function is overloaded and templated to handle these different types appropriately.
     */
    double makeReal( const double& value );
    double makeReal( const std::complex< double >& value );
    double makeReal( const autodiff::real& value );

    /** @brief Converts autodiff::dual numbers to double precision floating point numbers
     *  @tparam T Underlying type of the autodiff::dual number
     *  @param G Gradient type of the autodiff::dual number
     *  @param number Input autodiff::dual number
     *  @return Converted value as double
     *
     *  This function converts an autodiff::dual number to a double precision floating point number.
     *  It is templated to handle dual numbers with different underlying types.
     */
    template < typename T, typename G >
    double makeReal( const autodiff::detail::Dual< T, G >& number )
    {
      return double( number );
    }

    /** @brief Extracts the real part of a arbitrary scalartype-valued Matrix
     *  @tparam T scalar type
     *  @tparam Rest... parameter pack for additional matrix information, e.g, dimensions
     *  @param mat T-valued matrix
     *  @return double-valued matrix
     */
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

    /** @brief
     *
     */

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

    /** @brief Explicit Euler time integrator
     * @param yN Current value
     * @param dt Time step size
     * @param fRate Function that computes the rate of change
     * @param fRateArgs Additional arguments for the rate function
     * @return Updated value after time step
     *
     * This function computes one single time step using the explicit Euler method:
     * \f[ \boldsymbol{y}_{n+1} = \boldsymbol{y}_n + f(\boldsymbol{y}_n) * \Delta t \f]
     * where \f$ \boldsymbol{y}_n \f$ is the current value, \f$ \Delta t \f$ is the time step size,
     * and \f$ f(\boldsymbol{y}_n) \f$ is the rate of change.
     *
     */
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

      return yN + ( Iy - dt * fS ).colPivHouseholderQr.solve( Iy ) * dt * fRate( yN, fRateArgs... );
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
  } // namespace Math
} // namespace Marmot
