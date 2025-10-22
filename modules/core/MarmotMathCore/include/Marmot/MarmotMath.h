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

    /**
     * @brief Runtime factorial
     * @param n An unsigned integer
     * @return The factorial of n
     */
    unsigned long factorial( unsigned int n );

    /** @brief Extracts the exponent to the power of ten from a floating point number
     *  @param x Input floating point number
     *  @return Exponent to the power of ten
     *
     *  This function extracts the exponent to the power of ten from a given floating point number.
     *  For example, for an input of 3e5, the function will return 5.
     *  If the input number is very close to zero (between -1e-16 and 1e-16), the function returns 0.
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
     *  @return Positive part of scalar
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
     *  \text{sign}(x) =
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

    /** @brief Extracts the real part of a arbitrary scalartype-valued Vector
     *  @tparam T scalar type
     *  @param in T-valued vector
     *  @return double-valued vector
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
     * @brief Apply Macaulay function to a matrix
     * @tparam nRows Number of rows in the matrix
     * @tparam nCols Number of columns in the matrix
     * @param mat Input matrix
     * @return Matrix with Macaulay function applied element-wise
     *
     * Applies the Macaulay function element-wise to the input matrix.
     *
     * @todo: Can be replaced easily with Eigen's array() functionality ???
     */
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
     * @brief Implicit Euler time integrator for a linear vector-valued rate equation
     * @param yN Current value
     * @param dt Time step size
     * @param fRate Function that computes the rate of change
     * @param fRateArgs Additional arguments for the rate function
     * @return Updated value after time step
     *
     * This function computes one single time step using the implicit Euler method for a linear vector-valued rate
     * equation: \f[ \boldsymbol{y}_{n+1} = \boldsymbol{y}_n + f(\boldsymbol{y}_{n+1}) \Delta t \f] where \f$
     * \boldsymbol{y}_n \f$ is the current value, \f$ \Delta t \f$ is the time step size, and \f$ f(\boldsymbol{y}) \f$
     * is the linear rate of change evaluated at the next time step.
     *
     *
     * @todo: Is this really semi-implicit? It looks like a fully implicit Euler step for linear systems.
     * @todo: Use external central difference function?
     *
     * */
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

      return yN + ( Iy - dt * fS ).colPivHouseholderQr().solve( Iy ) * dt * fRate( yN, fRateArgs... );
    }

    /**
     * @brief Explicit Euler integration with error estimation based on Richardson extrapolation
     * @tparam functionType Type of the rate function
     * @tparam yType Type of the current value
     * @tparam Args Additional argument types for the rate function
     * @param yN Current value
     * @param dt Time step size
     * @param fRate Function that computes the rate of change
     * @param fRateArgs Additional arguments for the rate function
     * @return Updated value after time step
     *
     * This function computes one single time step using the explicit Euler method with Richardson extrapolation for
     * error estimation: \f[ \boldsymbol{y}_{n+1} = 2 \left( \boldsymbol{y}_n + f\left(\boldsymbol{y}_n
     * + \frac{f(\boldsymbol{y}_n) \Delta t}{2}\right) \frac{\Delta t}{2} \right) - \left( \boldsymbol{y}_n +
     * f(\boldsymbol{y}_n) \Delta t \right) \f] where \f$ \boldsymbol{y}_n \f$ is the current value, \f$ \Delta t \f$ is
     * the time step size, and \f$ f(\boldsymbol{y}) \f$ is the rate of change.
     *
     */
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
     * @brief Explicit Euler integration based on Richardson extrapolation with error estimation and time step
     * estimation
     * @tparam ySize Size of the state vector
     * @tparam yType Type of the current value
     * @tparam Args Additional argument types for the rate function
     * @param functionType Type of the rate function
     * @param Args Additional argument types for the rate function
     * @param yN Current value
     * @param dt Current time step size
     * @param TOL Desired tolerance for the error estimation
     * @param fRate Function that computes the rate of change
     * @param fRateArgs Additional arguments for the rate function
     * @return Tuple containing the updated value after time step and the new time step size
     *
     * This function computes one single time step using the explicit Euler method with Richardson extrapolation for
     * error estimation and adaptive time stepping: \f[ \boldsymbol{y}_{n+1} = 2 \left( \boldsymbol{y}_n +
     * f\left(\boldsymbol{y}_n
     * + \frac{f(\boldsymbol{y}_n) \Delta t}{2}\right) \frac{\Delta t}{2} \right) - \left( \boldsymbol{y}_n +
     * f(\boldsymbol{y}_n) \Delta t \right) \f] where \f$ \boldsymbol{y}_n \f$ is the current value, \f$ \Delta t \f$ is
     * the current time step size, and \f$ f(\boldsymbol{y}) \f$ is the rate of change.
     *
     * The function also estimates the error of the time step and adjusts the time step size for the next iteration
     * based on the desired tolerance \ref TOL. The new time step size is computed as: \f[ \Delta t_{\text{new}} =
     * \Delta t \cdot \min\left(2, \max\left(0.2, 0.9 \sqrt{\frac{TOL}{EST}}\right)\right) \f] where \f$ EST \f$ is the
     * estimated error.
     *
     * @todo: reuse explicitEulerRichardson function?
     */
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
     * @brief Computes the direction cosines between a given coordinate system and the global coordinate system
     * @param transformedCoordinateSystem 3x3 matrix representing the transformed coordinate system
     * @return 3x3 matrix of direction cosines
     *
     * The direction cosines are computed as the dot products between the basis vectors of the transformed coordinate
     * system and the global coordinate system.
     */
    Matrix3d directionCosines( const Matrix3d& transformedCoordinateSystem );

    /**
     * @brief Constructs a orthonormal coordinate system with a given normal vector as \f$x_1\f$-axis.
     * @param normalVector Input normal vector, which will be normalized and used as \f$x_1\f$-axis.
     * @return Orthonormal coordinate system as a 3x3 matrix.
     *
     * The orthonormal coordinate system is constructed such that:
     * - The first column corresponds to the normalized input normal vector (\f$x_1\f$-axis).
     * - The second column is a unit vector orthogonal to the first (\f$x_2\f$-axis).
     * - The third column is the cross product of the first two columns (\f$x_3\f$-axis).
     *
     * @note Do not use this function if you want to control the direction of the axes in the plane orthogonal to the
     * normal vector.
     * @todo: Maybe remove this function completely to avoid mistakes?
     */
    Matrix3d orthonormalCoordinateSystem( Vector3d& normalVector );

    /**
     * @brief Constructs an orthonormal coordinate system with two given normal vectors.
     * @param n1 First input vector.
     * @param n2 Second input vector.
     * @return Orthonormal coordinate system as a 3x3 matrix.
     *
     * @throws std::invalid_argument if the input normal vectors are not orthogonal.
     *
     * The orthonormal coordinate system is constructed such that:
     * - The first column corresponds to the normalized first input normal vector (\f$x_1\f$-axis).
     * - The second colum corresponds to the normalized second input normal vector (\f$x_2\f$-axis).
     * - The third column is the normalized cross product of the first two columns (\f$x_3\f$-axis).
     */
    Matrix3d orthonormalCoordinateSystem( const Vector3d& n1, const Vector3d& n2 );

    Matrix3d transformToLocalSystem( const Matrix3d& T, const Matrix3d& transformedCoordinateSystem );

    Matrix3d transformToGlobalSystem( const Matrix3d& T, const Matrix3d& transformedCoordinateSystem );

  } // namespace Math
} // namespace Marmot
