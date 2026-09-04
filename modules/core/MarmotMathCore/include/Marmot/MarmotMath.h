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
#include <autodiff/forward/dual/dual.hpp>
#include <autodiff/forward/real/real.hpp>
#include <complex>

namespace Marmot::Math {

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

  /** @brief Computes the exponential of value @p x with numerical limits check
   *  @param x Exponent to which e is raised
   *  @return exponential of x
   *
   * If x is larger than the maximum limit of double precision floating point numbers,
   * the maximum limit is returned. If @p x is smaller than the minimum limit, the minimum limit is returned.
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
   *  @tparam G Gradient type of the autodiff::dual number
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

  /**
   * @brief Transforms a second-order tensor from the global coordinate system to a local coordinate system.
   * @param T The tensor to be transformed, represented as a 3x3 matrix in the global coordinate system.
   * @param transformedCoordinateSystem The 3x3 transformation matrix representing the local coordinate system axes in
   * global coordinates.
   * @return The tensor represented in the local coordinate system as a 3x3 matrix.
   *
   * The transformation is performed as: \f$ T_{local} = Q^T \, T_{global} \, Q \f$, where \f$ Q \f$ is the
   * transformation matrix.
   */
  Matrix3d transformToLocalSystem( const Matrix3d& T, const Matrix3d& transformedCoordinateSystem );

  /**
   * @brief Transforms a second-order tensor from a local coordinate system to the global coordinate system.
   * @param T The tensor to be transformed, represented as a 3x3 matrix in the local coordinate system.
   * @param transformedCoordinateSystem The 3x3 transformation matrix representing the local coordinate system axes in
   * global coordinates.
   * @return The tensor represented in the global coordinate system as a 3x3 matrix.
   *
   * The transformation is performed as: \f$ T_{global} = Q \, T_{local} \, Q^T \f$, where \f$ Q \f$ is the
   * transformation matrix.
   */
  Matrix3d transformToGlobalSystem( const Matrix3d& T, const Matrix3d& transformedCoordinateSystem );

} // namespace Marmot::Math
