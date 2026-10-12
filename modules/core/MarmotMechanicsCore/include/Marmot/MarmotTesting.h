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

#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotMaterialPointSolverHypoElastic.h"
#include "autodiff/forward/dual/dual.hpp"
#include "unsupported/Eigen/CXX11/Tensor"
#include <Eigen/Core>
#include <Fastor/Fastor.h>
#include <iostream>

/**
 * @file MarmotTesting.h
 * @brief Lightweight unit-test helpers for Marmot material and numerical routines.
 *
 * Provides templated equality checks (@c checkIfEqual) for scalars, Eigen
 * matrices, Eigen tensors and Fastor tensors, as well as higher-level utilities
 * such as @c throwExceptionOnFailure, the Fibonacci hemisphere sampler
 * @c fibonacciLatticeHemisphere, and the @c spinTurbokreisel material-point test.
 */

namespace Marmot::Testing {

  /**
   * @brief Check whether two `double` values are equal within a tolerance.
   * @param a    First value.
   * @param b    Second value.
   * @param tol  Absolute tolerance (default: 1e-15).
   * @return `true` if `|a - b| <= tol`, `false` otherwise.
   */
  bool checkIfEqual( const double a, const double b, const double tol = 1e-15 );

  /**
   * @brief Check whether two `autodiff::dual` values are equal within a tolerance.
   * @param a    First dual number.
   * @param b    Second dual number.
   * @param tol  Absolute tolerance applied to the real part (default: 1e-15).
   * @return `true` if the values compare equal within @p tol, `false` otherwise.
   */
  bool checkIfEqual( const autodiff::dual a, const autodiff::dual b, const double tol = 1e-15 );

  /**
   * @brief Check whether two `std::complex<double>` values are equal within a tolerance.
   * @param a    First complex value.
   * @param b    Second complex value.
   * @param tol  Absolute tolerance (default: 1e-15).
   * @return `true` if both real and imaginary parts are within @p tol, `false` otherwise.
   */
  bool checkIfEqual( const std::complex< double > a, const std::complex< double > b, const double tol = 1e-15 );

  /**
   * @brief Convert a `double` to its string representation.
   * @param a  Value to convert.
   * @return   String representation of @p a.
   */
  std::string getString( const double a );

  /**
   * @brief Convert an `autodiff::dual` to its string representation.
   * @param a  Dual number to convert.
   * @return   String representation of @p a.
   */
  std::string getString( const autodiff::dual a );

  /**
   * @brief Check whether two dynamic Eigen matrices are element-wise equal within a tolerance.
   *
   * Prints a hint to `stdout` identifying the first mismatching entry before returning `false`.
   *
   * @tparam T   Scalar type of the matrices.
   * @param a    First matrix.
   * @param b    Second matrix.
   * @param tol  Absolute tolerance per element (default: 1e-15).
   * @return `true` if all elements satisfy `|a(i,j) - b(i,j)| <= tol`, `false` otherwise.
   */
  template < typename T >
  bool checkIfEqual( const Eigen::Matrix< T, -1, -1 >& a,
                     const Eigen::Matrix< T, -1, -1 >& b,
                     const double                      tol = 1e-15 )
  {
    if ( a.rows() != b.rows() || a.cols() != b.cols() ) {
      return false;
    }
    for ( int i = 0; i < a.rows(); i++ ) {
      for ( int j = 0; j < a.cols(); j++ ) {
        auto cond = checkIfEqual( a( i, j ), b( i, j ), tol );
        if ( !cond ) {
          std::cout << "  -> HINT:  a(" << i << "," << j << ") = " << getString( a( i, j ) ) << " !=  b(" << i << ","
                    << j << ") =" << getString( b( i, j ) ) << std::endl;
          return false;
        }
      }
    }
    return true;
  }

  /**
   * @brief Check whether two fixed-size Eigen tensors are element-wise equal within a tolerance.
   *
   * Prints a hint to `stdout` identifying the first mismatching entry before returning `false`.
   *
   * @tparam T     Scalar type of the tensors.
   * @tparam Rest  Compile-time dimension pack.
   * @param a    First tensor.
   * @param b    Second tensor.
   * @param tol  Absolute tolerance per element (default: 1e-15).
   * @return `true` if all elements satisfy `|a(i) - b(i)| <= tol`, `false` otherwise.
   */
  template < typename T, long int... Rest >
  bool checkIfEqual( const Eigen::TensorFixedSize< T, Eigen::Sizes< Rest... > >& a,
                     const Eigen::TensorFixedSize< T, Eigen::Sizes< Rest... > >& b,
                     const double                                                tol = 1e-15 )
  {
    const T* a_data = a.data();
    const T* b_data = b.data();

    for ( int i = 0; i < a.size(); i++ ) {
      auto cond = checkIfEqual( a_data[i], b_data[i], tol );
      if ( !cond ) {
        std::cout << "  -> HINT:  a(" << i << ") = " << getString( a_data[i] ) << " !=  b(" << i
                  << ") =" << getString( b_data[i] ) << std::endl;
        return false;
      }
    }
    return true;
  }

  /**
   * @brief Check whether two Fastor tensors are element-wise equal within a tolerance.
   *
   * Prints a hint to `stdout` identifying the first mismatching entry before returning `false`.
   *
   * @tparam T     Scalar type of the tensors.
   * @tparam Rest  Compile-time dimension pack.
   * @param a    First tensor.
   * @param b    Second tensor.
   * @param tol  Absolute tolerance per element (default: 1e-15).
   * @return `true` if all elements satisfy `|a(i) - b(i)| <= tol`, `false` otherwise.
   */
  template < typename T, size_t... Rest >
  bool checkIfEqual( const Fastor::Tensor< T, Rest... >& a,
                     const Fastor::Tensor< T, Rest... >& b,
                     const double                        tol = 1e-15 )
  {
    const T* a_data = a.data();
    const T* b_data = b.data();

    for ( size_t i = 0; i < a.size(); i++ ) {
      auto cond = checkIfEqual( a_data[i], b_data[i], tol );
      if ( !cond ) {
        std::cout << "  -> HINT:  a(" << i << ") = " << getString( a_data[i] ) << " !=  b(" << i
                  << ") =" << getString( b_data[i] ) << std::endl;
        return false;
      }
    }
    return true;
  }
  /**
   * @brief Throw a `std::runtime_error` when @p condition is `false`.
   * @param condition  The condition to check; an exception is thrown when `false`.
   * @param message    Optional message included in the exception (default: "").
   */
  void throwExceptionOnFailure( const bool condition, const std::string& message = "" );

  /**
   * @brief Execute a list of test functions and collect any exceptions they throw.
   *
   * Each function in @p testFunctions is called in order. Exceptions are caught,
   * their messages accumulated, and a single combined exception is re-thrown after
   * all tests have run so that every failure is reported.
   *
   * @param testFunctions  Vector of zero-argument callables to execute.
   */
  void executeTestsAndCollectExceptions( const std::vector< std::function< void() > >& testFunctions );

  /**
   * @brief Generate a Fibonacci-lattice sampling of the hemisphere.
   *
   * Returns @p N evenly-distributed (φ, θ) angle pairs on the upper hemisphere
   * using the golden-angle Fibonacci lattice.
   *
   * @tparam N  Number of sample points.
   * @return    An Nx2 matrix where each row is `(phi, theta)` in radians.
   */
  template < int N >
  Eigen::Matrix< double, N, 2 > fibonacciLatticeHemisphere()
  {
    double                        phi;
    double                        theta;
    Eigen::Matrix< double, 1, 2 > pt;
    Eigen::Matrix< double, N, 2 > pts;

    // i1 = 1, 2, ..., N
    for ( int i1 = 1; i1 <= N; i1++ ) {
      theta = acos( -( i1 - 1. ) / N );
      phi   = ( i1 - 1. ) * Constants::GoldenAngle;

      pts.row( i1 - 1 ) = Eigen::Vector2d( phi, theta );
    }

    return pts;
  };

  /**
   * @brief Stress/stiffness objectivity test ("spinning top" test) for hypo-elastic materials.
   *
   * Applies a sequence of large rigid-body rotations to a material point and verifies
   * that the Cauchy stress and algorithmic tangent remain objective to within the
   * specified tolerances.
   *
   * @param solver       The hypo-elastic material-point solver to test.
   * @param stressTol    Tolerance for the stress objectivity check (default: 1e-15).
   * @param stiffnessTol Tolerance for the stiffness objectivity check (default: 1e-15).
   * @return `true` if the solver passes both checks, `false` otherwise.
   */
  bool spinTurbokreisel( Marmot::Solvers::MarmotMaterialPointSolverHypoElastic& solver,
                         double                                                 stressTol    = 1e-15,
                         double                                                 stiffnessTol = 1e-15 );
} // namespace Marmot::Testing
