/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck
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
#include <cmath>

/**
 * @file MarmotBSpline.h
 * @brief B-spline basis functions of arbitrary degree and their derivatives (Cox--de Boor recursion).
 *
 * The functions evaluate the @f$ i @f$-th B-spline basis function @f$ N_{i,p}(u) @f$ of degree @f$ p @f$ over a
 * non-decreasing knot vector @f$ \{z_0, z_1, \dots\} @f$, and its derivative @f$ \mathrm{d}N_{i,p}/\mathrm{d}u @f$.
 * The degree is a template parameter, so the recursion is resolved at compile time. They are declared in the
 * global namespace.
 */

/**
 * @brief Evaluate the B-spline basis function @f$ N_{i,p}(u) @f$ by the Cox--de Boor recursion.
 *
 * @f[
 *   N_{i,p}(u) = \frac{u - z_i}{z_{i+p} - z_i} N_{i,p-1}(u)
 *              + \frac{z_{i+p+1} - u}{z_{i+p+1} - z_{i+1}} N_{i+1,p-1}(u),
 * @f]
 * where a term is omitted (@f$ 0/0 := 0 @f$) if its knot difference is below @f$ 10^{-14} @f$ (repeated knots).
 * The recursion ends at B< 0 >.
 *
 * @tparam p Degree @f$ p @f$ of the B-spline (order @f$ p + 1 @f$).
 * @param[in] u       Parameter value @f$ u @f$.
 * @param[in] knotVec Knot vector @f$ z @f$; entries @f$ z_i, \dots, z_{i+p+1} @f$ are accessed.
 * @param[in] i       Index @f$ i @f$ of the basis function.
 * @return @f$ N_{i,p}(u) @f$.
 */
template < int p >
double B( double u, const double* knotVec, int i )
{
  const auto& z = knotVec;
  return
    // clang-format off
      std::abs( ( z[p+i]   - z[i] )  >= 1e-14 ? ( u        - z[i] ) / (z[p+i] - z[i]     ) * B<p-1>(u, z, i)  : 0 )
      +
      std::abs( ( z[p+i+1] - z[i+1]) >= 1e-14 ? ( z[i+p+1] - u    ) / (z[p+i+1] - z[i+1] ) * B<p-1>(u, z, i+1) : 0 )
      ;
  // clang-format on
}

/**
 * @brief Evaluate the derivative @f$ \mathrm{d}N_{i,p}/\mathrm{d}u @f$ of a B-spline basis function.
 *
 * Differentiates the Cox--de Boor recursion of B() term by term,
 * @f[
 *   \frac{\mathrm{d}N_{i,p}}{\mathrm{d}u} = \frac{N_{i,p-1} + (u - z_i)\, N'_{i,p-1}}{z_{i+p} - z_i}
 *   + \frac{-N_{i+1,p-1} + (z_{i+p+1} - u)\, N'_{i+1,p-1}}{z_{i+p+1} - z_{i+1}},
 * @f]
 * omitting terms with a knot difference below @f$ 10^{-14} @f$. The recursion ends at dB_dU< 0 > @f$ = 0 @f$.
 *
 * @tparam p Degree @f$ p @f$ of the B-spline.
 * @param[in] u       Parameter value @f$ u @f$.
 * @param[in] knotVec Knot vector @f$ z @f$; entries @f$ z_i, \dots, z_{i+p+1} @f$ are accessed.
 * @param[in] i       Index @f$ i @f$ of the basis function.
 * @return @f$ \mathrm{d}N_{i,p}/\mathrm{d}u @f$.
 */
template < int p >
double dB_dU( double u, const double* knotVec, int i )
{
  const auto& z = knotVec;
  return
    // clang-format off
      ( std::abs( z[p+i]   - z[i] )  >= 1e-14 ?
        ( 1               ) / (z[p+i] - z[i]     ) * B<p-1>(u, z, i)  +
        ( u        - z[i] ) / (z[p+i] - z[i]     ) * dB_dU<p-1>(u, z, i)

        : 0 )
      +
      ( std::abs( z[p+i+1] - z[i+1]) >= 1e-14 ?
        (          - 1    ) / (z[p+i+1] - z[i+1] ) * B<p-1>(u, z, i+1) +
        ( z[i+p+1] - u    ) / (z[p+i+1] - z[i+1] ) * dB_dU<p-1>(u, z, i+1)
        : 0 )
      ;
  // clang-format on
}

/**
 * @brief Degree-0 B-spline, the indicator function of the half-open knot span @f$ [z_i, z_{i+1}) @f$.
 * @param[in] u       Parameter value @f$ u @f$.
 * @param[in] knotVec Knot vector @f$ z @f$.
 * @param[in] i       Index @f$ i @f$ of the basis function.
 * @return @f$ 1 @f$ if @f$ z_i \le u < z_{i+1} @f$, else @f$ 0 @f$.
 */
template <>
double inline B< 0 >( double u, const double* knotVec, int i )
{
  if ( knotVec[i] <= u && u < knotVec[i + 1] )
    return 1;
  return 0;
}

/**
 * @brief Derivative of the degree-0 B-spline, which is zero (almost everywhere).
 * @param[in] u       Parameter value @f$ u @f$ (unused).
 * @param[in] knotVec Knot vector (unused).
 * @param[in] i       Index of the basis function (unused).
 * @return @f$ 0 @f$.
 */
template <>
double inline dB_dU< 0 >( double u, const double* knotVec, int i )
{
  return 0;
}
