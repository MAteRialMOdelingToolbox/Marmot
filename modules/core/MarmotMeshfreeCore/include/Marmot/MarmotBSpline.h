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

template <>
double inline B< 0 >( double u, const double* knotVec, int i )
{
  if ( knotVec[i] <= u && u < knotVec[i + 1] )
    return 1;
  return 0;
}

template <>
double inline dB_dU< 0 >( double u, const double* knotVec, int i )
{
  return 0;
}
