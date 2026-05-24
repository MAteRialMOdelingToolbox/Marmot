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
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMath.h"
#include <algorithm>
#include <cmath>
#include <utility>

namespace Marmot {
  namespace Math {

    template < typename T = double >
    std::pair< FastorStandardTensors::Tensor3t< T >, FastorStandardTensors::Tensor33t< T > > computeEigenSystemJacobi(
      const FastorStandardTensors::Tensor33t< T >& A_in )
    {
      using namespace FastorStandardTensors;

      Tensor33t< T > A = A_in; // Copy to manipulate internally
      Tensor33t< T > V;

      // Initialize Eigenvectors to Identity Matrix
      V( 0, 0 ) = T( 1 );
      V( 0, 1 ) = T( 0 );
      V( 0, 2 ) = T( 0 );
      V( 1, 0 ) = T( 0 );
      V( 1, 1 ) = T( 1 );
      V( 1, 2 ) = T( 0 );
      V( 2, 0 ) = T( 0 );
      V( 2, 1 ) = T( 0 );
      V( 2, 2 ) = T( 1 );

      int max_sweeps = 50;

      for ( int sweep = 0; sweep < max_sweeps; ++sweep ) {
        // Evaluate primal sum of off-diagonals to check for convergence
        double sm = std::abs( Math::makeReal( A( 0, 1 ) ) ) + std::abs( Math::makeReal( A( 0, 2 ) ) ) +
                    std::abs( Math::makeReal( A( 1, 2 ) ) );

        if ( sm < 1e-13 )
          break; // Converged

        // Iterate over upper off-diagonal elements
        for ( int p = 0; p < 2; ++p ) {
          for ( int q = p + 1; q < 3; ++q ) {

            double apq_real = std::abs( Math::makeReal( A( p, q ) ) );

            // Only rotate if the off-diagonal is non-zero in the primal space
            if ( apq_real > 1e-14 ) {
              T      h = A( q, q ) - A( p, p );
              T      t;
              double h_real = std::abs( Math::makeReal( h ) );

              // Prevent division by zero and catastrophic cancellation
              if ( apq_real < 1e-6 * h_real ) {
                t = A( p, q ) / h;
              }
              else {
                T theta = T( 0.5 ) * h / A( p, q );
                T denom = sqrt( T( 1.0 ) + theta * theta );
                // Branch strictly on primal sign to maintain derivative continuity
                if ( Math::makeReal( theta ) < 0.0 ) {
                  t = T( -1.0 ) / ( -theta + denom );
                }
                else {
                  t = T( 1.0 ) / ( theta + denom );
                }
              }

              T c     = T( 1.0 ) / sqrt( T( 1.0 ) + t * t );
              T s     = t * c;
              T tau   = s / ( T( 1.0 ) + c );
              T h_val = t * A( p, q );

              // Apply rotations to A
              A( p, p ) -= h_val;
              A( q, q ) += h_val;
              A( p, q ) = T( 0.0 );
              A( q, p ) = T( 0.0 );

              for ( int r = 0; r < 3; ++r ) {
                if ( r != p && r != q ) {
                  T arp     = A( r, p );
                  T arq     = A( r, q );
                  A( r, p ) = arp - s * ( arq + arp * tau );
                  A( r, q ) = arq + s * ( arp - arq * tau );
                  A( p, r ) = A( r, p );
                  A( q, r ) = A( r, q );
                }
              }

              // Apply rotations to V (Eigenvectors)
              for ( int r = 0; r < 3; ++r ) {
                T vrp     = V( r, p );
                T vrq     = V( r, q );
                V( r, p ) = vrp - s * ( vrq + vrp * tau );
                V( r, q ) = vrq + s * ( vrp - vrq * tau );
              }
            }
          }
        }
      }

      Tensor3t< T > eigenvalues;
      eigenvalues( 0 ) = A( 0, 0 );
      eigenvalues( 1 ) = A( 1, 1 );
      eigenvalues( 2 ) = A( 2, 2 );

      // Sort descending (must swap eigenvectors to match)
      for ( int i = 0; i < 2; ++i ) {
        for ( int j = i + 1; j < 3; ++j ) {
          if ( Math::makeReal( eigenvalues( i ) ) < Math::makeReal( eigenvalues( j ) ) ) {
            std::swap( eigenvalues( i ), eigenvalues( j ) );
            for ( int r = 0; r < 3; ++r ) {
              std::swap( V( r, i ), V( r, j ) );
            }
          }
        }
      }

      return std::make_pair( eigenvalues, V );
    }

    template < typename T = double >
    std::pair< FastorStandardTensors::Tensor3t< T >, FastorStandardTensors::Tensor33t< T > > computeEigenSystemWithCardano(
      const FastorStandardTensors::Tensor33t< T >& A )
    {
      using namespace FastorStandardTensors;

      const T trA = trace( A );

      const T p1 = -trA;
      const T p2 = 0.5 * ( trA * trA - trace( A % A ) );
      const T p3 = -det( A );

      const T a = p2 - ( p1 * p1 ) / 3.0;
      const T b = ( 2.0 * p1 * p1 * p1 ) / 27.0 - ( p1 * p2 ) / 3.0 + p3;

      Tensor3t< T > eigenvalues;

      // 1. SAFE BRANCHING: Use the primal value to check for isotropic/repeated states
      // This prevents the derivative of sqrt(-a) from exploding when 'a' approaches 0.
      if ( std::abs( Math::makeReal( a ) ) < 1e-14 ) {
        T root           = -p1 / 3.0;
        eigenvalues( 0 ) = root;
        eigenvalues( 1 ) = root;
        eigenvalues( 2 ) = root;
      }
      else {
        T r = sqrt( -a / 3.0 );

        // Protect against division by zero in the dual part
        T aux = 3.0 * b / ( 2.0 * a ) * sqrt( -3.0 / a );

        // 2. SAFE CLAMPING: We must evaluate the primal part to clamp,
        // but if we hardcode `aux = 1.0`, we kill the derivative.
        // A common AD trick is to let the dual library handle the clamp,
        // or re-assign with a zeroed dual part if exactly at the boundary.
        double aux_real = Math::makeReal( aux );
        if ( aux_real <= -1.0 ) {
          aux = T( -1.0 ); // Kills gradient at exact boundary to prevent NaN in acos derivative
        }
        else if ( aux_real >= 1.0 ) {
          aux = T( 1.0 ); // Kills gradient at exact boundary to prevent NaN in acos derivative
        }

        T theta = acos( aux );

        for ( int k = 0; k < 3; ++k ) {
          // T(2.0) and T(3.0) ensure we don't accidentally cast down to double
          eigenvalues( k ) = T( 2.0 ) * r * cos( ( theta - T( 2.0 ) * T( M_PI ) * T( k ) ) / T( 3.0 ) ) - p1 / T( 3.0 );
        }
      }

      // Sort descending (using primal values for the logic operators)
      if ( Math::makeReal( eigenvalues( 0 ) ) < Math::makeReal( eigenvalues( 1 ) ) )
        std::swap( eigenvalues( 0 ), eigenvalues( 1 ) );
      if ( Math::makeReal( eigenvalues( 1 ) ) < Math::makeReal( eigenvalues( 2 ) ) )
        std::swap( eigenvalues( 1 ), eigenvalues( 2 ) );
      if ( Math::makeReal( eigenvalues( 0 ) ) < Math::makeReal( eigenvalues( 1 ) ) )
        std::swap( eigenvalues( 0 ), eigenvalues( 1 ) );

      // ------------------------------------------------------------------
      // EIGENVECTOR COMPUTATION
      // ------------------------------------------------------------------
      Tensor33t< T > eigenvectors;

      auto get_eigenvector = [&]( T lambda, Tensor3t< T >& v ) {
        T B[3][3] = { { A( 0, 0 ) - lambda, A( 0, 1 ), A( 0, 2 ) },
                      { A( 1, 0 ), A( 1, 1 ) - lambda, A( 1, 2 ) },
                      { A( 2, 0 ), A( 2, 1 ), A( 2, 2 ) - lambda } };

        T v1[3] = { B[1][1] * B[2][2] - B[1][2] * B[2][1],
                    B[1][2] * B[2][0] - B[1][0] * B[2][2],
                    B[1][0] * B[2][1] - B[1][1] * B[2][0] };
        T v2[3] = { B[2][1] * B[0][2] - B[2][2] * B[0][1],
                    B[2][2] * B[0][0] - B[2][0] * B[0][2],
                    B[2][0] * B[0][1] - B[2][1] * B[0][0] };
        T v3[3] = { B[0][1] * B[1][2] - B[0][2] * B[1][1],
                    B[0][2] * B[1][0] - B[0][0] * B[1][2],
                    B[0][0] * B[1][1] - B[0][1] * B[1][0] };

        // 3. USE PRIMAL FOR NORM COMPARISONS, BUT PRESERVE `T` FOR THE DIVISION
        double n1_real = Math::makeReal( T( sqrt( v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2] ) ) );
        double n2_real = Math::makeReal( T( sqrt( v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2] ) ) );
        double n3_real = Math::makeReal( T( sqrt( v3[0] * v3[0] + v3[1] * v3[1] + v3[2] * v3[2] ) ) );

        double max_norm_real = n1_real;
        T*     best_v        = v1;

        if ( n2_real > max_norm_real ) {
          max_norm_real = n2_real;
          best_v        = v2;
        }
        if ( n3_real > max_norm_real ) {
          max_norm_real = n3_real;
          best_v        = v3;
        }

        // Calculate the actual dual norm strictly for the winning vector
        T max_norm_dual = sqrt( best_v[0] * best_v[0] + best_v[1] * best_v[1] + best_v[2] * best_v[2] );

        if ( max_norm_real > 1e-10 ) {
          v( 0 ) = best_v[0] / max_norm_dual; // Divide by Dual to keep gradient!
          v( 1 ) = best_v[1] / max_norm_dual;
          v( 2 ) = best_v[2] / max_norm_dual;
        }
        else {
          v( 0 ) = T( 0.0 );
          v( 1 ) = T( 0.0 );
          v( 2 ) = T( 0.0 );
        }
      };

      auto make_orthogonal = []( const Tensor3t< T >& v, Tensor3t< T >& n ) {
        // Use real parts to find the smallest axis
        double abs_x = std::abs( Math::makeReal( v( 0 ) ) );
        double abs_y = std::abs( Math::makeReal( v( 1 ) ) );
        double abs_z = std::abs( Math::makeReal( v( 2 ) ) );

        Tensor3t< T > axis;
        axis( 0 ) = T( 0 );
        axis( 1 ) = T( 0 );
        axis( 2 ) = T( 0 );
        if ( abs_x <= abs_y && abs_x <= abs_z )
          axis( 0 ) = T( 1.0 );
        else if ( abs_y <= abs_x && abs_y <= abs_z )
          axis( 1 ) = T( 1.0 );
        else
          axis( 2 ) = T( 1.0 );

        n( 0 ) = v( 1 ) * axis( 2 ) - v( 2 ) * axis( 1 );
        n( 1 ) = v( 2 ) * axis( 0 ) - v( 0 ) * axis( 2 );
        n( 2 ) = v( 0 ) * axis( 1 ) - v( 1 ) * axis( 0 );

        // Normalize using dual math
        T len = sqrt( n( 0 ) * n( 0 ) + n( 1 ) * n( 1 ) + n( 2 ) * n( 2 ) );
        n( 0 ) /= len;
        n( 1 ) /= len;
        n( 2 ) /= len;
      };

      Tensor3t< T > e1, e2, e3;

      double max_eig = std::max( { std::abs( Math::makeReal( eigenvalues( 0 ) ) ),
                                   std::abs( Math::makeReal( eigenvalues( 1 ) ) ),
                                   std::abs( Math::makeReal( eigenvalues( 2 ) ) ) } );
      double tol     = 1e-6 * std::max( max_eig, 1.0 );

      bool root_01_repeated = std::abs( Math::makeReal( T( eigenvalues( 0 ) - eigenvalues( 1 ) ) ) ) < tol;
      bool root_12_repeated = std::abs( Math::makeReal( T( eigenvalues( 1 ) - eigenvalues( 2 ) ) ) ) < tol;

      if ( root_01_repeated && root_12_repeated ) {
        e1( 0 ) = T( 1.0 );
        e1( 1 ) = T( 0.0 );
        e1( 2 ) = T( 0.0 );
        e2( 0 ) = T( 0.0 );
        e2( 1 ) = T( 1.0 );
        e2( 2 ) = T( 0.0 );
        e3( 0 ) = T( 0.0 );
        e3( 1 ) = T( 0.0 );
        e3( 2 ) = T( 1.0 );
      }
      else if ( root_01_repeated ) {
        get_eigenvector( eigenvalues( 2 ), e3 );
        make_orthogonal( e3, e1 );
        e2( 0 ) = e3( 1 ) * e1( 2 ) - e3( 2 ) * e1( 1 );
        e2( 1 ) = e3( 2 ) * e1( 0 ) - e3( 0 ) * e1( 2 );
        e2( 2 ) = e3( 0 ) * e1( 1 ) - e3( 1 ) * e1( 0 );
      }
      else if ( root_12_repeated ) {
        get_eigenvector( eigenvalues( 0 ), e1 );
        make_orthogonal( e1, e2 );
        e3( 0 ) = e1( 1 ) * e2( 2 ) - e1( 2 ) * e2( 1 );
        e3( 1 ) = e1( 2 ) * e2( 0 ) - e1( 0 ) * e2( 2 );
        e3( 2 ) = e1( 0 ) * e2( 1 ) - e1( 1 ) * e2( 0 );
      }
      else {
        get_eigenvector( eigenvalues( 0 ), e1 );
        get_eigenvector( eigenvalues( 1 ), e2 );
        e3( 0 ) = e1( 1 ) * e2( 2 ) - e1( 2 ) * e2( 1 );
        e3( 1 ) = e1( 2 ) * e2( 0 ) - e1( 0 ) * e2( 2 );
        e3( 2 ) = e1( 0 ) * e2( 1 ) - e1( 1 ) * e2( 0 );
      }

      eigenvectors( 0, 0 ) = e1( 0 );
      eigenvectors( 0, 1 ) = e2( 0 );
      eigenvectors( 0, 2 ) = e3( 0 );
      eigenvectors( 1, 0 ) = e1( 1 );
      eigenvectors( 1, 1 ) = e2( 1 );
      eigenvectors( 1, 2 ) = e3( 1 );
      eigenvectors( 2, 0 ) = e1( 2 );
      eigenvectors( 2, 1 ) = e2( 2 );
      eigenvectors( 2, 2 ) = e3( 2 );

      return std::make_pair( eigenvalues, eigenvectors );
    }

    // template < typename T = double >
    // std::pair< FastorStandardTensors::Tensor3t< T >, FastorStandardTensors::Tensor33t< T > >
    // computeEigenSystemWithCardano( const FastorStandardTensors::Tensor33t< T >& A )
    // {
    //     using namespace FastorStandardTensors;

    //     // ------------------------------------------------------------------
    //     // 1. EIGENVALUE COMPUTATION
    //     // ------------------------------------------------------------------
    //     T p1 = -trace( A );
    //     T p2 = 0.5 * ( trace( A ) * trace( A ) - trace( A % A ) );
    //     T p3 = -det( A );

    //     T a = p2 - ( p1 * p1 ) / 3.0;
    //     T b = ( 2.0 * p1 * p1 * p1 ) / 27.0 - ( p1 * p2 ) / 3.0 + p3;

    //     Tensor3t< T > eigenvalues;

    //     // FIX: Catch the isotropic case to prevent division by zero
    //     if ( std::abs(Math::makeReal(a)) < 1e-14 ) {
    //         T root = -p1 / 3.0;
    //         eigenvalues(0) = root;
    //         eigenvalues(1) = root;
    //         eigenvalues(2) = root;
    //     }
    //     else {
    //         T r = sqrt( -a / 3.0 );
    //         T aux = 3.0 * b / ( 2.0 * a ) * sqrt( -3.0 / a );

    //         if ( Math::makeReal( aux ) <= -1.0 ) aux = -1.0;
    //         else if ( Math::makeReal( aux ) >= 1.0 ) aux = 1.0;

    //         T theta = acos( aux );

    //         for ( int k = 0; k < 3; ++k ) {
    //             eigenvalues( k ) = 2.0 * r * cos( ( theta - 2.0 * M_PI * k ) / 3.0 ) - p1 / 3.0;
    //         }
    //     }

    //     // Sort descending
    //     if (Math::makeReal(eigenvalues(0)) < Math::makeReal(eigenvalues(1))) std::swap(eigenvalues(0),
    //     eigenvalues(1)); if (Math::makeReal(eigenvalues(1)) < Math::makeReal(eigenvalues(2)))
    //     std::swap(eigenvalues(1), eigenvalues(2)); if (Math::makeReal(eigenvalues(0)) <
    //     Math::makeReal(eigenvalues(1))) std::swap(eigenvalues(0), eigenvalues(1));

    //     // ------------------------------------------------------------------
    //     // 2. EIGENVECTOR COMPUTATION (Rank-Deficiency Proofed)
    //     // ------------------------------------------------------------------
    //     Tensor33t< T > eigenvectors;

    //     // Returns norm. If norm < 1e-10, it implies a repeated root.
    //     auto get_eigenvector = [&](T lambda, Tensor3t<T>& v) -> double {
    //         T B[3][3] = {
    //             {A(0,0) - lambda, A(0,1),          A(0,2)},
    //             {A(1,0),          A(1,1) - lambda, A(1,2)},
    //             {A(2,0),          A(2,1),          A(2,2) - lambda}
    //         };

    //         // FIXED: Correctly computing cross products of the rows R1xR2, R2xR0, R0xR1
    //         T v1[3] = {B[1][1]*B[2][2] - B[1][2]*B[2][1], B[1][2]*B[2][0] - B[1][0]*B[2][2], B[1][0]*B[2][1] -
    //         B[1][1]*B[2][0]}; T v2[3] = {B[2][1]*B[0][2] - B[2][2]*B[0][1], B[2][2]*B[0][0] - B[2][0]*B[0][2],
    //         B[2][0]*B[0][1] - B[2][1]*B[0][0]}; T v3[3] = {B[0][1]*B[1][2] - B[0][2]*B[1][1], B[0][2]*B[1][0] -
    //         B[0][0]*B[1][2], B[0][0]*B[1][1] - B[0][1]*B[1][0]};

    //         double n1 = Math::makeReal(T(sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2])));
    //         double n2 = Math::makeReal(T(sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2])));
    //         double n3 = Math::makeReal(T(sqrt(v3[0]*v3[0] + v3[1]*v3[1] + v3[2]*v3[2])));

    //         double max_norm = n1;
    //         T* best_v = v1;
    //         if (n2 > max_norm) { max_norm = n2; best_v = v2; }
    //         if (n3 > max_norm) { max_norm = n3; best_v = v3; }

    //         if (max_norm > 1e-10) {
    //             v(0) = best_v[0] / max_norm;
    //             v(1) = best_v[1] / max_norm;
    //             v(2) = best_v[2] / max_norm;
    //         } else {
    //             v(0) = 0.0; v(1) = 0.0; v(2) = 0.0;
    //         }
    //         return max_norm;
    //     };

    //     // Safe orthogonal vector generation
    //     auto make_orthogonal = [](const Tensor3t<T>& v, Tensor3t<T>& n) {
    //         T abs_x = abs(v(0));
    //         T abs_y = abs(v(1));
    //         T abs_z = abs(v(2));

    //         // Find axis with smallest component to ensure non-zero cross product
    //         Tensor3t<T> axis;
    //         axis(0)=0; axis(1)=0; axis(2)=0;
    //         if (makeReal(abs_x) <= makeReal( abs_y ) && makeReal(abs_x) <= makeReal( abs_z )) axis(0) = 1.0;
    //         else if (makeReal( abs_y ) <= makeReal( abs_x ) && makeReal( abs_y ) <= makeReal( abs_z )) axis(1) = 1.0;
    //         else axis(2) = 1.0;

    //         n(0) = v(1)*axis(2) - v(2)*axis(1);
    //         n(1) = v(2)*axis(0) - v(0)*axis(2);
    //         n(2) = v(0)*axis(1) - v(1)*axis(0);

    //         T len = sqrt(n(0)*n(0) + n(1)*n(1) + n(2)*n(2));
    //         n(0) /= len; n(1) /= len; n(2) /= len;
    //     };

    //     Tensor3t<T> e1, e2, e3;

    //     // Determine eigenvalue multiplicity numerically
    //     // Use a relative tolerance (1e-6 is safe for double precision Cardano errors)
    //     double max_eig = std::max({ std::abs(Math::makeReal(eigenvalues(0))),
    //                                 std::abs(Math::makeReal(eigenvalues(1))),
    //                                 std::abs(Math::makeReal(eigenvalues(2))) });
    //     double tol = 1e-6 * std::max(max_eig, 1.0);

    //     bool root_01_repeated = std::abs(Math::makeReal(T(eigenvalues(0) - eigenvalues(1)))) < tol;
    //     bool root_12_repeated = std::abs(Math::makeReal(T(eigenvalues(1) - eigenvalues(2)))) < tol;

    //     if (root_01_repeated && root_12_repeated) {
    //         // Case D: All roots are repeated (Isotropic)
    //         e1(0)=1.0; e1(1)=0.0; e1(2)=0.0;
    //         e2(0)=0.0; e2(1)=1.0; e2(2)=0.0;
    //         e3(0)=0.0; e3(1)=0.0; e3(2)=1.0;
    //     }
    //     else if (root_01_repeated) {
    //         // Case B: lambda_0 and lambda_1 are repeated. lambda_2 is distinct.
    //         // Calculate e3 from distinct lambda_2, then span the orthogonal plane.
    //         get_eigenvector(eigenvalues(2), e3);
    //         make_orthogonal(e3, e1);
    //         // Ensure right-handed: e2 = e3 x e1
    //         e2(0) = e3(1)*e1(2) - e3(2)*e1(1);
    //         e2(1) = e3(2)*e1(0) - e3(0)*e1(2);
    //         e2(2) = e3(0)*e1(1) - e3(1)*e1(0);
    //     }
    //     else if (root_12_repeated) {
    //         // Case C: lambda_1 and lambda_2 are repeated. lambda_0 is distinct.
    //         // Calculate e1 from distinct lambda_0, then span the orthogonal plane.
    //         get_eigenvector(eigenvalues(0), e1);
    //         make_orthogonal(e1, e2);
    //         // Ensure right-handed: e3 = e1 x e2
    //         e3(0) = e1(1)*e2(2) - e1(2)*e2(1);
    //         e3(1) = e1(2)*e2(0) - e1(0)*e2(2);
    //         e3(2) = e1(0)*e2(1) - e1(1)*e2(0);
    //     }
    //     else {
    //         // Case A: All roots are distinct
    //         get_eigenvector(eigenvalues(0), e1);
    //         get_eigenvector(eigenvalues(1), e2);
    //         // Ensure right-handed: e3 = e1 x e2
    //         e3(0) = e1(1)*e2(2) - e1(2)*e2(1);
    //         e3(1) = e1(2)*e2(0) - e1(0)*e2(2);
    //         e3(2) = e1(0)*e2(1) - e1(1)*e2(0);
    //     }

    //     eigenvectors(0,0) = e1(0); eigenvectors(0,1) = e2(0); eigenvectors(0,2) = e3(0);
    //     eigenvectors(1,0) = e1(1); eigenvectors(1,1) = e2(1); eigenvectors(1,2) = e3(1);
    //     eigenvectors(2,0) = e1(2); eigenvectors(2,1) = e2(2); eigenvectors(2,2) = e3(2);

    //     return std::make_pair(eigenvalues, eigenvectors);
    // }

  } // namespace Math
} // namespace Marmot
