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
#include "Fastor/Fastor.h"
#include "Marmot/MarmotMath.h"
#include <iostream>

namespace Marmot {
  namespace ContinuumMechanics::TensorUtility {
    namespace TensorExponential {
      struct ExponentialMapFailed : std::exception {};

      /**
       * @brief Computes the exponential of a square second rank tensor
       *
       * @details Evaluates exp(@p theTensor) using a series expansion up to
       * @p maxIterations or until the terms converge below @p tolerance.
       *
       * @tparam T Element type of the tensor.
       * @tparam tensorSize Dimension of the square tensor.
       * @param theTensor Input tensor.
       * @param maxIterations Maximum number of series terms to compute.
       * @param tolerance Convergence threshold for the series terms.
       * @param alternativeTolerance Secondary threshold to detect failure.
       * @return Tensor representing the exponential of @p theTensor.
       * @throws ExponentialMapFailed If convergence is not achieved.
       */
      template < typename T, size_t tensorSize >
      Fastor::Tensor< T, tensorSize, tensorSize > computeTensorExponential(
        const Fastor::Tensor< T, tensorSize, tensorSize >& theTensor,
        int                                                maxIterations,
        double                                             tolerance,
        double                                             alternativeTolerance )
      {
        using TensorNN = Fastor::Tensor< T, tensorSize, tensorSize >;

        TensorNN theExponential;

        int n;
        int facN = 1;

        // term 1 ( n = 0 )
        theExponential.eye();

        // series 2 ( n = 1 )
        theExponential += theTensor;

        TensorNN tensorPower = theTensor;
        // loop starts at n = 2
        //
        double norm = 0.0;
        for ( n = 2; n < maxIterations; ++n ) {
          facN *= n;
          tensorPower           = tensorPower % theTensor;
          const TensorNN series = 1. / facN * tensorPower;
          theExponential += series;

          // workaround for autodiff::dual
          Fastor::Tensor< double, tensorSize, tensorSize > series_real;
          for ( size_t i = 0; i < tensorSize; i++ ) {
            for ( size_t j = 0; j < tensorSize; j++ ) {
              series_real( i, j ) = Math::makeReal( series( i, j ) );
            }
          }
          norm = Fastor::norm( series_real );

          if ( norm <= tolerance )
            break;
        }

        if ( n == maxIterations && norm > alternativeTolerance ) {
          throw ExponentialMapFailed();
        }

        return theExponential;
      }

      template < typename T, std::size_t tensorSize >
      struct TensorExponentialResult {
        Fastor::Tensor< T, tensorSize, tensorSize >                         theExponential;
        Fastor::Tensor< T, tensorSize, tensorSize, tensorSize, tensorSize > theExponentialDerivative;
      };

      namespace FirstOrderDerived {
        /// @brief Computes the exponential of a square second rank tensor and its first-order derivative.
        /// @details See the documentation of `computeTensorExponential` without derivative
        ///          for additional details.
        template < typename T, size_t tensorSize >
        TensorExponentialResult< T, tensorSize > computeTensorExponential(
          const Fastor::Tensor< T, tensorSize, tensorSize >& theTensor,
          int                                                maxIterations,
          double                                             tolerance )
        {
          // tensor exponential;
          using Tensor = Fastor::Tensor< T, tensorSize, tensorSize >;
          TensorExponentialResult< T, tensorSize > theResult;

          std::vector< Tensor > tensorPowers;

          int n;
          int facN = 1;

          // term 1 ( n = 0 )
          theResult.theExponential.eye();
          tensorPowers.push_back( theResult.theExponential ); // I

          // series 2 ( n = 1 )
          theResult.theExponential += theTensor;
          tensorPowers.push_back( theTensor ); // T

          // loop starts at n = 2
          for ( n = 2; n < maxIterations; ++n ) {
            facN *= n;
            const Tensor tensorPower = tensorPowers.back() % theTensor;
            const Tensor series      = 1. / facN * tensorPower;
            theResult.theExponential += series;

            tensorPowers.push_back( tensorPower );

            if ( Fastor::norm( series ) <= tolerance )
              break;
          }

          // derivative
          const int N = n;

          facN = 1;
          Fastor::Tensor< T, tensorSize, tensorSize, tensorSize, tensorSize > innerSeries;
          Fastor::Tensor< T, tensorSize, tensorSize, tensorSize, tensorSize > outerSeries;
          outerSeries.zeros();
          for ( n = 1; n <= N; ++n ) {
            innerSeries.zeros();
            for ( int m = 1; m <= n; ++m ) {
              innerSeries += Fastor::outer( tensorPowers[m - 1], tensorPowers[n - m] ); // ik_lj
            }
            facN *= n;
            outerSeries += 1. / facN * innerSeries;
          }
          theResult.theExponentialDerivative = Fastor::permute< Fastor::Index< 0, 3, 1, 2 > >(
            outerSeries ); // ik_lj -> ij_kl

          return theResult;
        }
      } // namespace FirstOrderDerived
    }   // namespace TensorExponential
  }     // namespace ContinuumMechanics::TensorUtility
} // namespace Marmot
