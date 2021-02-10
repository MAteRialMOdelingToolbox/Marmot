#pragma once
#include "Fastor/Fastor.h"
#include <iostream>

namespace Marmot {
    namespace ContinuumMechanics::TensorUtility {
        namespace TensorExponential {
            struct ExponentialMapFailed : std::exception {
            };

            template < typename T, size_t tensorSize >
            Fastor::Tensor< T, tensorSize, tensorSize > computeTensorExponentialOnly(
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
                    norm = Fastor::norm( series );
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

            template < typename T, size_t tensorSize >
            TensorExponentialResult< T, tensorSize > computeTensorExponentialAndDerivate(
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
        } // namespace TensorExponential
    }     // namespace ContinuumMechanics::TensorUtility
} // namespace Marmot
