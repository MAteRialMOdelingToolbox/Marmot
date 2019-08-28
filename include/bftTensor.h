#pragma once
#include "bftFunctions.h"
#include "bftVoigt.h"
#include "bftTypedefs.h"
#include <iostream>
#include <utility>

namespace bft {
    namespace TensorUtility {
        namespace IndexNotation {
            template <int nDim>
                constexpr std::pair<int, int> fromVoigt( int ij )
                {
                    if constexpr ( nDim == 1 )
                        return std::pair<int, int>( 0, 0 );

                    else if ( nDim == 2 )
                        switch ( ij ) {
                            case 0: return std::pair<int, int>( 0, 0 );
                            case 1: return std::pair<int, int>( 1, 1 );
                            case 2: return std::pair<int, int>( 0, 1 );
                        }

                    else if ( nDim == 3 ) {
                        switch ( ij ) {
                            case 0: return std::pair<int, int>( 0, 0 );
                            case 1: return std::pair<int, int>( 1, 1 );
                            case 2: return std::pair<int, int>( 2, 2 );
                            case 3: return std::pair<int, int>( 0, 1 );
                            case 4: return std::pair<int, int>( 0, 2 );
                            case 5: return std::pair<int, int>( 1, 2 );
                        }
                    }

                    throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid dimension / voigt index specified" );
                }

            template <int nDim>
                constexpr int toVoigt( int i, int j )
                {
                    if constexpr ( nDim == 1 )
                        return 0;
                    else if ( nDim == 2 )
                        return ( i == j ) ? ( i == 0 ? 0 : 1 ) : 2;

                    else if ( nDim == 3 ) {
                        constexpr int tensor2VgtIndicesMapping[3][3] = {{0, 3, 4}, {3, 1, 5}, {4, 5, 2}};
                        return tensor2VgtIndicesMapping[i][j];
                    }

                    throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid dimension specified" );
                }

            template <int nDim>
                Eigen::TensorFixedSize< double, Eigen::Sizes<  VOIGTFROMDIM(nDim), nDim, nDim> > voigtMap ()
                {
                    using namespace Eigen;
                    Eigen::TensorFixedSize< double, Eigen::Sizes<  VOIGTFROMDIM(nDim), nDim, nDim> >  result;
                    result.setZero();
                    for (int i = 0; i < nDim; i++)
                        for (int j = 0; j < nDim; j++)
                            result(toVoigt<nDim>(i,j), i,j) = 1;
                    return result;

                }

        } // namespace IndexNotation

        // namespace VoigtNotation

    } // namespace TensorUtility

} // namespace bft
