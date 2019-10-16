#pragma once
#include "bftFunctions.h"
#include "bftTypedefs.h"
#include "bftVoigt.h"
#include <iostream>
#include <utility>

namespace bft {
    namespace CommonTensors {
        extern const Tensor3333d I;
        extern const Tensor3333d Isym;
        extern const Tensor3333d Iskew;
        extern const Tensor3333d dDeviatoricStress_dStress;
    } // namespace CommonTensors

    namespace TensorUtility {

        constexpr int d( int a, int b ) { return a == b ? 1 : 0; }

        template <int x, int y, typename T, typename = std::enable_if<!std::is_const<std::remove_reference<T>>::value>>
        auto as( T& t )
        {
            return Eigen::Map<Eigen::Matrix<double, x, y>>( t.data() );
        }

        template <int x, int y, typename T, typename = void>
        auto as( const T& t )
        {
            return Eigen::Map<const Eigen::Matrix<double, x, y>>( t.data() );
        }

        template <typename Derived, typename = std::enable_if<!std::is_const<std::remove_reference<Derived>>::value>>
        auto flatten( Derived& t )
        {
            return Eigen::Map<Eigen::Matrix<double, Derived::RowsAtCompileTime * Derived::ColsAtCompileTime, 1>>(
                t.data() );
        }

        template <typename Derived, typename = void>
        auto flatten( const Derived& t )
        {
            return Eigen::Map<const Eigen::Matrix<double, Derived::RowsAtCompileTime * Derived::ColsAtCompileTime, 1>>(
                t.data() );
        }

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

                throw std::invalid_argument( MakeString()
                                             << __PRETTY_FUNCTION__ << ": invalid dimension / voigt index specified" );
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
            Eigen::TensorFixedSize<double, Eigen::Sizes<VOIGTFROMDIM( nDim ), nDim, nDim>> voigtMap()
            {
                using namespace Eigen;
                Eigen::TensorFixedSize<double, Eigen::Sizes<VOIGTFROMDIM( nDim ), nDim, nDim>> result;
                result.setZero();
                for ( int i = 0; i < nDim; i++ )
                    for ( int j = 0; j < nDim; j++ )
                        result( toVoigt<nDim>( i, j ), i, j ) = 1;
                return result;
            }

        } // namespace IndexNotation

        // namespace VoigtNotation

    } // namespace TensorUtility

} // namespace bft
