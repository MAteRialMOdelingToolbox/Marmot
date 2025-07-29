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
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include <utility>

namespace Marmot {
  namespace ContinuumMechanics::CommonTensors {

    EigenTensors::Tensor3333d              Initialize_I2xI2();
    inline const EigenTensors::Tensor3333d I2xI2 = Initialize_I2xI2();

    EigenTensors::Tensor3333d              Initialize_Isym();
    inline const EigenTensors::Tensor3333d Isym = Initialize_Isym();

    EigenTensors::Tensor3333d              Initialize_Iskew();
    inline const EigenTensors::Tensor3333d Iskew = Initialize_Iskew();

    EigenTensors::Tensor3333d              Initialize_IFourthOrder();
    inline const EigenTensors::Tensor3333d IFourthOrder = Initialize_IFourthOrder();

    EigenTensors::Tensor3333d              Initialize_IFourthOrderTranspose();
    inline const EigenTensors::Tensor3333d IFourthOrderTranspose = Initialize_IFourthOrderTranspose();

    EigenTensors::Tensor3333d              Initialize_dDeviatoricStress_dStress();
    inline const EigenTensors::Tensor3333d dDeviatoricStress_dStress = Initialize_dDeviatoricStress_dStress();
    ;

    EigenTensors::Tensor333d              Initialize_LeviCivita3D();
    inline const EigenTensors::Tensor333d LeviCivita3D = Initialize_LeviCivita3D();

    EigenTensors::Tensor122d              Initialize_LeviCivita2D();
    inline const EigenTensors::Tensor122d LeviCivita2D = Initialize_LeviCivita2D();

    EigenTensors::Tensor33d              Initialize_I2();
    inline const EigenTensors::Tensor33d I2 = Initialize_I2();

    constexpr int getNumberOfDofForRotation( int nDim )
    {
      if ( nDim == 2 )
        return 1;
      else
        return 3;
    }

    template < int sizeI, int sizeJ >
    Eigen::Matrix< double, sizeI * sizeJ, sizeI * sizeJ > makeIndexSwapTensor()
    {
      // Aux. Matrix, which helps to swap indices in Eigen::Matrices abused as higher order Tensors by
      // multiplication ,
      //
      // T_(ij)(kl) * IndexSwapTensor_(kl)(lk) = T_(ij)(lk)
      //
      // For instance:
      //
      // Matrix3d t;
      // t<<1,2,3,4,5,6,7,8,9;
      // const auto P  = makeIndexSwapTensor<3,3>();
      // std::cout << t.reshaped().transpose() << std::endl;
      // std::cout << t.reshaped().transpose() * P << std::endl;
      //
      // Output:
      //
      // 1,4,7,2,5,8,3,6,9
      // 1,2,3,4,5,6,7,8,9
      //

      Eigen::Matrix< double, sizeI * sizeJ, sizeI * sizeJ > P;

      auto d = []( int a, int b ) -> double { return a == b ? 1.0 : 0.0; };

      for ( int i = 0; i < sizeI; i++ )
        for ( int j = 0; j < sizeJ; j++ )
          for ( int l = 0; l < sizeI; l++ )
            for ( int k = 0; k < sizeJ; k++ )
              P( i + j * sizeI, l * sizeJ + k ) = d( i, l ) * d( k, j );

      return P;
    }

    template < int nDim >
    constexpr Eigen::TensorFixedSize< double, Eigen::Sizes< getNumberOfDofForRotation( nDim ), nDim, nDim > > getReferenceToCorrectLeviCivita()
    // template <int nDim>
    // template <int nDim=2>
    // constexpr Eigen::TensorFixedSize< double, Eigen::Sizes<getNumberOfDofForRotation (nDim), nDim, nDim> >
    // getReferenceToCorrectLeviCivita()
    {
      if constexpr ( nDim == 2 )
        return LeviCivita2D;
      else
        return LeviCivita3D;
    }

  } // namespace ContinuumMechanics::CommonTensors

  namespace ContinuumMechanics::TensorUtility {

    constexpr int d( int a, int b )
    {
      return a == b ? 1 : 0;
    }

    template < int x,
               int y,
               typename T,
               typename = std::enable_if< !std::is_const< std::remove_reference< T > >::value > >
    auto as( T& t )
    {
      return Eigen::Map< Eigen::Matrix< typename T::Scalar, x, y > >( t.data() );
    }

    template < int x, int y, typename T, typename = void >
    auto as( const T& t )
    {
      return Eigen::Map< const Eigen::Matrix< typename T::Scalar, x, y > >( t.data() );
    }

    template < typename Derived,
               typename = std::enable_if< !std::is_const< std::remove_reference< Derived > >::value > >
    auto flatten( Derived& t )
    {
      return Eigen::Map<
        Eigen::Matrix< typename Derived::Scalar, Derived::RowsAtCompileTime * Derived::ColsAtCompileTime, 1 > >(
        t.data() );
    }

    template < typename Derived, typename = void >
    auto flatten( const Derived& t )
    {
      return Eigen::Map<
        const Eigen::Matrix< typename Derived::Scalar, Derived::RowsAtCompileTime * Derived::ColsAtCompileTime, 1 > >(
        t.data() );
    }

    template < int... Pairs >
    constexpr auto contractionDims() // should be evaluated at compile time
    {
      static_assert( sizeof...( Pairs ) % 2 == 0, "Pairs must contain an even number of elements." );

      constexpr int numPairs = sizeof...( Pairs ) / 2;

      Eigen::array< Eigen::IndexPair< int >, numPairs > result{};

      if constexpr ( numPairs > 0 ) {          // return empty array if no pairs are given
        constexpr int values[] = { Pairs... }; // Pairs... expands the parameter pack

        // Fill the result array at compile-time
        for ( int i = 0; i < numPairs; ++i ) {
          result[i] = { values[2 * i], values[2 * i + 1] };
        }
      }
      return result;
    }

    Eigen::Matrix3d dyadicProduct( const Eigen::Vector3d& vector1, const Eigen::Vector3d& vector2 );

    namespace IndexNotation {
      template < int nDim >
      constexpr std::pair< int, int > fromVoigt( int ij )
      {
        if constexpr ( nDim == 1 )
          return std::pair< int, int >( 0, 0 );
        else if ( nDim == 2 )
          switch ( ij ) {
          case 0: return std::pair< int, int >( 0, 0 );
          case 1: return std::pair< int, int >( 1, 1 );
          case 2: return std::pair< int, int >( 0, 1 );
          }

        else if ( nDim == 3 ) {
          switch ( ij ) {
          case 0: return std::pair< int, int >( 0, 0 );
          case 1: return std::pair< int, int >( 1, 1 );
          case 2: return std::pair< int, int >( 2, 2 );
          case 3: return std::pair< int, int >( 0, 1 );
          case 4: return std::pair< int, int >( 0, 2 );
          case 5: return std::pair< int, int >( 1, 2 );
          }
        }

        throw std::invalid_argument( MakeString()
                                     << __PRETTY_FUNCTION__ << ": invalid dimension / voigt index specified" );
      }

      template < int nDim >
      constexpr int toVoigt( int i, int j )
      {
        if constexpr ( nDim == 1 )
          return 0;
        else if ( nDim == 2 )
          return ( i == j ) ? ( i == 0 ? 0 : 1 ) : 2;

        else if ( nDim == 3 ) {
          constexpr int tensor2VoigtNotationIndicesMapping[3][3] = { { 0, 3, 4 }, { 3, 1, 5 }, { 4, 5, 2 } };
          return tensor2VoigtNotationIndicesMapping[i][j];
        }

        throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid dimension specified" );
      }

      template < int nDim >
      Eigen::TensorFixedSize< double, Eigen::Sizes< VOIGTFROMDIM( nDim ), nDim, nDim > > voigtMap()
      {
        using namespace Eigen;
        Eigen::TensorFixedSize< double, Eigen::Sizes< VOIGTFROMDIM( nDim ), nDim, nDim > > result;
        result.setZero();
        for ( int i = 0; i < nDim; i++ )
          for ( int j = 0; j < nDim; j++ )
            result( toVoigt< nDim >( i, j ), i, j ) = 1;
        return result;
      }

    } // namespace IndexNotation

    // namespace ContinuumMechanics::VoigtNotation

  } // namespace ContinuumMechanics::TensorUtility

} // namespace Marmot
