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
#include "Eigen/Core"
#include "Fastor/Fastor.h"
#include "Marmot/MarmotTensor.h"
#include <autodiff/forward/dual/dual.hpp>

namespace Marmot {

  namespace FastorStandardTensors {

    using Tensor3d    = Fastor::Tensor< double, 3 >;
    using Tensor33d   = Fastor::Tensor< double, 3, 3 >;
    using Tensor333d  = Fastor::Tensor< double, 3, 3, 3 >;
    using Tensor3333d = Fastor::Tensor< double, 3, 3, 3, 3 >;

    template < typename T >
    using Tensor3t = Fastor::Tensor< T, 3 >;
    template < typename T >
    using Tensor33t = Fastor::Tensor< T, 3, 3 >;
    template < typename T >
    using Tensor333t = Fastor::Tensor< T, 3, 3, 3 >;
    template < typename T >
    using Tensor3333t = Fastor::Tensor< T, 3, 3, 3, 3 >;

    using TensorMap3d    = Fastor::TensorMap< double, 3 >;
    using TensorMap33d   = Fastor::TensorMap< double, 3, 3 >;
    using TensorMap333d  = Fastor::TensorMap< double, 3, 3, 3 >;
    using TensorMap3333d = Fastor::TensorMap< double, 3, 3, 3, 3 >;

    namespace Spatial3D {
      inline const Tensor33d I = Tensor33d( ( Eigen::Matrix3d() << Eigen::Matrix3d::Identity() ).finished().data(),
                                            Fastor::ColumnMajor );

      inline const Tensor333d LeviCivita = Tensor333d( Marmot::ContinuumMechanics::CommonTensors::LeviCivita3D.data(),
                                                       Fastor::ColumnMajor );

      inline const Tensor3333d IHyd = Tensor3333d( Marmot::ContinuumMechanics::CommonTensors::I2xI2.data(),
                                                   Fastor::ColumnMajor );

      inline const Tensor3333d ISymm = Tensor3333d( Marmot::ContinuumMechanics::CommonTensors::Isym.data(),
                                                    Fastor::ColumnMajor );

      inline const Tensor3333d ISkew = Tensor3333d( Marmot::ContinuumMechanics::CommonTensors::Iskew.data(),
                                                    Fastor::ColumnMajor );

      inline const Tensor3333d I4 = Tensor3333d( Marmot::ContinuumMechanics::CommonTensors::IFourthOrder.data(),
                                                 Fastor::ColumnMajor );

      inline const Tensor3333d
        ITranspose = Tensor3333d( Marmot::ContinuumMechanics::CommonTensors::IFourthOrderTranspose.data(),
                                  Fastor::ColumnMajor );

      inline const Tensor3333d Deviatoric = I4 - 1. / 3 * IHyd;

      inline const Tensor3333d DeviatoricTranspose = Fastor::transpose( DeviatoricTranspose );

      inline const Tensor3333d DeviatoricSymmetric = ISymm - 1. / 3 * IHyd;

    } // namespace Spatial3D

  } // namespace FastorStandardTensors

  namespace ContinuumMechanics::Micropolar::GeneralizedInvariants {

    FastorStandardTensors::Tensor3333d d2J2_dStress_dStress( double a1, double a2 );

  }

  namespace FastorIndices {
    enum { i_, j_, k_, l_, m_, n_, A_, B_, I_, J_, K_, L_, M_, N_, P_ };

    using A    = Fastor::Index< A_ >;
    using Ai   = Fastor::Index< A_, i_ >;
    using AB   = Fastor::Index< A_, B_ >;
    using B    = Fastor::Index< B_ >;
    using IJ   = Fastor::Index< I_, J_ >;
    using IJKL = Fastor::Index< I_, J_, K_, L_ >;
    using IJML = Fastor::Index< I_, J_, M_, L_ >;
    using IK   = Fastor::Index< I_, K_ >;
    using Ii   = Fastor::Index< I_, i_ >;
    using IikK = Fastor::Index< I_, i_, k_, K_ >;
    using Ik   = Fastor::Index< I_, k_ >;
    using IL   = Fastor::Index< I_, L_ >;
    using Im   = Fastor::Index< I_, m_ >;
    using JI   = Fastor::Index< J_, I_ >;
    using JK   = Fastor::Index< J_, K_ >;
    using JL   = Fastor::Index< J_, L_ >;
    using Jk   = Fastor::Index< J_, k_ >;
    using KI   = Fastor::Index< K_, I_ >;
    using KJ   = Fastor::Index< K_, J_ >;
    using KJN  = Fastor::Index< K_, J_, N_ >;
    using KL   = Fastor::Index< K_, L_ >;
    using KLMN = Fastor::Index< K_, L_, M_, N_ >;
    using KLMP = Fastor::Index< K_, L_, M_, P_ >;
    using KLNM = Fastor::Index< K_, L_, N_, M_ >;
    using KLPM = Fastor::Index< K_, L_, P_, M_ >;
    using KLm  = Fastor::Index< K_, L_, m_ >;
    using KMJ  = Fastor::Index< K_, M_, J_ >;
    using KMN  = Fastor::Index< K_, M_, N_ >;
    using Ki   = Fastor::Index< K_, i_ >;
    using Kk   = Fastor::Index< K_, k_ >;
    using L    = Fastor::Index< L_ >;
    using LI   = Fastor::Index< L_, I_ >;
    using LK   = Fastor::Index< L_, K_ >;
    using LN   = Fastor::Index< L_, N_ >;
    using Lm   = Fastor::Index< L_, m_ >;
    using LmN  = Fastor::Index< L_, m_, N_ >;
    using MJKL = Fastor::Index< M_, J_, K_, L_ >;
    using MK   = Fastor::Index< M_, K_ >;
    using ML   = Fastor::Index< M_, L_ >;
    using MNL  = Fastor::Index< M_, N_, L_ >;
    using MPm  = Fastor::Index< M_, P_, m_ >;
    using Mi   = Fastor::Index< M_, i_ >;
    using NLJl = Fastor::Index< N_, L_, J_, l_ >;
    using Nm   = Fastor::Index< N_, m_ >;
    using Pm   = Fastor::Index< P_, m_ >;
    using i    = Fastor::Index< i_ >;
    using iA   = Fastor::Index< i_, A_ >;
    using iAkB = Fastor::Index< i_, A_, k_, B_ >;
    using iB   = Fastor::Index< i_, B_ >;
    using iI   = Fastor::Index< i_, I_ >;
    using iIKL = Fastor::Index< i_, I_, K_, L_ >;
    using iIjJ = Fastor::Index< i_, I_, j_, J_ >;
    using iIkK = Fastor::Index< i_, I_, k_, K_ >;
    using iIkL = Fastor::Index< i_, I_, k_, L_ >;
    using iImn = Fastor::Index< i_, I_, m_, n_ >;
    using iJ   = Fastor::Index< i_, J_ >;
    using iJKL = Fastor::Index< i_, J_, K_, L_ >;
    using iJLl = Fastor::Index< i_, J_, L_, l_ >;
    using iJkL = Fastor::Index< i_, J_, k_, L_ >;
    using iJl  = Fastor::Index< i_, J_, l_ >;
    using iK   = Fastor::Index< i_, K_ >;
    using iKjL = Fastor::Index< i_, K_, j_, L_ >;
    using iL   = Fastor::Index< i_, L_ >;
    using iM   = Fastor::Index< i_, M_ >;
    using iN   = Fastor::Index< i_, N_ >;
    using iNL  = Fastor::Index< i_, N_, L_ >;
    using ij   = Fastor::Index< i_, j_ >;
    using ijB  = Fastor::Index< i_, j_, B_ >;
    using ijKJ = Fastor::Index< i_, j_, K_, J_ >;
    using ijKL = Fastor::Index< i_, j_, K_, L_ >;
    using ijL  = Fastor::Index< i_, j_, L_ >;
    using ijLm = Fastor::Index< i_, j_, L_, m_ >;
    using ijk  = Fastor::Index< i_, j_, k_ >;
    using ijkB = Fastor::Index< i_, j_, k_, B_ >;
    using ijkK = Fastor::Index< i_, j_, k_, K_ >;
    using ijkl = Fastor::Index< i_, j_, k_, l_ >;
    using ijl  = Fastor::Index< i_, j_, l_ >;
    using ijm  = Fastor::Index< i_, j_, m_ >;
    using ijmn = Fastor::Index< i_, j_, m_, n_ >;
    using ijnk = Fastor::Index< i_, j_, n_, k_ >;
    using ijnm = Fastor::Index< i_, j_, n_, m_ >;
    using ik   = Fastor::Index< i_, k_ >;
    using im   = Fastor::Index< i_, m_ >;
    using imk  = Fastor::Index< i_, m_, k_ >;
    using imkl = Fastor::Index< i_, m_, k_, l_ >;
    using imL  = Fastor::Index< i_, m_, L_ >;
    using imLk = Fastor::Index< i_, m_, L_, k_ >;
    using in   = Fastor::Index< i_, n_ >;
    using inB  = Fastor::Index< i_, n_, B_ >;
    using inkB = Fastor::Index< i_, n_, k_, B_ >;
    using j    = Fastor::Index< j_ >;
    using jA   = Fastor::Index< j_, A_ >;
    using jB   = Fastor::Index< j_, B_ >;
    using jJ   = Fastor::Index< j_, J_ >;
    using jL   = Fastor::Index< j_, L_ >;
    using jLm  = Fastor::Index< j_, L_, m_ >;
    using ji   = Fastor::Index< j_, i_ >;
    using jin  = Fastor::Index< j_, i_, n_ >;
    using jk   = Fastor::Index< j_, k_ >;
    using jK   = Fastor::Index< j_, K_ >;
    using jkB  = Fastor::Index< j_, k_, B_ >;
    using jkl  = Fastor::Index< j_, k_, l_ >;
    using jl   = Fastor::Index< j_, l_ >;
    using k    = Fastor::Index< k_ >;
    using kA   = Fastor::Index< k_, A_ >;
    using kB   = Fastor::Index< k_, B_ >;
    using kI   = Fastor::Index< k_, I_ >;
    using kK   = Fastor::Index< k_, K_ >;
    using kL   = Fastor::Index< k_, L_ >;
    using kM   = Fastor::Index< k_, M_ >;
    using kNL  = Fastor::Index< k_, N_, L_ >;
    using kj   = Fastor::Index< k_, j_ >;
    using kJ   = Fastor::Index< k_, J_ >;
    using kl   = Fastor::Index< k_, l_ >;
    using km   = Fastor::Index< k_, m_ >;
    using l    = Fastor::Index< l_ >;
    using lB   = Fastor::Index< l_, B_ >;
    using lm   = Fastor::Index< l_, m_ >;
    using m    = Fastor::Index< m_ >;
    using mK   = Fastor::Index< m_, K_ >;
    using mLl  = Fastor::Index< m_, L_, l_ >;
    using mj   = Fastor::Index< m_, j_ >;
    using mjL  = Fastor::Index< m_, j_, L_ >;
    using mn   = Fastor::Index< m_, n_ >;
    using mnKL = Fastor::Index< m_, n_, K_, L_ >;
    using mnij = Fastor::Index< m_, n_, i_, j_ >;
    using mnkB = Fastor::Index< m_, n_, k_, B_ >;
    using mnkL = Fastor::Index< m_, n_, k_, L_ >;
    using nB   = Fastor::Index< n_, B_ >;

    using to_IJKL = Fastor::OIndex< I_, J_, K_, L_ >;
    using to_IJkK = Fastor::OIndex< I_, J_, k_, K_ >;
    using to_IJkL = Fastor::OIndex< I_, J_, k_, L_ >;
    using to_IikK = Fastor::OIndex< I_, i_, k_, K_ >;
    using to_IjkK = Fastor::OIndex< I_, j_, k_, K_ >;
    using to_Ii   = Fastor::OIndex< I_, i_ >;
    using to_NLJl = Fastor::OIndex< N_, L_, J_, l_ >;
    using to_iIKL = Fastor::OIndex< i_, I_, K_, L_ >;
    using to_iIjJ = Fastor::OIndex< i_, I_, j_, J_ >;
    using to_iImn = Fastor::OIndex< i_, I_, m_, n_ >;
    using to_ij   = Fastor::OIndex< i_, j_ >;
    using to_ijIJ = Fastor::OIndex< i_, j_, I_, J_ >;
    using to_ijKL = Fastor::OIndex< i_, j_, K_, L_ >;
    using to_ijL  = Fastor::OIndex< i_, j_, L_ >;
    using to_ijLk = Fastor::OIndex< i_, j_, L_, k_ >;
    using to_ijLm = Fastor::OIndex< i_, j_, L_, m_ >;
    using to_ijk  = Fastor::OIndex< i_, j_, k_ >;
    using to_ijkK = Fastor::OIndex< i_, j_, k_, K_ >;
    using to_ijkL = Fastor::OIndex< i_, j_, k_, L_ >;
    using to_ijKl = Fastor::OIndex< i_, j_, K_, l_ >;
    using to_ijkl = Fastor::OIndex< i_, j_, k_, l_ >;
    using to_ijm  = Fastor::OIndex< i_, j_, m_ >;
    using to_ijmM = Fastor::OIndex< i_, j_, m_, M_ >;
    using to_jAB  = Fastor::OIndex< j_, A_, B_ >;
    using to_jAkB = Fastor::OIndex< j_, A_, k_, B_ >;
    using to_ji   = Fastor::OIndex< j_, i_ >;
    using to_jikL = Fastor::OIndex< j_, i_, k_, L_ >;
    using to_jikl = Fastor::OIndex< j_, i_, k_, l_ >;
    using to_jkiB = Fastor::OIndex< j_, k_, i_, B_ >;
    using to_kK   = Fastor::OIndex< k_, K_ >;
    using to_kL   = Fastor::OIndex< k_, L_ >;
  } // namespace FastorIndices

  template < typename T,
             size_t nRows,
             size_t nCols,
             typename = std::enable_if< !std::is_const< std::remove_reference< T > >::value > >
  auto inline mapEigenToFastor( const Fastor::Tensor< T, nRows, nCols >& fastor )
  {
    return Eigen::Map< const Eigen::Matrix< T, nRows, nCols, Eigen::RowMajor > >( fastor.data() );
  }

  template < typename T, size_t nRows, size_t nCols, typename = void >
  auto inline mapEigenToFastor( Fastor::Tensor< T, nRows, nCols >& fastor )
  {
    return Eigen::Map< Eigen::Matrix< T, nRows, nCols, Eigen::RowMajor > >( fastor.data() );
  }

  template < typename T, size_t nRows, typename = void >
  auto inline mapEigenToFastor( const Fastor::Tensor< T, nRows >& fastor )
  {
    return Eigen::Map< Eigen::Matrix< T, nRows, 1 > >( fastor.data() );
  }

  template < typename T, size_t nRows, typename = void >
  auto inline mapEigenToFastor( const Fastor::TensorMap< T, nRows >& fastor )
  {
    return Eigen::Map< Eigen::Matrix< T, nRows, 1 > >( fastor.data() );
  }

  template < typename T, size_t nRows, size_t nCols, typename = void >
  auto inline mapEigenToFastor( const Fastor::TensorMap< T, nRows, nCols >& fastor )
  {
    return Eigen::Map< Eigen::Matrix< T, nRows, nCols, Eigen::RowMajor > >( fastor.data() );
  }

  template < template < typename, size_t... > class TensorType, typename T, size_t... Rest >
  void inline copyFastorToColumnMajor( T* target, const TensorType< T, Rest... >& source )
  {
    ( Fastor::TensorMap< T, Rest... >( target ) ) = Fastor::torowmajor( source );
  }

  enum DimensionType { U, W };
  template < DimensionType... dims, typename T, size_t... dims3D >
  auto inline reduceTo2D( const Fastor::Tensor< T, dims3D... >& theTensor3D )
  {

    static_assert( sizeof...( dims ) == sizeof...( dims3D ),
                   "CONVERSION PACK (1) AND TENSOR-ORDER PACK (2) MUST HAVE SAME LENGTH" );

    return theTensor3D( Fastor::fseq < dims == U ? 0 : 2, dims == U ? 2 : 3 > ()... );
  }

  template < typename Derived, size_t order >
  auto inline reduceTo2D( const Fastor::AbstractTensor< Derived, order >& theTensor3D )
  {
    using result_type = typename Derived::result_type;
    return reduceTo2D( result_type( theTensor3D ) );
  }

  constexpr int const3( size_t x )
  {
    return 3;
  }

  template < typename T, size_t... dims2D >
  auto inline expandTo3D( const Fastor::Tensor< T, dims2D... >& theTensor2D )
  {
    static_assert( ( ( ( dims2D > 0 ) && ( dims2D < 3 ) ) && ... ),
                   "INPUT TENSOR IS NOT A VALID 2D MIXED DISPLACEMENT/ROTATION TENSOR" );

    Fastor::Tensor< double, const3( dims2D )... > theTensor3D( 0.0 );

    theTensor3D( Fastor::fseq < dims2D == 2 ? 0 : 2, dims2D == 2 ? 2 : 3 > ()... ) = theTensor2D;

    return theTensor3D;
  }

  template < typename Derived, size_t order >
  auto inline expandTo3D( const Fastor::AbstractTensor< Derived, order >& theTensor2D )
  {
    using result_type = typename Derived::result_type;
    return expandTo3D( result_type( theTensor2D ) );
  }

  template < typename T >
  FastorStandardTensors::Tensor33t< T > multiplyFastorTensor33WithScalar( FastorStandardTensors::Tensor33t< T > tensor,
                                                                          T                                     scalar )
  {
    // workaround for fastor bug (issue #149)
    FastorStandardTensors::Tensor33t< T > res;

    for ( int i = 0; i < 3; i++ ) {
      for ( int j = 0; j < 3; j++ ) {
        res( i, j ) = tensor( i, j ) * scalar;
      }
    }
    return res;
  }

  template < typename T, size_t... Rest >
  Fastor::Tensor< T, Rest... > multiplyFastorTensorWithScalar( Fastor::Tensor< T, Rest... > tensor, T scalar )
  {
    // workaround for fastor bug (issue #149)
    Fastor::Tensor< T, Rest... > out;
    T*                           out_data = out.data();
    T*                           in_data  = tensor.data();

    for ( Fastor::FASTOR_INDEX i = 0; i < tensor.size(); ++i ) {
      out_data[out.get_mem_index( i )] = in_data[tensor.get_mem_index( i )] * scalar;
    }
    return out;
  }

  template < typename T >
  T einsum_ij_ij_hardcoded( const FastorStandardTensors::Tensor33t< T >& A,
                            const FastorStandardTensors::Tensor33t< T >& B )
  {
    T result( 0.0 );

    for ( int i = 0; i < 3; i++ ) {
      for ( int j = 0; j < 3; j++ ) {
        result += A( i, j ) * B( i, j );
      }
    }

    return result;
  }

  template < typename T, size_t... Rest >
  Fastor::Tensor< T, Rest... > fastorTensorFromDoubleTensor( const Fastor::Tensor< double, Rest... >& in )
  {
    // workaround for lack of casting of pointer types
    // e.g. no cast available for autoduff::dual* and double*

    Fastor::Tensor< T, Rest... > out;
    T*                           out_data = out.data();
    double*                      in_data  = in.data();

    for ( Fastor::FASTOR_INDEX i = 0; i < in.size(); ++i ) {
      out_data[out.get_mem_index( i )] = static_cast< T >( in_data[in.get_mem_index( i )] );
    }
    return out;
  }

  template < typename T, size_t... Rest >
  Fastor::Tensor< T, Rest... > fastorTensorFromDoubleTensorMap( const Fastor::TensorMap< double, Rest... >& in )
  {
    // workaround for lack of casting of pointer types
    // e.g. no cast available for autoduff::dual* and double*

    Fastor::Tensor< T, Rest... > out;
    T*                           out_data = out.data();
    double*                      in_data  = in.data();

    for ( Fastor::FASTOR_INDEX i = 0; i < in.size(); ++i ) {
      out_data[out.get_mem_index( i )] = static_cast< T >( in_data[in.get_mem_index( i )] );
    }
    return out;
  }

  template < typename T, size_t dim = 3 >
  Fastor::Tensor< T, dim, dim > secondRankTensorFromSecondRankDoubleTensor(
    const Fastor::Tensor< double, dim, dim >& in )
  {
    return fastorTensorFromDoubleTensor< T >( in );
  }

  template < typename T, size_t... Rest >
  Fastor::Tensor< double, Rest... > makeReal( const Fastor::Tensor< T, Rest... >& in )
  {

    Fastor::Tensor< double, Rest... > out;
    double*                           out_data = out.data();
    T*                                in_data  = in.data();

    for ( Fastor::FASTOR_INDEX i = 0; i < in.size(); ++i ) {
      out_data[out.get_mem_index( i )] = static_cast< double >( in_data[in.get_mem_index( i )] );
    }
    return out;
  }

  template < typename T, size_t... Rest >
  Fastor::Tensor< autodiff::dual, Rest... > makeDual( const Fastor::Tensor< T, Rest... >& in )
  {

    Fastor::Tensor< autodiff::dual, Rest... > out;
    autodiff::dual*                           out_data = out.data();
    T*                                        in_data  = in.data();

    for ( Fastor::FASTOR_INDEX i = 0; i < in.size(); ++i ) {
      out_data[out.get_mem_index( i )] = autodiff::dual( in_data[in.get_mem_index( i )] );
    }
    return out;
  }

  template < size_t order, size_t... Rest >
  Fastor::Tensor< autodiff::HigherOrderDual< order, double >, Rest... > makeHigherOrderDual(
    const Fastor::Tensor< double, Rest... >& in )
  {

    Fastor::Tensor< autodiff::HigherOrderDual< order, double >, Rest... > out;
    autodiff::HigherOrderDual< order, double >*                           out_data = out.data();
    double*                                                               in_data  = in.data();

    for ( Fastor::FASTOR_INDEX i = 0; i < in.size(); ++i ) {
      out_data[out.get_mem_index( i )] = autodiff::HigherOrderDual< order, double >( in_data[in.get_mem_index( i )] );
    }
    return out;
  }

} // namespace Marmot
