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
#include "Fastor/Fastor.h"
#include "Marmot/MarmotAutomaticDifferentiation.h"
#include "Marmot/MarmotMicromorphicTensorBasics.h"
#include "Marmot/MarmotTensor.h"
#include "autodiff/forward/dual.hpp"
#include "autodiff/forward/dual/eigen.hpp"
#include <autodiff/forward/dual/dual.hpp>
#include <functional>

namespace Marmot {

  namespace AutomaticDifferentiation {

    using namespace autodiff;
    using namespace Eigen;

    template < size_t dim >
    Fastor::Tensor< dual2nd, dim, dim > secondRankDual2ndFromDualTensorWithShift(
      const Fastor::Tensor< dual, dim, dim >& in )
    {
      Fastor::Tensor< dual2nd, dim, dim > out( 0.0 );
      for ( size_t i = 0; i < dim; i++ ) {
        for ( size_t j = 0; j < dim; j++ ) {
          // out( i, j ) = shiftTo2ndOrderDual( in( i, j ) );
          out( i, j ).val = in( i, j );
          seed< 2 >( out( i, j ), (double)in( i, j ).grad );
        }
      }
      return out;
    }

    template < size_t order, size_t... Rest >
    Fastor::Tensor< HigherOrderDual< order + 1, double >, Rest... > increaseDualOrderWithShift(
      const Fastor::Tensor< HigherOrderDual< order, double >, Rest... >& in )
    {
      using in_scalar_type  = HigherOrderDual< order, double >;
      using out_scalar_type = HigherOrderDual< order + 1, double >;

      Fastor::Tensor< out_scalar_type, Rest... > out( 0.0 );
      out_scalar_type*                           out_data = out.data();
      in_scalar_type*                            in_data  = in.data();

      for ( Fastor::FASTOR_INDEX i = 0; i < in.size(); ++i ) {
        out_data[out.get_mem_index( i )] = increaseDualOrderWithShift< order >( in_data[in.get_mem_index( i )] );
      }
      return out;
    }

    template < size_t dim >
    using tensor_to_scalar_function_type = std::function< dual( const Fastor::Tensor< dual, dim, dim >& T ) >;

    template < size_t dim >
    Fastor::Tensor< double, dim, dim > df_dTensor( const tensor_to_scalar_function_type< dim >& f,
                                                   const Fastor::Tensor< double, dim, dim >&    tensor )
    {
      Fastor::Tensor< double, dim, dim > df_dT;
      Fastor::Tensor< dual, dim, dim >   tensorDual = secondRankTensorFromSecondRankDoubleTensor< dual, dim >( tensor );

      for ( size_t i = 0; i < dim; i++ ) {
        for ( size_t j = 0; j < dim; j++ ) {

          seed< 1 >( tensorDual( i, j ), 1.0 );

          df_dT( i, j ) = f( tensorDual ).grad;

          seed< 1 >( tensorDual( i, j ), 0.0 );
        }
      }

      return df_dT;
    }

    template < size_t dim >
    using tensor_to_scalar_function_type_2nd = std::function< dual2nd( const Fastor::Tensor< dual2nd, dim, dim >& T ) >;

    template < size_t dim >
    Fastor::Tensor< dual, dim, dim > df_dTensor( const tensor_to_scalar_function_type_2nd< dim >& f,
                                                 const Fastor::Tensor< dual, dim, dim >&          T )
    {
      Fastor::Tensor< dual, dim, dim > df_dT_( 0.0 );

      // shift gradient part to 2nd order part
      Fastor::Tensor< dual2nd, dim, dim > T_right = secondRankDual2ndFromDualTensorWithShift< dim >( T );

      for ( size_t i = 0; i < dim; i++ ) {
        for ( size_t j = 0; j < dim; j++ ) {

          seed< 1 >( T_right( i, j ), 1.0 );
          const dual2nd f_right = f( T_right );
          df_dT_( i, j ).val    = derivative< 1 >( f_right );
          df_dT_( i, j ).grad   = derivative< 2 >( f_right );

          seed< 1 >( T_right( i, j ), 0.0 );
        }
      }

      return df_dT_;
    }

    // also return f
    template < size_t dim >
    std::tuple< dual, Fastor::Tensor< dual, dim, dim > > df_dTensor_(
      const tensor_to_scalar_function_type_2nd< dim >& f,
      const Fastor::Tensor< dual, dim, dim >&          T )
    {
      dual                             f_;
      Fastor::Tensor< dual, dim, dim > df_dT_( 0.0 );

      // shift gradient part to 2nd order part
      Fastor::Tensor< dual2nd, dim, dim > T_right = secondRankDual2ndFromDualTensorWithShift< dim >( T );

      for ( size_t i = 0; i < dim; i++ ) {
        for ( size_t j = 0; j < dim; j++ ) {

          seed< 1 >( T_right( i, j ), 1.0 );
          const dual2nd f_right = f( T_right );
          df_dT_( i, j ).val    = derivative< 1 >( f_right );
          df_dT_( i, j ).grad   = derivative< 2 >( f_right );

          seed< 1 >( T_right( i, j ), 0.0 );
          f_.val  = double( f_right );
          f_.grad = f_right.val.grad;
        }
      }
      return { f_, df_dT_ };
    }

    template < size_t order, size_t... Rest >
    using tensor_to_scalar_function_type_arbitrary_dual_order = std::function< HigherOrderDual< order, double >(
      const Fastor::Tensor< HigherOrderDual< order, double >, Rest... >& T ) >;

    template < size_t order, size_t... Rest >
    Fastor::Tensor< HigherOrderDual< order, double >, Rest... > df_dT(
      tensor_to_scalar_function_type_arbitrary_dual_order< order + 1, Rest... >& f,
      const Fastor::Tensor< HigherOrderDual< order, double >, Rest... >&         T )
    {

      using scalartype            = HigherOrderDual< order, double >;
      using higherOrderScalartype = HigherOrderDual< order + 1, double >;

      higherOrderScalartype                            f_;
      Fastor::Tensor< scalartype, Rest... >            df_dT_( 0.0 );
      Fastor::Tensor< higherOrderScalartype, Rest... > T_right = increaseDualOrderWithShift< order >( T );

      higherOrderScalartype* T_right_data = T_right.data();
      scalartype*            df_dT_data   = df_dT_.data();

      for ( Fastor::FASTOR_INDEX i = 0; i < T.size(); ++i ) {
        seed< 1 >( T_right_data[T_right.get_mem_index( i )], 1.0 );
        f_ = f( T_right );
        /* std::cout <<  f_.val << f_.grad << std::endl; */

        df_dT_data[df_dT_.get_mem_index( i )] = decreaseDualOrderWithShift< order + 1 >( f_ );
        seed< 1 >( T_right_data[T_right.get_mem_index( i )], 0.0 );
      }

      return df_dT_;
    }

    template < size_t... RestF, size_t... RestT >
    std::pair< Fastor::Tensor< double, RestF... >, Fastor::Tensor< double, RestF..., RestT... > > dF_dT(
      std::function< Fastor::Tensor< dual, RestF... >( const Fastor::Tensor< dual, RestT... >& ) >& F,
      const Fastor::Tensor< double, RestT... >&                                                     T )
    {

      Fastor::Tensor< double, RestF... >           F_( 0.0 );
      Fastor::Tensor< double, RestF..., RestT... > dF_dT_( 0.0 );
      Fastor::Tensor< dual, RestT... >             T_right = makeDual( T );

      Fastor::Tensor< dual, RestF... > F_at_T_right( 0.0 );

      double* dF_dT_data        = dF_dT_.data();
      dual*   T_right_data      = T_right.data();
      dual*   F_at_T_right_data = F_at_T_right.data();

      for ( Fastor::FASTOR_INDEX i = 0; i < T.size(); ++i ) {
        const int T_right_mem_idx = T_right.get_mem_index( i );
        T_right_data[T_right_mem_idx].grad += 1.0;
        F_at_T_right = F( T_right );

        for ( Fastor::FASTOR_INDEX j = 0; j < F_at_T_right.size(); ++j ) {
          dF_dT_data[dF_dT_.get_mem_index( j * T.size() + i )] = F_at_T_right_data[F_at_T_right.get_mem_index( j )]
                                                                   .grad;
        }
        T_right_data[T_right_mem_idx].grad -= 1.0;
      }

      F_ = makeReal( F_at_T_right );

      return { F_, dF_dT_ };
    }

    namespace SecondOrder {

      template < size_t dim >
      using tensor_to_scalar_function_type = std::function< dual2nd( const Fastor::Tensor< dual2nd, dim, dim >& T ) >;

      template < size_t dim >
      std::tuple< double, Fastor::Tensor< double, dim, dim >, Fastor::Tensor< double, dim, dim, dim, dim > > d2f_dTensor_dTensor(
        const tensor_to_scalar_function_type< dim >& F,
        const Fastor::Tensor< double, dim, dim >&    T )
      {
        double                                       F_;
        dual2nd                                      F_right;
        Fastor::Tensor< double, dim, dim >           dF_dT_;
        Fastor::Tensor< double, dim, dim, dim, dim > d2F_dT2;
        Fastor::Tensor< dual2nd, dim, dim > T_right = secondRankTensorFromSecondRankDoubleTensor< dual2nd, dim >( T );

        for ( size_t i = 0; i < dim; i++ ) {
          for ( size_t j = 0; j < dim; j++ ) {
            seed< 1 >( T_right( i, j ), 1.0 );

            for ( size_t k = 0; k < dim; k++ ) {
              for ( size_t l = 0; l < dim; l++ ) {

                seed< 2 >( T_right( k, l ), 1.0 );
                F_right               = F( T_right );
                d2F_dT2( i, j, k, l ) = derivative< 2 >( F_right );

                seed< 2 >( T_right( k, l ), 0.0 );
              }
            }
            dF_dT_( i, j ) = derivative< 1 >( F_right );
            F_             = double( F_right );
            seed< 1 >( T_right( i, j ), 0.0 );
          }
        }

        return { F_, dF_dT_, d2F_dT2 };
      }

      template < size_t dim >
      using tensor_and_scalar_to_scalar_function_type = std::function<
        dual2nd( const Fastor::Tensor< dual2nd, dim, dim >& T, const dual2nd scalar ) >;

      template < size_t dim >
      Fastor::Tensor< double, dim, dim > d2f_dTensor_dScalar( const tensor_and_scalar_to_scalar_function_type< dim >& F,
                                                              const Fastor::Tensor< double, dim, dim >&               T,
                                                              const double scalar )
      {
        Fastor::Tensor< double, dim, dim >  d2F_dTdScalar;
        Fastor::Tensor< dual2nd, dim, dim > T_right = secondRankTensorFromSecondRankDoubleTensor< dual2nd, dim >( T );

        dual2nd scalar_right( scalar );
        seed< 2 >( scalar_right, 1.0 );

        for ( size_t i = 0; i < dim; i++ ) {
          for ( size_t j = 0; j < dim; j++ ) {

            seed< 1 >( T_right( i, j ), 1.0 );

            d2F_dTdScalar( i, j ) = derivative< 2 >( F( T_right, scalar_right ) );

            seed< 1 >( T_right( i, j ), 0.0 );
          }
        }

        return d2F_dTdScalar;
      }
    } // namespace SecondOrder

  } // namespace AutomaticDifferentiation

} // namespace Marmot
