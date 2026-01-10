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
#include "Marmot/MarmotFastorTensorBasics.h"
#include "autodiff/forward/dual.hpp"
#include "autodiff/forward/dual/eigen.hpp"
#include <autodiff/forward/dual/dual.hpp>
#include <functional>

namespace Marmot {

  namespace AutomaticDifferentiation {

    using namespace autodiff;

    /** @brief Increase the order of an arbitray order dual-valued Fastor tensor by one and shift the derivatives
     *  @tparam order current order of the dual numers in the tensor
     *  @tparam Rest dimensions of the tensor
     *  @param in input tensor
     *  @return output tensor with increased order dual numbers
     */
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

    /** @typedef tensor_to_scalar_function_type
     *  @brief Alias for a function mapping a tensor to a scalar with dual numbers
     *  @tparam Rest dimensions of the tensor
     * */
    template < size_t... Rest >
    using tensor_to_scalar_function_type = std::function< dual( const Fastor::Tensor< dual, Rest... >& T ) >;

    /** @brief Compute the gradient of a tensor-to-scalar function
     *  @tparam Rest dimensions of the tensor
     *  @param f function mapping a tensor to a scalar
     *  @param T input tensor at which the gradient is evaluated
     *  @return gradient of f with respect to T, same shape as T
     */
    template < size_t... Rest >
    Fastor::Tensor< double, Rest... > df_dT( const tensor_to_scalar_function_type< Rest... >& f,
                                             const Fastor::Tensor< double, Rest... >&         T )
    {
      Fastor::Tensor< double, Rest... > df_dT( 0.0 );
      Fastor::Tensor< dual, Rest... >   T_right = makeDual( T );

      double* df_dT_data   = df_dT.data();
      dual*   T_right_data = T_right.data();

      for ( Fastor::FASTOR_INDEX i = 0; i < T.size(); ++i ) {
        const int T_right_mem_idx = T_right.get_mem_index( i );
        seed< 1 >( T_right_data[T_right_mem_idx], 1.0 );
        df_dT_data[df_dT.get_mem_index( i )] = derivative< 1 >( f( T_right ) );
        seed< 1 >( T_right_data[T_right_mem_idx], 0.0 );
      }

      return df_dT;
    }

    /** @typedef tensor_to_scalar_function_type_arbitrary_dual_order
     *  @brief A type alias for a function mapping a tensor to a scalar with dual numbers of arbitrary order
     *  @tparam order current order of the dual numbers in the tensor
     *  @tparam Rest dimensions of the tensor
     * */
    template < size_t order, size_t... Rest >
    using tensor_to_scalar_function_type_arbitrary_dual_order = std::function< HigherOrderDual< order, double >(
      const Fastor::Tensor< HigherOrderDual< order, double >, Rest... >& T ) >;

    /** @brief Compute the gradient of a tensor-to-scalar function where the with dual numbers of arbitrary order
     *  @tparam order current order of the dual numers in the tensor
     *  @tparam Rest dimensions of the tensor
     *  @param f function mapping a tensor to a scalar
     *  @param T input tensor at which the gradient is evaluated
     *  @return pair of function value and gradient of f with respect to T, same shape as T
     */
    template < size_t order, size_t... Rest >
    std::pair< HigherOrderDual< order, double >, Fastor::Tensor< HigherOrderDual< order, double >, Rest... > > df_dT(
      const tensor_to_scalar_function_type_arbitrary_dual_order< order + 1, Rest... >& f,
      const Fastor::Tensor< HigherOrderDual< order, double >, Rest... >&               T )
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
        f_                                    = f( T_right );
        df_dT_data[df_dT_.get_mem_index( i )] = decreaseDualOrderWithShift< order + 1 >( f_ );
        seed< 1 >( T_right_data[T_right.get_mem_index( i )], 0.0 );
      }
      f_ = f( T_right );
      return { decreaseDualOrder< order + 1 >( f_ ), df_dT_ };
    }

    /** @brief Compute the gradient of a tensor-to-tensor function
     *  @tparam RestF dimensions of the output tensor
     *  @tparam RestT dimensions of the input tensor
     *  @param F function mapping a tensor to a tensor
     *  @param T input tensor at which the gradient is evaluated
     *  @return pair of function value and gradient of F with respect to T, gradient has shape (RestF..., RestT...)
     */
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

      /** @typedef tensor_to_scalar_function_type
       *  @brief Alias for a function mapping a tensor to a scalar with second order dual numbers
       *  @tparam dim dimension of the input tensor (dim x dim)
       *
       *  @note The implementation is currently limited to second rank (dim x dim) tensors
       */
      template < size_t dim >
      using tensor_to_scalar_function_type = std::function< dual2nd( const Fastor::Tensor< dual2nd, dim, dim >& T ) >;

      /** @brief Computes the second derivative of a function that maps a tensor to a scalar
       *  @tparam dim dimension of the input tensor (dim x dim)
       *  @param F function mapping a (dim x dim) tensor to a scalar
       *  @param T input tensor at which the derivative is evaluated
       *  @return tuple of function value, first derivative and second derivative of F with respect to T
       *          first derivative has shape (dim, dim)
       *          second derivative has shape (dim, dim, dim, dim)
       *
       *  @note The implementation is currently limited to second rank (dim x dim) tensors
       */
      template < size_t dim >
      std::tuple< double, Fastor::Tensor< double, dim, dim >, Fastor::Tensor< double, dim, dim, dim, dim > > d2f_dT2(
        const tensor_to_scalar_function_type< dim >& F,
        const Fastor::Tensor< double, dim, dim >&    T )
      {
        double                                       F_;
        dual2nd                                      F_right;
        Fastor::Tensor< double, dim, dim >           dF_dT_;
        Fastor::Tensor< double, dim, dim, dim, dim > d2F_dT2;
        Fastor::Tensor< dual2nd, dim, dim >          T_right = makeHigherOrderDual< 2 >( T );

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

      /** @typedef tensor_and_scalar_to_scalar_function_type
       *  @brief Alias for a function mapping a tensor and a scalar to a scalar with second order dual numbers
       *  @tparam dim dimension of the input tensor (dim x dim)
       *
       *  @note The implementation is currently limited to second rank (dim x dim) tensors
       */
      template < size_t dim >
      using tensor_and_scalar_to_scalar_function_type = std::function<
        dual2nd( const Fastor::Tensor< dual2nd, dim, dim >& T, const dual2nd scalar ) >;

      /** @brief Computes the mixed second derivative of a function that maps a tensor and a scalar to a scalar
       *  @tparam dim dimension of the input tensor (dim x dim)
       *  @param F function mapping a (dim x dim) tensor and a scalar to a scalar
       *  @param T input tensor at which the derivative is evaluated
       *  @param scalar input scalar at which the derivative is evaluated
       *  @return mixed second derivative of F with respect to T and the scalar, has shape (dim, dim)
       *
       *  @note The implementation is currently limited to second rank (dim x dim) tensors
       */
      template < size_t dim >
      Fastor::Tensor< double, dim, dim > d2f_dTensor_dScalar( const tensor_and_scalar_to_scalar_function_type< dim >& F,
                                                              const Fastor::Tensor< double, dim, dim >&               T,
                                                              const double scalar )
      {
        Fastor::Tensor< double, dim, dim >  d2F_dTdScalar;
        Fastor::Tensor< dual2nd, dim, dim > T_right = makeHigherOrderDual< 2 >( T );

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
    namespace ThirdOrder {

      /** @typedef tensor_to_scalar_function_type
       *  @brief Alias for a function mapping a tensor to a scalar with third order dual numbers
       *  @tparam dim dimension of the input tensor (dim x dim)
       *
       *  @note The implementation is currently limited to second rank (dim x dim) tensors
       */
      template < size_t dim >
      using tensor_to_scalar_function_type = std::function< dual3rd( const Fastor::Tensor< dual3rd, dim, dim >& T ) >;

      /** @brief Computes the third derivative of a function that maps a tensor to a scalar
       *  @tparam dim dimension of the input tensor (dim x dim)
       *  @param F function mapping a (dim x dim) tensor to a scalar
       *  @param T input tensor at which the derivative is evaluated
       *  @return tuple of function value, first derivative, second derivative and third derivative of F with respect to T
       *          first derivative has shape (dim, dim)
       *          second derivative has shape (dim, dim, dim, dim)
       *          third derivative has shape (dim, dim, dim, dim, dim, dim)
       *
       *  @note The implementation is currently limited to second rank (dim x dim) tensors
       */
      template < size_t dim >
      std::tuple< double,
                  Fastor::Tensor< double, dim, dim >,
                  Fastor::Tensor< double, dim, dim, dim, dim >,
                  Fastor::Tensor< double, dim, dim, dim, dim, dim, dim > >
      d3f_dT3( const tensor_to_scalar_function_type< dim >& F, const Fastor::Tensor< double, dim, dim >& T )
      {
        double                                                 F_;
        dual3rd                                                F_right;
        Fastor::Tensor< double, dim, dim >                     dF_dT_;
        Fastor::Tensor< double, dim, dim, dim, dim >           d2F_dT2;
        Fastor::Tensor< double, dim, dim, dim, dim, dim, dim > d3F_dT3;
        Fastor::Tensor< dual3rd, dim, dim >                    T_right = makeHigherOrderDual< 3 >( T );

        for ( size_t i = 0; i < dim; i++ ) {
          for ( size_t j = 0; j < dim; j++ ) {
            seed< 1 >( T_right( i, j ), 1.0 );

            for ( size_t k = 0; k < dim; k++ ) {
              for ( size_t l = 0; l < dim; l++ ) {

                seed< 2 >( T_right( k, l ), 1.0 );

                for ( size_t m = 0; m < dim; m++ ) {
                  for ( size_t n = 0; n < dim; n++ ) {

                    seed< 3 >( T_right( m, n ), 1.0 );
                    F_right                     = F( T_right );
                    d3F_dT3( i, j, k, l, m, n ) = derivative< 3 >( F_right );

                    seed< 3 >( T_right( m, n ), 0.0 );
                  }
                }

                d2F_dT2( i, j, k, l ) = derivative< 2 >( F_right );

                seed< 2 >( T_right( k, l ), 0.0 );
              }
            }
            dF_dT_( i, j ) = derivative< 1 >( F_right );
            F_             = double( F_right );
            seed< 1 >( T_right( i, j ), 0.0 );
          }
        }

        return { F_, dF_dT_, d2F_dT2, d3F_dT3 };
      }

    } // namespace ThirdOrder

  }   // namespace AutomaticDifferentiation

} // namespace Marmot
