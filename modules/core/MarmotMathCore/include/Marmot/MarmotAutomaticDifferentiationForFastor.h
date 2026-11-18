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
#include <algorithm>
/* Fastor/Fastor.h ends with `#pragma warning( default : ... )`; scope it so it neither overrides the build's warning
 * flags nor leaks into consumers. */
#ifdef _MSC_VER
#  pragma warning( push )
#endif
#include "Fastor/Fastor.h"
#ifdef _MSC_VER
#  pragma warning( pop )
#endif
#include "Marmot/MarmotAutomaticDifferentiation.h"
#include "Marmot/MarmotAutomaticDifferentiationTypes.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include <functional>

namespace Marmot::Math::Differentiation::Automatic {

  /** @brief Increase the order of an arbitray order dual-valued Fastor tensor by one and shift the derivatives
   *  @tparam order current order of the dual numers in the tensor
   *  @tparam Rest dimensions of the tensor
   *  @param in input tensor
   *  @return output tensor with increased order dual numbers
   */
  template < size_t order, size_t... Rest >
  Fastor::Tensor< HyperDual< order + 1 >, Rest... > increaseDualOrderWithShift(
    const Fastor::Tensor< HyperDual< order >, Rest... >& in )
  {
    using in_scalar_type  = HyperDual< order >;
    using out_scalar_type = HyperDual< order + 1 >;

    Fastor::Tensor< out_scalar_type, Rest... > out( 0.0 );
    out_scalar_type*                           out_data = out.data();
    in_scalar_type*                            in_data  = in.data();

    for ( Fastor::FASTOR_INDEX i = 0; i < in.size(); ++i ) {
      out_data[out.get_mem_index( i )] = increaseDualOrderWithShift< order >( in_data[in.get_mem_index( i )] );
    }
    return out;
  }

  /** @typedef tensor_to_scalar_function_type
   *  @brief Alias for a function mapping a tensor to a scalar with JacobianDual numbers, one direction per
   *  tensor component
   *  @tparam Rest dimensions of the tensor
   * */
  template < size_t... Rest >
  using tensor_to_scalar_function_type = std::function< JacobianDual< ( Rest * ... ) >(
    const Fastor::Tensor< JacobianDual< ( Rest * ... ) >, Rest... >& T ) >;

  /** @brief Compute the gradient of a tensor-to-scalar function
   *  @tparam Rest dimensions of the tensor
   *  @param f function mapping a tensor to a scalar
   *  @param T input tensor at which the gradient is evaluated
   *  @param isSymmetric if true and T is a square rank-2 tensor, T (and thus the resulting gradient) is assumed to
   *  be symmetric, e.g. the right Cauchy-Green tensor. Only the dim*(dim+1)/2 upper-triangle entries of T are
   *  seeded, and the result is mirrored to the lower triangle, instead of all dim^2 entries -- a ~1.5x reduction
   *  in function evaluations for dim=3; the fraction of evaluations needed, (dim+1)/(2*dim), approaches but never
   *  reaches 1/2 as dim grows, so this never fully halves the evaluation count. Ignored for tensors that are not
   *  square rank-2.
   *  @return gradient of f with respect to T, same shape as T
   */
  template < size_t... Rest >
  Fastor::Tensor< double, Rest... > df_dT( const tensor_to_scalar_function_type< Rest... >& f,
                                           const Fastor::Tensor< double, Rest... >&         T,
                                           bool                                             isSymmetric = false )
  {
    using scalar = JacobianDual< ( Rest * ... ) >;

    Fastor::Tensor< double, Rest... > df_dT( 0.0 );
    Fastor::Tensor< scalar, Rest... > T_right = Tensors::makeOtherScalarType< scalar >( T );

    double* df_dT_data   = df_dT.data();
    scalar* T_right_data = T_right.data();

    if constexpr ( isVectorDual< scalar > ) {
      // all directions in one evaluation; isSymmetric would not save anything
      for ( Fastor::FASTOR_INDEX i = 0; i < T.size(); ++i )
        seedDirection( T_right_data[T_right.get_mem_index( i )], i );
      const scalar f_right = f( T_right );
      for ( Fastor::FASTOR_INDEX i = 0; i < T.size(); ++i )
        df_dT_data[df_dT.get_mem_index( i )] = partialDerivative( f_right, i );
    }
    else {

      if constexpr ( Marmot::Math::Tensors::IsSquareRank2Tensor< Rest... >::value ) {
        if ( isSymmetric ) {
          constexpr size_t dim = Marmot::Math::Tensors::IsSquareRank2Tensor< Rest... >::dim;

          for ( size_t row = 0; row < dim; ++row ) {
            for ( size_t col = row; col < dim; ++col ) {
              const Fastor::FASTOR_INDEX lin_ij  = row * dim + col;
              const Fastor::FASTOR_INDEX lin_ji  = col * dim + row;
              const int                  mem_idx = static_cast< int >( T_right.get_mem_index( lin_ij ) );

              seed< 1 >( T_right_data[mem_idx], 1.0 );
              const double val = derivative< 1 >( f( T_right ) );
              seed< 1 >( T_right_data[mem_idx], 0.0 );

              df_dT_data[df_dT.get_mem_index( lin_ij )] = val;
              if ( col != row )
                df_dT_data[df_dT.get_mem_index( lin_ji )] = val;
            }
          }

          return df_dT;
        }
      }

      for ( Fastor::FASTOR_INDEX i = 0; i < T.size(); ++i ) {
        const int T_right_mem_idx = static_cast< int >( T_right.get_mem_index( i ) );
        seed< 1 >( T_right_data[T_right_mem_idx], 1.0 );
        df_dT_data[df_dT.get_mem_index( i )] = derivative< 1 >( f( T_right ) );
        seed< 1 >( T_right_data[T_right_mem_idx], 0.0 );
      }
    }

    return df_dT;
  }

  /** @typedef tensor_to_scalar_function_type_arbitrary_dual_order
   *  @brief A type alias for a function mapping a tensor to a scalar with dual numbers of arbitrary order
   *  @tparam order current order of the dual numbers in the tensor
   *  @tparam Rest dimensions of the tensor
   * */
  template < size_t order, size_t... Rest >
  using tensor_to_scalar_function_type_arbitrary_dual_order = std::function< HyperDual< order >(
    const Fastor::Tensor< HyperDual< order >, Rest... >& T ) >;

  /** @brief Compute the gradient of a tensor-to-scalar function where the with dual numbers of arbitrary order
   *  @tparam order current order of the dual numers in the tensor
   *  @tparam Rest dimensions of the tensor
   *  @param f function mapping a tensor to a scalar
   *  @param T input tensor at which the gradient is evaluated
   *  @param isSymmetric if true and T is a square rank-2 tensor, T (and thus the resulting gradient) is assumed to
   *  be symmetric, e.g. the right Cauchy-Green tensor. Only the dim*(dim+1)/2 upper-triangle entries of T are
   *  seeded, and the result is mirrored to the lower triangle, instead of all dim^2 entries -- a ~1.5x reduction
   *  in function evaluations for dim=3; the fraction of evaluations needed, (dim+1)/(2*dim), approaches but never
   *  reaches 1/2 as dim grows, so this never fully halves the evaluation count. Ignored for tensors that are not
   *  square rank-2.
   *  @return pair of function value and gradient of f with respect to T, same shape as T
   */
  template < size_t order, size_t... Rest >
  std::pair< HyperDual< order >, Fastor::Tensor< HyperDual< order >, Rest... > > df_dT(
    const tensor_to_scalar_function_type_arbitrary_dual_order< order + 1, Rest... >& f,
    const Fastor::Tensor< HyperDual< order >, Rest... >&                             T,
    bool                                                                             isSymmetric = false )
  {

    using scalartype            = HyperDual< order >;
    using higherOrderScalartype = HyperDual< order + 1 >;

    higherOrderScalartype                            f_;
    Fastor::Tensor< scalartype, Rest... >            df_dT_( 0.0 );
    Fastor::Tensor< higherOrderScalartype, Rest... > T_right = increaseDualOrderWithShift< order >( T );

    higherOrderScalartype* T_right_data = T_right.data();
    scalartype*            df_dT_data   = df_dT_.data();

    if constexpr ( Marmot::Math::Tensors::IsSquareRank2Tensor< Rest... >::value ) {
      if ( isSymmetric ) {
        constexpr size_t dim = Marmot::Math::Tensors::IsSquareRank2Tensor< Rest... >::dim;

        for ( size_t row = 0; row < dim; ++row ) {
          for ( size_t col = row; col < dim; ++col ) {
            const Fastor::FASTOR_INDEX lin_ij = row * dim + col;
            const Fastor::FASTOR_INDEX lin_ji = col * dim + row;

            seed< 1 >( T_right_data[T_right.get_mem_index( lin_ij )], 1.0 );
            f_                   = f( T_right );
            const scalartype val = decreaseDualOrderWithShift< order + 1 >( f_ );
            seed< 1 >( T_right_data[T_right.get_mem_index( lin_ij )], 0.0 );

            df_dT_data[df_dT_.get_mem_index( lin_ij )] = val;
            if ( col != row )
              df_dT_data[df_dT_.get_mem_index( lin_ji )] = val;
          }
        }

        f_ = f( T_right );
        return { decreaseDualOrder< order + 1 >( f_ ), df_dT_ };
      }
    }

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
   *  @param isSymmetric if true, T is a square rank-2 tensor and F's output is also a square rank-2 tensor of the
   *  same dimension (e.g. the second Piola-Kirchhoff stress as a function of the right Cauchy-Green tensor), T is
   *  assumed to be symmetric and F is assumed to be transpose-equivariant (F(A^T) = F(A)^T for general, not
   *  necessarily symmetric, A), which holds for any tensor function built from tensor invariants/products. This
   *  reduces the number of function evaluations from dim^2 to dim*(dim+1)/2 -- a ~1.5x reduction for dim=3,
   *  approaching (but never reaching) a 2x reduction as dim grows -- by only seeding the upper triangle of T and,
   *  for each evaluation, filling both the direct entry and its transpose-mirrored counterpart
   *  \f$ \partial F_{ab}/\partial T_{kl} = \partial F_{ba}/\partial T_{lk} \f$. Ignored otherwise.
   *  @return pair of function value and gradient of F with respect to T, gradient has shape (RestF..., RestT...)
   */
  template < size_t... RestF, size_t... RestT >
  std::pair< Fastor::Tensor< double, RestF... >, Fastor::Tensor< double, RestF..., RestT... > > dF_dT(
    std::function< Fastor::Tensor< JacobianDual< ( RestT * ... ) >, RestF... >(
      const Fastor::Tensor< JacobianDual< ( RestT * ... ) >, RestT... >& ) >& F,
    const Fastor::Tensor< double, RestT... >&                                 T,
    bool                                                                      isSymmetric = false )
  {
    using scalar = JacobianDual< ( RestT * ... ) >;

    Fastor::Tensor< double, RestF... >           F_( 0.0 );
    Fastor::Tensor< double, RestF..., RestT... > dF_dT_( 0.0 );
    Fastor::Tensor< scalar, RestT... >           T_right = Tensors::makeOtherScalarType< scalar >( T );

    Fastor::Tensor< scalar, RestF... > F_at_T_right( 0.0 );

    double* dF_dT_data        = dF_dT_.data();
    scalar* T_right_data      = T_right.data();
    scalar* F_at_T_right_data = F_at_T_right.data();

    if constexpr ( isVectorDual< scalar > ) {
      // all directions in one evaluation; isSymmetric would not save anything
      for ( Fastor::FASTOR_INDEX i = 0; i < T.size(); ++i )
        seedDirection( T_right_data[T_right.get_mem_index( i )], i );
      F_at_T_right = F( T_right );
      for ( Fastor::FASTOR_INDEX j = 0; j < F_at_T_right.size(); ++j ) {
        const scalar& F_j                = F_at_T_right_data[F_at_T_right.get_mem_index( j )];
        F_.data()[F_.get_mem_index( j )] = value( F_j );
        for ( Fastor::FASTOR_INDEX i = 0; i < T.size(); ++i )
          dF_dT_data[dF_dT_.get_mem_index( j * T.size() + i )] = partialDerivative( F_j, i );
      }
    }
    else {

      if constexpr ( Marmot::Math::Tensors::IsSquareRank2Tensor< RestT... >::value &&
                     Marmot::Math::Tensors::IsSquareRank2Tensor< RestF... >::value &&
                     Marmot::Math::Tensors::IsSquareRank2Tensor< RestT... >::dim ==
                       Marmot::Math::Tensors::IsSquareRank2Tensor< RestF... >::dim ) {
        if ( isSymmetric ) {
          constexpr size_t dim = Marmot::Math::Tensors::IsSquareRank2Tensor< RestT... >::dim;

          for ( size_t row = 0; row < dim; ++row ) {
            for ( size_t col = row; col < dim; ++col ) {
              const Fastor::FASTOR_INDEX lin_ij          = row * dim + col;
              const Fastor::FASTOR_INDEX lin_ji          = col * dim + row;
              const int                  T_right_mem_idx = static_cast< int >( T_right.get_mem_index( lin_ij ) );
              seed< 1 >( T_right_data[T_right_mem_idx], 1.0 );
              F_at_T_right = F( T_right );

              for ( size_t outRow = 0; outRow < dim; ++outRow ) {
                for ( size_t outCol = 0; outCol < dim; ++outCol ) {
                  const Fastor::FASTOR_INDEX out_ab = outRow * dim + outCol;
                  const Fastor::FASTOR_INDEX out_ba = outCol * dim + outRow;
                  const double val = derivative< 1 >( F_at_T_right_data[F_at_T_right.get_mem_index( out_ab )] );
                  dF_dT_data[dF_dT_.get_mem_index( out_ab * T.size() + lin_ij )] = val;
                  if ( col != row )
                    dF_dT_data[dF_dT_.get_mem_index( out_ba * T.size() + lin_ji )] = val;
                }
              }
              seed< 1 >( T_right_data[T_right_mem_idx], 0.0 );
            }
          }

          F_ = Tensors::makeReal( F_at_T_right );
          return { F_, dF_dT_ };
        }
      }

      for ( Fastor::FASTOR_INDEX i = 0; i < T.size(); ++i ) {
        const int T_right_mem_idx = static_cast< int >( T_right.get_mem_index( i ) );
        seed< 1 >( T_right_data[T_right_mem_idx], 1.0 );
        F_at_T_right = F( T_right );

        for ( Fastor::FASTOR_INDEX j = 0; j < F_at_T_right.size(); ++j ) {
          dF_dT_data[dF_dT_.get_mem_index( j * T.size() + i )] = derivative< 1 >(
            F_at_T_right_data[F_at_T_right.get_mem_index( j )] );
        }
        seed< 1 >( T_right_data[T_right_mem_idx], 0.0 );
      }

      F_ = Tensors::makeReal( F_at_T_right );
    }

    return { F_, dF_dT_ };
  }

  namespace SecondOrder {

    /** @typedef tensor_to_scalar_function_type
     *  @brief Alias for a function mapping a tensor to a scalar with HessianDual numbers, whose dim * dim directions
     *  are the components of the tensor
     *  @tparam dim dimension of the input tensor (dim x dim)
     *
     *  @note The implementation is currently limited to second rank (dim x dim) tensors
     */
    template < size_t dim >
    using tensor_to_scalar_function_type = std::function< HessianDual< dim * dim >(
      const Fastor::Tensor< HessianDual< dim * dim >, dim, dim >& T ) >;

    /** @brief Computes the second derivative of a function that maps a tensor to a scalar
     *  @tparam dim dimension of the input tensor (dim x dim)
     *  @param F function mapping a (dim x dim) tensor to a scalar
     *  @param T input tensor at which the derivative is evaluated
     *  @return tuple of function value, first derivative and second derivative of F with respect to T
     *          first derivative has shape (dim, dim)
     *          second derivative has shape (dim, dim, dim, dim)
     *
     *  @note The implementation is currently limited to second rank (dim x dim) tensors
     *  @param isSymmetric if true, T is assumed to be symmetric, e.g. the right Cauchy-Green tensor, and F is
     *  assumed to be a transpose-invariant scalar function of T (as is the case for any scalar function built
     *  from tensor invariants). This reduces the number of outer-pair seed directions from dim^2 to dim*(dim+1)/2
     *  -- a ~1.5x reduction for dim=3, approaching (but never reaching) a 2x reduction as dim grows -- by
     *  exploiting the identity
     *  \f$ \partial^2 F / \partial T_{ij}\partial T_{kl} = \partial^2 F/\partial T_{ji}\partial T_{lk} \f$, which
     *  holds at a symmetric T for such functions (note: this mirrors both index pairs simultaneously; mirroring
     *  either pair independently is in general NOT valid)
     *
     *  With the vector-mode backend (isHessianDual), F is evaluated once and isSymmetric is irrelevant; otherwise F is
     *  evaluated once per pair of directions.
     */
    template < size_t dim >
    std::tuple< double, Fastor::Tensor< double, dim, dim >, Fastor::Tensor< double, dim, dim, dim, dim > > d2f_dT2(
      const tensor_to_scalar_function_type< dim >& F,
      const Fastor::Tensor< double, dim, dim >&    T,
      bool                                         isSymmetric = false )
    {
      using scalar = HessianDual< dim * dim >;

      double                                       F_;
      Fastor::Tensor< double, dim, dim >           dF_dT_;
      Fastor::Tensor< double, dim, dim, dim, dim > d2F_dT2;

      if constexpr ( isHessianDual< scalar > ) {
        Fastor::Tensor< scalar, dim, dim > T_right;
        for ( size_t i = 0; i < dim; i++ )
          for ( size_t j = 0; j < dim; j++ ) {
            T_right( i, j ).val() = T( i, j );
            seedDirection( T_right( i, j ), dim * i + j );
          }
        const scalar F_right = F( T_right );
        F_                   = value( F_right );
        for ( size_t i = 0; i < dim; i++ )
          for ( size_t j = 0; j < dim; j++ ) {
            dF_dT_( i, j ) = partialDerivative( F_right, dim * i + j );
            for ( size_t k = 0; k < dim; k++ )
              for ( size_t l = 0; l < dim; l++ )
                d2F_dT2( i, j, k, l ) = secondPartialDerivative( F_right, dim * i + j, dim * k + l );
          }
        return { F_, dF_dT_, d2F_dT2 };
      }
      else {
        scalar                             F_right;
        Fastor::Tensor< scalar, dim, dim > T_right = Tensors::makeHigherOrderDual< 2 >( T );

        for ( size_t i = 0; i < dim; i++ ) {
          for ( size_t j = isSymmetric ? i : 0; j < dim; j++ ) {
            seed< 1 >( T_right( i, j ), 1.0 );

            for ( size_t k = 0; k < dim; k++ ) {
              for ( size_t l = 0; l < dim; l++ ) {

                seed< 2 >( T_right( k, l ), 1.0 );
                F_right          = F( T_right );
                const double val = derivative< 2 >( F_right );

                d2F_dT2( i, j, k, l ) = val;
                if ( isSymmetric && j != i )
                  d2F_dT2( j, i, l, k ) = val;

                seed< 2 >( T_right( k, l ), 0.0 );
              }
            }
            dF_dT_( i, j ) = derivative< 1 >( F_right );
            F_             = value( F_right );
            if ( isSymmetric && j != i )
              dF_dT_( j, i ) = dF_dT_( i, j );
            seed< 1 >( T_right( i, j ), 0.0 );
          }
        }

        return { F_, dF_dT_, d2F_dT2 };
      }
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
     *  @param isSymmetric if true, T (and thus the resulting derivative) is assumed to be symmetric, e.g. the
     *  right Cauchy-Green tensor. Only the dim*(dim+1)/2 upper-triangle entries of T are seeded, and the result is
     *  mirrored to the lower triangle, instead of all dim^2 entries -- a ~1.5x reduction in function evaluations
     *  for dim=3; the fraction of evaluations needed, (dim+1)/(2*dim), approaches but never reaches 1/2 as dim
     *  grows, so this never fully halves the evaluation count
     */
    template < size_t dim >
    Fastor::Tensor< double, dim, dim > d2f_dTensor_dScalar( const tensor_and_scalar_to_scalar_function_type< dim >& F,
                                                            const Fastor::Tensor< double, dim, dim >&               T,
                                                            const double scalar,
                                                            bool         isSymmetric = false )
    {
      Fastor::Tensor< double, dim, dim >  d2F_dTdScalar;
      Fastor::Tensor< dual2nd, dim, dim > T_right = Tensors::makeHigherOrderDual< 2 >( T );

      dual2nd scalar_right( scalar );
      seed< 2 >( scalar_right, 1.0 );

      for ( size_t i = 0; i < dim; i++ ) {
        for ( size_t j = isSymmetric ? i : 0; j < dim; j++ ) {

          seed< 1 >( T_right( i, j ), 1.0 );

          d2F_dTdScalar( i, j ) = derivative< 2 >( F( T_right, scalar_right ) );

          seed< 1 >( T_right( i, j ), 0.0 );

          if ( isSymmetric && j != i )
            d2F_dTdScalar( j, i ) = d2F_dTdScalar( i, j );
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
     *  @return tuple of function value, first derivative, second derivative and third derivative of F with respect to
     * T first derivative has shape (dim, dim) second derivative has shape (dim, dim, dim, dim) third derivative has
     * shape (dim, dim, dim, dim, dim, dim)
     *
     *  @note The implementation is currently limited to second rank (dim x dim) tensors
     *
     *  F is evaluated once per unordered triple of the dim^2 directions (165 for dim = 3 instead of 729), since mixed
     *  partial derivatives do not depend on the order of differentiation.
     *  @param isSymmetric if true, T is assumed to be symmetric, e.g. the right Cauchy-Green tensor, and F is
     *  assumed to be a transpose-invariant scalar function of T (as is the case for any scalar function built
     *  from tensor invariants). This further skips the triples whose transpose is evaluated anyway, exploiting the
     *  identity
     *  \f$ \partial^3 F/\partial T_{ij}\partial T_{kl}\partial T_{mn} = \partial^3
     *  F/\partial T_{ji}\partial T_{lk}\partial T_{nm} \f$, which holds at a symmetric T for such functions (note:
     *  this mirrors all three index pairs simultaneously; mirroring any pair independently is in general NOT
     *  valid)
     */
    template < size_t dim >
    std::tuple< double,
                Fastor::Tensor< double, dim, dim >,
                Fastor::Tensor< double, dim, dim, dim, dim >,
                Fastor::Tensor< double, dim, dim, dim, dim, dim, dim > >
    d3f_dT3( const tensor_to_scalar_function_type< dim >& F,
             const Fastor::Tensor< double, dim, dim >&    T,
             bool                                         isSymmetric = false )
    {
      // Mixed partial derivatives do not depend on the order of differentiation: F is evaluated only for the sorted
      // triples a <= b <= c of directions (a = dim * i + j for T_ij), seeded on the levels 1, 2 and 3, and the third
      // derivative is copied to all permutations. With isSymmetric, a triple is skipped if its transpose (all three
      // directions T_ij -> T_ji) sorts lower, and its derivatives are copied from the transpose. The first and second
      // derivatives come out of the same evaluations; the few that none of them provides are evaluated separately.
      constexpr size_t n          = dim * dim;
      auto             transposed = []( size_t a ) { return ( a % dim ) * dim + a / dim; };

      double                                                 F_ = 0.0;
      Fastor::Tensor< double, dim, dim >                     dF_dT_;
      Fastor::Tensor< double, dim, dim, dim, dim >           d2F_dT2;
      Fastor::Tensor< double, dim, dim, dim, dim, dim, dim > d3F_dT3;
      Fastor::Tensor< dual3rd, dim, dim >                    T_right = Tensors::makeHigherOrderDual< 3 >( T );
      dual3rd*                                               t       = T_right.data();
      std::array< bool, n >                                  hasD1{};
      std::array< std::array< bool, n >, n >                 hasD2{};

      auto setD3 = [&]( size_t a, size_t b, size_t c, double v ) {
        const std::array< size_t, 3 >                   d     = { a, b, c };
        constexpr std::array< std::array< int, 3 >, 6 > perms = {
          { { 0, 1, 2 }, { 0, 2, 1 }, { 1, 0, 2 }, { 1, 2, 0 }, { 2, 0, 1 }, { 2, 1, 0 } } };
        for ( const auto& p : perms )
          d3F_dT3( d[p[0]] / dim, d[p[0]] % dim, d[p[1]] / dim, d[p[1]] % dim, d[p[2]] / dim, d[p[2]] % dim ) = v;
      };
      auto setD2 = [&]( size_t a, size_t b, double v ) {
        d2F_dT2( a / dim, a % dim, b / dim, b % dim ) = v;
        d2F_dT2( b / dim, b % dim, a / dim, a % dim ) = v;
        hasD2[a][b] = hasD2[b][a] = true;
      };
      auto setD1 = [&]( size_t a, double v ) {
        dF_dT_( a / dim, a % dim ) = v;
        hasD1[a]                   = true;
      };

      for ( size_t a = 0; a < n; a++ )
        for ( size_t b = a; b < n; b++ )
          for ( size_t c = b; c < n; c++ ) {
            std::array< size_t, 3 > dT = { transposed( a ), transposed( b ), transposed( c ) };
            std::sort( dT.begin(), dT.end() );
            if ( isSymmetric && dT < std::array< size_t, 3 >{ a, b, c } )
              continue;

            seed< 1 >( t[a], 1.0 );
            seed< 2 >( t[b], 1.0 );
            seed< 3 >( t[c], 1.0 );
            const dual3rd F_right = F( T_right );
            seed< 3 >( t[c], 0.0 );
            seed< 2 >( t[b], 0.0 );
            seed< 1 >( t[a], 0.0 );

            F_ = value( F_right );
            setD3( a, b, c, derivative< 3 >( F_right ) );
            setD2( a, b, derivative< 2 >( F_right ) );
            setD1( a, derivative< 1 >( F_right ) );
            if ( isSymmetric ) {
              setD3( transposed( a ), transposed( b ), transposed( c ), derivative< 3 >( F_right ) );
              setD2( transposed( a ), transposed( b ), derivative< 2 >( F_right ) );
              setD1( transposed( a ), derivative< 1 >( F_right ) );
            }
          }

      // second derivatives not yet provided (only possible with isSymmetric)
      for ( size_t a = 0; a < n; a++ )
        for ( size_t b = a; b < n; b++ ) {
          if ( hasD2[a][b] )
            continue;
          seed< 1 >( t[a], 1.0 );
          seed< 2 >( t[b], 1.0 );
          const dual3rd F_right = F( T_right );
          seed< 2 >( t[b], 0.0 );
          seed< 1 >( t[a], 0.0 );
          setD2( a, b, derivative< 2 >( F_right ) );
          setD1( a, derivative< 1 >( F_right ) );
        }

      // first derivatives not yet provided (only possible with isSymmetric)
      for ( size_t a = 0; a < n; a++ ) {
        if ( hasD1[a] )
          continue;
        if ( isSymmetric && hasD1[transposed( a )] ) {
          setD1( a, dF_dT_( transposed( a ) / dim, transposed( a ) % dim ) );
          continue;
        }
        seed< 1 >( t[a], 1.0 );
        const dual3rd F_right = F( T_right );
        seed< 1 >( t[a], 0.0 );
        setD1( a, derivative< 1 >( F_right ) );
      }

      return { F_, dF_dT_, d2F_dT2, d3F_dT3 };
    }

  } // namespace ThirdOrder

} // namespace Marmot::Math::Differentiation::Automatic
