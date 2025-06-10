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
#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotNumericalDifferentiation.h"
#include "Marmot/MarmotTypedefs.h"
#include <Fastor/config/macros.h>
#include <Fastor/tensor/Tensor.h>
#include <functional>

namespace Marmot {
  namespace NumericalAlgorithms::Differentiation {

    namespace ScalarToTensor {

      template < size_t... Rest >
      Fastor::Tensor< double, Rest... > forwardDifference(
        const std::function< Fastor::Tensor< double, Rest... >( const double ) >& F,
        const double                                                              x )
      {

        double volatile h                       = std::max( 1.0, std::abs( x ) ) * Marmot::Constants::SquareRootEps;
        Fastor::Tensor< double, Rest... > dF    = F( x + h ) - F( x );
        Fastor::Tensor< double, Rest... > dF_dx = multiplyFastorTensorWithScalar( dF, 1. / h );
        return dF_dx;
      }

      template < size_t... Rest >
      Fastor::Tensor< double, Rest... > centralDifference(
        const std::function< Fastor::Tensor< double, Rest... >( const double ) >& F,
        const double                                                              x )
      {

        double volatile h                       = std::max( 1.0, std::abs( x ) ) * Marmot::Constants::CubicRootEps;
        Fastor::Tensor< double, Rest... > dF    = F( x + h ) - F( x - h );
        Fastor::Tensor< double, Rest... > dF_dx = multiplyFastorTensorWithScalar( dF, 1. / ( 2. * h ) );
        return dF_dx;
      }
    } // namespace ScalarToTensor

    namespace TensorToScalar {

      template < size_t... Rest >
      using tensor_to_scalar_function_type = std::function< double( const Fastor::Tensor< double, Rest... >& T ) >;

      /* template < size_t... Rest > */
      /* Fastor::Tensor< double, Rest... > forwardDifference( const tensor_to_scalar_function_type< Rest... >& F, */
      /*                                                      const Fastor::Tensor< double, Rest... >&         T ) */
      /* { */

      /*   Fastor::Tensor< double, Rest... >       dF_dT; */
      /*   const Fastor::Tensor< double, Rest... > F_      = F( T ); */
      /*   Fastor::Tensor< double, Rest... >       T_right = T; */

      /*   double* dF_dT_data   = dF_dT.data(); */
      /*   double* T_right_data = T_right.data(); */

      /*   for ( Fastor::FASTOR_INDEX i = 0; i < T.size(); ++i ) { */
      /*     const int dF_dT_mem_idx   = dF_dT.get_mem_index( i ); */
      /*     const int T_right_mem_idx = T_right.get_mem_index( i ); */
      /*     double volatile h         = std::max( 1.0, std::abs( double( T_right_data[T_right_mem_idx] ) ) ) * */
      /*                         Marmot::Constants::SquareRootEps; */
      /*     T_right_data[T_right_mem_idx] += h; */
      /*     dF_dT_data[dF_dT_mem_idx] = ( F( T_right ) - F_ ) / ( 1. * h ); */
      /*     T_right_data[T_right_mem_idx] -= h; */
      /*   } */

      template < size_t dim >
      Fastor::Tensor< double, dim, dim > forwardDifference( const tensor_to_scalar_function_type< dim, dim >& f,
                                                            const Fastor::Tensor< double, dim, dim >&         T )
      {
        Fastor::Tensor< double, dim, dim > T_right( T );
        Fastor::Tensor< double, dim, dim > dF_dT( 0.0 );
        const double                       f_ = f( T );

        for ( size_t i = 0; i < dim; i++ ) {
          for ( size_t j = 0; j < dim; j++ ) {
            double volatile h = std::max( 1.0, std::abs( T( i, j ) ) ) * Marmot::Constants::SquareRootEps;

            T_right( i, j ) += h;
            dF_dT( i, j ) = ( f( T_right ) - f_ ) / ( 1. * h );
            T_right( i, j ) -= h;
          }
        }

        return dF_dT;
      }

      template < size_t dim >
      Fastor::Tensor< double, dim, dim > centralDifference( const tensor_to_scalar_function_type< dim, dim >& F,
                                                            const Fastor::Tensor< double, dim, dim >&         T )
      {

        Fastor::Tensor< double, dim, dim > dF_dT;
        Fastor::Tensor< double, dim, dim > T_right( T );
        Fastor::Tensor< double, dim, dim > T_left( T );

        for ( size_t i = 0; i < dim; i++ ) {
          for ( size_t j = 0; j < dim; j++ ) {
            double volatile h = std::max( 1.0, std::abs( T( i, j ) ) ) * Marmot::Constants::CubicRootEps;

            T_right = T;
            T_left  = T;
            T_right( i, j ) += h;
            T_left( i, j ) -= h;

            dF_dT( i, j ) = ( F( T_right ) - F( T_left ) ) / ( 2. * h );
          }
        }

        return dF_dT;
      }
    } // namespace TensorToScalar

    namespace TensorToTensor {

      template < size_t... RestF, size_t... RestT >
      Fastor::Tensor< double, RestF..., RestT... > forwardDifference(
        const std::function< Fastor::Tensor< double, RestF... >( const Fastor::Tensor< double, RestT... >& ) >& F,
        const Fastor::Tensor< double, RestT... >&                                                               T )
      {

        Fastor::Tensor< double, RestF..., RestT... > dF_dT( 0.0 );
        Fastor::Tensor< double, RestT... >           T_right( T );

        Fastor::Tensor< double, RestF... > F_at_T = F( T );
        Fastor::Tensor< double, RestF... > F_at_T_right( 0.0 );

        double* dF_dT_data        = dF_dT.data();
        double* T_right_data      = T_right.data();
        double* F_at_T_data       = F_at_T.data();
        double* F_at_T_right_data = F_at_T_right.data();

        for ( Fastor::FASTOR_INDEX i = 0; i < T.size(); ++i ) {
          const int T_right_mem_idx = T_right.get_mem_index( i );
          double volatile h         = std::max( 1.0, std::abs( double( T.data()[T_right_mem_idx] ) ) ) *
                              Marmot::Constants::SquareRootEps;
          T_right = T;
          T_right_data[T_right_mem_idx] += h;
          F_at_T_right = F( T_right );

          for ( Fastor::FASTOR_INDEX j = 0; j < F_at_T_right.size(); ++j ) {
            dF_dT_data[dF_dT.get_mem_index( j * T.size() + i )] = ( F_at_T_right_data[F_at_T_right.get_mem_index( j )] -
                                                                    F_at_T_data[F_at_T.get_mem_index( j )] ) /
                                                                  ( 1. * h );
          }
        }

        return dF_dT;
      }

      template < size_t... Rest1, size_t... Rest2 >
      Fastor::Tensor< double, Rest1..., Rest2... > centralDifference(
        const std::function< Fastor::Tensor< double, Rest1... >( const Fastor::Tensor< double, Rest2... >& ) >& F,
        const Fastor::Tensor< double, Rest2... >&                                                               T )
      {

        Fastor::Tensor< double, Rest1..., Rest2... > dF_dT( 0.0 );
        Fastor::Tensor< double, Rest2... >           T_right( T );
        Fastor::Tensor< double, Rest2... >           T_left( T );

        Fastor::Tensor< double, Rest1... > F_at_T_right( 0.0 );
        Fastor::Tensor< double, Rest1... > F_at_T_left( 0.0 );

        double* dF_dT_data        = dF_dT.data();
        double* T_right_data      = T_right.data();
        double* T_left_data       = T_left.data();
        double* F_at_T_right_data = F_at_T_right.data();
        double* F_at_T_left_data  = F_at_T_left.data();

        for ( Fastor::FASTOR_INDEX i = 0; i < T.size(); ++i ) {
          const int T_right_mem_idx = T_right.get_mem_index( i );
          const int T_left_mem_idx  = T_left.get_mem_index( i );
          double volatile h         = std::max( 1.0, std::abs( double( T.data()[T_right_mem_idx] ) ) ) *
                              Marmot::Constants::CubicRootEps;
          T_left = T;
          T_left_data[T_left_mem_idx] -= h;
          F_at_T_left = F( T_left );

          T_right = T;
          T_right_data[T_right_mem_idx] += h;
          F_at_T_right = F( T_right );

          for ( Fastor::FASTOR_INDEX j = 0; j < F_at_T_right.size(); ++j ) {
            dF_dT_data[dF_dT.get_mem_index( j * T.size() + i )] = ( F_at_T_right_data[F_at_T_right.get_mem_index( j )] -
                                                                    F_at_T_left_data[F_at_T_left.get_mem_index( j )] ) /
                                                                  ( 2. * h );
          }
        }

        return dF_dT;
      }
    } // namespace TensorToTensor

    namespace Complex {

      using complexDouble                                  = std::complex< double >;
      static const double        imaginaryPerturbationSize = 1e-20;
      static const complexDouble imaginaryPerturbation     = imaginaryPerturbationSize * imaginaryUnit;

      template < size_t dim >
      using tensor_to_scalar_function_type = std::function< std::complex< double >(
        const Fastor::Tensor< std::complex< double >, dim, dim >& T ) >;

      template < size_t dim >
      Fastor::Tensor< double, dim, dim > forwardDifference( const tensor_to_scalar_function_type< dim >& F,
                                                            const Fastor::Tensor< double, dim, dim >&    T )
      {
        Fastor::Tensor< double, dim, dim > dF_dT;
        Fastor::Tensor< std::complex< double >, dim, dim >
          T_right = fastorTensorFromDoubleTensor< std::complex< double >, dim >( T );

        for ( size_t i = 0; i < dim; i++ ) {
          for ( size_t j = 0; j < dim; j++ ) {
            T_right( i, j ) += imaginaryPerturbation;

            dF_dT( i, j ) = F( T_right ).imag() / imaginaryPerturbationSize;

            T_right( i, j ) -= imaginaryPerturbation;
          }
        }

        return dF_dT;
      }
      namespace ScalarToTensor {

        template < size_t... Rest >
        std::tuple< Fastor::Tensor< double, Rest... >, Fastor::Tensor< double, Rest... > > forwardDifference(
          const std::function< Fastor::Tensor< complexDouble, Rest... >( const complexDouble ) >& F,
          const double                                                                            x )
        {
          const complexDouble                      xRight = x + imaginaryPerturbation;
          Fastor::Tensor< complexDouble, Rest... > F_     = F( xRight );
          Fastor::Tensor< double, Rest... >        dF_dx( 0.0 );
          Fastor::Tensor< double, Rest... >        FReal( 0.0 );
          complexDouble*                           F_data     = F_.data();
          double*                                  dF_dx_data = dF_dx.data();
          double*                                  FReal_data = FReal.data();

          for ( Fastor::FASTOR_INDEX j = 0; j < F_.size(); ++j ) {
            FReal_data[FReal.get_mem_index( j )] = F_data[F_.get_mem_index( j )].real();
            dF_dx_data[dF_dx.get_mem_index( j )] = F_data[F_.get_mem_index( j )].imag() / imaginaryPerturbationSize;
          }

          return { FReal, dF_dx };
        }
      } // namespace ScalarToTensor

      namespace TensorToScalar {

        template < size_t... Rest >
        using tensor_to_scalar_function_type = std::function< complexDouble(
          const Fastor::Tensor< complexDouble, Rest... >& T ) >;

        template < size_t dim >
        Fastor::Tensor< double, dim, dim > forwardDifference( const tensor_to_scalar_function_type< dim, dim >& f,
                                                              const Fastor::Tensor< double, dim, dim >&         T )
        {
          Fastor::Tensor< complexDouble, dim, dim > T_right = fastorTensorFromDoubleTensor< complexDouble >( T );
          Fastor::Tensor< double, dim, dim >        dF_dT( 0.0 );

          for ( size_t i = 0; i < dim; i++ ) {
            for ( size_t j = 0; j < dim; j++ ) {

              T_right( i, j ) += imaginaryPerturbation;
              dF_dT( i, j ) = f( T_right ).imag() / imaginaryPerturbationSize;
              T_right( i, j ) -= imaginaryPerturbation;
            }
          }

          return dF_dT;
          /* Fastor::Tensor< double, Rest... >       dF_dT; */
          /* Fastor::Tensor< complexDouble, Rest... >       T_right = fastorTensorFromDoubleTensor<complexDouble>(T); */

          /* double* dF_dT_data   = dF_dT.data(); */
          /* complexDouble* T_right_data = T_right.data(); */

          /* for ( Fastor::FASTOR_INDEX i = 0; i < T.size(); ++i ) { */
          /*   const int dF_dT_mem_idx   = dF_dT.get_mem_index( i ); */
          /*   const int T_right_mem_idx = T_right.get_mem_index( i ); */
          /*   T_right_data[T_right_mem_idx] += imaginaryPerturbation; */
          /*   dF_dT_data[dF_dT_mem_idx] = f( T_right )/ imaginaryPerturbationSize; */
          /*   T_right_data[T_right_mem_idx] -= imaginaryPerturbation; */
          /* } */

          /* return dF_dT; */
        }
      } // namespace TensorToScalar
      namespace TensorToTensor {
        template < size_t... RestF, size_t... RestT >
        Fastor::Tensor< double, RestF..., RestT... > forwardDifference(
          std::function<
            Fastor::Tensor< complexDouble, RestF... >( const Fastor::Tensor< complexDouble, RestT... >& ) >& F,
          const Fastor::Tensor< double, RestT... >&                                                          T )
        {

          Fastor::Tensor< double, RestF..., RestT... > dF_dT( 0.0 );
          Fastor::Tensor< complexDouble, RestT... >    T_right = fastorTensorFromDoubleTensor< complexDouble >( T );

          Fastor::Tensor< complexDouble, RestF... > F_at_T_right( 0.0 );

          double*        dF_dT_data        = dF_dT.data();
          complexDouble* T_right_data      = T_right.data();
          complexDouble* F_at_T_right_data = F_at_T_right.data();

          for ( Fastor::FASTOR_INDEX i = 0; i < T.size(); ++i ) {
            const int T_right_mem_idx = T_right.get_mem_index( i );

            T_right_data[T_right_mem_idx] += imaginaryPerturbation;
            F_at_T_right = F( T_right );

            for ( Fastor::FASTOR_INDEX j = 0; j < F_at_T_right.size(); ++j ) {
              dF_dT_data[dF_dT.get_mem_index( j * T.size() +
                                              i )] = ( F_at_T_right_data[F_at_T_right.get_mem_index( j )] ).imag() /
                                                     imaginaryPerturbationSize;
            }
            T_right_data[T_right_mem_idx] -= imaginaryPerturbation;
          }

          return dF_dT;
        }
      } // namespace TensorToTensor
    }   // namespace Complex
  }     // namespace NumericalAlgorithms::Differentiation
} // namespace Marmot
