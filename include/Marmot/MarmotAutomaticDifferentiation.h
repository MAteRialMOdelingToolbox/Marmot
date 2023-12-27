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
#include "Marmot/MarmotTensor.h"
#include "autodiff/forward/dual.hpp"
#include "autodiff/forward/dual/eigen.hpp"
#include <autodiff/forward/dual/dual.hpp>
#include <functional>

namespace Marmot {

  namespace AutomaticDifferentiation {

    using namespace autodiff;
    using namespace Eigen;

    dual2nd        shiftTo2ndOrderDual( const dual& x );
    VectorXdual2nd shiftTo2ndOrderDual( const VectorXdual& X );

    template < size_t order, typename T, typename G >
    auto& valnode( const Dual< T, G >& dual )
    {
      constexpr auto N = detail::Order< Dual< T, G > >;
      static_assert( order <= N );
      if constexpr ( order == 0 )
        return dual.val;
      else if constexpr ( order == 1 )
        return dual.val;
      else
        return valnode< order - 1 >( dual.val );
    }
    template < size_t order, typename T, typename G >
    auto& valnode( Dual< T, G >& dual )
    {
      constexpr auto N = detail::Order< Dual< T, G > >;
      static_assert( order <= N );
      if constexpr ( order == 0 )
        return dual.val;
      else if constexpr ( order == 1 )
        return dual.val;
      else
        return valnode< order - 1 >( dual.val );
    }

    template < size_t order >
    autodiff::HigherOrderDual< order + 1, double > increaseDualOrderWithShift(
      const autodiff::HigherOrderDual< order, double >& in )
    {
      using out_scalar_type = autodiff::HigherOrderDual< order + 1, double >;
      using namespace autodiff::detail;

      out_scalar_type out( 0.0 );
      const double*   in_point  = &valnode< order >( in );
      double*         out_point = &valnode< order + 1 >( out );

      for ( size_t i = 0; i < size_t( std::pow( 2, order ) ); i++ ) {
        *( out_point + i ) = *( in_point + i );
      }

      return out;
    }
    template < size_t order >
    autodiff::HigherOrderDual< order - 1, double > decreaseDualOrderWithShift(
      autodiff::HigherOrderDual< order, double >& in )
    {
      using out_scalar_type = autodiff::HigherOrderDual< order - 1, double >;
      using namespace autodiff::detail;

      out_scalar_type out( 0.0 );
      double*         in_point  = &valnode< order >( in );
      double*         out_point = &valnode< order - 1 >( out );

      for ( size_t i = 0; i < size_t( std::pow( 2, order - 1 ) ); i++ ) {
        *( out_point + i ) = *( in_point + i + size_t( std::pow( 2, order - 1 ) ) );
      }

      return out;
    }

    template < size_t order >
    Vector< HigherOrderDual< order + 1, double >, -1 > increaseDualOrderWithShift(
      const Vector< HigherOrderDual< order, double >, -1 >& in )
    {
      using in_scalar_type  = HigherOrderDual< order, double >;
      using out_scalar_type = HigherOrderDual< order + 1, double >;

      Vector< out_scalar_type, -1 > out      = Vector< out_scalar_type, -1 >( in.size() );
      out_scalar_type*              out_data = out.data();
      const in_scalar_type*         in_data  = in.data();

      for ( int i = 0; i < in.size(); i++ ) {
        out_data[i] = increaseDualOrderWithShift< order >( in_data[i] );
      }
      return out;
    }

    using scalar_to_scalar_function_type = std::function< dual( const dual& ) >;
    double df_dx( const scalar_to_scalar_function_type& f, const double& x );

    using scalar_to_scalar_function_type_2nd = std::function< dual2nd( const dual2nd& ) >;
    dual df_dx( const scalar_to_scalar_function_type_2nd& f, const dual& x );

    using vector_to_vector_function_type = std::function< VectorXdual( const VectorXdual& X ) >;
    MatrixXd forwardMode( const vector_to_vector_function_type& F, const VectorXd& X );

    using vector_to_vector_function_type_dual = std::function< VectorXdual( const VectorXdual& X ) >;
    std::pair< VectorXd, MatrixXd > jacobian( const vector_to_vector_function_type_dual& F, const VectorXd& X );

    using vector_to_vector_function_type_dual2nd = std::function< VectorXdual2nd( const VectorXdual2nd& X ) >;
    std::pair< VectorXdual, MatrixXdual > jacobian2nd( const vector_to_vector_function_type_dual2nd& F,
                                                       const VectorXdual&                            X );
  } // namespace AutomaticDifferentiation

} // namespace Marmot
