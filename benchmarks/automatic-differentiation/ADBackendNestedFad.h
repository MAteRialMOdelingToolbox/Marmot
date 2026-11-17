#pragma once
// Shared implementation of the adapter vocabulary for "Fad-like" number types that expose a value node and an
// indexable derivative node (Sacado::Fad::SFad, ceres::Jet). `Access` provides val( x ) and dx( x, i ).
//
// Scalar mode (one direction per level, same semantics as autodiff):
//   HyperDual<N> = Fad< Fad< ... Fad< double, 1 > ... , 1 >, 1 >   (N levels)
//   seed<level>( x, v )      sets the epsilon of nesting level `level` (1 = outermost)
//   derivative<level>( x )   mixed derivative along levels 1..level
//
// Vector mode (all directions in one pass):
//   VecDual<N, depth> = Fad< Fad< ... Fad< double, N > ... , N >, N >   (depth levels)
//   seedDirection( x, i )    marks x as independent variable i on every level
//   vderivative( x, i, j, ... ) = d^k x / dx_i dx_j ...
#include <cstddef>
#include <type_traits>
#include <utility>

namespace ADBench {

  template < typename Access >
  struct NestedFadOps {

    template < typename S >
    static double value( const S& x )
    {
      if constexpr ( std::is_arithmetic_v< S > )
        return x;
      else
        return value( Access::val( x ) );
    }

    template < std::size_t level, typename S >
    static void seed( S& x, double v )
    {
      if constexpr ( level == 1 )
        Access::dx( x, 0 ) = Access::template constant< std::decay_t< decltype( Access::dx( x, 0 ) ) > >( v );
      else
        seed< level - 1 >( Access::val( x ), v );
    }

    template < std::size_t level, typename S >
    static double derivative( const S& x )
    {
      if constexpr ( level == 1 )
        return value( Access::dx( x, 0 ) );
      else
        return derivative< level - 1 >( Access::dx( x, 0 ) );
    }

    template < typename S >
    static void seedDirection( S& x, std::size_t i )
    {
      if constexpr ( !std::is_arithmetic_v< S > ) {
        Access::dx( x, i ) = Access::template constant< std::decay_t< decltype( Access::dx( x, i ) ) > >( 1.0 );
        seedDirection( Access::val( x ), i );
      }
    }

    template < typename S >
    static double vderivative( const S& x, std::size_t i )
    {
      return value( Access::dx( x, i ) );
    }
    template < typename S, typename... I >
    static double vderivative( const S& x, std::size_t i, std::size_t j, I... rest )
    {
      return vderivative( Access::dx( x, i ), j, rest... );
    }
  };

} // namespace ADBench
