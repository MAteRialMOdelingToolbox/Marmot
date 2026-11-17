#pragma once
// Backend adapter: autodiff (https://github.com/autodiff/autodiff), the library Marmot currently uses.
//
// The adapter exposes the minimal, backend-neutral vocabulary Marmot needs:
//   HyperDual<N>             N-th order nested hyper-dual number (dual, dual2nd, dual3rd, ...)
//   Taylor<N>                univariate truncated Taylor number of order N (autodiff::Real<N>)
//   seed<level>( x, v )      set the unit-epsilon of nesting level `level` (1 = outermost), autodiff semantics
//   derivative<level>( x )   mixed derivative along levels 1..level (autodiff::derivative<level> semantics)
//   value( x )               primal value
//   taylorDerivatives<N>( f, x ) all derivatives 0..N of a univariate function
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <autodiff/forward/real.hpp>
#include <array>

namespace ADBench {

  struct AutodiffBackend {

    static constexpr const char* name = "autodiff";
    static constexpr bool        hasTaylor     = true;
    static constexpr bool        hasVectorMode = false;
    static constexpr int         maxOrder      = 3;

    template < typename S >
    static S constant( double v )
    {
      return S( v );
    }


  template < std::size_t N >
  using HyperDual = autodiff::HigherOrderDual< N, double >;

  template < std::size_t N >
  using Taylor = autodiff::Real< N, double >;

  template < std::size_t level, typename S >
  static void seed( S& x, double v )
  {
    autodiff::seed< level >( x, v );
  }

  template < std::size_t level, typename S >
  static double derivative( const S& x )
  {
    return autodiff::derivative< level >( x );
  }

  template < typename S >
  static double value( const S& x )
  {
    return autodiff::detail::val( x );
  }

  template < std::size_t N, typename F >
  static std::array< double, N + 1 > taylorDerivatives( const F& f, double x )
  {
    Taylor< N > x_( x );
    const auto  d = autodiff::derivatives( f, autodiff::along( 1. ), autodiff::at( x_ ) );
    std::array< double, N + 1 > out;
    for ( std::size_t i = 0; i <= N; ++i )
      out[i] = d[i];
    return out;
  }

  };

} // namespace ADBench
