#pragma once
// Backend adapter: Marmot's own AD vocabulary (Marmot/MarmotAutomaticDifferentiationTypes.h), currently backed by
// autodiff (https://github.com/autodiff/autodiff). This is the reference every other backend is compared against.
//
// The adapters expose the vocabulary of the Marmot types header:
//   HyperDual<N>                 N-th order nested hyper-dual number (dual, dual2nd, dual3rd, ...)
//   Taylor<N>                    univariate truncated Taylor number of order N
//   seed<level>( x, v )          set the epsilon of nesting level `level` (1 = outermost)
//   derivative<level>( x )       mixed derivative along levels 1..level
//   value( x )                   primal value
//   taylorDerivatives<N>( f, x ) all derivatives 0..N of a univariate function
#include "Marmot/MarmotAutomaticDifferentiationTypes.h"
#include <array>

namespace ADBench {

  struct AutodiffBackend {
    static constexpr const char* name           = "autodiff";
    static constexpr bool        hasReverseMode = false;
    static constexpr bool        hasTaylor      = true;
    static constexpr int         maxScalarOrder = 3;
    static constexpr bool        hasVectorMode  = false;
    static constexpr int         maxVectorOrder = 0;

    template < std::size_t N >
    using HyperDual = Marmot::Math::Differentiation::Automatic::HyperDual< N >;

    template < std::size_t N >
    using Taylor = Marmot::Math::Differentiation::Automatic::Taylor< N >;

    template < typename S >
    static S constant( double v )
    {
      return S( v );
    }

    template < std::size_t level, typename S >
    static void seed( S& x, double v )
    {
      Marmot::Math::Differentiation::Automatic::seed< level >( x, v );
    }

    template < std::size_t level, typename S >
    static double derivative( const S& x )
    {
      return Marmot::Math::Differentiation::Automatic::derivative< level >( x );
    }

    template < typename S >
    static double value( const S& x )
    {
      return Marmot::Math::Differentiation::Automatic::value( x );
    }

    template < std::size_t N, typename F >
    static std::array< double, N + 1 > taylorDerivatives( const F& f, double x )
    {
      return Marmot::Math::Differentiation::Automatic::taylorDerivatives< N >( f, x );
    }
  };

} // namespace ADBench
