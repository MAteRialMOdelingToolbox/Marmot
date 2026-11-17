#pragma once
// Backend adapter: Stan Math (https://github.com/stan-dev/math), forward mode stan::math::fvar< T >.
// fvar carries a single direction (val_, d_), so only scalar mode is available; it has no Taylor mode.
#include "ADBackendNestedFad.h"
#include <stan/math/fwd.hpp>

namespace ADBench {

  namespace StanDetail {
    template < std::size_t N >
    struct HyperDualImpl {
      using type = stan::math::fvar< typename HyperDualImpl< N - 1 >::type >;
    };
    template <>
    struct HyperDualImpl< 1 > {
      using type = stan::math::fvar< double >;
    };

    struct Access {
      template < typename S >
      static auto& val( S& x )
      {
        return x.val_;
      }
      template < typename S >
      static const auto& val( const S& x )
      {
        return x.val_;
      }
      template < typename S >
      static auto& dx( S& x, std::size_t )
      {
        return x.d_;
      }
      template < typename S >
      static const auto& dx( const S& x, std::size_t )
      {
        return x.d_;
      }
      template < typename S >
      static S constant( double v )
      {
        return S( v );
      }
    };
  } // namespace StanDetail

  struct StanBackend : NestedFadOps< StanDetail::Access > {
    static constexpr const char* name           = "stan";
    static constexpr bool        hasReverseMode = false;
    static constexpr bool        hasTaylor      = false;
    static constexpr bool        hasVectorMode  = false;
    static constexpr int         maxScalarOrder = 3;
    static constexpr int         maxVectorOrder = 0;

    template < std::size_t N >
    using HyperDual = typename StanDetail::HyperDualImpl< N >::type;

    template < typename S >
    static S constant( double v )
    {
      return S( v );
    }
  };

} // namespace ADBench
