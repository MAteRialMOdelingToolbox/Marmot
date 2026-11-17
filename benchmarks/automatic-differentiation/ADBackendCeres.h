#pragma once
// Backend adapter: ceres::Jet<T, N> (header-only use of ceres/jet.h). Jet has no Taylor mode.
#include "ADBackendNestedFad.h"
#include "ceres/jet.h"
#include <array>

namespace ADBench {

  namespace CeresDetail {
    template < std::size_t N, std::size_t depth >
    struct VecDualImpl {
      using type = ceres::Jet< typename VecDualImpl< N, depth - 1 >::type, int( N ) >;
    };
    template < std::size_t N >
    struct VecDualImpl< N, 1 > {
      using type = ceres::Jet< double, int( N ) >;
    };

    // Jet's value constructor is explicit, so nested Jets are built level by level
    template < typename S >
    S constant( double v )
    {
      if constexpr ( std::is_arithmetic_v< S > )
        return v;
      else
        return S( constant< std::decay_t< decltype( std::declval< S >().a ) > >( v ) );
    }

    struct Access {
      template < typename S >
      static auto& val( S& x )
      {
        return x.a;
      }
      template < typename S >
      static const auto& val( const S& x )
      {
        return x.a;
      }
      template < typename S >
      static auto& dx( S& x, std::size_t i )
      {
        return x.v[i];
      }
      template < typename S >
      static const auto& dx( const S& x, std::size_t i )
      {
        return x.v[i];
      }
      template < typename S >
      static S constant( double v )
      {
        return CeresDetail::constant< S >( v );
      }
    };
  } // namespace CeresDetail

  struct CeresBackend : NestedFadOps< CeresDetail::Access > {
    static constexpr const char* name          = "ceres";
    static constexpr bool        hasTaylor     = false;
    static constexpr bool        hasVectorMode = true;
    static constexpr int         maxOrder      = 3; // needs the patched jet.h, see ceres-2.2.0-jet-fastor.patch

    template < std::size_t N >
    using HyperDual = typename CeresDetail::VecDualImpl< 1, N >::type;

    template < std::size_t N, std::size_t depth >
    using VecDual = typename CeresDetail::VecDualImpl< N, depth >::type;

    template < typename S >
    static S constant( double v )
    {
      return CeresDetail::constant< S >( v );
    }

    template < std::size_t N >
    using Taylor = void;
  };

} // namespace ADBench

// ceres::Jet<T, N> only defines scalar arithmetic with its own T. For nested Jets (T = Jet<...>) that rules out
// expressions like `K / 8. * J` that Marmot's templates use throughout, so add the double overloads (found by ADL).
namespace ceres {
#define MARMOT_JET_DOUBLE_OP( OP )                                                                                    \
  template < typename T, int N, std::enable_if_t< !std::is_arithmetic_v< T >, int > = 0 >                             \
  inline Jet< T, N > operator OP( const Jet< T, N >& f, double s )                                                    \
  {                                                                                                                   \
    return f OP ADBench::CeresDetail::constant< T >( s );                                                             \
  }                                                                                                                   \
  template < typename T, int N, std::enable_if_t< !std::is_arithmetic_v< T >, int > = 0 >                             \
  inline Jet< T, N > operator OP( double s, const Jet< T, N >& f )                                                    \
  {                                                                                                                   \
    return ADBench::CeresDetail::constant< T >( s ) OP f;                                                             \
  }
  MARMOT_JET_DOUBLE_OP( + )
  MARMOT_JET_DOUBLE_OP( - )
  MARMOT_JET_DOUBLE_OP( * )
  MARMOT_JET_DOUBLE_OP( / )
#undef MARMOT_JET_DOUBLE_OP
} // namespace ceres
