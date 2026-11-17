#pragma once
// Backend adapter: Sacado (Trilinos), header-only use of Sacado::Fad::SFad (static, expression templates)
// and Sacado::Tay::Taylor (dynamic degree).
#include "ADBackendNestedFad.h"
#include "Sacado.hpp"
#include <Eigen/Core>
#include <array>

namespace ADBench {

  namespace SacadoDetail {
    template < std::size_t N, std::size_t depth >
    struct VecDualImpl {
      using type = Sacado::Fad::SFad< typename VecDualImpl< N, depth - 1 >::type, N >;
    };
    template < std::size_t N >
    struct VecDualImpl< N, 1 > {
      using type = Sacado::Fad::SFad< double, N >;
    };

    struct Access {
      template < typename S >
      static auto& val( S& x )
      {
        return x.val();
      }
      template < typename S >
      static const auto& val( const S& x )
      {
        return x.val();
      }
      template < typename S >
      static auto& dx( S& x, std::size_t i )
      {
        return x.fastAccessDx( i );
      }
      template < typename S >
      static const auto& dx( const S& x, std::size_t i )
      {
        return x.fastAccessDx( i );
      }
      template < typename S >
      static S constant( double v )
      {
        return S( v );
      }
    };
  } // namespace SacadoDetail

  struct SacadoBackend : NestedFadOps< SacadoDetail::Access > {
    static constexpr const char* name           = "sacado";
    static constexpr bool        hasReverseMode = false;
    static constexpr bool        hasTaylor      = true;
    static constexpr bool        hasVectorMode  = true;
    static constexpr int         maxScalarOrder = 3;
    static constexpr int         maxVectorOrder = 3;

    template < std::size_t N >
    using HyperDual = typename SacadoDetail::VecDualImpl< 1, N >::type;

    template < std::size_t N, std::size_t depth >
    using VecDual = typename SacadoDetail::VecDualImpl< N, depth >::type;

    template < typename S >
    static S constant( double v )
    {
      return S( v );
    }

    // Sacado::Tay::Taylor has a run-time degree; the alias keeps the vocabulary uniform
    template < std::size_t N >
    using Taylor = Sacado::Tay::Taylor< double >;

    template < std::size_t N, typename F >
    static std::array< double, N + 1 > taylorDerivatives( const F& f, double x )
    {
      Sacado::Tay::Taylor< double > x_{ int( N ), x };
      x_.fastAccessCoeff( 1 )               = 1.0;
      const Sacado::Tay::Taylor< double > y = f( x_ );
      std::array< double, N + 1 >         out;
      double                              fact = 1.;
      for ( std::size_t i = 0; i <= N; ++i ) {
        if ( i > 1 )
          fact *= double( i );
        out[i] = y.coeff( int( i ) ) * fact;
      }
      return out;
    }
  };

} // namespace ADBench

namespace Eigen {
  template < typename T, int N >
  struct NumTraits< Sacado::Fad::SFad< T, N > > : NumTraits< double > {
    using Real       = Sacado::Fad::SFad< T, N >;
    using NonInteger = Real;
    using Literal    = double;
    using Nested     = Real;
    enum {
      IsComplex             = 0,
      IsInteger             = 0,
      IsSigned              = 1,
      RequireInitialization = 1,
      ReadCost              = 1,
      AddCost               = 3,
      MulCost               = 3
    };
  };
  template < typename T, int N, typename BinOp >
  struct ScalarBinaryOpTraits< Sacado::Fad::SFad< T, N >, double, BinOp > {
    using ReturnType = Sacado::Fad::SFad< T, N >;
  };
  template < typename T, int N, typename BinOp >
  struct ScalarBinaryOpTraits< double, Sacado::Fad::SFad< T, N >, BinOp > {
    using ReturnType = Sacado::Fad::SFad< T, N >;
  };
} // namespace Eigen
