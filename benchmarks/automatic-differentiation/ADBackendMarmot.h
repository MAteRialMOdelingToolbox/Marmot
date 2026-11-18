#pragma once
// Backend adapter: Marmot's own AD types, independent of the backend the types header is configured with.
//   scalar mode: ScalarMode::Dual (nested hyper-dual numbers, eager) and ScalarMode::Taylor
//   vector mode: VectorMode::Dual (static, expression templates with eagerly evaluated values and partial
//                derivatives); first derivatives only, it does not nest.
#include "ADBackendNestedFad.h"
#include "Marmot/MarmotAutomaticDifferentiationScalarMode.h"
#include "Marmot/MarmotAutomaticDifferentiationVectorMode.h"
#include <array>
#include <cstddef>

namespace ADBench {

  namespace MarmotDetail {
    struct Access {
      template < typename S >
      static auto& val( S& x )
      {
        return x.val();
      }
      template < typename S >
      static auto val( const S& x )
      {
        return x.val();
      }
      template < typename S >
      static auto& dx( S& x, std::size_t i )
      {
        return x.dx( i );
      }
      template < typename S >
      static auto dx( const S& x, std::size_t i )
      {
        return x.dx( i );
      }
      template < typename S >
      static S constant( double v )
      {
        return S( v );
      }
    };
  } // namespace MarmotDetail

  struct MarmotBackend : NestedFadOps< MarmotDetail::Access > {
    static constexpr const char* name           = "marmot";
    static constexpr bool        hasReverseMode = false;
    static constexpr bool        hasTaylor      = true;
    static constexpr bool        hasVectorMode  = true;
    static constexpr int         maxScalarOrder = 3;
    static constexpr int         maxVectorOrder = 1;

    template < std::size_t N >
    using HyperDual = Marmot::Math::Differentiation::Automatic::ScalarMode::HyperDual< N >;

    template < std::size_t N >
    using Taylor = Marmot::Math::Differentiation::Automatic::ScalarMode::Taylor< N >;

    template < std::size_t level, typename S >
    static void seed( S& x, double v )
    {
      if constexpr ( level == 1 )
        x.grad() = v;
      else
        seed< level - 1 >( x.val(), v );
    }

    template < std::size_t level, typename S >
    static double derivative( const S& x )
    {
      if constexpr ( level == 0 )
        return Marmot::Math::Differentiation::Automatic::ScalarMode::primalValue( x );
      else
        return derivative< level - 1 >( x.grad() );
    }

    template < typename S >
    static double value( const S& x )
    {
      if constexpr ( Marmot::Math::Differentiation::Automatic::VectorMode::isDual< S > )
        return x.val();
      else
        return Marmot::Math::Differentiation::Automatic::ScalarMode::primalValue( x );
    }

    template < std::size_t N, typename F >
    static std::array< double, N + 1 > taylorDerivatives( const F& f, double x )
    {
      Taylor< N > x_( x );
      x_[1]                         = 1.0;
      const Taylor< N >           y = f( x_ );
      std::array< double, N + 1 > out;
      double                      factorial = 1.0;
      for ( std::size_t k = 0; k <= N; ++k ) {
        factorial *= k ? double( k ) : 1.0;
        out[k] = y[k] * factorial;
      }
      return out;
    }

    template < std::size_t N, std::size_t depth >
      requires( depth == 1 )
    using VecDual = Marmot::Math::Differentiation::Automatic::VectorMode::Dual< N >;

    template < typename S >
    static S constant( double v )
    {
      return S( v );
    }
  };

} // namespace ADBench
