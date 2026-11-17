#pragma once
// Backend adapter: Boost.Math forward-mode autodiff (boost/math/differentiation/autodiff.hpp).
//
// Same vocabulary as ADBackendAutodiff.h, so the drivers and kernels are shared verbatim.
//
// Mapping:
//   HyperDual<1> = fvar<double, 1>                       (autodiff::dual)
//   HyperDual<N> = fvar<HyperDual<N-1>, 1>               (autodiff::HigherOrderDual<N>)
//   Taylor<N>    = fvar<double, N>                       (autodiff::Real<N>)
//
// fvar keeps its coefficients in a private std::array and offers no mutable accessor, so seeding uses the
// coefficient storage directly. The nested layout is identical to autodiff's nested Dual<T, G>
// (2^N doubles, level k epsilon at index 2^(N-k)); this is guarded by static_asserts below.
#include <array>
#include <boost/math/differentiation/autodiff.hpp>
#include <cstddef>
#include <type_traits>

#include <Eigen/Core>

namespace ADBench::BoostDetail {

  namespace bmd = boost::math::differentiation;

  namespace detail {
    template < std::size_t N >
    struct HyperDualImpl {
      using type = bmd::detail::fvar< typename HyperDualImpl< N - 1 >::type, 1 >;
    };
    template <>
    struct HyperDualImpl< 1 > {
      using type = bmd::detail::fvar< double, 1 >;
    };

    template < typename S >
    struct Depth : std::integral_constant< std::size_t, 0 > {};
    template < typename R, std::size_t O >
    struct Depth< bmd::detail::fvar< R, O > > : std::integral_constant< std::size_t, 1 + Depth< R >::value > {};

    template < typename S >
    inline double* coefficients( S& x )
    {
      static_assert( std::is_standard_layout_v< S > );
      static_assert( sizeof( S ) == ( std::size_t( 1 ) << Depth< S >::value ) * sizeof( double ),
                     "fvar layout differs from the expected dense nested hyper-dual layout" );
      return reinterpret_cast< double* >( &x );
    }
    template < typename S >
    inline const double* coefficients( const S& x )
    {
      return coefficients( const_cast< S& >( x ) );
    }
  } // namespace detail
} // namespace ADBench::BoostDetail

namespace ADBench {

  struct BoostBackend {
    static constexpr const char* name = "boost";
    static constexpr bool        hasTaylor     = true;
    static constexpr bool        hasVectorMode = false;
    static constexpr int         maxOrder      = 3;

    template < typename S >
    static S constant( double v )
    {
      return S( v );
    }


  template < std::size_t N >
  using HyperDual = typename BoostDetail::detail::HyperDualImpl< N >::type;

  template < std::size_t N >
  using Taylor = boost::math::differentiation::detail::fvar< double, N >;

  template < std::size_t level, typename S >
  static void seed( S& x, double v )
  {
    constexpr std::size_t N = BoostDetail::detail::Depth< S >::value;
    static_assert( level >= 1 && level <= N );
    BoostDetail::detail::coefficients( x )[std::size_t( 1 ) << ( N - level )] = v;
  }

  template < std::size_t level, typename S >
  static double derivative( const S& x )
  {
    constexpr std::size_t N = BoostDetail::detail::Depth< S >::value;
    static_assert( level >= 1 && level <= N );
    std::size_t idx = 0;
    for ( std::size_t i = 1; i <= level; ++i )
      idx += std::size_t( 1 ) << ( N - i );
    return BoostDetail::detail::coefficients( x )[idx];
  }

  template < typename S >
  static double value( const S& x )
  {
    return BoostDetail::detail::coefficients( x )[0];
  }

  template < std::size_t N, typename F >
  static std::array< double, N + 1 > taylorDerivatives( const F& f, double x )
  {
    const Taylor< N >           y = f( boost::math::differentiation::make_fvar< double, N >( x ) );
    std::array< double, N + 1 > out;
    for ( std::size_t i = 0; i <= N; ++i )
      out[i] = y.derivative( i );
    return out;
  }

  };

} // namespace ADBench

// Eigen needs NumTraits for a custom scalar, and ScalarBinaryOpTraits to mix it with double
// (autodiff ships the equivalent in autodiff/forward/dual/eigen.hpp).
namespace Eigen {
  template < typename R, std::size_t O >
  struct NumTraits< boost::math::differentiation::detail::fvar< R, O > > : NumTraits< double > {
    using Real       = boost::math::differentiation::detail::fvar< R, O >;
    using NonInteger = Real;
    using Literal    = double; // not fvar: its unconstrained converting ctor would make Eigen treat any matrix as a scalar
    using Nested     = Real;
    enum { IsComplex = 0, IsInteger = 0, IsSigned = 1, RequireInitialization = 1, ReadCost = 1, AddCost = 3, MulCost = 3 };
  };

  template < typename R, std::size_t O, typename BinOp >
  struct ScalarBinaryOpTraits< boost::math::differentiation::detail::fvar< R, O >, double, BinOp > {
    using ReturnType = boost::math::differentiation::detail::fvar< R, O >;
  };
  template < typename R, std::size_t O, typename BinOp >
  struct ScalarBinaryOpTraits< double, boost::math::differentiation::detail::fvar< R, O >, BinOp > {
    using ReturnType = boost::math::differentiation::detail::fvar< R, O >;
  };

  namespace internal {
    // fvar's converting constructor accepts any type, so Eigen's scalar-promotion fallback (is_convertible<T, Scalar>)
    // would treat a whole matrix as a scalar and make `A * B` ambiguous. Only arithmetic types promote to fvar.
    template < typename S, bool >
    struct PromoteToFvarIf {};
    template < typename S >
    struct PromoteToFvarIf< S, true > {
      using type = S;
    };
    template < typename R, std::size_t O, typename T >
    struct promote_scalar_arg< boost::math::differentiation::detail::fvar< R, O >, T, false >
      : PromoteToFvarIf< boost::math::differentiation::detail::fvar< R, O >, std::is_arithmetic_v< T > > {};
  } // namespace internal
} // namespace Eigen
