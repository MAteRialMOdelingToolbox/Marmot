#pragma once
// Backend adapter: TinyAD (https://github.com/patr-schm/TinyAD, MIT). TinyAD::Scalar< k, double, with_hessian >
// carries the dense gradient (and optionally the dense Hessian) with respect to k variables, i.e. it is a vector
// mode of first or second order. It does not nest, so scalar mode is first order only (no mixed hyper-duals), there
// are no third derivatives and no Taylor mode.
#include <Eigen/Core>
#include <TinyAD/Scalar.hh>
#include <cstddef>

// TinyAD::Scalar has no unary plus, which Fastor's inverse() uses (found by ADL)
namespace TinyAD {
  template < int k, typename PassiveT, bool with_hessian, int hrs, int hcs, int hr, int hc >
  Scalar< k, PassiveT, with_hessian, hrs, hcs, hr, hc > operator+(
    const Scalar< k, PassiveT, with_hessian, hrs, hcs, hr, hc >& x )
  {
    return x;
  }
} // namespace TinyAD

namespace ADBench {

  struct TinyADBackend {
    static constexpr const char* name           = "tinyad";
    static constexpr bool        hasReverseMode = false;
    static constexpr bool        hasTaylor      = false;
    static constexpr int         maxScalarOrder = 1;
    static constexpr bool        hasVectorMode  = true;
    static constexpr int         maxVectorOrder = 2;

    // only N = 1 is used (maxScalarOrder = 1)
    template < std::size_t N >
    using HyperDual = TinyAD::Scalar< 1, double, false >;

    template < std::size_t N, std::size_t depth >
    using VecDual = TinyAD::Scalar< int( N ), double, ( depth >= 2 ) >;

    template < typename S >
    static S constant( double v )
    {
      return S( v );
    }

    template < std::size_t level, typename S >
    static void seed( S& x, double v )
    {
      static_assert( level == 1, "TinyAD::Scalar does not nest" );
      x.grad( 0 ) = v;
    }

    template < std::size_t level, typename S >
    static double derivative( const S& x )
    {
      static_assert( level == 1, "TinyAD::Scalar does not nest" );
      return x.grad( 0 );
    }

    template < typename S >
    static double value( const S& x )
    {
      return x.val;
    }

    template < typename S >
    static void seedDirection( S& x, std::size_t i )
    {
      x.grad( Eigen::Index( i ) ) = 1.0;
    }

    template < typename S >
    static double vderivative( const S& x, std::size_t i )
    {
      return x.grad( Eigen::Index( i ) );
    }

    template < typename S >
    static double vderivative( const S& x, std::size_t i, std::size_t j )
    {
      return x.Hess( Eigen::Index( i ), Eigen::Index( j ) );
    }
  };

} // namespace ADBench
