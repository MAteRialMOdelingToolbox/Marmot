// Per-operation cost of the AD number types (ns per op), to locate where the backends differ.
#include "ADBackendAutodiff.h"
#include "ADBackendBoost.h"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>

using namespace ADBench;

namespace {

  template < typename S >
  double valueOf( const S& x )
  {
    return *reinterpret_cast< const double* >( &x );
  }

  // time `op` applied in a dependent chain of length n, so the result of one step feeds the next
  template < typename S, typename Op >
  double timeOp( Op op, S x0, long n = 2000000 )
  {
    S                  x  = x0;
    volatile double    sink;
    auto               t0 = std::chrono::steady_clock::now();
    for ( long i = 0; i < n; i++ ) {
      x = op( x );
      asm volatile( "" : "+m"( x ) : : "memory" ); // DoNotOptimize: materialize every coefficient each step
    }
    auto t1 = std::chrono::steady_clock::now();
    sink    = valueOf( x );
    (void)sink;
    return std::chrono::duration< double, std::nano >( t1 - t0 ).count() / double( n );
  }

  template < typename S >
  void row( const char* type, const char* backend, S x0 )
  {
    S c( 1.0000001 );
    asm volatile( "" : "+m"( c ) );
    double  tMul  = timeOp< S >( [&]( const S& x ) -> S { return x * c; }, x0 );
    double  tDiv  = timeOp< S >( [&]( const S& x ) -> S { return x / c; }, x0 );
    double  tSqrt = timeOp< S >( [&]( const S& x ) -> S { return sqrt( x * x ); }, x0 );
    double  tLog  = timeOp< S >( [&]( const S& x ) -> S { return exp( log( x ) ); }, x0 );
    double  tPow  = timeOp< S >( [&]( const S& x ) -> S { return pow( pow( x, 0.75 ), 1. / 0.75 ); }, x0 );
    std::printf( "%-10s %-9s %8.2f %8.2f %10.2f %10.2f %10.2f\n", type, backend, tMul, tDiv, tSqrt, tLog, tPow );
  }

  template < std::size_t N >
  void hyperDualRows()
  {
    const std::string t = "dual" + std::to_string( N );
    AutodiffBackend::HyperDual< N > a( 1.3 );
    BoostBackend::HyperDual< N >    b( 1.3 );
    AutodiffBackend::seed< 1 >( a, 1.0 );
    BoostBackend::seed< 1 >( b, 1.0 );
    if constexpr ( N > 1 ) {
      AutodiffBackend::seed< N >( a, 1.0 );
      BoostBackend::seed< N >( b, 1.0 );
    }
    row( t.c_str(), "autodiff", a );
    row( t.c_str(), "boost", b );
  }

  template < std::size_t N >
  void taylorRows()
  {
    const std::string                    t = "taylor" + std::to_string( N );
    AutodiffBackend::Taylor< N >         a( 1.3 );
    a[1] = 1.0;
    auto b = boost::math::differentiation::make_fvar< double, N >( 1.3 );
    row( t.c_str(), "autodiff", a );
    row( t.c_str(), "boost", b );
  }

} // namespace

int main()
{
  std::printf( "%-10s %-9s %8s %8s %10s %10s %10s   [ns per chained op]\n", "type", "backend", "mul", "div",
               "sqrt(x*x)", "exp(log)", "pow(pow)" );
  hyperDualRows< 1 >();
  hyperDualRows< 2 >();
  hyperDualRows< 3 >();
  taylorRows< 2 >();
  taylorRows< 5 >();
  taylorRows< 7 >();
}
