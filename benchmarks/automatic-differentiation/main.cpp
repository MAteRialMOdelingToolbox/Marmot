// Compares forward-mode AD backends (autodiff, Boost.Math, Sacado, Ceres Jet) on Marmot's AD code paths:
// first checks that every backend reproduces autodiff's derivatives, then times each kernel.
#include "ADSuite.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>

using namespace ADBench;

namespace {

  volatile double sink = 0;

  // median over `batches` of the mean time per call (ns), each batch running at least minBatchSeconds
  double timeCase( const Case& c, int batches, double minBatchSeconds )
  {
    using clock = std::chrono::steady_clock;
    long n      = 1;
    for ( ;; ) {
      auto t0 = clock::now();
      for ( long i = 0; i < n; i++ )
        sink = sink + c.evaluate()[0];
      if ( std::chrono::duration< double >( clock::now() - t0 ).count() > minBatchSeconds )
        break;
      n *= 2;
    }
    std::vector< double > samples;
    for ( int b = 0; b < batches; b++ ) {
      auto t0 = clock::now();
      for ( long i = 0; i < n; i++ )
        sink = sink + c.evaluate()[0];
      samples.push_back( std::chrono::duration< double, std::nano >( clock::now() - t0 ).count() / double( n ) );
    }
    std::sort( samples.begin(), samples.end() );
    return samples[samples.size() / 2];
  }

  double maxRelativeDifference( const std::vector< double >& a, const std::vector< double >& b )
  {
    if ( a.size() != b.size() )
      return INFINITY;
    double scale = 0, diff = 0;
    for ( double x : a )
      scale = std::max( scale, std::abs( x ) );
    for ( std::size_t i = 0; i < a.size(); i++ )
      diff = std::max( diff, std::abs( a[i] - b[i] ) );
    return diff / std::max( scale, 1e-300 );
  }

} // namespace

int main( int argc, char** argv )
{
  const int    batches  = argc > 1 ? std::atoi( argv[1] ) : 15;
  const double minBatch = argc > 2 ? std::atof( argv[2] ) : 0.05;

  std::vector< Suite > suites;
  for ( auto make : { makeAutodiffSuites, makeBoostSuites, makeSacadoSuites, makeCeresSuites } )
    for ( auto& s : make() )
      suites.push_back( std::move( s ) );

  const Suite& reference = suites.front(); // autodiff, scalar mode: the current Marmot implementation

  // id -> label -> (ns, relerr)
  std::map< std::string, std::map< std::string, std::pair< double, double > > > results;
  bool                                                                           ok = true;
  for ( const Case& ref : reference.cases ) {
    const std::vector< double > expected = ref.evaluate();
    for ( const Suite& s : suites )
      for ( const Case& c : s.cases )
        if ( c.id == ref.id ) {
          const double err = maxRelativeDifference( expected, c.evaluate() );
          ok               = ok && err < 1e-10;
          results[ref.id][s.label] = { timeCase( c, batches, minBatch ), err };
        }
  }

  auto printTable = [&]( const char* title, auto cell ) {
    std::printf( "\n%s\n%-42s", title, "kernel" );
    for ( const Suite& s : suites )
      std::printf( " %11s", s.label.c_str() );
    std::printf( "\n" );
    for ( const Case& ref : reference.cases ) {
      std::printf( "%-42s", ref.name.c_str() );
      for ( const Suite& s : suites ) {
        auto it = results[ref.id].find( s.label );
        if ( it == results[ref.id].end() )
          std::printf( " %11s", "n/a" );
        else
          std::printf( " %11s", cell( it->second, results[ref.id][reference.label] ).c_str() );
      }
      std::printf( "\n" );
    }
  };

  char buf[64];
  printTable( "Median time per kernel call [ns]", [&]( auto r, auto ) {
    std::snprintf( buf, sizeof buf, "%.0f", r.first );
    return std::string( buf );
  } );
  printTable( "Time relative to autodiff (current Marmot; < 1 is faster)", [&]( auto r, auto ref ) {
    std::snprintf( buf, sizeof buf, "%.2f", r.first / ref.first );
    return std::string( buf );
  } );
  printTable( "Max relative deviation from autodiff", [&]( auto r, auto ) {
    std::snprintf( buf, sizeof buf, "%.1e", r.second );
    return std::string( buf );
  } );
  std::printf( "\nResults %s.\n", ok ? "agree (all max relerr < 1e-10)" : "DISAGREE" );
  return ok ? 0 : 1;
}
