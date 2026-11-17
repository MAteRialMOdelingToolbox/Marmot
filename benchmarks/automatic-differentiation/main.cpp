// Compares the autodiff and Boost.Math forward-mode AD backends on Marmot's AD code paths:
// first checks that both produce the same derivatives, then times each kernel.
#include "ADSuite.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

using namespace ADBench;

namespace {

  volatile double sink = 0;

  // median over `batches` of the mean time per call (ns), each batch running long enough to be measurable
  double timeCase( const Case& c, int batches, double minBatchSeconds )
  {
    using clock = std::chrono::steady_clock;
    // calibrate the number of calls per batch
    long   n = 1;
    double t = 0;
    for ( ;; ) {
      auto t0 = clock::now();
      for ( long i = 0; i < n; i++ )
        sink = sink + c.evaluate()[0];
      t = std::chrono::duration< double >( clock::now() - t0 ).count();
      if ( t > minBatchSeconds )
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
    double scale = 0;
    for ( double x : a )
      scale = std::max( scale, std::abs( x ) );
    double diff = 0;
    for ( std::size_t i = 0; i < a.size(); i++ )
      diff = std::max( diff, std::abs( a[i] - b[i] ) );
    return diff / std::max( scale, 1e-300 );
  }

} // namespace

int main( int argc, char** argv )
{
  const int    batches = argc > 1 ? std::atoi( argv[1] ) : 15;
  const double minBatch = argc > 2 ? std::atof( argv[2] ) : 0.05;

  const Suite ad = makeAutodiffSuite();
  const Suite bm = makeBoostSuite();

  bool ok = true;
  std::printf( "%-46s %12s %12s %8s %10s\n", "kernel", "autodiff[ns]", "boost[ns]", "ratio", "max relerr" );
  for ( std::size_t i = 0; i < ad.cases.size(); i++ ) {
    const Case& a   = ad.cases[i];
    const Case& b   = bm.cases[i];
    const double err = maxRelativeDifference( a.evaluate(), b.evaluate() );
    ok               = ok && err < 1e-10;
    const double ta  = timeCase( a, batches, minBatch );
    const double tb  = timeCase( b, batches, minBatch );
    std::printf( "%-46s %12.1f %12.1f %8.2f %10.1e\n", a.name.c_str(), ta, tb, tb / ta, err );
  }
  std::printf( "\nratio = boost / autodiff (lower is better for boost). Results %s.\n",
               ok ? "agree (max relerr < 1e-10)" : "DISAGREE" );
  return ok ? 0 : 1;
}
