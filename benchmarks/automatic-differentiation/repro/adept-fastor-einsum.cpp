// Minimal reproduction: Adept-2 records Fastor einsum contractions with correct values but wrong Jacobians.
// Build like the benchmark (see README), linking libadept and MarmotTensor.cpp.
#include "ADBackendAdept.h"
#include "ADBackendAutodiff.h"
#include "ADSuiteImpl.h"
#include <cstdio>
using namespace ADBench;
using namespace Fastor;
enum { i_, j_, I_, J_ };
using iI    = Index< i_, I_ >;
using IJ    = Index< I_, J_ >;
using jJ    = Index< j_, J_ >;
using iJ    = Index< i_, J_ >;
using to_ij = OIndex< i_, j_ >;
template < typename S >
Tensor< S, 3, 3 > f( const Tensor< S, 3, 3 >& F_, int k )
{
  Tensor< S, 3, 3 > P = F_ * F_;                     // elementwise, some function of F
  if ( k == 0 )
    return einsum< iI, IJ, to_ij >( F_, P );         // two operands
  if ( k == 1 )
    return einsum< iI, IJ, jJ, to_ij >( F_, P, F_ ); // three operands, F twice (Marmot's KirchhoffStressFromPK2)
  Tensor< S, 3, 3 > G = F_ + 0.;                     // three operands, distinct objects
  return einsum< iI, IJ, jJ, to_ij >( F_, P, G );
}
int main()
{
  auto        F       = testDeformationGradient();
  const char* label[] = { "einsum(F,P) 2 operands",
                          "einsum(F,P,F) 3 operands, F twice",
                          "einsum(F,P,G) 3 operands, G=F+0" };
  for ( int k = 0; k < 3; k++ ) {
    std::function< Tensor< AutodiffBackend::HyperDual< 1 >, 3, 3 >(
      const Tensor< AutodiffBackend::HyperDual< 1 >, 3, 3 >& ) >
      g                  = [k]( const auto& X ) { return f( X, k ); };
    auto [v1, J1]        = dF_dT< AutodiffBackend >( g, F );
    using S              = adept::adouble;
    auto&             st = AdeptBackend::stack();
    Tensor< S, 3, 3 > F_;
    for ( int a = 0; a < 9; a++ )
      F_.data()[a] = F.data()[a];
    st.new_recording();
    Tensor< S, 3, 3 > o = f( F_, k );
    st.independent( F_.data(), 9 );
    st.dependent( o.data(), 9 );
    double jac[81];
    st.jacobian( jac );
    double dj = 0, dv = 0;
    for ( int j = 0; j < 9; j++ ) {
      dv = std::max( dv, std::abs( v1.data()[j] - adept::value( o.data()[j] ) ) );
      for ( int a = 0; a < 9; a++ )
        dj = std::max( dj, std::abs( J1.data()[j * 9 + a] - jac[j + 9 * a] ) );
    }
    std::printf( "%-36s max |dvalue| = %.1e  max |dJ| = %.3g\n", label[k], dv, dj );
  }
}
