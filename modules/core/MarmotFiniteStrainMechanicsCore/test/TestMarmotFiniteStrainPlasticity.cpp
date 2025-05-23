#include "Fastor/Fastor.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotFiniteStrainPlasticity.h"
#include "Marmot/MarmotTesting.h"

using namespace Marmot::Testing;
using namespace Marmot::FastorStandardTensors;

Tensor33d make_dGp()
{
  // create example matrix
  double    a1 = std::log( 1.05 );
  double    a2 = std::log( 1.1 );
  Tensor33d dGp;
  dGp( 0, 0 ) = a1;
  dGp( 0, 1 ) = 0.;
  dGp( 0, 2 ) = 0.;
  dGp( 1, 0 ) = 0.;
  dGp( 1, 1 ) = a2;
  dGp( 1, 2 ) = a2;
  dGp( 2, 0 ) = 0.;
  dGp( 2, 1 ) = 0.;
  dGp( 2, 2 ) = a2;
  return dGp;
}

auto testExponentialMap()
{
  using namespace Marmot::ContinuumMechanics::FiniteStrain::Plasticity::FlowIntegration;
  Tensor33d dGp = make_dGp();
  // compute exponential map
  Tensor33d exMap = exponentialMap( dGp );
  // create reference solution
  double    a23 = 0.10484119778475742;
  Tensor33d exMapexpect;
  exMapexpect( 0, 0 ) = 1.05;
  exMapexpect( 0, 1 ) = 0;
  exMapexpect( 0, 2 ) = 0;
  exMapexpect( 1, 0 ) = 0;
  exMapexpect( 1, 1 ) = 1.1;
  exMapexpect( 1, 2 ) = 0.;
  exMapexpect( 2, 0 ) = 0;
  exMapexpect( 2, 1 ) = a23;
  exMapexpect( 2, 2 ) = 1.1;
  // check if expected and calculated values are the same
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ ) {
      throwExceptionOnFailure( checkIfEqual( exMap( i, j ), exMapexpect( i, j ), 1e-15 ),
                               MakeString() << __PRETTY_FUNCTION__ << " exponential mapping failed." );
    }
}

auto testExplicitIntegration()
{
  using namespace Marmot::ContinuumMechanics::FiniteStrain::Plasticity::FlowIntegration::FirstOrderDerived;
  Tensor33d dGp = make_dGp();
  // make explicit integration
  std::pair< Tensor33d, Tensor3333d > res     = explicitIntegration( dGp );
  Tensor33d                           dFp     = res.first;
  Tensor3333d                         dFp_dGp = res.second;
  // check if expected and calculated values are the same
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ ) {
      double dFpexpect = ( i == j ? 1 : 0 ) + dGp( j, i );
      throwExceptionOnFailure( checkIfEqual( dFp( i, j ), dFpexpect, 1e-15 ),
                               MakeString() << __PRETTY_FUNCTION__ << " explicit integration dFp failed." );
      for ( int k = 0; k < 3; k++ )
        for ( int l = 0; l < 3; l++ ) {
          double dFp_dGpexpect = ( i == j ? 1 : 0 ) * ( k == l ? 1 : 0 );
          throwExceptionOnFailure( checkIfEqual( dFp_dGp( i, k, l, j ), dFp_dGpexpect, 1e-15 ),
                                   MakeString() << __PRETTY_FUNCTION__ << " explicit integration dFp_dGp failed." );
        }
    }
}

auto testExponentialMapAndDerivative()
{
  using namespace Marmot::ContinuumMechanics::FiniteStrain::Plasticity::FlowIntegration::FirstOrderDerived;
  Tensor33d dGp = make_dGp();
  // compute exponential map
  std::pair< Tensor33d, Tensor3333d > exMapTot = exponentialMap( dGp );
  Tensor33d                           exMap    = exMapTot.first;
  Tensor3333d                         DexMap   = exMapTot.second;
  // create reference solution for exponential map
  double    a23 = 0.10484119778475742;
  Tensor33d exMapexpect;
  exMapexpect( 0, 0 ) = 1.05;
  exMapexpect( 0, 1 ) = 0;
  exMapexpect( 0, 2 ) = 0;
  exMapexpect( 1, 0 ) = 0;
  exMapexpect( 1, 1 ) = 1.1;
  exMapexpect( 1, 2 ) = 0.;
  exMapexpect( 2, 0 ) = 0;
  exMapexpect( 2, 1 ) = a23;
  exMapexpect( 2, 2 ) = 1.1;
  // create reference solution for exponantial map derivative
  double      a1 = 1.074806173592447;
  double      a2 = 0.051617096256127606;
  double      a3 = 0.05242059889237864;
  Tensor3333d DexMapexpect;
  DexMapexpect( 0, 0, 0, 0 ) = 1.05;
  DexMapexpect( 0, 1, 0, 1 ) = a1;
  DexMapexpect( 1, 0, 1, 0 ) = a1;
  DexMapexpect( 0, 2, 0, 2 ) = a1;
  DexMapexpect( 2, 0, 2, 0 ) = a1;
  DexMapexpect( 0, 1, 0, 2 ) = a2;
  DexMapexpect( 2, 0, 1, 0 ) = a2;
  DexMapexpect( 1, 1, 1, 1 ) = 1.1;
  DexMapexpect( 1, 2, 1, 2 ) = 1.1;
  DexMapexpect( 2, 1, 2, 1 ) = 1.1;
  DexMapexpect( 2, 2, 2, 2 ) = 1.1;
  DexMapexpect( 1, 1, 1, 2 ) = a3;
  DexMapexpect( 2, 1, 1, 1 ) = a3;
  DexMapexpect( 2, 1, 2, 2 ) = a3;
  DexMapexpect( 2, 2, 1, 2 ) = a3;
  DexMapexpect( 2, 1, 1, 2 ) = 0.0016654055686274117;
  // check if expected and calculated values are the same
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ ) {
      throwExceptionOnFailure( checkIfEqual( exMap( i, j ), exMapexpect( i, j ), 1e-15 ),
                               MakeString() << __PRETTY_FUNCTION__ << " exponential mapping failed." );
      for ( int k = 0; k < 3; k++ )
        for ( int l = 0; l < 3; l++ ) {
          // check derivative
          DexMapexpect( i, j, l, k ) = ( DexMapexpect( i, j, l, k ) == 1 ) ? 0.0 : DexMapexpect( i, j, l, k );
          throwExceptionOnFailure( checkIfEqual( DexMap( i, j, k, l ), DexMapexpect( i, j, l, k ), 1e-15 ),
                                   MakeString()
                                     << __PRETTY_FUNCTION__ << " derivative of exponential mapping failed at " << j << i
                                     << k << l );
        }
    }
}

int main()
{
  auto tests = std::vector< std::function< void() > >{ testExponentialMap,
                                                       testExplicitIntegration,
                                                       testExponentialMapAndDerivative };

  executeTestsAndCollectExceptions( tests );
  return 0;
}
