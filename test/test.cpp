#include "Marmot/GradientEnhancedFiniteStrainMaterialPoint.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotMPMLibrary.h"
#include "Marmot/MarmotTesting.h"
#include <cmath>
#include <functional>
#include <memory>
#include <utility>
#include <vector>

using namespace Marmot;
using namespace Marmot::MaterialPoints;
using namespace Marmot::Testing;

namespace {

  // GRADIENTENHANCEDCOMPRESSIBLENEOHOOKEDAMAGE: K, G, kappa0, kappaF, l, rho
  const std::vector< double > matProps = { 3500., 1500., 1e-3, 1e-2, 0.3, 2.0 };

  template < int nDim >
  struct Setup {
    std::unique_ptr< MarmotMaterialPoint > mp;
    std::vector< double >                  stateVars;

    Setup()
    {
      const std::string name = nDim == 2 ? "GradientEnhancedFiniteStrain/PlaneStrain"
                                         : "GradientEnhancedFiniteStrain/3D";
      const double      x[3] = { 0.1, 0.2, 0.3 };
      mp.reset( MarmotLibrary::MarmotMaterialPointFactory::createMaterialPoint( name, 1, x, nDim, 0.5 ) );
      mp->assignMaterial(
        MarmotMaterialSection( "GRADIENTENHANCEDCOMPRESSIBLENEOHOOKEDAMAGE", matProps.data(), matProps.size() ) );
      stateVars.assign( mp->getNumberOfRequiredStateVars(), 0.0 );
      mp->assignStateVars( stateVars.data(), stateVars.size() );
      mp->initializeYourself();
    }

    GradientEnhancedFiniteStrainMaterialPoint< nDim >& point()
    {
      return *dynamic_cast< GradientEnhancedFiniteStrainMaterialPoint< nDim >* >( mp.get() );
    }

    // one increment with displacement gradient increment H and nonlocal increment dn, not committed
    void compute( const Fastor::Tensor< double, nDim, nDim >& H, double dn )
    {
      mp->prepareYourself( 1.0, 1.0 );
      point().incrementDeformation( Fastor::Tensor< double, nDim >( 0.0 ), H, dn );
      mp->computeYourself( 1.0, 1.0 );
    }

    auto trial( const Fastor::Tensor< double, nDim, nDim >& H, double dn )
    {
      const auto backup = stateVars;
      compute( H, dn );
      auto r    = point().response;
      auto t    = point().tangents;
      stateVars = backup;
      return std::make_pair( r, t );
    }
  };

  template < int nDim >
  Fastor::Tensor< double, nDim, nDim > gradient( double scale )
  {
    Fastor::Tensor< double, nDim, nDim > H;
    for ( int i = 0; i < nDim; i++ )
      for ( int j = 0; j < nDim; j++ )
        H( i, j ) = scale * std::sin( 1.3 * ( nDim * i + j + 1 ) );
    return H;
  }

  // the tangents w.r.t. the INCREMENT of the deformation gradient, DeltaF = I + H, against central differences
  template < int nDim >
  void checkTangents( Setup< nDim >& s, const Fastor::Tensor< double, nDim, nDim >& H0, double dn0 )
  {
    const auto [r0, t0] = s.trial( H0, dn0 );

    const double h = 1e-7;
    for ( int k = 0; k < nDim; k++ )
      for ( int l = 0; l < nDim; l++ ) {
        auto Hp = H0, Hm = H0;
        Hp( k, l ) += h;
        Hm( k, l ) -= h;
        const auto rp = s.trial( Hp, dn0 ).first;
        const auto rm = s.trial( Hm, dn0 ).first;
        for ( int i = 0; i < nDim; i++ )
          for ( int j = 0; j < nDim; j++ ) {
            const double num = ( rp.S( i, j ) - rm.S( i, j ) ) / ( 2 * h );
            throwExceptionOnFailure( std::abs( t0.dS_dDeltaF( i, j, k, l ) - num ) < 1e-5 * ( 1. + std::abs( num ) ),
                                     MakeString() << nDim << "D: dS_dDeltaF inconsistent" );
          }
        const double numL = ( rp.dL - rm.dL ) / ( 2 * h );
        throwExceptionOnFailure( std::abs( t0.dL_dDeltaF( k, l ) - numL ) < 1e-6 * ( 1. + std::abs( numL ) ),
                                 MakeString() << nDim << "D: dL_dDeltaF inconsistent" );
      }

    const double hN = 1e-9;
    const auto   rp = s.trial( H0, dn0 + hN ).first;
    const auto   rm = s.trial( H0, dn0 - hN ).first;
    for ( int i = 0; i < nDim; i++ )
      for ( int j = 0; j < nDim; j++ ) {
        const double num = ( rp.S( i, j ) - rm.S( i, j ) ) / ( 2 * hN );
        throwExceptionOnFailure( std::abs( t0.dS_dN( i, j ) - num ) < 1e-4 * ( 1. + std::abs( num ) ),
                                 MakeString() << nDim << "D: dS_dN inconsistent" );
      }
  }

  template < int nDim >
  void checkMaterialPoint()
  {
    // from the undeformed state
    {
      Setup< nDim > s;
      throwExceptionOnFailure( s.mp->getDensityUndeformed() == matProps[5],
                               "density must be known after initialization" );
      checkTangents< nDim >( s, gradient< nDim >( 0.03 ), 3e-3 );
    }

    // after an accepted increment: F = DeltaF . F_n, and the tangents chain through F_n
    {
      Setup< nDim > s;
      const auto    H1 = gradient< nDim >( 0.03 );
      s.compute( H1, 3e-3 );
      s.mp->acceptStateAndPosition();

      const auto F = s.mp->getStateView( "deformation gradient" ); // Fastor, row major
      for ( int i = 0; i < nDim; i++ )
        for ( int j = 0; j < nDim; j++ )
          throwExceptionOnFailure( std::abs( F.stateLocation[3 * i + j] - ( ( i == j ) + H1( i, j ) ) ) < 1e-14,
                                   MakeString() << nDim << "D: accepted deformation gradient" );

      throwExceptionOnFailure( std::abs( *s.mp->getStateView( "nonlocal damage" ).stateLocation - 3e-3 ) < 1e-15,
                               "the nonlocal field must accumulate" );

      checkTangents< nDim >( s, gradient< nDim >( 0.01 ), 1e-3 );
    }
  }

} // namespace

void testPlaneStrainMaterialPoint()
{
  checkMaterialPoint< 2 >();
}

void test3DMaterialPoint()
{
  checkMaterialPoint< 3 >();
}

int main()
{
  auto testFunctions = std::vector< std::function< void() > >{ testPlaneStrainMaterialPoint, test3DMaterialPoint };
  executeTestsAndCollectExceptions( testFunctions );
  return 0;
}
