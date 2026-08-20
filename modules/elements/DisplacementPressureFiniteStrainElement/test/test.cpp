#include "Marmot/DisplacementPressureFiniteStrainElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotTesting.h"
#include <array>
#include <cmath>
#include <vector>

using namespace Marmot::Testing;
using namespace Marmot::Elements;

using Element = DisplacementPressureFiniteStrainElement< 10, 4 >;

namespace {

  // Tetra10 node coordinates chosen identical to the natural coordinates (isoparametric identity
  // mapping), i.e. a valid, non-degenerate reference tetrahedron with corner nodes 0-3 followed by
  // the 6 edge-midpoint nodes, matching Marmot's Tetra10 natural-coordinate node numbering.
  std::vector< double > tetra10Coordinates()
  {
    return {
      0., 1., 0., // node 0 (corner)
      0., 0., 1., // node 1 (corner)
      0., 0., 0., // node 2 (corner)
      1., 0., 0., // node 3 (corner)
      0., .5, .5, // node 4
      0., 0., .5, // node 5
      0., .5, 0., // node 6
      .5, .5, 0., // node 7
      .5, 0., .5, // node 8
      .5, 0., 0., // node 9
    };
  }

  Element buildElement( std::vector< double >& stateVarsStorage )
  {
    using namespace Marmot::FiniteElement::Quadrature;

    Element element( 1, FullIntegration );

    static std::vector< double > coords = tetra10Coordinates();
    element.assignNodeCoordinates( coords.data() );
    element.initializeYourself();

    // 2-term Ogden material: nTerms=2, mu={1500,300}, alpha={2,-2}
    static std::array< double, 5 > matProps = { 2., 1500., 300., 2., -2. };
    MarmotMaterialSection          section( "OGDEN", matProps.data(), matProps.size() );
    element.assignProperty( section );

    const int nStateVars = element.getNumberOfRequiredStateVars();
    stateVarsStorage.assign( nStateVars, 0.0 );
    element.assignStateVars( stateVarsStorage.data(), nStateVars );

    return element;
  }

} // namespace

// Test E-1: at the undeformed configuration with zero pressure, the residual must vanish
// everywhere (zero isochoric stress, and J=1 exactly satisfies the incompressibility residual).
void testZeroStateResponse()
{
  std::vector< double > stateVars;
  Element               element = buildElement( stateVars );

  constexpr int n = Element::sizeLoadVector;

  std::vector< double > QTotal( n, 0.0 );
  std::vector< double > dQ( n, 0.0 );
  std::vector< double > Pint( n, 0.0 );
  std::vector< double > Ke( n * n, 0.0 );

  element.computeKernels( QTotal.data(), dQ.data(), Pint.data(), Ke.data(), 0.0, 1.0 );

  for ( int i = 0; i < n; i++ )
    throwExceptionOnFailure( std::abs( Pint[i] ) < 1e-10,
                             "E-1: Residual must vanish at the undeformed, zero-pressure state (dof " +
                               std::to_string( i ) + ")" );
}

// Test E-2: the analytic tangent Ke returned by computeKernels must match a central-difference
// approximation of the residual Pint(QTotal), for an arbitrary deformed and pressurized state.
// This exercises the full u-u, u-p and p-u coupling blocks derived for the mixed formulation.
void testFiniteDifferenceTangent()
{
  std::vector< double > stateVars;
  Element               element = buildElement( stateVars );

  constexpr int n = Element::sizeLoadVector;

  std::vector< double > QTotal( n, 0.0 );
  for ( int i = 0; i < Element::sizeDoFU; i++ )
    QTotal[i] = 0.02 * std::sin( 0.3 * i + 1.0 );
  for ( int i = 0; i < Element::sizeDoFP; i++ )
    QTotal[Element::sizeDoFU + i] = 50.0 * std::cos( 0.2 * i );

  std::vector< double > dQ( n, 0.0 );
  std::vector< double > Pint( n, 0.0 );
  std::vector< double > Ke( n * n, 0.0 );

  element.computeKernels( QTotal.data(), dQ.data(), Pint.data(), Ke.data(), 0.0, 1.0 );

  const Eigen::Map< const Eigen::Matrix< double, n, n > > KeMap( Ke.data() );

  const double h = 1e-6;
  for ( int col = 0; col < n; col++ ) {
    std::vector< double > QP = QTotal;
    std::vector< double > QM = QTotal;
    QP[col] += h;
    QM[col] -= h;

    std::vector< double > PintP( n, 0.0 ), PintM( n, 0.0 );
    std::vector< double > KeDummy( n * n, 0.0 );

    element.computeKernels( QP.data(), dQ.data(), PintP.data(), KeDummy.data(), 0.0, 1.0 );
    element.computeKernels( QM.data(), dQ.data(), PintM.data(), KeDummy.data(), 0.0, 1.0 );

    for ( int row = 0; row < n; row++ ) {
      const double dPint_dQ_numeric = ( PintP[row] - PintM[row] ) / ( 2 * h );

      throwExceptionOnFailure( std::abs( dPint_dQ_numeric - KeMap( row, col ) ) < 1e-2,
                               "E-2: Tangent mismatch at row=" + std::to_string( row ) +
                                 " col=" + std::to_string( col ) + " (numeric=" + std::to_string( dPint_dQ_numeric ) +
                                 ", analytic=" + std::to_string( KeMap( row, col ) ) + ")" );
    }
  }
}

int main()
{
  auto tests = std::vector< std::function< void() > >{ testZeroStateResponse, testFiniteDifferenceTangent };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
