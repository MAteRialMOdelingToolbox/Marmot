#include "Marmot/DisplacementFiniteStrainFBarElement.h"
#include "Marmot/DisplacementFiniteStrainULElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotTesting.h"
#include <cmath>
#include <vector>

using namespace Marmot::Testing;
using namespace Marmot::Elements;
using namespace Marmot::FiniteElement::Quadrature;

namespace {

  // Unit right tetrahedron, corner nodes matching Marmot's Tetra4 natural-coordinate node
  // numbering ( N1:(0,0,0), N2:(1,0,0), N3:(0,1,0), N4:(0,0,1) ).
  std::vector< double > tetra4Coordinates()
  {
    return {
      0.,
      0.,
      0., // node 0
      1.,
      0.,
      0., // node 1
      0.,
      1.,
      0., // node 2
      0.,
      0.,
      1., // node 3
    };
  }

  // Unit cube, corner nodes matching Marmot's Hexa8 natural-coordinate node numbering (bottom
  // face 0-1-2-3, top face 4-5-6-7, node i+4 directly above node i).
  std::vector< double > hexa8Coordinates()
  {
    return {
      0., 0., 0., // node 0
      1., 0., 0., // node 1
      1., 1., 0., // node 2
      0., 1., 0., // node 3
      0., 0., 1., // node 4
      1., 0., 1., // node 5
      1., 1., 1., // node 6
      0., 1., 1., // node 7
    };
  }

  // CompressibleNeoHooke: K (bulk modulus), G (shear modulus). A large K/G ratio is exactly the
  // near-incompressible regime the F-bar method targets.
  const std::vector< double > compressibleNeoHookeProps = { 1.0e4, 1.0 };

  template < class Element >
  Element buildElement( std::vector< double >& coordinates, std::vector< double >& stateVarsStorage )
  {
    Element element( 1, FullIntegration, Element::SectionType::Solid );

    element.assignNodeCoordinates( coordinates.data() );
    element.initializeYourself();

    MarmotMaterialSection section( "COMPRESSIBLENEOHOOKE",
                                   compressibleNeoHookeProps.data(),
                                   compressibleNeoHookeProps.size() );
    element.assignProperty( section );

    const int nStateVars = element.getNumberOfRequiredStateVars();
    stateVarsStorage.assign( nStateVars, 0.0 );
    element.assignStateVars( stateVarsStorage.data(), nStateVars );

    return element;
  }

} // namespace

// Test E-1: a fully-integrated Tetra4 has a single quadrature point coincident with the element
// centroid, so F_qp == F_0 identically (Tetra4's shape functions are linear, hence dNdX is
// constant throughout the element): alpha == 1 and the F-bar correction term vanishes exactly.
// The F-bar element must therefore reproduce the plain UL element's residual and tangent bit for
// bit (up to floating point round-off), for an arbitrary deformed state.
void testTetra4ReducesToStandardElement()
{
  using FBarElement = DisplacementFiniteStrainFBarElement< 3, 4 >;
  using ULElement   = DisplacementFiniteStrainULElement< 3, 4 >;

  auto coordsFBar = tetra4Coordinates();
  auto coordsUL   = tetra4Coordinates();

  std::vector< double > stateVarsFBar, stateVarsUL;
  FBarElement           fBarElement = buildElement< FBarElement >( coordsFBar, stateVarsFBar );
  ULElement             ulElement   = buildElement< ULElement >( coordsUL, stateVarsUL );

  constexpr int n = FBarElement::sizeLoadVector;

  std::vector< double > QTotal( n, 0.0 );
  for ( int i = 0; i < n; i++ )
    QTotal[i] = 0.05 * std::sin( 0.7 * i + 0.3 );

  std::vector< double > dQ( n, 0.0 );

  std::vector< double > PintFBar( n, 0.0 ), PintUL( n, 0.0 );
  std::vector< double > KeFBar( n * n, 0.0 ), KeUL( n * n, 0.0 );

  fBarElement.computeKernels( QTotal.data(), dQ.data(), PintFBar.data(), KeFBar.data(), 0.0, 1.0 );
  ulElement.computeKernels( QTotal.data(), dQ.data(), PintUL.data(), KeUL.data(), 0.0, 1.0 );

  for ( int i = 0; i < n; i++ )
    throwExceptionOnFailure( std::abs( PintFBar[i] - PintUL[i] ) < 1e-10,
                             "E-1: F-bar Tetra4 residual must match the standard UL Tetra4 residual (dof " +
                               std::to_string( i ) + ")" );

  for ( int i = 0; i < n * n; i++ )
    throwExceptionOnFailure( std::abs( KeFBar[i] - KeUL[i] ) < 1e-10,
                             "E-1: F-bar Tetra4 tangent must match the standard UL Tetra4 tangent (entry " +
                               std::to_string( i ) + ")" );
}

// Test E-2: the analytic tangent Ke returned by computeKernels for a Hexa8 (8 quadrature points,
// genuinely distinct from the centroid) must match a central-difference approximation of the
// residual Pint(QTotal), exercising the full F-bar correction term.
void testHexa8FiniteDifferenceTangent()
{
  using Element = DisplacementFiniteStrainFBarElement< 3, 8 >;

  auto coords = hexa8Coordinates();

  std::vector< double > stateVars;
  Element               element = buildElement< Element >( coords, stateVars );

  constexpr int n = Element::sizeLoadVector;

  std::vector< double > QTotal( n, 0.0 );
  for ( int i = 0; i < n; i++ )
    QTotal[i] = 0.03 * std::sin( 0.37 * i + 0.7 ) + 0.01 * std::cos( 0.11 * i );

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
  auto tests = std::vector< std::function< void() > >{ testTetra4ReducesToStandardElement,
                                                       testHexa8FiniteDifferenceTangent };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
