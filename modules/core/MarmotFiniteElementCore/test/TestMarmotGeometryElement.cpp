#include "Marmot/MarmotGeometryElement.h"
#include "Marmot/MarmotNumericalDifferentiation.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>

using namespace Marmot;
using namespace Marmot::Testing;
namespace NumDiff = Marmot::NumericalAlgorithms::Differentiation;

// ---------------------------------------------------------------------------------------------
// Tetra4 / Tetra10: neither is used by any element tested elsewhere in this codebase, so their
// shape functions have never been exercised. Rather than trust a hand-transcribed formula (this
// file's Tetra10 dNdXi in particular is a long hand-expanded product-rule derivative, exactly the
// kind of code most prone to a transcription slip), correctness is checked two ways that need no
// prior knowledge of nodal coordinates: N must partition unity (and reduce to the Kronecker delta
// at each of Tetra4's easily-derived nodal points), and dNdXi must equal the numerical derivative
// of N with respect to the same natural coordinates.
// ---------------------------------------------------------------------------------------------

void testTetra4PartitionOfUnityAndKroneckerDeltaAtNodes()
{
  using namespace Marmot::FiniteElement::Spatial3D;

  // Derived directly from N: N0=1-xi0-xi1-xi2, N1=xi0, N2=xi1, N3=xi2.
  const std::vector< Eigen::Vector3d > nodalXi = { { 0.0, 0.0, 0.0 },
                                                   { 1.0, 0.0, 0.0 },
                                                   { 0.0, 1.0, 0.0 },
                                                   { 0.0, 0.0, 1.0 } };

  for ( int i = 0; i < 4; i++ ) {
    const auto N = Tetra4::N( nodalXi[i] );
    throwExceptionOnFailure( checkIfEqual( N.sum(), 1.0, 1e-12 ), "Tetra4::N() does not partition unity." );
    for ( int j = 0; j < 4; j++ ) {
      const double expected = ( i == j ) ? 1.0 : 0.0;
      throwExceptionOnFailure( checkIfEqual( N( j ), expected, 1e-12 ),
                               "Tetra4::N() does not reduce to the Kronecker delta at its own nodes." );
    }
  }
}

void testTetra4DNdXiMatchesNumericalDifferentiationOfN()
{
  using namespace Marmot::FiniteElement::Spatial3D;

  const Eigen::Vector3d xi( 0.2, 0.3, 0.1 );

  NumDiff::vector_to_vector_function_type Nfunc = []( const Eigen::VectorXd& xiVec ) -> Eigen::VectorXd {
    return Tetra4::N( Eigen::Vector3d( xiVec ) );
  };

  const Eigen::MatrixXd numJacobian = NumDiff::centralDifference( Nfunc, xi ); // (dN_i/dXi_j), 4x3

  const auto analyticDNdXi = Tetra4::dNdXi( xi ); // (dN_i/dXi_j) stored as (nDim x nNodes) = dNdXi(j,i)

  throwExceptionOnFailure( checkIfEqual( numJacobian, Eigen::MatrixXd( analyticDNdXi.transpose() ), 1e-6 ),
                           "Tetra4::dNdXi() does not match the numerical derivative of Tetra4::N()." );
}

void testTetra10PartitionOfUnityAtSamplePoints()
{
  using namespace Marmot::FiniteElement::Spatial3D;

  const std::vector< Eigen::Vector3d > samples = { { 0.0, 0.0, 0.0 },
                                                   { 1.0, 0.0, 0.0 },
                                                   { 0.0, 1.0, 0.0 },
                                                   { 0.0, 0.0, 1.0 },
                                                   { 0.2, 0.3, 0.1 },
                                                   { 0.25, 0.25, 0.25 } };

  for ( const auto& xi : samples ) {
    const auto N = Tetra10::N( xi );
    throwExceptionOnFailure( checkIfEqual( N.sum(), 1.0, 1e-10 ), "Tetra10::N() does not partition unity." );
  }
}

void testTetra10DNdXiMatchesNumericalDifferentiationOfN()
{
  using namespace Marmot::FiniteElement::Spatial3D;

  const Eigen::Vector3d xi( 0.2, 0.3, 0.1 );

  NumDiff::vector_to_vector_function_type Nfunc = []( const Eigen::VectorXd& xiVec ) -> Eigen::VectorXd {
    return Tetra10::N( Eigen::Vector3d( xiVec ) );
  };

  const Eigen::MatrixXd numJacobian   = NumDiff::centralDifference( Nfunc, xi );
  const auto            analyticDNdXi = Tetra10::dNdXi( xi );

  throwExceptionOnFailure( checkIfEqual( numJacobian, Eigen::MatrixXd( analyticDNdXi.transpose() ), 1e-6 ),
                           "Tetra10::dNdXi() does not match the numerical derivative of Tetra10::N() -- likely "
                           "a transcription error in the hand-expanded derivative." );
}

// ---------------------------------------------------------------------------------------------
// BGreen(dNdX, F): at F = Identity, the Green-Lagrange strain operator reduces algebraically to
// the standard linear B-operator (every cross term multiplies a zero off-diagonal entry of F, and
// every retained term multiplies F's unit diagonal) -- an exact identity, not an approximation,
// and one that exercises MarmotGeometryElement's BGreen() specialization (which has no callers
// anywhere in the codebase) for every shape.
// ---------------------------------------------------------------------------------------------

template < int nDim, int nNodes >
void checkBGreenAtIdentityMatchesB( const std::vector< double >&            nodeCoordsVec,
                                    const Eigen::Matrix< double, nDim, 1 >& xi,
                                    const std::string&                      label )
{
  MarmotGeometryElement< nDim, nNodes > geo;
  geo.assignNodeCoordinates( nodeCoordsVec.data() );

  const auto dNdXi = geo.dNdXi( xi );
  const auto J     = geo.Jacobian( dNdXi );
  const auto dNdX  = geo.dNdX( dNdXi, J.inverse() );

  const auto B      = geo.B( dNdX );
  const auto BGreen = geo.BGreen( dNdX, Eigen::Matrix< double, nDim, nDim >::Identity() );

  throwExceptionOnFailure( checkIfEqual( Eigen::MatrixXd( BGreen ), Eigen::MatrixXd( B ), 1e-10 ),
                           label + ": BGreen(dNdX, Identity) must equal B(dNdX)." );
}

void testBGreenAtIdentityMatchesBForEveryShape()
{
  const std::vector< double > quad4Coords = { 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0 };
  checkBGreenAtIdentityMatchesB< 2, 4 >( quad4Coords, Eigen::Vector2d( 0.2, 0.3 ), "Quad4" );

  const std::vector< double > quad8Coords =
    { 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.5, 0.0, 1.0, 0.5, 0.5, 1.0, 0.0, 0.5 };
  checkBGreenAtIdentityMatchesB< 2, 8 >( quad8Coords, Eigen::Vector2d( 0.2, 0.3 ), "Quad8" );

  const std::vector< double > hexa8Coords = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0,
                                              0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0 };
  checkBGreenAtIdentityMatchesB< 3, 8 >( hexa8Coords, Eigen::Vector3d( 0.2, 0.3, 0.1 ), "Hexa8" );

  const std::vector< double > hexa20Coords = {
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, // bottom corners
    0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, // top corners
    0.5, 0.0, 0.0, 1.0, 0.5, 0.0, 0.5, 1.0, 0.0, 0.0, 0.5, 0.0, // bottom-face edge midsides
    0.5, 0.0, 1.0, 1.0, 0.5, 1.0, 0.5, 1.0, 1.0, 0.0, 0.5, 1.0, // top-face edge midsides
    0.0, 0.0, 0.5, 1.0, 0.0, 0.5, 1.0, 1.0, 0.5, 0.0, 1.0, 0.5  // vertical edge midsides
  };
  checkBGreenAtIdentityMatchesB< 3, 20 >( hexa20Coords, Eigen::Vector3d( 0.2, 0.3, 0.1 ), "Hexa20" );

  const std::vector< double > tetra4Coords = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };
  checkBGreenAtIdentityMatchesB< 3, 4 >( tetra4Coords, Eigen::Vector3d( 0.2, 0.2, 0.2 ), "Tetra4" );
}

// ---------------------------------------------------------------------------------------------
// B_bar(dNdX, dNdX0): when evaluated at the SAME point (dNdX0 == dNdX), the B-bar correction terms
// (B4, B6, B8) vanish identically and B_bar algebraically collapses to the standard B operator --
// again an exact identity, not an approximation.
// ---------------------------------------------------------------------------------------------

void testBBarAtEqualPointsMatchesBForHexa8()
{
  const std::vector< double > hexa8Coords = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0,
                                              0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0 };

  MarmotGeometryElement< 3, 8 > geo;
  geo.assignNodeCoordinates( hexa8Coords.data() );

  const Eigen::Vector3d xi( 0.2, 0.3, 0.1 );
  const auto            dNdXi = geo.dNdXi( xi );
  const auto            J     = geo.Jacobian( dNdXi );
  const auto            dNdX  = geo.dNdX( dNdXi, J.inverse() );

  const auto B     = geo.B( dNdX );
  const auto B_bar = geo.B_bar( dNdX, dNdX );

  throwExceptionOnFailure( checkIfEqual( Eigen::MatrixXd( B_bar ), Eigen::MatrixXd( B ), 1e-10 ),
                           "B_bar(dNdX, dNdX) must equal B(dNdX) when evaluated at the same point." );
}

// ---------------------------------------------------------------------------------------------
// B_axisymmetric: no algebraic identity collapses this to the standard B operator (it introduces
// a genuinely new N/r hoop-strain term), so it is checked against hand-computed values instead.
// ---------------------------------------------------------------------------------------------

void testBAxisymmetricMatchesHandDerivedValuesQuad4()
{
  MarmotGeometryElement< 2, 4 >::dNdXiSized dNdX;
  dNdX << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0; // arbitrary but fixed 2x4 values

  MarmotGeometryElement< 2, 4 >::NSized N;
  N << 0.1, 0.2, 0.3, 0.4;

  const Eigen::Vector2d x_gauss( 2.0, 0.0 ); // r = 2

  MarmotGeometryElement< 2, 4 > geo;
  const auto                    B = geo.B_axisymmetric( dNdX, N, x_gauss );

  throwExceptionOnFailure( B.rows() == 4 && B.cols() == 8, "B_axisymmetric() returned the wrong size." );

  for ( int i = 0; i < 4; i++ ) {
    throwExceptionOnFailure( checkIfEqual( B( 0, 2 * i ), dNdX( 0, i ), 1e-12 ), "B_axisymmetric row 0 mismatch." );
    throwExceptionOnFailure( checkIfEqual( B( 0, 2 * i + 1 ), 0.0, 1e-12 ), "B_axisymmetric row 0 mismatch." );
    throwExceptionOnFailure( checkIfEqual( B( 1, 2 * i + 1 ), dNdX( 1, i ), 1e-12 ), "B_axisymmetric row 1 mismatch." );
    throwExceptionOnFailure( checkIfEqual( B( 1, 2 * i ), 0.0, 1e-12 ), "B_axisymmetric row 1 mismatch." );
    throwExceptionOnFailure( checkIfEqual( B( 2, 2 * i ), N( i ) / x_gauss( 0 ), 1e-12 ),
                             "B_axisymmetric row 2 (hoop strain N/r) mismatch." );
    throwExceptionOnFailure( checkIfEqual( B( 3, 2 * i ), dNdX( 1, i ), 1e-12 ), "B_axisymmetric row 3 mismatch." );
    throwExceptionOnFailure( checkIfEqual( B( 3, 2 * i + 1 ), dNdX( 0, i ), 1e-12 ), "B_axisymmetric row 3 mismatch." );
  }
}

int main()
{
  auto tests = std::vector< std::function< void() > >{
    testTetra4PartitionOfUnityAndKroneckerDeltaAtNodes,
    testTetra4DNdXiMatchesNumericalDifferentiationOfN,
    testTetra10PartitionOfUnityAtSamplePoints,
    testTetra10DNdXiMatchesNumericalDifferentiationOfN,
    testBGreenAtIdentityMatchesBForEveryShape,
    testBBarAtEqualPointsMatchesBForHexa8,
    testBAxisymmetricMatchesHandDerivedValuesQuad4,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
