#include "Marmot/MarmotConsistentMass.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryElement.h"
#include "Marmot/MarmotMassLumping.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>
#include <cmath>
#include <tuple>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::Discretizations::Elements;

// ---------------------------------------------------------------------------------------------
// getElementShapeByMetric()
// ---------------------------------------------------------------------------------------------

void testGetElementShapeByMetricAllValidCombinations()
{
  throwExceptionOnFailure( getElementShapeByMetric( 1, 2 ) == Bar2, "nDim=1, nNodes=2 must be Bar2." );
  throwExceptionOnFailure( getElementShapeByMetric( 1, 3 ) == Bar3, "nDim=1, nNodes=3 must be Bar3." );
  throwExceptionOnFailure( getElementShapeByMetric( 2, 3 ) == Tri3, "nDim=2, nNodes=3 must be Tri3." );
  throwExceptionOnFailure( getElementShapeByMetric( 2, 4 ) == Quad4, "nDim=2, nNodes=4 must be Quad4." );
  throwExceptionOnFailure( getElementShapeByMetric( 2, 8 ) == Quad8, "nDim=2, nNodes=8 must be Quad8." );
  throwExceptionOnFailure( getElementShapeByMetric( 2, 9 ) == Quad9, "nDim=2, nNodes=9 must be Quad9." );
  throwExceptionOnFailure( getElementShapeByMetric( 2, 16 ) == Quad16, "nDim=2, nNodes=16 must be Quad16." );
  throwExceptionOnFailure( getElementShapeByMetric( 3, 4 ) == Tetra4, "nDim=3, nNodes=4 must be Tetra4." );
  throwExceptionOnFailure( getElementShapeByMetric( 3, 10 ) == Tetra10, "nDim=3, nNodes=10 must be Tetra10." );
  throwExceptionOnFailure( getElementShapeByMetric( 3, 8 ) == Hexa8, "nDim=3, nNodes=8 must be Hexa8." );
  throwExceptionOnFailure( getElementShapeByMetric( 3, 20 ) == Hexa20, "nDim=3, nNodes=20 must be Hexa20." );
  throwExceptionOnFailure( getElementShapeByMetric( 3, 27 ) == Hexa27, "nDim=3, nNodes=27 must be Hexa27." );
  throwExceptionOnFailure( getElementShapeByMetric( 3, 64 ) == Hexa64, "nDim=3, nNodes=64 must be Hexa64." );
}

void testGetElementShapeByMetricThrowsForInvalidCombinations()
{
  const std::vector< std::pair< int, int > > invalidCombinations = {
    { 1, 4 }, // no 4-node 1D shape
    { 2, 5 }, // no 5-node 2D shape
    { 3, 5 }, // no 5-node 3D shape
    { 4, 8 }, // no 4-D shapes at all
  };

  for ( const auto& [nDim, nNodes] : invalidCombinations ) {
    bool threw = false;
    try {
      getElementShapeByMetric( nDim, nNodes );
    }
    catch ( const std::invalid_argument& ) {
      threw = true;
    }
    throwExceptionOnFailure( threw,
                             "getElementShapeByMetric() must throw for nDim=" + std::to_string( nDim ) +
                               ", nNodes=" + std::to_string( nNodes ) + "." );
  }
}

// ---------------------------------------------------------------------------------------------
// NB() / Jacobian(): dynamic (runtime-sized) versions have no callers anywhere in the codebase --
// every element uses the templated (compile-time-sized) overloads declared alongside them in
// MarmotFiniteElement.h instead. Cross-validated against those already-trusted templated versions
// for a known Quad4 configuration, rather than re-deriving expected values independently.
// ---------------------------------------------------------------------------------------------

void testDynamicNBMatchesTemplatedNBForQuad4()
{
  Marmot::Discretizations::Elements::GeometryElement< 2, 4 > geo;
  const auto N = geo.N( Eigen::Vector2d( 0.2, 0.3 ) ); // NSized = Matrix<double,1,4>

  const auto templatedResult = Marmot::Discretizations::Elements::NB< 2, 4 >( N );
  const auto dynamicResult   = Marmot::Discretizations::Elements::NB( Eigen::VectorXd( N.transpose() ), 2 );

  throwExceptionOnFailure( checkIfEqual( Eigen::MatrixXd( dynamicResult ), Eigen::MatrixXd( templatedResult ), 1e-12 ),
                           "The dynamic NB() must match the templated NB<nDim,nNodes>() for the same input." );
}

void testDynamicJacobianMatchesTemplatedJacobianForQuad4()
{
  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 2.0, 0.0, 3.0, 2.0, 0.0, 1.0 }; // general quadrilateral

  Marmot::Discretizations::Elements::GeometryElement< 2, 4 > geo;
  geo.assignNodeCoordinates( nodeCoordsVec.data() );

  const Eigen::Vector2d xi( 0.2, 0.3 );
  const auto            dNdXi = geo.dNdXi( xi );

  const auto templatedJ = geo.Jacobian( dNdXi ); // uses the templated Jacobian<2,4> internally

  const Eigen::MatrixXd dynamicDNdXi = dNdXi;
  const Eigen::VectorXd dynamicCoords( Eigen::Map< const Eigen::VectorXd >( nodeCoordsVec.data(), 8 ) );
  const Eigen::MatrixXd dynamicJ = Marmot::Discretizations::Elements::Jacobian( dynamicDNdXi, dynamicCoords );

  throwExceptionOnFailure( checkIfEqual( dynamicJ, Eigen::MatrixXd( templatedJ ), 1e-12 ),
                           "The dynamic Jacobian() must match the templated Jacobian<nDim,nNodes>() for the same "
                           "input." );
}

// ---------------------------------------------------------------------------------------------
// getReferenceDomain() / getInterpolationOrder()
// ---------------------------------------------------------------------------------------------

void testReferenceDomainAndInterpolationOrderOfAllShapes()
{
  const std::vector< std::tuple< ElementShapes, ReferenceDomain, int > > expected = {
    { Bar2, ReferenceDomain::Line, 1 },
    { Bar3, ReferenceDomain::Line, 2 },
    { Tri3, ReferenceDomain::Triangle, 1 },
    { Quad4, ReferenceDomain::Quadrilateral, 1 },
    { Quad8, ReferenceDomain::Quadrilateral, 2 },
    { Quad9, ReferenceDomain::Quadrilateral, 2 },
    { Quad16, ReferenceDomain::Quadrilateral, 3 },
    { Tetra4, ReferenceDomain::Tetrahedron, 1 },
    { Tetra10, ReferenceDomain::Tetrahedron, 2 },
    { Hexa8, ReferenceDomain::Hexahedron, 1 },
    { Hexa20, ReferenceDomain::Hexahedron, 2 },
    { Hexa27, ReferenceDomain::Hexahedron, 2 },
    { Hexa64, ReferenceDomain::Hexahedron, 3 },
  };
  for ( const auto& [shape, domain, order] : expected ) {
    throwExceptionOnFailure( getReferenceDomain( shape ) == domain,
                             "getReferenceDomain() is wrong for shape " + std::to_string( shape ) + "." );
    throwExceptionOnFailure( getInterpolationOrder( shape ) == order,
                             "getInterpolationOrder() is wrong for shape " + std::to_string( shape ) + "." );
  }
}

void testNumberOfCornerNodesOfAllReferenceDomains()
{
  const std::vector< std::pair< ReferenceDomain, int > > expected = { { ReferenceDomain::Line, 2 },
                                                                      { ReferenceDomain::Triangle, 3 },
                                                                      { ReferenceDomain::Quadrilateral, 4 },
                                                                      { ReferenceDomain::Tetrahedron, 4 },
                                                                      { ReferenceDomain::Hexahedron, 8 } };
  for ( const auto& [domain, nCorners] : expected )
    throwExceptionOnFailure( getNumberOfCornerNodes( domain ) == nCorners,
                             "getNumberOfCornerNodes() is wrong for reference domain " +
                               std::to_string( static_cast< int >( domain ) ) + "." );
}

// ---------------------------------------------------------------------------------------------
// MassLumping::nLinearNodes(): the linear element blended into the manifold-based lumping has the corner nodes of
// the element; checked at compile time, since the callers use it as a template argument
// ---------------------------------------------------------------------------------------------

namespace {
  using Marmot::Discretizations::Elements::MassLumping::nLinearNodes;

  static_assert( nLinearNodes< 1, 2 >() == 2, "Bar2 is its own linear element." );
  static_assert( nLinearNodes< 1, 3 >() == 2, "The linear element of Bar3 is Bar2." );
  static_assert( nLinearNodes< 2, 3 >() == 3, "Tri3 is its own linear element." );
  static_assert( nLinearNodes< 2, 4 >() == 4, "Quad4 is its own linear element." );
  static_assert( nLinearNodes< 2, 8 >() == 4, "The linear element of Quad8 is Quad4." );
  static_assert( nLinearNodes< 2, 9 >() == 4, "The linear element of Quad9 is Quad4." );
  static_assert( nLinearNodes< 2, 16 >() == 4, "The linear element of Quad16 is Quad4." );
  static_assert( nLinearNodes< 3, 4 >() == 4, "Tetra4 is its own linear element." );
  static_assert( nLinearNodes< 3, 10 >() == 4, "The linear element of Tetra10 is Tetra4, not Hexa8." );
  static_assert( nLinearNodes< 3, 8 >() == 8, "Hexa8 is its own linear element." );
  static_assert( nLinearNodes< 3, 20 >() == 8, "The linear element of Hexa20 is Hexa8." );
  static_assert( nLinearNodes< 3, 27 >() == 8, "The linear element of Hexa27 is Hexa8." );
  static_assert( nLinearNodes< 3, 64 >() == 8, "The linear element of Hexa64 is Hexa8." );
} // namespace

// ---------------------------------------------------------------------------------------------
// Quadrature: the catalogue quadratureRules() and the lookup getQuadraturePoints()
// ---------------------------------------------------------------------------------------------

namespace {
  using namespace Marmot::Discretizations::Elements::Quadrature;

  double factorial( int n )
  {
    return n <= 1 ? 1.0 : n * factorial( n - 1 );
  }

  /// exact integral of the monomial xi^e[0] * eta^e[1] * zeta^e[2] over a reference domain
  double exactIntegralOfMonomial( ReferenceDomain domain, const std::vector< int >& e )
  {
    switch ( domain ) {
    case ReferenceDomain::Line:
    case ReferenceDomain::Quadrilateral:
    case ReferenceDomain::Hexahedron: {
      double result = 1;
      for ( const int a : e )
        result *= a % 2 == 0 ? 2.0 / ( a + 1 ) : 0.0;
      return result;
    }
    case ReferenceDomain::Triangle:
    case ReferenceDomain::Tetrahedron: {
      // a! b! (c!) / ( a + b (+ c) + nDim )!
      double numerator = 1;
      int    sum       = 0;
      for ( const int a : e ) {
        numerator *= factorial( a );
        sum += a;
      }
      return numerator / factorial( sum + static_cast< int >( e.size() ) );
    }
    }
    throw std::invalid_argument( "invalid reference domain" );
  }

  /// all exponents of the monomials of the polynomial space a rule claims: Q_degree on tensor-product domains,
  /// P_degree on simplices
  std::vector< std::vector< int > > exponentsOfPolynomialSpace( const QuadratureRuleData& rule )
  {
    const int  nDim    = static_cast< int >( rule.points.front().xi.size() );
    const int  degree  = rule.degree;
    const bool simplex = isSimplex( rule.domain );

    std::vector< std::vector< int > > exponents;
    std::vector< int >                e( nDim, 0 );
    while ( true ) {
      int sum = 0;
      for ( const int a : e )
        sum += a;
      if ( !simplex || sum <= degree )
        exponents.push_back( e );

      // next combination with every exponent in [0, degree]
      int i = 0;
      while ( i < nDim && e[i] == degree )
        e[i++] = 0;
      if ( i == nDim )
        break;
      e[i]++;
    }
    return exponents;
  }

  int dimension( ReferenceDomain domain )
  {
    switch ( domain ) {
    case ReferenceDomain::Line: return 1;
    case ReferenceDomain::Triangle:
    case ReferenceDomain::Quadrilateral: return 2;
    case ReferenceDomain::Tetrahedron:
    case ReferenceDomain::Hexahedron: return 3;
    }
    throw std::invalid_argument( "invalid reference domain" );
  }

  /// the linear element of a reference domain, whose nodes are the corners of the domain
  ElementShapes linearShape( ReferenceDomain domain )
  {
    return getElementShapeByMetric( dimension( domain ), getNumberOfCornerNodes( domain ) );
  }

  /// a fixed-size row of shape function values as a dynamic column vector
  /// @details Copied value by value: assigning a fixed-size row of fewer than four values to a dynamic vector with
  /// Eigen makes GCC 14 at x86-64-v3 warn about the AVX packet loop (-Warray-bounds), although it never runs.
  template < typename Derived >
  Eigen::VectorXd asColumn( const Eigen::MatrixBase< Derived >& N )
  {
    Eigen::VectorXd column( N.size() );
    for ( Eigen::Index i = 0; i < N.size(); i++ )
      column( i ) = N( i );
    return column;
  }

  /// the shape functions of every shape that implements them, as a column vector
  Eigen::VectorXd shapeFunctions( ElementShapes shape, const Eigen::VectorXd& xi )
  {
    namespace Elements = Marmot::Discretizations::Elements;
    switch ( shape ) {
    case Bar2: return asColumn( Elements::Spatial1D::Bar2::N( xi( 0 ) ) );
    case Bar3: return asColumn( Elements::Spatial1D::Bar3::N( xi( 0 ) ) );
    case Tri3: return asColumn( Elements::Spatial2D::Tri3::N( xi ) );
    case Quad4: return asColumn( Elements::Spatial2D::Quad4::N( xi ) );
    case Quad8: return asColumn( Elements::Spatial2D::Quad8::N( xi ) );
    case Tetra4: return asColumn( Elements::Spatial3D::Tetra4::N( xi ) );
    case Tetra10: return asColumn( Elements::Spatial3D::Tetra10::N( xi ) );
    case Hexa8: return asColumn( Elements::Spatial3D::Hexa8::N( xi ) );
    case Hexa20: return asColumn( Elements::Spatial3D::Hexa20::N( xi ) );
    default: throw std::invalid_argument( "no shape functions for this shape" );
    }
  }

  /// the corners of a reference domain in corner order, written out independently of the catalogue
  std::vector< Eigen::VectorXd > corners( ReferenceDomain domain )
  {
    auto v = []( std::initializer_list< double > c ) {
      Eigen::VectorXd x( static_cast< Eigen::Index >( c.size() ) );
      Eigen::Index    i = 0;
      for ( const double a : c )
        x( i++ ) = a;
      return x;
    };
    switch ( domain ) {
    case ReferenceDomain::Line: return { v( { -1 } ), v( { 1 } ) };
    case ReferenceDomain::Triangle: return { v( { 0, 0 } ), v( { 1, 0 } ), v( { 0, 1 } ) };
    case ReferenceDomain::Quadrilateral: return { v( { -1, -1 } ), v( { 1, -1 } ), v( { 1, 1 } ), v( { -1, 1 } ) };
    case ReferenceDomain::Tetrahedron:
      return { v( { 0, 0, 0 } ), v( { 1, 0, 0 } ), v( { 0, 1, 0 } ), v( { 0, 0, 1 } ) };
    case ReferenceDomain::Hexahedron:
      return { v( { -1, -1, -1 } ),
               v( { 1, -1, -1 } ),
               v( { 1, 1, -1 } ),
               v( { -1, 1, -1 } ),
               v( { -1, -1, 1 } ),
               v( { 1, -1, 1 } ),
               v( { 1, 1, 1 } ),
               v( { -1, 1, 1 } ) };
    }
    throw std::invalid_argument( "invalid reference domain" );
  }

  template < typename Call >
  bool throwsInvalidArgument( Call call )
  {
    try {
      call();
    }
    catch ( const std::invalid_argument& ) {
      return true;
    }
    return false;
  }
} // namespace

void testEveryRuleOfTheCatalogueIntegratesItsDegreeExactly()
{
  for ( const auto& rule : quadratureRules() ) {
    for ( const auto& e : exponentsOfPolynomialSpace( rule ) ) {
      double integral = 0;
      for ( const auto& qp : rule.points ) {
        double value = qp.weight;
        for ( size_t i = 0; i < e.size(); i++ )
          value *= std::pow( qp.xi( i ), e[i] );
        integral += value;
      }
      // the 3-point Gauss abscissa and the 3x3x3 weights are tabulated to about 15 digits
      throwExceptionOnFailure( checkIfEqual( integral, exactIntegralOfMonomial( rule.domain, e ), 1e-13 ),
                               "A rule of degree " + std::to_string( rule.degree ) + " with " +
                                 std::to_string( rule.points.size() ) +
                                 " points does not integrate a monomial of its polynomial space exactly." );
    }
  }
}

void testGaussRulesHavePositiveWeightsAndInteriorPoints()
{
  // a negative weight can make the consistent mass indefinite, and a point outside the domain breaks the evaluation
  // of the material there; the Newton-Cotes rules have their points on the boundary by construction
  for ( const auto& rule : quadratureRules() ) {
    if ( rule.family != QuadratureFamily::Gauss )
      continue;
    for ( const auto& qp : rule.points ) {
      const bool inside = isSimplex( rule.domain ) ? qp.xi.minCoeff() > 0 && qp.xi.sum() < 1
                                                   : qp.xi.cwiseAbs().maxCoeff() < 1;
      throwExceptionOnFailure( qp.weight > 0 && inside,
                               "The Gauss rule of degree " + std::to_string( rule.degree ) + " with " +
                                 std::to_string( rule.points.size() ) +
                                 " points has a non-positive weight or a point outside its reference domain." );
    }
  }
}

void testCatalogueIsSortedByNumberOfPointsAndDegree()
{
  const auto& rules = quadratureRules();
  for ( size_t i = 0; i < rules.size(); i++ )
    for ( size_t j = i + 1; j < rules.size(); j++ )
      if ( rules[i].domain == rules[j].domain && rules[i].family == rules[j].family )
        throwExceptionOnFailure( rules[i].points.size() < rules[j].points.size() && rules[i].degree < rules[j].degree,
                                 "The rules of one domain and family must be sorted by increasing number of points "
                                 "and degree, so that the lookup finds the cheapest rule first." );
}

void testGetQuadraturePointsReturnsTheCheapestRule()
{
  // domain, degree, expected number of points; together with the sorted catalogue this covers the lookup
  const std::vector< std::tuple< ReferenceDomain, int, size_t > > expected = {
    { ReferenceDomain::Line, 0, 1 },          // degree 0 takes the first rule
    { ReferenceDomain::Line, 2, 2 },          // a degree between two rules is rounded up
    { ReferenceDomain::Quadrilateral, 2, 4 }, // degree in each coordinate (Q_2): 2x2 points ...
    { ReferenceDomain::Triangle, 2, 3 },      // ... total degree (P_2) on simplices: 3 points
    { ReferenceDomain::Tetrahedron, 2, 4 },   // the consistent mass of C3D4
    { ReferenceDomain::Tetrahedron, 3, 14 },  // a gap in the catalogue: the next rule has degree 5
  };
  for ( const auto& [domain, degree, nPoints] : expected )
    throwExceptionOnFailure( getQuadraturePoints( domain, gauss( degree ) ).size() == nPoints,
                             "getQuadraturePoints() must return the Gauss rule with the fewest points for degree " +
                               std::to_string( degree ) + "." );

  // the shape overload looks up the reference domain of the shape: all shapes on it get the same rule, Newton-Cotes
  // rules included
  throwExceptionOnFailure( &getQuadraturePoints( Hexa8, gauss( 3 ) ) ==
                               &getQuadraturePoints( ReferenceDomain::Hexahedron, gauss( 3 ) ) &&
                             &getQuadraturePoints( Hexa20, gauss( 3 ) ) ==
                               &getQuadraturePoints( ReferenceDomain::Hexahedron, gauss( 3 ) ),
                           "Hexa8 and Hexa20 must get the Gauss rule of the hexahedron." );
  throwExceptionOnFailure( &getQuadraturePoints( Quad4, newtonCotes( 1 ) ) ==
                               &getQuadraturePoints( ReferenceDomain::Quadrilateral, newtonCotes( 1 ) ) &&
                             &getQuadraturePoints( Quad8, newtonCotes( 1 ) ) ==
                               &getQuadraturePoints( ReferenceDomain::Quadrilateral, newtonCotes( 1 ) ),
                           "Quad4 and Quad8 must get the Newton-Cotes rule of the quadrilateral." );
}

void testGetQuadraturePointsThrowsWithoutAMatchingRule()
{
  throwExceptionOnFailure( throwsInvalidArgument( [] { getQuadraturePoints( Quad4, gauss( 6 ) ); } ),
                           "There is no Gauss rule of degree 6 on a quadrilateral." );
  throwExceptionOnFailure( throwsInvalidArgument( [] { getQuadraturePoints( Tetra10, gauss( 6 ) ); } ),
                           "There is no Gauss rule of degree 6 on a tetrahedron." );
  throwExceptionOnFailure( throwsInvalidArgument( [] { getQuadraturePoints( Hexa8, newtonCotes( 1 ) ); } ),
                           "There is no Newton-Cotes rule on a hexahedron." );
  throwExceptionOnFailure( throwsInvalidArgument( [] { getQuadraturePoints( Quad4, newtonCotes( 2 ) ); } ),
                           "There is no Newton-Cotes rule of degree 2 on a quadrilateral." );
}

void testEveryNewtonCotesRuleIsTheNodalQuadratureOfItsLinearElement()
{
  // the points of a Newton-Cotes rule of degree 1 are the nodes of the linear element of its domain, in node order
  // (N_i(xi_j) = delta_ij), its weights are the integrals of their shape functions, which must be positive, and
  // newtonCotes( 1 ) finds it for that element
  for ( const auto& rule : quadratureRules() ) {
    if ( rule.family != QuadratureFamily::NewtonCotes )
      continue;
    const ElementShapes shape = linearShape( rule.domain );
    const std::string   name  = "The Newton-Cotes rule of degree " + std::to_string( rule.degree ) + " of shape " +
                             std::to_string( shape );
    const auto& points = rule.points;

    throwExceptionOnFailure( &getQuadraturePoints( shape, newtonCotes( rule.degree ) ) == &points,
                             name + " must be the one newtonCotes() returns for the shape." );

    for ( size_t j = 0; j < points.size(); j++ ) {
      const Eigen::VectorXd Nj = shapeFunctions( shape, points[j].xi );
      throwExceptionOnFailure( Nj.size() == static_cast< Eigen::Index >( points.size() ),
                               name + " must have one point per node." );
      for ( size_t i = 0; i < points.size(); i++ )
        throwExceptionOnFailure( checkIfEqual( Nj( i ), i == j ? 1.0 : 0.0, 1e-15 ),
                                 name + " must have point " + std::to_string( j ) + " at node " + std::to_string( j ) +
                                   "." );
    }

    Eigen::VectorXd integralOfN = Eigen::VectorXd::Zero( static_cast< Eigen::Index >( points.size() ) );
    for ( const auto& qp : getQuadraturePoints( shape, gauss( 2 * getInterpolationOrder( shape ) ) ) )
      integralOfN += qp.weight * shapeFunctions( shape, qp.xi );

    for ( size_t i = 0; i < points.size(); i++ ) {
      throwExceptionOnFailure( checkIfEqual( points[i].weight, integralOfN( i ), 1e-15 ),
                               name + " must have the integral of N_" + std::to_string( i ) + " as weight " +
                                 std::to_string( i ) + "." );
      throwExceptionOnFailure( points[i].weight > 0, name + " must have positive weights." );
    }
  }
}

void testEveryShapeListsTheCornersOfItsDomainFirstInCornerOrder()
{
  // Node i < nCorners of every shape is corner i of its reference domain (the node order of the linear element):
  // N_i(corner_j) = delta_ij for all nodes i. The manifold-based mass lumping blends the linear element into the
  // first nodes, and the Newton-Cotes rules of degree 1 have the corners as points, so both rely on this.
  //
  // Tetra10 is left out: its shape functions number the corners (0, 1, 0), (0, 0, 1), (0, 0, 0), (1, 0, 0), a
  // permutation of the corner order. No registered element uses Tetra10; it must be renumbered before one does.
  for ( const ElementShapes shape : { Bar2, Bar3, Tri3, Quad4, Quad8, Tetra4, Hexa8, Hexa20 } ) {
    const auto corners_ = corners( getReferenceDomain( shape ) );
    for ( size_t j = 0; j < corners_.size(); j++ ) {
      const Eigen::VectorXd N = shapeFunctions( shape, corners_[j] );
      for ( Eigen::Index i = 0; i < N.size(); i++ )
        throwExceptionOnFailure( checkIfEqual( N( i ), i == static_cast< Eigen::Index >( j ) ? 1.0 : 0.0, 1e-15 ),
                                 "Node " + std::to_string( j ) + " of shape " + std::to_string( shape ) +
                                   " must be corner " + std::to_string( j ) + " of its reference domain." );
    }
  }

  // the Newton-Cotes rules of degree 1 have the corners as points, in corner order
  for ( const auto& rule : quadratureRules() ) {
    if ( rule.family != QuadratureFamily::NewtonCotes || rule.degree != 1 )
      continue;
    const auto corners_ = corners( rule.domain );
    throwExceptionOnFailure( rule.points.size() == corners_.size(),
                             "A Newton-Cotes rule of degree 1 must have one point per corner." );
    for ( size_t j = 0; j < corners_.size(); j++ )
      throwExceptionOnFailure( rule.points[j].xi == corners_[j],
                               "Point " + std::to_string( j ) + " of a Newton-Cotes rule of degree 1 must be corner " +
                                 std::to_string( j ) + "." );
  }
}

void testConsistentMassRuleIsOfDegree2p()
{
  // the consistent mass N^T N has degree 2p, whatever rule the element uses otherwise
  using Marmot::Discretizations::Elements::ConsistentMass::integrationRule;
  const std::vector< std::pair< ElementShapes, size_t > > expected = {
    { Bar3, 3 },     // degree 4 on the line: 3 points
    { Quad4, 4 },    // degree 2 in each coordinate: 2x2
    { Quad8, 9 },    // degree 4 in each coordinate: 3x3
    { Tetra4, 4 },   // total degree 2: 4 points, where the element rule of C3D4 has 1
    { Tetra10, 14 }, // total degree 4: the 14-point rule of degree 5
    { Hexa20, 27 },  // 3x3x3
  };
  for ( const auto& [shape, nPoints] : expected )
    throwExceptionOnFailure( integrationRule( shape ).size() == nPoints,
                             "The consistent mass of shape " + std::to_string( shape ) + " must be integrated with " +
                               std::to_string( nPoints ) + " points." );

  throwExceptionOnFailure( throwsInvalidArgument( [] { integrationRule( Hexa64 ); } ),
                           "There is no Gauss rule of degree 6 for the consistent mass of the cubic Hexa64." );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{
    testGetElementShapeByMetricAllValidCombinations,
    testGetElementShapeByMetricThrowsForInvalidCombinations,
    testDynamicNBMatchesTemplatedNBForQuad4,
    testDynamicJacobianMatchesTemplatedJacobianForQuad4,
    testReferenceDomainAndInterpolationOrderOfAllShapes,
    testNumberOfCornerNodesOfAllReferenceDomains,
    testEveryRuleOfTheCatalogueIntegratesItsDegreeExactly,
    testGaussRulesHavePositiveWeightsAndInteriorPoints,
    testCatalogueIsSortedByNumberOfPointsAndDegree,
    testGetQuadraturePointsReturnsTheCheapestRule,
    testGetQuadraturePointsThrowsWithoutAMatchingRule,
    testEveryNewtonCotesRuleIsTheNodalQuadratureOfItsLinearElement,
    testEveryShapeListsTheCornersOfItsDomainFirstInCornerOrder,
    testConsistentMassRuleIsOfDegree2p,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
