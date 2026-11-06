#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTypedefs.h"
#include <cmath>
#include <stdexcept>

using namespace Eigen;
namespace Marmot::Discretizations::Elements {
  // create a large N matrix dependent on the degrees of freedom per node
  MatrixXd NB( const VectorXd& N, const int nDoFPerNode )
  {

    MatrixXd N_( nDoFPerNode, N.size() * nDoFPerNode );
    N_ = MatrixXd::Zero( nDoFPerNode, N.size() * nDoFPerNode );

    for ( int i = 0; i < N.size(); i++ ) {
      for ( int j = 0; j < nDoFPerNode; j++ ) {
        N_( j, nDoFPerNode * i + j ) = N( i );
      }
    }

    return N_;
  }

  MatrixXd Jacobian( const MatrixXd& dNdXi, const VectorXd& coordinates )
  {

    /* Dynamic version of jacobian, and not necessarily square!
     *
     * Notation:
     *
     * /                          \
     * | x1,xi1,  x1,xi2,  x1,xi3 |
     * |                          |
     * | x2,xi1,  x2,xi2,  x2,xi3 |
     * |                          |
     * | x3,xi1,  x3,xi2,  x3,xi3 |
     * \                          /
     *
     * */

    const int nDimXi = static_cast< int >( dNdXi.rows() );
    const int nNodes = static_cast< int >( dNdXi.cols() );
    const int nDimX  = static_cast< int >( coordinates.size() ) / nNodes;

    MatrixXd J_ = MatrixXd::Zero( nDimX, nDimXi );

    for ( int i = 0; i < nDimX; i++ )      // loop over global dimensions
      for ( int j = 0; j < nDimXi; j++ )   // loop over local dimensions
        for ( int k = 0; k < nNodes; k++ ) // Loop over nodes
          J_( i, j ) += dNdXi( j, k ) * coordinates( i + k * nDimX );
    return J_;
  }

  VectorXi expandNodeIndicesToCoordinateIndices( const VectorXi& nodeIndices, int nDim )
  {
    VectorXi indices( nDim * nodeIndices.size() );

    for ( int i = 0; i < nodeIndices.size(); i++ )
      for ( int k = 0; k < nDim; k++ )
        indices( ( i * nDim ) + k ) = ( nodeIndices( i ) * nDim ) + k;

    return indices;
  }

  namespace Quadrature {

    const std::vector< QuadratureRuleData >& quadratureRules()
    {
      using Domain = ReferenceDomain;
      using Family = QuadratureFamily;

      // coordinates of the 4-point tetrahedron rule
      const double a = ( 5 - std::sqrt( 5 ) ) / 20;
      const double b = ( 5 + 3 * std::sqrt( 5 ) ) / 20;

      // the 7-point triangle rule of degree 5 (Radon 1948): the centroid and two orbits (t, t, 1 - 2t)
      const double s15 = std::sqrt( 15. );
      const double t1  = ( 6 - s15 ) / 21;
      const double t2  = ( 6 + s15 ) / 21;
      const double v1  = ( 155 - s15 ) / 2400;
      const double v2  = ( 155 + s15 ) / 2400;

      // the 14-point tetrahedron rule of degree 5 (Walkington 2000): two orbits (a, a, a, 1 - 3a) and one orbit
      // (c, c, 1/2 - c, 1/2 - c), all weights positive; solved to 50 digits from the moment equations
      const double a1 = 0.0927352503108912264023239, w1 = 0.0122488405193936582572850;
      const double a2 = 0.310885919263300609797346, w2 = 0.0187813209530026417998643;
      const double c3 = 0.0455037041256496494918805, w3 = 0.00709100346284691107301157;
      const double r1 = 1 - 3 * a1, r2 = 1 - 3 * a2, r3 = 0.5 - c3;

      // Within each domain and family, the rules MUST be sorted by increasing number of points and degree:
      // getQuadraturePoints() returns the first rule that reaches the requested degree, so a rule placed too early
      // would be returned for requests a cheaper rule would serve, e.g. 3x3x3 points for gauss( 1 ). A new rule goes
      // to its place in this order (checked by testCatalogueIsSortedByNumberOfPointsAndDegree). The points of a
      // Newton-Cotes rule of degree 1 are the corners of the domain in corner order (see ReferenceDomain).
      // clang-format off
      static const std::vector< QuadratureRuleData > rules = {
        // ------------------------------------------------------------------------------------------------ Line
        { Domain::Line, Family::Gauss, 1, {
            { ( VectorXd( 1 ) << 0 ).finished(),                    2.0 } } },
        { Domain::Line, Family::Gauss, 3, {
            { ( VectorXd( 1 ) << -gp2 ).finished(),                 1.0 },
            { ( VectorXd( 1 ) << +gp2 ).finished(),                 1.0 } } },
        { Domain::Line, Family::Gauss, 5, {
            { ( VectorXd( 1 ) << -gp3 ).finished(),                 5./9 },
            { ( VectorXd( 1 ) << 0. ).finished(),                   8./9 },
            { ( VectorXd( 1 ) << +gp3 ).finished(),                 5./9 } } },
        // trapezoidal rule = 2-point Gauss-Lobatto rule
        { Domain::Line, Family::NewtonCotes, 1, {
            { ( VectorXd( 1 ) << -1. ).finished(),                  1.0 },
            { ( VectorXd( 1 ) << +1. ).finished(),                  1.0 } } },

        // ------------------------------------------------------------------------------------------------ Triangle
        { Domain::Triangle, Family::Gauss, 1, {
            { ( Vector2d() << 1./3, 1./3 ).finished(),              1./2 } } },
        { Domain::Triangle, Family::Gauss, 2, {
            { ( Vector2d() << 1./6, 1./6 ).finished(),              1./6 },
            { ( Vector2d() << 2./3, 1./6 ).finished(),              1./6 },
            { ( Vector2d() << 1./6, 2./3 ).finished(),              1./6 } } },
        { Domain::Triangle, Family::Gauss, 5, {
            { ( Vector2d() << 1./3, 1./3 ).finished(),              9./80 },
            { ( Vector2d() << t1, t1 ).finished(),                  v1 },
            { ( Vector2d() << 1 - 2 * t1, t1 ).finished(),          v1 },
            { ( Vector2d() << t1, 1 - 2 * t1 ).finished(),          v1 },
            { ( Vector2d() << t2, t2 ).finished(),                  v2 },
            { ( Vector2d() << 1 - 2 * t2, t2 ).finished(),          v2 },
            { ( Vector2d() << t2, 1 - 2 * t2 ).finished(),          v2 } } },
        // vertex rule; there is no Lobatto counterpart on triangles
        { Domain::Triangle, Family::NewtonCotes, 1, {
            { ( Vector2d() << 0., 0. ).finished(),                  1./6 },
            { ( Vector2d() << 1., 0. ).finished(),                  1./6 },
            { ( Vector2d() << 0., 1. ).finished(),                  1./6 } } },

        // ------------------------------------------------------------------------------------------------ Quadrilateral
        { Domain::Quadrilateral, Family::Gauss, 1, {
            { Vector2d::Zero(),                                     4. } } },
        { Domain::Quadrilateral, Family::Gauss, 3, {
            { ( Vector2d() << +gp2, +gp2 ).finished(),              1.0 },
            { ( Vector2d() << -gp2, +gp2 ).finished(),              1.0 },
            { ( Vector2d() << -gp2, -gp2 ).finished(),              1.0 },
            { ( Vector2d() << +gp2, -gp2 ).finished(),              1.0 } } },
        { Domain::Quadrilateral, Family::Gauss, 5, {
            { ( Vector2d() << 0, 0. ).finished(),                   64./81 },
            { ( Vector2d() << -gp3, -gp3 ).finished(),              25./81 },
            { ( Vector2d() << +gp3, -gp3 ).finished(),              25./81 },
            { ( Vector2d() << +gp3, +gp3 ).finished(),              25./81 },
            { ( Vector2d() << -gp3, +gp3 ).finished(),              25./81 },
            { ( Vector2d() << 0, -gp3 ).finished(),                 40./81 },
            { ( Vector2d() << gp3, 0. ).finished(),                 40./81 },
            { ( Vector2d() << 0, +gp3 ).finished(),                 40./81 },
            { ( Vector2d() << -gp3, 0. ).finished(),                40./81 } } },
        // tensor product of the trapezoidal rule
        { Domain::Quadrilateral, Family::NewtonCotes, 1, {
            { ( Vector2d() << -1., -1. ).finished(),                1.0 },
            { ( Vector2d() << +1., -1. ).finished(),                1.0 },
            { ( Vector2d() << +1., +1. ).finished(),                1.0 },
            { ( Vector2d() << -1., +1. ).finished(),                1.0 } } },

        // ------------------------------------------------------------------------------------------------ Tetrahedron
        { Domain::Tetrahedron, Family::Gauss, 1, {
            { ( Vector3d() << 1./4, 1./4, 1./4 ).finished(),        1./6 } } },
        { Domain::Tetrahedron, Family::Gauss, 2, {
            { ( Vector3d() << a, a, a ).finished(),                 1./24 },
            { ( Vector3d() << a, a, b ).finished(),                 1./24 },
            { ( Vector3d() << a, b, a ).finished(),                 1./24 },
            { ( Vector3d() << b, a, a ).finished(),                 1./24 } } },
        { Domain::Tetrahedron, Family::Gauss, 5, {
            { ( Vector3d() << a1, a1, a1 ).finished(),              w1 },
            { ( Vector3d() << a1, a1, r1 ).finished(),              w1 },
            { ( Vector3d() << a1, r1, a1 ).finished(),              w1 },
            { ( Vector3d() << r1, a1, a1 ).finished(),              w1 },
            { ( Vector3d() << a2, a2, a2 ).finished(),              w2 },
            { ( Vector3d() << a2, a2, r2 ).finished(),              w2 },
            { ( Vector3d() << a2, r2, a2 ).finished(),              w2 },
            { ( Vector3d() << r2, a2, a2 ).finished(),              w2 },
            { ( Vector3d() << c3, c3, r3 ).finished(),              w3 },
            { ( Vector3d() << c3, r3, c3 ).finished(),              w3 },
            { ( Vector3d() << r3, c3, c3 ).finished(),              w3 },
            { ( Vector3d() << c3, r3, r3 ).finished(),              w3 },
            { ( Vector3d() << r3, c3, r3 ).finished(),              w3 },
            { ( Vector3d() << r3, r3, c3 ).finished(),              w3 } } },

        // ------------------------------------------------------------------------------------------------ Hexahedron
        { Domain::Hexahedron, Family::Gauss, 1, {
            { Vector3d::Zero(),                                     8.0 } } },
        { Domain::Hexahedron, Family::Gauss, 3, {
            { ( Vector3d() << -gp2, -gp2, -gp2 ).finished(),        1.0 },
            { ( Vector3d() << +gp2, -gp2, -gp2 ).finished(),        1.0 },
            { ( Vector3d() << +gp2, +gp2, -gp2 ).finished(),        1.0 },
            { ( Vector3d() << -gp2, +gp2, -gp2 ).finished(),        1.0 },
            { ( Vector3d() << -gp2, -gp2, +gp2 ).finished(),        1.0 },
            { ( Vector3d() << +gp2, -gp2, +gp2 ).finished(),        1.0 },
            { ( Vector3d() << +gp2, +gp2, +gp2 ).finished(),        1.0 },
            { ( Vector3d() << -gp2, +gp2, +gp2 ).finished(),        1.0 } } },
        { Domain::Hexahedron, Family::Gauss, 5, {
            { ( Vector3d() << -gp3, -gp3, -gp3 ).finished(),        0.171467764060357 },
            { ( Vector3d() << 0, -gp3, -gp3 ).finished(),           0.274348422496571 },
            { ( Vector3d() << +gp3, -gp3, -gp3 ).finished(),        0.171467764060357 },
            { ( Vector3d() << -gp3, 0, -gp3 ).finished(),           0.274348422496571 },
            { ( Vector3d() << 0, 0, -gp3 ).finished(),              0.438957475994513 },
            { ( Vector3d() << gp3, 0, -gp3 ).finished(),            0.274348422496571 },
            { ( Vector3d() << -gp3, +gp3, -gp3 ).finished(),        0.171467764060357 },
            { ( Vector3d() << 0, +gp3, -gp3 ).finished(),           0.274348422496571 },
            { ( Vector3d() << +gp3, +gp3, -gp3 ).finished(),        0.171467764060357 },

            { ( Vector3d() << -gp3, -gp3, 0 ).finished(),           0.274348422496571 },
            { ( Vector3d() << 0, -gp3, 0 ).finished(),              0.438957475994513 },
            { ( Vector3d() << +gp3, -gp3, 0 ).finished(),           0.274348422496571 },
            { ( Vector3d() << -gp3, 0, 0 ).finished(),              0.438957475994513 },
            { ( Vector3d() << 0, 0, 0 ).finished(),                 0.702331961591221 },
            { ( Vector3d() << gp3, 0, 0 ).finished(),               0.438957475994513 },
            { ( Vector3d() << -gp3, +gp3, 0 ).finished(),           0.274348422496571 },
            { ( Vector3d() << 0, +gp3, 0 ).finished(),              0.438957475994513 },
            { ( Vector3d() << +gp3, +gp3, 0 ).finished(),           0.274348422496571 },

            { ( Vector3d() << -gp3, -gp3, +gp3 ).finished(),        0.171467764060357 },
            { ( Vector3d() << 0, -gp3, +gp3 ).finished(),           0.274348422496571 },
            { ( Vector3d() << +gp3, -gp3, +gp3 ).finished(),        0.171467764060357 },
            { ( Vector3d() << -gp3, 0, +gp3 ).finished(),           0.274348422496571 },
            { ( Vector3d() << 0, 0, +gp3 ).finished(),              0.438957475994513 },
            { ( Vector3d() << gp3, 0, +gp3 ).finished(),            0.274348422496571 },
            { ( Vector3d() << -gp3, +gp3, +gp3 ).finished(),        0.171467764060357 },
            { ( Vector3d() << 0, +gp3, +gp3 ).finished(),           0.274348422496571 },
            { ( Vector3d() << +gp3, +gp3, +gp3 ).finished(),        0.171467764060357 } } },
      };
      // clang-format on

      return rules;
    }

    const std::vector< QuadraturePointInfo >& getQuadraturePoints( ReferenceDomain domain, QuadratureRule rule )
    {
      // the first match is the cheapest, since the catalogue is sorted by number of points within domain and family
      for ( const auto& entry : quadratureRules() )
        if ( entry.domain == domain && entry.family == rule.family && entry.degree >= rule.degree )
          return entry.points;

      static const char* domainNames[] = { "Line", "Triangle", "Quadrilateral", "Tetrahedron", "Hexahedron" };
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << "getQuadraturePoints(): there is no "
                                   << ( rule.family == QuadratureFamily::Gauss ? "Gauss" : "Newton-Cotes" )
                                   << " rule of degree " << rule.degree << " on the reference domain "
                                   << domainNames[static_cast< int >( domain )] << "." );
    }

  } // namespace Quadrature
} // namespace Marmot::Discretizations::Elements
