#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>
#include <cmath>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Marmot::Testing;
using Marmot::Discretizations::ElementProperties;
using Marmot::Discretizations::MaterialSection;
using Marmot::Discretizations::Elements::Element;
using Marmot::Utilities::MakeString;

namespace {

  const std::vector< double > linearBond = { 100., 1000. }; // Kt, Kn

  /// unit host elements (Abaqus node orders): corners, then midsides
  const std::vector< double > quad4 = { 0., 0., 1., 0., 1., 1., 0., 1. };
  const std::vector< double > quad8 = { 0., 0., 1., 0., 1., 1., 0., 1., 0.5, 0., 1., 0.5, 0.5, 1., 0., 0.5 };
  const std::vector< double > tet4  = { 0., 0., 0., 1., 0., 0., 0., 1., 0., 0., 0., 1. };
  const std::vector< double > tet10 = { 0.,  0.,  0., 1., 0.,  0., 0., 1., 0.,  0.,  0., 1.,  0.5, 0.,  0.,
                                        0.5, 0.5, 0., 0., 0.5, 0., 0., 0., 0.5, 0.5, 0., 0.5, 0.,  0.5, 0.5 };
  const std::vector< double > hex8  = { 0., 0., 0., 1., 0., 0., 1., 1., 0., 0., 1., 0.,
                                        0., 0., 1., 1., 0., 1., 1., 1., 1., 0., 1., 1. };
  std::vector< double >       hex20()
  {
    std::vector< double > X = hex8;
    // midsides of the edges 1-2, 2-3, 3-4, 4-1, 5-6, 6-7, 7-8, 8-5, 1-5, 2-6, 3-7, 4-8
    const int e[12][2] = { { 0, 1 },
                           { 1, 2 },
                           { 2, 3 },
                           { 3, 0 },
                           { 4, 5 },
                           { 5, 6 },
                           { 6, 7 },
                           { 7, 4 },
                           { 0, 4 },
                           { 1, 5 },
                           { 2, 6 },
                           { 3, 7 } };
    for ( const auto& edge : e )
      for ( int i = 0; i < 3; i++ )
        X.push_back( 0.5 * ( hex8[edge[0] * 3 + i] + hex8[edge[1] * 3 + i] ) );
    return X;
  }

  /// the element types: name, dimension, bar nodes, host coordinates, bar (end, end[, mid]) inside the host
  struct Type {
    std::string           name;
    int                   nDim, nBarNodes;
    std::vector< double > host, bar;
  };
  const std::vector< Type > types = {
    { "EB2D2Q4", 2, 2, quad4, { 0.2, 0.3, 0.8, 0.6 } },
    { "EB2D3Q8", 2, 3, quad8, { 0.2, 0.3, 0.8, 0.6, 0.5, 0.45 } },
    { "EB3D2T4", 3, 2, tet4, { 0.1, 0.1, 0.1, 0.4, 0.2, 0.2 } },
    { "EB3D2T10", 3, 2, tet10, { 0.1, 0.1, 0.1, 0.4, 0.2, 0.2 } },
    { "EB3D2H8", 3, 2, hex8, { 0.2, 0.3, 0.4, 0.8, 0.6, 0.7 } },
    { "EB3D3H20", 3, 3, hex20(), { 0.2, 0.3, 0.4, 0.8, 0.6, 0.7, 0.5, 0.45, 0.55 } },
  };

  /// an element with all its data (the element keeps pointers to them)
  struct Bond {
    std::unique_ptr< Element > el;
    std::vector< double >      X, props, law, stateVars;

    Bond( const Type& t, std::vector< double > props_ )
      : el( Marmot::Registry::ElementFactory::createElement( t.name, 1 ) ),
        props( std::move( props_ ) ),
        law( linearBond )
    {
      X = t.bar;
      X.insert( X.end(), t.host.begin(), t.host.end() );
      el->assignNodeCoordinates( X.data() );
      el->assignProperty( ElementProperties( props.data(), static_cast< int >( props.size() ) ) );
      el->assignProperty( MaterialSection( "LINEARELASTICBONDSLIP", law.data(), static_cast< int >( law.size() ) ) );
      el->initializeYourself();
      stateVars.assign( el->getNumberOfRequiredStateVars(), 0.0 );
      el->assignStateVars( stateVars.data(), static_cast< int >( stateVars.size() ) );
      el->setInitialConditions( Element::MarmotMaterialInitialization, nullptr );
    }
  };

  template < class Exception >
  bool throws( const std::function< void() >& f )
  {
    try {
      f();
    }
    catch ( const Exception& ) {
      return true;
    }
    return false;
  }

} // namespace

/// the generic element interface, as used by a finite element framework, for bond elements with all host shapes
void testInterface()
{
  for ( const auto& t : types ) {
    const double perimeter = 0.3;
    Bond         bond( t, { perimeter, -1., 1. } );
    Element&     el     = *bond.el;
    const int    nHost  = static_cast< int >( t.host.size() ) / t.nDim;
    const int    nNodes = t.nBarNodes + nHost, nDof = t.nDim * nNodes;

    throwExceptionOnFailure( el.getNNodes() == nNodes && el.getNSpatialDimensions() == t.nDim &&
                               el.getNDofPerElement() == nDof,
                             t.name + ": nodes, dimension, dofs" );
    throwExceptionOnFailure( el.getElementShape() == ( t.nBarNodes == 2 ? "bar2" : "bar3" ), t.name + ": shape" );
    const auto fields = el.getNodeFields();
    throwExceptionOnFailure( static_cast< int >( fields.size() ) == nNodes &&
                               fields[0] == std::vector< std::string >{ "displacement" },
                             t.name + ": node fields" );
    const auto pattern = el.getDofIndicesPermutationPattern();
    for ( int i = 0; i < nDof; i++ )
      throwExceptionOnFailure( pattern[i] == i, t.name + ": dof pattern" );
    throwExceptionOnFailure( el.getNumberOfQuadraturePoints() == ( t.nBarNodes == 2 ? 2 : 3 ),
                             t.name + ": default number of integration points" );

    // coordinates: the center and the integration points lie on the bar
    Eigen::Vector3d a = Eigen::Vector3d::Zero(), d = Eigen::Vector3d::Zero();
    for ( int i = 0; i < t.nDim; i++ )
      a( i ) = t.bar[i], d( i ) = t.bar[t.nDim + i] - t.bar[i];
    auto onBar = [&]( const std::vector< double >& x ) {
      Eigen::Vector3d p = Eigen::Vector3d::Zero();
      for ( int i = 0; i < t.nDim; i++ )
        p( i ) = x[i];
      return ( p - a ).cross( d ).norm() < 1e-12 * d.squaredNorm();
    };
    throwExceptionOnFailure( onBar( el.getCoordinatesAtCenter() ), t.name + ": center on the bar" );
    for ( const auto& x : el.getCoordinatesAtQuadraturePoints() )
      throwExceptionOnFailure( onBar( x ), t.name + ": integration point on the bar" );

    // the bar slides relative to the host: explicit = implicit residual, elastic energy = P . U / 2 (linear law)
    Eigen::VectorXd U = Eigen::VectorXd::Zero( nDof ), P = Eigen::VectorXd::Zero( nDof ),
                    PExplicit = Eigen::VectorXd::Zero( nDof );
    Eigen::MatrixXd K         = Eigen::MatrixXd::Zero( nDof, nDof );
    for ( int n = 0; n < t.nBarNodes; n++ )
      for ( int i = 0; i < t.nDim; i++ )
        U( n * t.nDim + i ) = 1e-3 * ( i + 1 );
    const auto backup = bond.stateVars;
    el.computeKernelsExplicit( U.data(), U.data(), PExplicit.data(), 0.0, 1.0 );
    bond.stateVars = backup;
    el.computeKernels( U.data(), U.data(), P.data(), K.data(), 0.0, 1.0 );
    throwExceptionOnFailure( ( P - PExplicit ).norm() < 1e-12 * P.norm() && P.norm() > 0,
                             t.name + ": explicit kernels" );
    double energy = 0;
    el.computeInternalEnergy( energy );
    throwExceptionOnFailure( std::abs( energy - 0.5 * P.dot( U ) ) < 1e-6 * std::abs( energy ),
                             MakeString() << t.name << ": internal energy " << energy << " vs " << 0.5 * P.dot( U ) );

    // bond elements carry no mass and no body force
    Eigen::VectorXd Pb = Eigen::VectorXd::Zero( nDof ), m = Eigen::VectorXd::Zero( nDof ),
                    load = Eigen::VectorXd::Ones( t.nDim );
    Eigen::MatrixXd M    = Eigen::MatrixXd::Zero( nDof, nDof );
    el.computeBodyForce( Pb.data(), nullptr, load.data(), nullptr, 0.0, 1.0 );
    el.computeLumpedInertia( m.data() );
    el.computeConsistentInertia( M.data() );
    throwExceptionOnFailure( Pb.norm() == 0 && m.norm() == 0 && M.norm() == 0, t.name + ": no body force, no mass" );

    // states
    for ( const std::string name : { "slip", "bond stress", "elastic energy", "dissipation" } )
      throwExceptionOnFailure( el.getStateView( name, 0 ).stateSize >= 1, t.name + ": state " + name );

    // not supported or invalid
    throwExceptionOnFailure( throws< std::invalid_argument >( [&]() {
                               el.computeDistributedLoad( Element::Pressure,
                                                          P.data(),
                                                          K.data(),
                                                          0,
                                                          load.data(),
                                                          U.data(),
                                                          0.0,
                                                          1.0 );
                             } ),
                             t.name + ": distributed loads" );
    throwExceptionOnFailure( throws< std::invalid_argument >(
                               [&]() { el.setInitialConditions( Element::GeostaticStress, nullptr ); } ),
                             t.name + ": invalid initial condition" );
  }
}

/// invalid properties are rejected; 2 to 5 Gauss-Lobatto points are accepted
void testProperties()
{
  const Type& t = types[0];
  for ( int n = 2; n <= 5; n++ ) {
    Bond bond( t, { 0.3, -1., 1., double( n ) } );
    throwExceptionOnFailure( bond.el->getNumberOfQuadraturePoints() == n, MakeString() << n << " integration points" );
  }
  throwExceptionOnFailure( throws< std::invalid_argument >( [&]() {
                             Bond bond( t, { 0.3, -1., 1., 6. } );
                           } ),
                           "6 integration points must be rejected" );
  throwExceptionOnFailure( throws< std::invalid_argument >( [&]() {
                             Bond bond( t, { 0.3, -1. } );
                           } ),
                           "missing properties must be rejected" );
}

int main()
{
  executeTestsAndCollectExceptions( { testInterface, testProperties } );
  return 0;
}
