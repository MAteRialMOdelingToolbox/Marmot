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

  /// the element types with their spatial dimension, number of nodes, material and its properties [.., density]
  struct Type {
    std::string           name;
    int                   nDim, nNodes;
    std::string           material;
    std::vector< double > materialProperties;
  };
  const std::vector< double > linearElastic = { 210000., 0.3, 7.85e-9 };
  const std::vector< double > neoHooke      = { 3000., 1000., 1e-9 };
  const std::vector< Type >   types         = { { "T2D2", 2, 2, "LINEARELASTIC", linearElastic },
                                                { "T2D3", 2, 3, "LINEARELASTIC", linearElastic },
                                                { "T3D2", 3, 2, "LINEARELASTIC", linearElastic },
                                                { "T3D3", 3, 3, "LINEARELASTIC", linearElastic },
                                                { "T2D2CR", 2, 2, "COMPRESSIBLENEOHOOKE", neoHooke },
                                                { "T2D3CR", 2, 3, "COMPRESSIBLENEOHOOKE", neoHooke },
                                                { "T3D2CR", 3, 2, "COMPRESSIBLENEOHOOKE", neoHooke },
                                                { "T3D3CR", 3, 3, "COMPRESSIBLENEOHOOKE", neoHooke } };

  /// a straight bar from X1 to X2 (node order end, end, mid)
  std::vector< double > nodes( const Type& t, const std::vector< double >& X1, const std::vector< double >& X2 )
  {
    std::vector< double > X = X1;
    X.insert( X.end(), X2.begin(), X2.end() );
    if ( t.nNodes == 3 )
      for ( int i = 0; i < t.nDim; i++ )
        X.push_back( 0.5 * ( X1[i] + X2[i] ) );
    return X;
  }

  /// an element with all its data (the element keeps pointers to them)
  struct Truss {
    std::unique_ptr< Element > el;
    std::vector< double >      X, props, mat, stateVars;

    Truss( const Type& t, std::vector< double > X_, double area )
      : el( Marmot::Registry::ElementFactory::createElement( t.name, 1 ) ),
        X( std::move( X_ ) ),
        props( { area } ),
        mat( t.materialProperties )
    {
      el->assignNodeCoordinates( X.data() );
      el->assignProperty( ElementProperties( props.data(), 1 ) );
      el->initializeYourself();
      el->assignProperty( MaterialSection( t.material, mat.data(), static_cast< int >( mat.size() ) ) );
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

/// the generic element interface, as used by a finite element framework, for all truss types
void testInterface()
{
  for ( const auto& t : types ) {
    const std::vector< double > X1 = t.nDim == 2 ? std::vector< double >{ 1., 2. }
                                                 : std::vector< double >{ 1., 2., 3. };
    const std::vector< double > X2 = t.nDim == 2 ? std::vector< double >{ 4., 6. }
                                                 : std::vector< double >{ 4., 6., 3. };
    const double                A = 2.0, L = 5.0;
    Truss                       truss( t, nodes( t, X1, X2 ), A );
    Element&                    el   = *truss.el;
    const int                   nDof = t.nDim * t.nNodes;

    throwExceptionOnFailure( el.getNNodes() == t.nNodes && el.getNSpatialDimensions() == t.nDim &&
                               el.getNDofPerElement() == nDof,
                             t.name + ": nodes, dimension, dofs" );
    throwExceptionOnFailure( el.getElementShape() == ( t.nNodes == 2 ? "bar2" : "bar3" ), t.name + ": shape" );
    const auto fields = el.getNodeFields();
    throwExceptionOnFailure( static_cast< int >( fields.size() ) == t.nNodes &&
                               fields[0] == std::vector< std::string >{ "displacement" },
                             t.name + ": node fields" );
    const auto pattern = el.getDofIndicesPermutationPattern();
    for ( int i = 0; i < nDof; i++ )
      throwExceptionOnFailure( pattern[i] == i, t.name + ": dof pattern" );

    // coordinates: the center and the quadrature points lie on the bar
    const auto center = el.getCoordinatesAtCenter();
    for ( int i = 0; i < t.nDim; i++ )
      throwExceptionOnFailure( std::abs( center[i] - 0.5 * ( X1[i] + X2[i] ) ) < 1e-12, t.name + ": center" );
    const auto qps = el.getCoordinatesAtQuadraturePoints();
    throwExceptionOnFailure( static_cast< int >( qps.size() ) == el.getNumberOfQuadraturePoints(),
                             t.name + ": number of quadrature points" );
    for ( const auto& x : qps ) {
      Eigen::Vector3d a = Eigen::Vector3d::Zero(), d = Eigen::Vector3d::Zero(), p = Eigen::Vector3d::Zero();
      for ( int i = 0; i < t.nDim; i++ )
        a( i ) = X1[i], d( i ) = X2[i] - X1[i], p( i ) = x[i];
      throwExceptionOnFailure( ( p - a ).cross( d ).norm() < 1e-12 * d.squaredNorm(), t.name + ": qp on the bar" );
    }

    // mass: the consistent and the lumped mass matrix carry rho A L per direction
    const double    rho = t.materialProperties.back();
    Eigen::MatrixXd M   = Eigen::MatrixXd::Zero( nDof, nDof );
    el.computeConsistentInertia( M.data() );
    Eigen::VectorXd m = Eigen::VectorXd::Zero( nDof );
    el.computeLumpedInertia( m.data() );
    throwExceptionOnFailure( std::abs( M.sum() / t.nDim - rho * A * L ) < 1e-12 * rho * A * L &&
                               std::abs( m.sum() / t.nDim - rho * A * L ) < 1e-12 * rho * A * L && m.minCoeff() > 0,
                             MakeString()
                               << t.name << ": mass " << M.sum() / t.nDim << ", lumped " << m.sum() / t.nDim );

    // body force: the nodal forces add up to b A L
    Eigen::VectorXd Pb = Eigen::VectorXd::Zero( nDof ), load = Eigen::VectorXd::Zero( t.nDim );
    load( 1 ) = -9.81;
    el.computeBodyForce( Pb.data(), nullptr, load.data(), nullptr, 0.0, 1.0 );
    double sumY = 0;
    for ( int a = 0; a < t.nNodes; a++ )
      sumY += Pb( a * t.nDim + 1 );
    throwExceptionOnFailure( std::abs( sumY - load( 1 ) * A * L ) < 1e-12 * std::abs( load( 1 ) * A * L ),
                             MakeString() << t.name << ": body force resultant " << sumY );

    // a small axial stretch: explicit kernels give the implicit residual, and the elastic energy is P . U / 2
    Eigen::VectorXd U = Eigen::VectorXd::Zero( nDof ), P = Eigen::VectorXd::Zero( nDof ),
                    PExplicit = Eigen::VectorXd::Zero( nDof );
    Eigen::MatrixXd K         = Eigen::MatrixXd::Zero( nDof, nDof );
    for ( int i = 0; i < t.nDim; i++ ) {
      U( t.nDim + i ) = 1e-4 * ( X2[i] - X1[i] );
      if ( t.nNodes == 3 )
        U( 2 * t.nDim + i ) = 0.5e-4 * ( X2[i] - X1[i] );
    }
    const auto backup = truss.stateVars;
    el.computeKernelsExplicit( U.data(), U.data(), PExplicit.data(), 0.0, 1.0 );
    truss.stateVars = backup;
    el.computeKernels( U.data(), U.data(), P.data(), K.data(), 0.0, 1.0 );
    throwExceptionOnFailure( ( P - PExplicit ).norm() < 1e-10 * P.norm(), t.name + ": explicit kernels" );
    double energy = 0;
    el.computeInternalEnergy( energy );
    throwExceptionOnFailure( std::abs( energy - 0.5 * P.dot( U ) ) < 1e-3 * std::abs( energy ),
                             MakeString() << t.name << ": internal energy " << energy << " vs " << 0.5 * P.dot( U ) );

    // states of the element
    for ( const std::string name : { "stress", "strain", "normal force", "elastic energy", "dissipation" } )
      throwExceptionOnFailure( el.getStateView( name, 0 ).stateSize == 1, t.name + ": state " + name );

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
    throwExceptionOnFailure( throws< std::invalid_argument >( [&]() {
                               auto other = std::unique_ptr< Element >(
                                 Marmot::Registry::ElementFactory::createElement( t.name, 2 ) );
                               other->assignProperty( ElementProperties( nullptr, 0 ) );
                             } ),
                             t.name + ": missing cross section area" );
  }
}

int main()
{
  executeTestsAndCollectExceptions( { testInterface } );
  return 0;
}
