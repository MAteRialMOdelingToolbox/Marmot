#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>
#include <cmath>
#include <functional>
#include <memory>
#include <numbers>
#include <string>
#include <vector>

using namespace Marmot::Testing;
using Marmot::Discretizations::ElementProperties;
using Marmot::Discretizations::MaterialSection;
using Marmot::Discretizations::Elements::Element;
using Marmot::Utilities::MakeString;

namespace {

  const std::vector< double > linearElastic = { 200000., 0.3 };

  /// the element types with their spatial dimension and number of nodes
  struct Type {
    std::string name;
    int         nDim, nNodes;
  };
  const std::vector< Type > types =
    { { "B2D2", 2, 2 }, { "B2D3", 2, 3 }, { "B3D2", 3, 2 }, { "B3D3", 3, 3 }, { "B2D2CR", 2, 2 }, { "B3D2CR", 3, 2 } };

  /// a straight beam from X1 to X2 (mid node at the middle for 3 nodes)
  std::vector< double > nodes( const Type& t, const std::vector< double >& X1, const std::vector< double >& X2 )
  {
    std::vector< double > X = X1;
    X.insert( X.end(), X2.begin(), X2.end() );
    if ( t.nNodes == 3 )
      for ( int i = 0; i < t.nDim; i++ )
        X.push_back( 0.5 * ( X1[i] + X2[i] ) );
    return X;
  }

  /// a 2 x 2 point rectangle 20 x 10 (y along the local y axis)
  std::vector< double > properties( int nDim, double J = 4000. )
  {
    if ( nDim == 2 )
      return { 2., -5., 100., 5., 100. };
    return { 0., 0., 1., J, 4., -5., -2.5, 50., 5., -2.5, 50., -5., 2.5, 50., 5., 2.5, 50. };
  }

  /// an element with all its data (the element keeps pointers to them)
  struct Beam {
    std::unique_ptr< Element > el;
    std::vector< double >      X, props, mat, stateVars;

    Beam( const Type& t, std::vector< double > X_, std::vector< double > props_ )
      : el( Marmot::Registry::ElementFactory::createElement( t.name, 1 ) ),
        X( std::move( X_ ) ),
        props( std::move( props_ ) ),
        mat( linearElastic )
    {
      el->assignNodeCoordinates( X.data() );
      el->assignProperty( ElementProperties( props.data(), static_cast< int >( props.size() ) ) );
      el->initializeYourself();
      el->assignProperty( MaterialSection( "LINEARELASTIC", mat.data(), static_cast< int >( mat.size() ) ) );
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

/// the generic element interface, as used by a finite element framework, for all beam types
void testInterface()
{
  for ( const auto& t : types ) {
    const std::vector< double > X1 = t.nDim == 2 ? std::vector< double >{ 1., 2. }
                                                 : std::vector< double >{ 1., 2., 3. };
    const std::vector< double > X2 = t.nDim == 2 ? std::vector< double >{ 61., 82. }
                                                 : std::vector< double >{ 61., 82., 3. };
    Beam                        b( t, nodes( t, X1, X2 ), properties( t.nDim ) );
    Element&                    el          = *b.el;
    const int                   nDofPerNode = t.nDim == 2 ? 3 : 6, nDof = nDofPerNode * t.nNodes;

    throwExceptionOnFailure( el.getNNodes() == t.nNodes && el.getNSpatialDimensions() == t.nDim &&
                               el.getNDofPerElement() == nDof,
                             t.name + ": nodes, dimension, dofs" );
    throwExceptionOnFailure( el.getElementShape() == ( t.nNodes == 2 ? "bar2" : "bar3" ), t.name + ": shape" );
    const auto fields = el.getNodeFields();
    throwExceptionOnFailure( static_cast< int >( fields.size() ) == t.nNodes &&
                               fields[0] == std::vector< std::string >{ "displacement", "rotation" },
                             t.name + ": node fields" );
    const auto pattern = el.getDofIndicesPermutationPattern();
    for ( int i = 0; i < nDof; i++ )
      throwExceptionOnFailure( pattern[i] == i, t.name + ": dofs ordered node by node" );

    // coordinates: center and quadrature points on the axis
    const auto center = el.getCoordinatesAtCenter();
    for ( int i = 0; i < t.nDim; i++ )
      throwExceptionOnFailure( std::abs( center[i] - 0.5 * ( X1[i] + X2[i] ) ) < 1e-12, t.name + ": center" );
    const auto qps = el.getCoordinatesAtQuadraturePoints();
    throwExceptionOnFailure( static_cast< int >( qps.size() ) == el.getNumberOfQuadraturePoints() &&
                               el.getNumberOfQuadraturePoints() == ( t.nNodes == 2 ? 3 : 4 ),
                             t.name + ": quadrature points" );
    for ( const auto& x : qps ) {
      // on the segment X1 X2: the distance to the line vanishes
      Eigen::Vector3d a = Eigen::Vector3d::Zero(), d = Eigen::Vector3d::Zero(), p = Eigen::Vector3d::Zero();
      for ( int i = 0; i < t.nDim; i++ )
        a( i ) = X1[i], d( i ) = X2[i] - X1[i], p( i ) = x[i];
      throwExceptionOnFailure( ( p - a ).cross( d ).norm() < 1e-9 * d.squaredNorm(), t.name + ": qp on the axis" );
    }

    // a small axial stretch: explicit kernels give the implicit residual, and the elastic energy is P . U / 2
    Eigen::VectorXd U = Eigen::VectorXd::Zero( nDof ), dU, P = Eigen::VectorXd::Zero( nDof ),
                    PExplicit = Eigen::VectorXd::Zero( nDof );
    Eigen::MatrixXd K         = Eigen::MatrixXd::Zero( nDof, nDof );
    for ( int i = 0; i < t.nDim; i++ )
      U( nDofPerNode + i ) = 1e-3 * ( X2[i] - X1[i] );
    if ( t.nNodes == 3 )
      for ( int i = 0; i < t.nDim; i++ )
        U( 2 * nDofPerNode + i ) = 0.5e-3 * ( X2[i] - X1[i] );
    const auto backup = b.stateVars;
    dU                = U;
    el.computeKernelsExplicit( U.data(), dU.data(), PExplicit.data(), 0.0, 1.0 );
    b.stateVars = backup;
    el.computeKernels( U.data(), dU.data(), P.data(), K.data(), 0.0, 1.0 );
    throwExceptionOnFailure( ( P - PExplicit ).norm() < 1e-10 * P.norm(), t.name + ": explicit kernels" );
    double energy = 0;
    el.computeInternalEnergy( energy );
    throwExceptionOnFailure( std::abs( energy - 0.5 * P.dot( U ) ) < 1e-6 * std::abs( energy ),
                             MakeString() << t.name << ": internal energy " << energy << " vs " << 0.5 * P.dot( U ) );

    // states: all named states, the section forces / strains, a fiber
    const std::vector< std::string > named = t.nDim == 2 ? std::vector< std::string >{ "normal force",
                                                                                       "bending moment",
                                                                                       "axial strain",
                                                                                       "curvature",
                                                                                       "elastic energy",
                                                                                       "dissipation" }
                                                         : std::vector< std::string >{ "normal force",
                                                                                       "bending moments",
                                                                                       "torque",
                                                                                       "axial strain",
                                                                                       "curvatures",
                                                                                       "twist rate",
                                                                                       "elastic energy",
                                                                                       "dissipation" };
    for ( const auto& name : named )
      throwExceptionOnFailure( el.getStateView( name, 0 ).stateSize >= 1, t.name + ": state " + name );
    throwExceptionOnFailure( el.getStateView( "section forces", 0 ).stateSize == ( t.nDim == 2 ? 2 : 4 ),
                             t.name + ": section forces" );
    throwExceptionOnFailure( el.getStateView( "section strains", 0 ).stateSize == ( t.nDim == 2 ? 2 : 4 ),
                             t.name + ": section strains" );
    throwExceptionOnFailure( el.getStateView( "fiber 0 stress", 0 ).stateSize == 3, t.name + ": fiber stress" );
    throwExceptionOnFailure( throws< std::invalid_argument >( [&]() { el.getStateView( "no such state", 0 ); } ),
                             t.name + ": unknown state" );
    throwExceptionOnFailure( throws< std::invalid_argument >( [&]() { el.getStateView( "fiber 99 stress", 0 ); } ),
                             t.name + ": unknown fiber" );

    // body force: the nodal forces add up to b A L
    Eigen::VectorXd Pb = Eigen::VectorXd::Zero( nDof ), load( t.nDim );
    load.setZero();
    load( 1 ) = 7.85e-9 * 9810.;
    el.computeBodyForce( Pb.data(), nullptr, load.data(), nullptr, 0.0, 1.0 );
    double sumY = 0;
    for ( int a = 0; a < t.nNodes; a++ )
      sumY += Pb( a * nDofPerNode + 1 );
    const double L = std::sqrt( 60. * 60. + 80. * 80. ), A = 200.;
    throwExceptionOnFailure( std::abs( sumY - load( 1 ) * A * L ) < 1e-10 * load( 1 ) * A * L,
                             MakeString() << t.name << ": body force resultant " << sumY );

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
    throwExceptionOnFailure( throws< std::invalid_argument >( [&]() { el.assignStateVars( b.stateVars.data(), 1 ); } ),
                             t.name + ": too few state variables" );
  }
}

/// invalid geometry and section input is rejected
void testInvalidInput()
{
  for ( const auto& t : types ) {
    const std::vector< double > X1 = t.nDim == 2 ? std::vector< double >{ 0., 0. }
                                                 : std::vector< double >{ 0., 0., 0. };
    const std::vector< double > X2 = t.nDim == 2 ? std::vector< double >{ 10., 0. }
                                                 : std::vector< double >{ 10., 0., 0. };
    auto make = [&]( std::vector< double > X, std::vector< double > props ) { Beam b( t, X, props ); };

    // wrong number of properties
    auto wrong = properties( t.nDim );
    wrong.pop_back();
    throwExceptionOnFailure( throws< std::invalid_argument >( [&]() { make( nodes( t, X1, X2 ), wrong ); } ),
                             t.name + ": wrong number of properties" );
    // a section point without area
    auto noArea   = properties( t.nDim );
    noArea.back() = 0.;
    throwExceptionOnFailure( throws< std::invalid_argument >( [&]() { make( nodes( t, X1, X2 ), noArea ); } ),
                             t.name + ": section point without area" );
    // zero length
    throwExceptionOnFailure( throws< std::invalid_argument >(
                               [&]() { make( nodes( t, X1, X1 ), properties( t.nDim ) ); } ),
                             t.name + ": zero length" );
    if ( t.nDim == 3 ) // no torsion constant
      throwExceptionOnFailure( throws< std::invalid_argument >(
                                 [&]() { make( nodes( t, X1, X2 ), properties( t.nDim, 0. ) ); } ),
                               t.name + ": torsion constant" );
    if ( t.nNodes == 3 ) { // mid node off the middle
      auto X = nodes( t, X1, X2 );
      X.back() += 1.;
      throwExceptionOnFailure( throws< std::invalid_argument >( [&]() { make( X, properties( t.nDim ) ); } ),
                               t.name + ": mid node off the middle" );
    }
  }
}

/// co-rotational beams request a cutback (StressUpdateFailed) for degenerate deformations
void testCorotationalCutbacks()
{
  {
    const Type      t = types[4]; // B2D2CR: the second node moved onto the first one
    Beam            b( t, { 0., 0., 10., 0. }, properties( 2 ) );
    Eigen::VectorXd U = Eigen::VectorXd::Zero( 6 ), P = Eigen::VectorXd::Zero( 6 );
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero( 6, 6 );
    U( 3 )            = -10.;
    throwExceptionOnFailure( throws< Marmot::Utilities::StressUpdateFailed >(
                               [&]() { b.el->computeKernels( U.data(), U.data(), P.data(), K.data(), 0.0, 1.0 ); } ),
                             "B2D2CR: collapsed chord" );
  }
  {
    const Type t = types[5]; // B3D2CR
    for ( int component : { 0, 2 } ) {
      // node 2 rotated by pi about x (twist: the co-rotated frame is undefined) or about z, the local y axis of this
      // section (bending: the local rotation reaches pi)
      Beam            b( t, { 0., 0., 0., 10., 0., 0. }, properties( 3 ) );
      Eigen::VectorXd U = Eigen::VectorXd::Zero( 12 ), P = Eigen::VectorXd::Zero( 12 );
      Eigen::MatrixXd K  = Eigen::MatrixXd::Zero( 12, 12 );
      U( 9 + component ) = std::numbers::pi;
      throwExceptionOnFailure( throws< Marmot::Utilities::StressUpdateFailed >(
                                 [&]() { b.el->computeKernels( U.data(), U.data(), P.data(), K.data(), 0.0, 1.0 ); } ),
                               MakeString() << "B3D2CR: rotation by pi about axis " << component );
    }
  }
}

int main()
{
  executeTestsAndCollectExceptions( { testInterface, testInvalidInput, testCorotationalCutbacks } );
  return 0;
}
