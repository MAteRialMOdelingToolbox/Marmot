#include "Marmot/DisplacementCell.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotMPMLibrary.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>
#include <cmath>
#include <functional>
#include <memory>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;

namespace {

  // COMPRESSIBLENEOHOOKE: K, G, rho
  const std::vector< double > matProps = { 3500., 1500., 2.0 };

  constexpr int nDim   = 2;
  constexpr int nNodes = 4;
  constexpr int nDof   = nDim * nNodes;

  const std::vector< double > cellCoordinates = { 0.0, 0.0, 2.0, 0.0, 2.0, 1.0, 0.0, 1.0 };

  // one cell with 2x2 material points, as the host (EdelweissMeshfree) drives it
  struct Setup {
    std::unique_ptr< MarmotCell >                         cell;
    std::vector< std::unique_ptr< MarmotMaterialPoint > > mps;
    std::vector< std::vector< double > >                  stateVars;
    MarmotMaterialSection section{ "COMPRESSIBLENEOHOOKE", matProps.data(), (int)matProps.size() };

    Setup()
      : cell( MarmotLibrary::MarmotCellFactory::createCell( "Displacement/Quad4",
                                                            1,
                                                            cellCoordinates.data(),
                                                            cellCoordinates.size() ) )
    {
      int label = 1;
      for ( double y : { 0.25, 0.75 } )
        for ( double x : { 0.5, 1.5 } ) {
          const double xy[2] = { x, y };
          mps.emplace_back( MarmotLibrary::MarmotMaterialPointFactory::createMaterialPoint( "Displacement/PlaneStrain",
                                                                                            label++,
                                                                                            xy,
                                                                                            2,
                                                                                            0.5 ) );
        }

      for ( auto& mp : mps ) {
        mp->assignMaterial( section );
        stateVars.emplace_back( mp->getNumberOfRequiredStateVars(), 0.0 );
        mp->assignStateVars( stateVars.back().data(), stateVars.back().size() );
        mp->initializeYourself();
      }

      std::vector< MarmotMaterialPoint* > raw;
      for ( auto& mp : mps )
        raw.push_back( mp.get() );
      cell->assignMaterialPoints( raw );
    }

    // one MPM increment: reset the increment, interpolate, compute the points, assemble
    std::pair< Eigen::VectorXd, Eigen::MatrixXd > kernels( const Eigen::VectorXd& dQ )
    {
      for ( auto& mp : mps )
        mp->prepareYourself( 1.0, 1.0 );
      cell->interpolateFieldsToMaterialPoints( dQ.data() );
      for ( auto& mp : mps )
        mp->computeYourself( 1.0, 1.0 );

      Eigen::VectorXd P = Eigen::VectorXd::Zero( nDof );
      Eigen::MatrixXd K = Eigen::MatrixXd::Zero( nDof, nDof );
      cell->computeMaterialPointKernels( dQ.data(), P.data(), K.data(), 1.0, 1.0 );
      return { P, K };
    }

    // evaluate without committing: the material history is restored afterwards
    std::pair< Eigen::VectorXd, Eigen::MatrixXd > trial( const Eigen::VectorXd& dQ )
    {
      const auto backup = stateVars;
      auto       result = kernels( dQ );
      stateVars         = backup;
      for ( size_t i = 0; i < mps.size(); i++ )
        mps[i]->assignStateVars( stateVars[i].data(), stateVars[i].size() );
      return result;
    }

    void accept( const Eigen::VectorXd& dQ )
    {
      kernels( dQ );
      for ( auto& mp : mps )
        mp->acceptStateAndPosition();
    }
  };

  Eigen::VectorXd increment( double scaleU )
  {
    Eigen::VectorXd dQ( nDof );
    for ( int i = 0; i < nDof; i++ )
      dQ[i] = scaleU * std::sin( 1.7 * ( i + 1 ) );
    return dQ;
  }

  Eigen::MatrixXd numericalTangent( const std::function< Eigen::VectorXd( const Eigen::VectorXd& ) >& f,
                                    const Eigen::VectorXd&                                            Q0 )
  {
    const double    h = 1e-7;
    Eigen::MatrixXd numK( f( Q0 ).size(), Q0.size() );
    for ( int j = 0; j < Q0.size(); j++ ) {
      Eigen::VectorXd Qp = Q0, Qm = Q0;
      Qp[j] += h;
      Qm[j] -= h;
      numK.col( j ) = ( f( Qp ) - f( Qm ) ) / ( 2 * h );
    }
    return numK;
  }

} // namespace

void testLayout()
{
  Setup s;
  throwExceptionOnFailure( s.cell->getNDofPerCell() == nDof, "2 dofs per node expected" );
  throwExceptionOnFailure( s.cell->getNNodes() == nNodes, "Quad4 cell" );
  for ( const auto& f : s.cell->getNodeFields() )
    throwExceptionOnFailure( f == std::vector< std::string >{ "displacement" }, "node fields" );
}

void testConsistentTangentFirstStep()
{
  Setup                 s;
  const Eigen::VectorXd dQ = increment( 0.05 );
  const auto [P, K]        = s.trial( dQ );
  const auto numK          = numericalTangent( [&]( const Eigen::VectorXd& q ) { return s.trial( q ).first; }, dQ );

  const double err = ( K - numK ).norm() / numK.norm();
  throwExceptionOnFailure( err < 1e-7, MakeString() << "first step: tangent inconsistent, relative error " << err );
}

void testConsistentTangentAfterAcceptedStep()
{
  // the second increment starts from an accepted finite deformation (dY_dX != I)
  Setup s;
  s.accept( increment( 0.05 ) );

  const Eigen::VectorXd dQ = increment( 0.02 );
  const auto [P, K]        = s.trial( dQ );
  const auto numK          = numericalTangent( [&]( const Eigen::VectorXd& q ) { return s.trial( q ).first; }, dQ );

  const double err = ( K - numK ).norm() / numK.norm();
  throwExceptionOnFailure( err < 1e-7, MakeString() << "second step: tangent inconsistent, relative error " << err );
}

void testRigidRotationIsStressFree()
{
  Setup           s;
  const double    phi = 0.5;
  Eigen::Matrix2d R;
  R << std::cos( phi ), -std::sin( phi ), std::sin( phi ), std::cos( phi );

  Eigen::VectorXd dQ = Eigen::VectorXd::Zero( nDof );
  for ( int A = 0; A < nNodes; A++ )
    dQ.segment< 2 >( 2 * A ) = ( R - Eigen::Matrix2d::Identity() ) * Eigen::Vector2d( &cellCoordinates[2 * A] );

  const auto P = s.trial( dQ ).first;
  throwExceptionOnFailure( P.norm() < 1e-10,
                           MakeString() << "rigid rotation must not load the cell, |P| = " << P.norm() );
}

void testBodyForceAndInertia()
{
  Setup           s;
  const double    b[2] = { 0.3, -1.1 };
  Eigen::VectorXd fExt = Eigen::VectorXd::Zero( nDof );
  Eigen::MatrixXd K    = Eigen::MatrixXd::Zero( nDof, nDof );
  s.cell->computeBodyLoad( s.cell->getSupportedBodyLoadTypes().at( "BODYFORCE" ), b, fExt.data(), K.data(), 0., 1. );

  const double V = 4 * 0.5; // four points of volume 0.5
  for ( int i = 0; i < 2; i++ ) {
    double sum = 0;
    for ( int A = 0; A < nNodes; A++ )
      sum += fExt[2 * A + i];
    throwExceptionOnFailure( checkIfEqual( sum, -b[i] * V, 1e-12 ), "body force must add up to -b V" );
  }

  s.trial( increment( 0.0 ) ); // the points learn their density on the first computation
  Eigen::VectorXd I = Eigen::VectorXd::Zero( nDof );
  s.cell->computeLumpedInertia( I.data() );
  throwExceptionOnFailure( checkIfEqual( I.sum(), nDim * matProps[2] * V, 1e-12 ),
                           "lumped mass must add up to nDim rho V" );
}

void testPressureLoadTangent()
{
  Setup                 s;
  const Eigen::VectorXd dQ   = increment( 0.05 );
  const double          p[2] = { 0.0, -2.0 };

  auto load = [&]( const Eigen::VectorXd& q ) {
    s.trial( q ); // put the material points into the deformed state
    for ( auto& mp : s.mps )
      mp->prepareYourself( 1.0, 1.0 );
    s.cell->interpolateFieldsToMaterialPoints( q.data() );
    Eigen::VectorXd fExt = Eigen::VectorXd::Zero( nDof );
    Eigen::MatrixXd K    = Eigen::MatrixXd::Zero( nDof, nDof );
    s.cell->computeDistributedLoad( 0, 1, 3, p, fExt.data(), K.data(), 1.0, 1.0 );
    return std::make_pair( fExt, K );
  };

  const auto [f0, K0] = load( dQ );
  const auto numK     = numericalTangent( [&]( const Eigen::VectorXd& q ) { return load( q ).first; }, dQ );

  throwExceptionOnFailure( f0.norm() > 0, "pressure must load the cell" );
  throwExceptionOnFailure( ( K0 - numK ).norm() < 1e-7 * numK.norm(), "pressure load tangent inconsistent" );
}

int main()
{
  auto testFunctions = std::vector< std::function< void() > >{ testLayout,
                                                               testConsistentTangentFirstStep,
                                                               testConsistentTangentAfterAcceptedStep,
                                                               testRigidRotationIsStressFree,
                                                               testBodyForceAndInertia,
                                                               testPressureLoadTangent };
  executeTestsAndCollectExceptions( testFunctions );
  return 0;
}
