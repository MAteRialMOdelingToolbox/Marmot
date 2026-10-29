#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotSmallStrainKinematics.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>
#include <cmath>
#include <functional>
#include <string>
#include <vector>

using namespace Marmot::Testing;
namespace SmallStrain = Marmot::FiniteElement::SmallStrain;

namespace {

  // The B-operator of nDim spatial dimensions: the dense matrix the helpers replace.
  template < int nDim, int nNodes >
  Eigen::Matrix< double, SmallStrain::Voigt< nDim >::RowsAtCompileTime, nNodes * nDim > denseB(
    const Eigen::Matrix< double, nDim, nNodes >& dNdX )
  {
    if constexpr ( nDim == 1 )
      return dNdX;
    else if constexpr ( nDim == 2 )
      return Marmot::FiniteElement::Spatial2D::B< nNodes >( dNdX );
    else
      return Marmot::FiniteElement::Spatial3D::B< nNodes >( dNdX );
  }

  // Shape function gradients without any symmetry, so that a mixed-up index shows.
  template < int nDim, int nNodes >
  Eigen::Matrix< double, nDim, nNodes > someGradients()
  {
    Eigen::Matrix< double, nDim, nNodes > dNdX;
    for ( int i = 0; i < nDim; i++ )
      for ( int j = 0; j < nNodes; j++ )
        dNdX( i, j ) = std::sin( 1.7 * ( i + 1 ) + 0.9 * ( j + 1 ) );
    return dNdX;
  }

  template < int nDim, int nNodes >
  void checkAgainstTheDenseBOperator( const std::string& label )
  {
    const auto dNdX = someGradients< nDim, nNodes >();
    const auto B    = denseB< nDim, nNodes >( dNdX );

    Eigen::Matrix< double, nDim * nNodes, 1 > dQU;
    for ( int i = 0; i < nDim * nNodes; i++ )
      dQU( i ) = std::cos( 0.7 * ( i + 1 ) );

    const SmallStrain::Voigt< nDim > dE = SmallStrain::strainIncrement< nDim, nNodes >( dQU, dNdX );
    throwExceptionOnFailure( checkIfEqual( Eigen::MatrixXd( dE ), Eigen::MatrixXd( B * dQU ), 1e-13 ),
                             label + ": strainIncrement() does not equal B * dQU." );

    SmallStrain::Voigt< nDim > stress;
    for ( int i = 0; i < stress.size(); i++ )
      stress( i ) = std::sin( 2.9 * ( i + 1 ) );
    const double J0xW = 0.37;

    // starting from a non-zero force: the helper must add, not overwrite
    Eigen::Matrix< double, nDim * nNodes, 1 >       fU = Eigen::Matrix< double, nDim * nNodes, 1 >::Constant( 0.25 );
    const Eigen::Matrix< double, nDim * nNodes, 1 > expected = fU + B.transpose() * stress * J0xW;

    SmallStrain::addInternalForce< nDim, nNodes >( fU, dNdX, stress, J0xW );
    throwExceptionOnFailure( checkIfEqual( Eigen::MatrixXd( fU ), Eigen::MatrixXd( expected ), 1e-13 ),
                             label + ": addInternalForce() does not equal fU + B^T * sigma * J0xW." );
  }

} // namespace

void testHexa20()
{
  checkAgainstTheDenseBOperator< 3, 20 >( "3D, 20 nodes" );
}

void testHexa8()
{
  checkAgainstTheDenseBOperator< 3, 8 >( "3D, 8 nodes" );
}

void testQuad8()
{
  checkAgainstTheDenseBOperator< 2, 8 >( "2D, 8 nodes" );
}

void testQuad4()
{
  checkAgainstTheDenseBOperator< 2, 4 >( "2D, 4 nodes" );
}

void testBar3()
{
  checkAgainstTheDenseBOperator< 1, 3 >( "1D, 3 nodes" );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{ testHexa20, testHexa8, testQuad8, testQuad4, testBar3 };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
