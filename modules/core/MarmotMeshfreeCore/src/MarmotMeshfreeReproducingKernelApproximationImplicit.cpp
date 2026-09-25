#include "Marmot/MarmotMeshfreeReproducingKernelApproximationImplicit.h"
#include "Marmot/MarmotMonomialBasisFunctions.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/src/Core/Matrix.h>
#include <cmath>

namespace Marmot::Meshfree {
  MarmotMeshfreeReproducingKernelApproximationImplicit::MarmotMeshfreeReproducingKernelApproximationImplicit(
    int dim,
    int completenessOrder )
    : MarmotMeshfreeReproducingKernelApproximation( dim, completenessOrder )
  {
  }

  void MarmotMeshfreeReproducingKernelApproximationImplicit::computeShapeFunctionGradients(
    const double*                                             coord,
    const std::vector< const MarmotMeshfreeKernelFunction* >& kernelFunctions,
    double*                                                   shapeFunctionValueGradients ) const
  {
    throw std::runtime_error( "Not implemented" );
  }

  void MarmotMeshfreeReproducingKernelApproximationImplicit::computeShapeFunctionsAndGradients(
    const double*                                             coord,
    const std::vector< const MarmotMeshfreeKernelFunction* >& kernelFunctionCandidates,
    double*                                                   shapeFunctionValues_,
    double*                                                   shapeFunctionValueGradients_ ) const
  {

    const auto coveringKernelFunctionIndices = findCoveringKernelFunctionIndices( coord, kernelFunctionCandidates );

    std::vector< const MarmotMeshfreeKernelFunction* > coveringKernelFunctions;
    for ( const auto& idx : coveringKernelFunctionIndices )
      coveringKernelFunctions.push_back( kernelFunctionCandidates[idx] );

    const auto correctedCompletenessOrder = getCorrectedCompletenessOrder( coveringKernelFunctionIndices.size() );

    const auto sizeH = Math::computeSizeOfMonomialBasisVector( correctedCompletenessOrder, _dim );

    if ( sizeH < 1 ) {
      throw std::runtime_error( "Size of H vector is less than 1" );
    }

    const Eigen::Map< const Eigen::VectorXd > coordVec( coord, _dim );

    Eigen::Map< Eigen::VectorXd > shapeFunctionValues( shapeFunctionValues_, kernelFunctionCandidates.size() );
    shapeFunctionValues.setZero();

    Eigen::Map< Eigen::MatrixXd > shapeFunctionValueGradients( shapeFunctionValueGradients_,
                                                               _dim,
                                                               kernelFunctionCandidates.size() );
    shapeFunctionValueGradients.setZero();

    const auto M = computeMMatrix( coordVec, kernelFunctionCandidates, correctedCompletenessOrder );

    // solve for b(x)
    // b = M^-1 * H0
    const auto H0 = H0Vector( M.rows() );

    using namespace Eigen;

    // first column is for the base function
    // the rest of the columns are for the gradients
    MatrixXd H0Mat( M.rows(), 1 + _dim );
    H0Mat.setZero();

    H0Mat.col( 0 ) = H0;
    H0Mat.block( 1, 1, _dim, _dim ).diagonal().setConstant( -1.0 );

    const Eigen::MatrixXd bMat = M.colPivHouseholderQr().solve( H0Mat );

    const VectorXd b0    = bMat.col( 0 );
    const MatrixXd bGrad = bMat.block( 0, 1, M.rows(), _dim );

    for ( const auto& A : coveringKernelFunctionIndices ) {
      const auto H = computeHVector( coordVec - Eigen::Map< const Eigen::VectorXd >( kernelFunctionCandidates[A]
                                                                                       ->getCenterCoordinates(),
                                                                                     _dim ),
                                     coveringKernelFunctions,
                                     correctedCompletenessOrder );

      const auto phi_A = kernelFunctionCandidates[A]->computeKernelFunction( coord );

      shapeFunctionValues[A]               = b0.dot( H ) * phi_A;
      shapeFunctionValueGradients.col( A ) = bGrad.transpose() * H * phi_A;
    }
  }

}; // namespace Marmot::Meshfree
