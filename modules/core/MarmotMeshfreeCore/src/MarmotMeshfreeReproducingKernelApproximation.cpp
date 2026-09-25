#include "Marmot/MarmotMeshfreeReproducingKernelApproximation.h"
#include "Marmot/MarmotMeshfreeKernelFunction.h"
#include "Marmot/MarmotMonomialBasisFunctions.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/src/Core/Matrix.h>
#include <cmath>
#include <vector>

namespace Marmot::Meshfree {
  MarmotMeshfreeReproducingKernelApproximation::MarmotMeshfreeReproducingKernelApproximation( int dim,
                                                                                              int completenessOrder )
    : _dim( dim ), _desiredCompletenessOrder( completenessOrder )
  {
  }

  Eigen::VectorXd MarmotMeshfreeReproducingKernelApproximation::computeHVector(
    const Eigen::VectorXd&                                    x_minus_xI,
    const std::vector< const MarmotMeshfreeKernelFunction* >& coveringShapeFunctions,
    const int                                                 completenessOrder )
  {
    const auto _sizeH = Math::computeSizeOfMonomialBasisVector( completenessOrder, x_minus_xI.size() );

    Eigen::VectorXd res = Eigen::VectorXd::Ones( _sizeH );

    Math::computeMonomialBasis( completenessOrder, x_minus_xI, res );

    return res;
  }

  Eigen::MatrixXd MarmotMeshfreeReproducingKernelApproximation::computeHVectorGradient(
    const Eigen::VectorXd&                                    x_minus_xI,
    const std::vector< const MarmotMeshfreeKernelFunction* >& coveringShapeFunctions,
    const int                                                 completenessOrder )
  {
    const auto _sizeH = Math::computeSizeOfMonomialBasisVector( completenessOrder, x_minus_xI.size() );

    Eigen::MatrixXd res = Eigen::MatrixXd::Ones( _sizeH, x_minus_xI.size() );

    Math::computeMonomialBasisGradient( completenessOrder, x_minus_xI, res );

    return res;
  }

  Eigen::VectorXd MarmotMeshfreeReproducingKernelApproximation::H0Vector( int sizeHVector )
  {
    Eigen::VectorXd H0 = Eigen::VectorXd::Zero( sizeHVector );
    H0( 0 )            = 1;
    return H0;
  }

  Eigen::MatrixXd MarmotMeshfreeReproducingKernelApproximation::computeMMatrix(
    const Eigen::VectorXd&                                    globalCoord,
    const std::vector< const MarmotMeshfreeKernelFunction* >& coveringShapeFunctions,
    int                                                       completenessOrder )
  {
    const auto _sizeH = Math::computeSizeOfMonomialBasisVector( completenessOrder, globalCoord.size() );

    Eigen::MatrixXd M = Eigen::MatrixXd::Zero( _sizeH, _sizeH );

    for ( const auto& coveringShapeFunction : coveringShapeFunctions ) {
      const Eigen::Map< const Eigen::VectorXd > center( coveringShapeFunction->getCenterCoordinates(),
                                                        globalCoord.size() );

      const Eigen::VectorXd Hx = computeHVector( globalCoord - center, coveringShapeFunctions, completenessOrder );

      M += Hx * Hx.transpose() * coveringShapeFunction->computeKernelFunction( globalCoord.data() );
    }

    return M;
  }

  std::pair< Eigen::MatrixXd, std::vector< Eigen::MatrixXd > > MarmotMeshfreeReproducingKernelApproximation::
    computeMMatrixAndGradient( const Eigen::VectorXd&                                    globalCoord,
                               const std::vector< const MarmotMeshfreeKernelFunction* >& coveringShapeFunctions,
                               int                                                       completenessOrder )
  {

    const int  _dim   = globalCoord.size();
    const auto _sizeH = Math::computeSizeOfMonomialBasisVector( completenessOrder, _dim );

    Eigen::MatrixXd                M = Eigen::MatrixXd::Zero( _sizeH, _sizeH );
    std::vector< Eigen::MatrixXd > MGradients( _dim );
    for ( auto& MGradient : MGradients )
      MGradient = Eigen::MatrixXd::Zero( _sizeH, _sizeH );

    for ( const auto& coveringShapeFunction : coveringShapeFunctions ) {
      const Eigen::Map< const Eigen::VectorXd > center( coveringShapeFunction->getCenterCoordinates(),
                                                        globalCoord.size() );

      const Eigen::VectorXd Hx = computeHVector( globalCoord - center, coveringShapeFunctions, completenessOrder );
      const Eigen::MatrixXd HxGradient = computeHVectorGradient( globalCoord - center,
                                                                 coveringShapeFunctions,
                                                                 completenessOrder );

      const double    phi         = coveringShapeFunction->computeKernelFunction( globalCoord.data() );
      Eigen::VectorXd phiGradient = Eigen::VectorXd::Zero( _dim );
      coveringShapeFunction->computeKernelFunctionGradient( globalCoord.data(), phiGradient.data() );

      const Eigen::MatrixXd HxHT = Hx * Hx.transpose();
      M += HxHT * phi;

      for ( int i = 0; i < _dim; i++ ) {
        const Eigen::MatrixXd HGrad_i_xHT = HxGradient.col( i ) * Hx.transpose();
        MGradients[i] += ( HGrad_i_xHT + HGrad_i_xHT.transpose() ) * phi + HxHT * phiGradient( i );
      }
    }

    return { M, MGradients };
  }

  void MarmotMeshfreeReproducingKernelApproximation::computeShapeFunctions(
    const double*                                             coord,
    const std::vector< const MarmotMeshfreeKernelFunction* >& kernelFunctionCandidates,
    double*                                                   shapeFunctionValues ) const
  {
    const auto coveringKernelFunctionsIndices = findCoveringKernelFunctionIndices( coord, kernelFunctionCandidates );

    const auto correctedCompletenessOrder = getCorrectedCompletenessOrder( coveringKernelFunctionsIndices.size() );

    std::vector< const MarmotMeshfreeKernelFunction* > coveringKernelFunctions;
    for ( const auto& idx : coveringKernelFunctionsIndices )
      coveringKernelFunctions.push_back( kernelFunctionCandidates[idx] );

    /* // MAp: */
    const Eigen::Map< const Eigen::VectorXd > coordVec( coord, _dim );

    const auto M = computeMMatrix( coordVec, coveringKernelFunctions, correctedCompletenessOrder );

    // solve for b(x)
    // b = M^-1 * H0
    const auto H0 = H0Vector( M.rows() );

    const Eigen::VectorXd b = M.colPivHouseholderQr().solve( H0 );

    // compute the shape function values

    for ( int A = 0; A < static_cast< int >( kernelFunctionCandidates.size() ); A++ )
      shapeFunctionValues[A] = 0;

    for ( const auto& A : coveringKernelFunctionsIndices ) {

      const auto phi_A = kernelFunctionCandidates[A]->computeKernelFunction( coord );

      const auto H = computeHVector( coordVec - Eigen::Map< const Eigen::VectorXd >( kernelFunctionCandidates[A]
                                                                                       ->getCenterCoordinates(),
                                                                                     _dim ),
                                     coveringKernelFunctions,
                                     correctedCompletenessOrder );

      shapeFunctionValues[A] = b.dot( H ) * phi_A;
    }
  }

  void MarmotMeshfreeReproducingKernelApproximation::computeShapeFunctionGradients(
    const double*                                             coord,
    const std::vector< const MarmotMeshfreeKernelFunction* >& kernelFunctions,
    double*                                                   shapeFunctionValueGradients ) const
  {
    throw std::runtime_error( "Not implemented" );
  }

  const std::vector< int > MarmotMeshfreeReproducingKernelApproximation::findCoveringKernelFunctionIndices(
    const double*                                             coord,
    const std::vector< const MarmotMeshfreeKernelFunction* >& kernelFunctions ) const

  {
    std::vector< int > coveringFunctionsIndices;

    for ( int i = 0; i < static_cast< int >( kernelFunctions.size() ); i++ )
      if ( std::abs( kernelFunctions[i]->computeKernelFunction( coord ) ) > 1e-14 )
        coveringFunctionsIndices.push_back( i );

    return coveringFunctionsIndices;
  }

  void MarmotMeshfreeReproducingKernelApproximation::computeShapeFunctionsAndGradients(
    const double*                                             coord,
    const std::vector< const MarmotMeshfreeKernelFunction* >& kernelFunctionCandidates,
    double*                                                   shapeFunctionValues,
    double*                                                   shapeFunctionValueGradients_ ) const
  {
    const auto coveringKernelFunctionIndices = findCoveringKernelFunctionIndices( coord, kernelFunctionCandidates );

    std::vector< const MarmotMeshfreeKernelFunction* > coveringKernelFunctions;
    for ( const auto& idx : coveringKernelFunctionIndices )
      coveringKernelFunctions.push_back( kernelFunctionCandidates[idx] );

    const auto correctedCompletenessOrder = getCorrectedCompletenessOrder( coveringKernelFunctionIndices.size() );

    const Eigen::Map< const Eigen::VectorXd > coordVec( coord, _dim );
    const auto sizeH = Math::computeSizeOfMonomialBasisVector( correctedCompletenessOrder, _dim );

    if ( sizeH < 1 ) {
      throw std::runtime_error( "Size of H vector is less than 1" );
    }

    Eigen::Map< Eigen::MatrixXd > shapeFunctionValueGradients( shapeFunctionValueGradients_,
                                                               _dim,
                                                               kernelFunctionCandidates.size() );

    // Per-covering-kernel quantities are computed once here and reused for both the moment-matrix
    // assembly and the final shape-function/gradient evaluation, instead of being recomputed twice.
    const size_t                   nCovering = coveringKernelFunctionIndices.size();
    std::vector< Eigen::VectorXd > H_cache( nCovering );
    std::vector< Eigen::MatrixXd > HGradient_cache( nCovering );
    std::vector< double >          phi_cache( nCovering );
    std::vector< Eigen::VectorXd > phiGradient_cache( nCovering );

    Eigen::MatrixXd                M = Eigen::MatrixXd::Zero( sizeH, sizeH );
    std::vector< Eigen::MatrixXd > MGradients( _dim, Eigen::MatrixXd::Zero( sizeH, sizeH ) );

    for ( size_t k = 0; k < nCovering; k++ ) {
      const auto* kf = coveringKernelFunctions[k];

      const Eigen::VectorXd x_minus_center = coordVec -
                                             Eigen::Map< const Eigen::VectorXd >( kf->getCenterCoordinates(), _dim );

      const double    phi         = kf->computeKernelFunction( coord );
      Eigen::VectorXd phiGradient = Eigen::VectorXd::Zero( _dim );
      kf->computeKernelFunctionGradient( coord, phiGradient.data() );

      Eigen::VectorXd H         = computeHVector( x_minus_center, coveringKernelFunctions, correctedCompletenessOrder );
      Eigen::MatrixXd HGradient = computeHVectorGradient( x_minus_center,
                                                          coveringKernelFunctions,
                                                          correctedCompletenessOrder );

      const Eigen::MatrixXd HxHT = H * H.transpose();
      M += HxHT * phi;
      for ( int i = 0; i < _dim; i++ ) {
        const Eigen::MatrixXd HGrad_i_xHT = HGradient.col( i ) * H.transpose();
        MGradients[i] += ( HGrad_i_xHT + HGrad_i_xHT.transpose() ) * phi + HxHT * phiGradient( i );
      }

      phi_cache[k]         = phi;
      phiGradient_cache[k] = std::move( phiGradient );
      H_cache[k]           = std::move( H );
      HGradient_cache[k]   = std::move( HGradient );
    }

    // solve for b(x)
    // b = M^-1 * H0
    const auto H0 = H0Vector( sizeH );

    const auto MHr = M.colPivHouseholderQr();

    const Eigen::VectorXd b = MHr.solve( H0 );

    // set all shape function values to zero
    Eigen::Map< Eigen::VectorXd > shapeFunctionValues_( shapeFunctionValues, kernelFunctionCandidates.size() );
    shapeFunctionValues_.setZero();

    shapeFunctionValueGradients.setZero();

    // The gradient of b(x) is independent of the individual kernel function A, so it is
    // computed once here instead of redundantly inside the loop over covering kernels.
    //
    // dPsiA_dxi = H0_J * ( InvM_JK * H_K * phi_A ),xi
    //
    // dPsiA_dxi = H0_J * ( InvM_JK,xi * H_K    * phi_A +
    //                      InvM_JK    * H_K,xi * phi_A +
    //                      InvM_JK    * H_K    * phi_A,xi )
    //
    // dPsiA_dxi = H0_J * ( - [InvM_JA * M_AB,xi * InvM_BK]    * H_K    * phi_A +
    //                      InvM_JK                            * H_K,xi * phi_A +
    //                      InvM_JK                            * H_K    * phi_A,xi )
    Eigen::MatrixXd bGradient = Eigen::MatrixXd::Zero( sizeH, _dim );
    for ( int i = 0; i < _dim; i++ ) {
      // b_{,i} = -M^-1 ( M_{,i} b ). Solve once with a vector right-hand side rather than solving
      // for the full sizeH x sizeH matrix M^-1 M_{,i} and then multiplying by b.
      bGradient.col( i ) = -MHr.solve( MGradients[i] * b );
    }
    // Transpose view (no copy); bGradient outlives the loop below.
    const auto bGradientTransposed = bGradient.transpose();

    for ( size_t k = 0; k < nCovering; k++ ) {
      const int              A             = coveringKernelFunctionIndices[k];
      const double           phi_A         = phi_cache[k];
      const Eigen::VectorXd& phiGradient_A = phiGradient_cache[k];
      const Eigen::VectorXd& H             = H_cache[k];
      const Eigen::MatrixXd& HGradient     = HGradient_cache[k];

      shapeFunctionValues[A] = b.dot( H ) * phi_A;

      shapeFunctionValueGradients.col( A ) += bGradientTransposed * H * phi_A;
      shapeFunctionValueGradients.col( A ) += b.dot( H ) * phiGradient_A;
      // Column-vector form of ( b^T H_{,} ), keeping the result a dim x 1 column.
      shapeFunctionValueGradients.col( A ) += ( HGradient.transpose() * b ) * phi_A;
    }
  }

}; // namespace Marmot::Meshfree
