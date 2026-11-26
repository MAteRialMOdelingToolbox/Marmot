/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
 *
 * This file is part of the MAteRialMOdellingToolbox (marmot).
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of marmot.
 * ---------------------------------------------------------------------
 */
#include "Marmot/MarmotMeshfreeReproducingKernelApproximation.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMeshfreeKernelFunction.h"
#include "Marmot/MarmotMonomialBasisFunctions.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace Marmot::Discretizations::Meshfree {

  ReproducingKernelApproximation::ReproducingKernelApproximation( int    dim,
                                                                  int    completenessOrder,
                                                                  double momentMatrixConditionLimit )
    : _dim( dim ),
      _desiredCompletenessOrder( completenessOrder ),
      _momentMatrixConditionLimit( momentMatrixConditionLimit )
  {
    if ( !( momentMatrixConditionLimit >= 1.0 ) )
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__ << ": the limit on the condition of the moment matrix is "
                                   << momentMatrixConditionLimit
                                   << ", but a condition number is at least 1 (infinity for no limit)" );
    if ( dim < 1 || dim > maxDim || completenessOrder < 0 ||
         Math::computeSizeOfMonomialBasisVector( completenessOrder, dim ) > maxSizeOfBasis )
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__ << ": dimension " << dim << " and completeness order "
                                   << completenessOrder << " are not supported: the dimension is at most " << maxDim
                                   << ", the order at least 0, and the monomial basis at most " << maxSizeOfBasis
                                   << " entries (cubic in 3D)" );
    for ( int order = 0; order <= completenessOrder; order++ )
      _basisOfOrder.emplace_back( order, dim );
  }

  int ReproducingKernelApproximation::evaluateKernelFunctions(
    const double*                               coord,
    const std::vector< const KernelFunction* >& kernelFunctions,
    double*                                     phi,
    double*                                     dPhi_dx ) const
  {
    int nCovering = 0;
    for ( size_t A = 0; A < kernelFunctions.size(); A++ ) {
      double* dPhi_A_dx = dPhi_dx ? dPhi_dx + _dim * A : nullptr;
      phi[A]            = dPhi_A_dx ? kernelFunctions[A]->computeKernelFunctionAndGradient( coord, dPhi_A_dx )
                                    : kernelFunctions[A]->computeKernelFunction( coord );

      if ( std::abs( phi[A] ) > 1e-14 )
        nCovering++;
      else {
        phi[A] = 0.0;
        if ( dPhi_A_dx )
          for ( int i = 0; i < _dim; i++ )
            dPhi_A_dx[i] = 0.0;
      }
    }
    return nCovering;
  }

  void ReproducingKernelApproximation::evaluateBasisOfKernels(
    const double*                               coord,
    const std::vector< const KernelFunction* >& kernelFunctions,
    const double*                               phi,
    int                                         completenessOrder,
    Eigen::MatrixXd&                            H,
    Eigen::MatrixXd*                            dH_dx ) const
  {
    const Math::MonomialBasis& basis = _basisOfOrder[completenessOrder];
    const Eigen::Index         nC    = static_cast< Eigen::Index >( kernelFunctions.size() );
    H.resize( basis.size(), nC );
    if ( dH_dx )
      dH_dx->resize( basis.size(), _dim * nC );

    Coordinates x_minus_xA( _dim );
    for ( Eigen::Index A = 0; A < nC; A++ ) {
      if ( phi[A] == 0.0 )
        continue;
      const double* xA = kernelFunctions[A]->getCenterCoordinates();
      for ( int i = 0; i < _dim; i++ )
        x_minus_xA( i ) = coord[i] - xA[i];

      basis.evaluate( x_minus_xA, H.col( A ) );
      if ( dH_dx )
        basis.evaluateGradient( x_minus_xA, dH_dx->middleCols( _dim * A, _dim ) );
    }
  }

  void ReproducingKernelApproximation::assembleMomentMatrix( const Eigen::MatrixXd& H,
                                                             const Eigen::MatrixXd* dH_dx,
                                                             const double*          phi,
                                                             const double*          dPhi_dx,
                                                             MomentMatrix&          M,
                                                             MomentMatrixGradient*  dM_dx ) const
  {
    // the sums are formed in local matrices: the compiler then knows that adding to them does not change the bases,
    // and keeps the entries of a basis in registers
    const Eigen::Index   sizeH = H.rows();
    MomentMatrix         sumM  = MomentMatrix::Zero( sizeH, sizeH );
    MomentMatrixGradient sumdM_dx;
    if ( dM_dx )
      for ( int i = 0; i < _dim; i++ )
        sumdM_dx[i].setZero( sizeH, sizeH );

    for ( Eigen::Index B = 0; B < H.cols(); B++ ) {
      if ( phi[B] == 0.0 )
        continue;

      // M += H H^T phi_B, and M_,i += ( H_,i H^T + H H_,i^T ) phi_B + H H^T phi_B,i, entry by entry: for the small
      // matrices of a basis, plain loops are faster than products of matrices of a runtime size
      const double phi_B = phi[B];
      const auto   H_B   = H.col( B );
      for ( Eigen::Index b = 0; b < sizeH; b++ )
        for ( Eigen::Index a = 0; a < sizeH; a++ )
          sumM( a, b ) += H_B( a ) * H_B( b ) * phi_B;
      if ( dM_dx )
        for ( int i = 0; i < _dim; i++ ) {
          const double  dPhi_B_dxi = dPhi_dx[_dim * B + i];
          const auto    dH_B_dxi   = dH_dx->col( _dim * B + i );
          MomentMatrix& sumdM_dxi  = sumdM_dx[i];
          for ( Eigen::Index b = 0; b < sizeH; b++ )
            for ( Eigen::Index a = 0; a < sizeH; a++ )
              sumdM_dxi( a, b ) += ( dH_B_dxi( a ) * H_B( b ) + H_B( a ) * dH_B_dxi( b ) ) * phi_B +
                                   H_B( a ) * H_B( b ) * dPhi_B_dxi;
        }
    }

    M = sumM;
    if ( dM_dx )
      *dM_dx = sumdM_dx;
  }

  Eigen::ColPivHouseholderQR< ReproducingKernelApproximation::MomentMatrix > ReproducingKernelApproximation::
    factorizeMomentMatrix( const MomentMatrix& M, const double* coord, int nCoveringKernels ) const
  {
    Eigen::ColPivHouseholderQR< MomentMatrix > MQr( M );
    if ( MQr.rank() < M.rows() )
      throw std::
        runtime_error( Marmot::Utilities::MakeString()
                       << __PRETTY_FUNCTION__ << ": singular moment matrix (rank " << MQr.rank() << " < " << M.rows()
                       << ") at " << Eigen::Map< const Eigen::VectorXd >( coord, _dim ).transpose() << ", covered by "
                       << nCoveringKernels
                       << " kernels: no kernel covers the point, or the kernel centres are degenerate for the "
                          "completeness order (e.g. all coplanar in 3D, as for a single layer of particles through the "
                          "thickness). Use at least two layers of kernels in every direction, larger supports, or a "
                          "lower completeness order" );
    return MQr;
  }

  int ReproducingKernelApproximation::chooseCompletenessOrderAndAssembleMomentMatrix(
    const double*                               coord,
    const std::vector< const KernelFunction* >& kernelFunctions,
    const double*                               phi,
    const double*                               dPhi_dx,
    int                                         nCovering,
    Eigen::MatrixXd&                            H,
    Eigen::MatrixXd*                            dH_dx,
    MomentMatrix&                               M,
    MomentMatrixGradient*                       dM_dx ) const
  {
    int completenessOrder = getCorrectedCompletenessOrder( nCovering );
    evaluateBasisOfKernels( coord, kernelFunctions, phi, completenessOrder, H, dH_dx );
    assembleMomentMatrix( H, dH_dx, phi, dPhi_dx, M, dM_dx );
    if ( std::isinf( _momentMatrixConditionLimit ) )
      return completenessOrder;

    // nearly singular: one order lower, down to the Shepard functions of the covering kernels (order 0)
    while ( completenessOrder > 0 && conditionOfMomentMatrix( M ) > _momentMatrixConditionLimit ) {
      completenessOrder--;
      evaluateBasisOfKernels( coord, kernelFunctions, phi, completenessOrder, H, dH_dx );
      assembleMomentMatrix( H, dH_dx, phi, dPhi_dx, M, dM_dx );
    }
    return completenessOrder;
  }

  double ReproducingKernelApproximation::conditionOfMomentMatrix( const MomentMatrix& M )
  {
    const Eigen::Index n = M.rows();
    BasisVector        D( n );
    for ( Eigen::Index k = 0; k < n; k++ ) {
      if ( !( M( k, k ) > 0.0 ) )
        return std::numeric_limits< double >::infinity();
      D( k ) = 1.0 / std::sqrt( M( k, k ) );
    }
    const MomentMatrix scaledM = D.asDiagonal() * M * D.asDiagonal();

    // The scaled matrix has a unit diagonal, so its pivots are of order 1 when it is healthy and never below its
    // smallest eigenvalue. A pivot at the rounding level marks it singular: Eigen's LDLT then still succeeds and
    // solves with the pseudo-inverse of D, so its condition estimate would be finite.
    const Eigen::LDLT< MomentMatrix > scaledMLdlt( scaledM );
    if ( scaledMLdlt.info() != Eigen::Success ||
         !( scaledMLdlt.vectorD().minCoeff() > n * std::numeric_limits< double >::epsilon() ) )
      return std::numeric_limits< double >::infinity();
    const double reciprocalCondition = scaledMLdlt.rcond();
    return reciprocalCondition > 0.0 ? 1.0 / reciprocalCondition : std::numeric_limits< double >::infinity();
  }

  void ReproducingKernelApproximation::computeShapeFunctions(
    const double*                               coord,
    const std::vector< const KernelFunction* >& kernelFunctions,
    double*                                     shapeFunctionValues ) const
  {
    evaluateShapeFunctions( coord, kernelFunctions, shapeFunctionValues, nullptr );
  }

  void ReproducingKernelApproximation::computeShapeFunctionsAndGradients(
    const double*                               coord,
    const std::vector< const KernelFunction* >& kernelFunctions,
    double*                                     shapeFunctionValues,
    double*                                     shapeFunctionValueGradients_ColMajor ) const
  {
    evaluateShapeFunctions( coord, kernelFunctions, shapeFunctionValues, shapeFunctionValueGradients_ColMajor );
  }

  void ReproducingKernelApproximation::evaluateShapeFunctions(
    const double*                               coord,
    const std::vector< const KernelFunction* >& kernelFunctions,
    double*                                     Psi,
    double*                                     dPsi_dx ) const
  {
    // the kernel values and gradients are kept in the output arrays until the shape functions overwrite them
    double*   phi       = Psi;
    double*   dPhi_dx   = dPsi_dx;
    const int nCovering = evaluateKernelFunctions( coord, kernelFunctions, phi, dPhi_dx );

    Eigen::MatrixXd      H, dH_dx;
    MomentMatrix         M;
    MomentMatrixGradient dM_dx;
    const int            completenessOrder = chooseCompletenessOrderAndAssembleMomentMatrix( coord,
                                                                                  kernelFunctions,
                                                                                  phi,
                                                                                  dPhi_dx,
                                                                                  nCovering,
                                                                                  H,
                                                                                  dPsi_dx ? &dH_dx : nullptr,
                                                                                  M,
                                                                                  dPsi_dx ? &dM_dx : nullptr );
    const int            sizeH             = _basisOfOrder[completenessOrder].size();

    // b = M^-1 H0, with H0 = H(0) = [1, 0, ..., 0]^T
    const auto        MQr = factorizeMomentMatrix( M, coord, nCovering );
    const BasisVector b   = MQr.solve( BasisVector::Unit( sizeH, 0 ) );

    // b_,i = -M^-1 ( M_,i b ): the same for every kernel, one solve per dimension
    BasisGradient db_dx( sizeH, _dim );
    if ( dPsi_dx )
      for ( int i = 0; i < _dim; i++ )
        db_dx.col( i ) = -MQr.solve( dM_dx[i] * b );

    // Psi_A = b^T H phi_A, and Psi_A,i = ( b_,i^T H + b^T H_,i ) phi_A + b^T H phi_A,i
    for ( size_t A = 0; A < kernelFunctions.size(); A++ ) {
      const double phi_A = phi[A];
      if ( phi_A == 0.0 )
        continue; // not covering: zero value and gradient, already written
      const auto H_A = H.col( A );

      const double bH = b.dot( H_A );
      Psi[A]          = bH * phi_A;
      if ( dPsi_dx )
        for ( int i = 0; i < _dim; i++ ) {
          const double dPhi_A_dxi = dPhi_dx[_dim * A + i];
          dPsi_dx[_dim * A + i]   = ( db_dx.col( i ).dot( H_A ) + b.dot( dH_dx.col( _dim * A + i ) ) ) * phi_A +
                                  bH * dPhi_A_dxi;
        }
    }
  }

} // namespace Marmot::Discretizations::Meshfree
