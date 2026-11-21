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
#include <stdexcept>
#include <vector>

namespace Marmot::Discretizations::Meshfree {

  ReproducingKernelApproximation::ReproducingKernelApproximation( int dim, int completenessOrder )
    : _dim( dim ), _desiredCompletenessOrder( completenessOrder )
  {
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
      phi[A]              = kernelFunctions[A]->computeKernelFunction( coord );
      const bool covering = std::abs( phi[A] ) > 1e-14;
      if ( !covering )
        phi[A] = 0.0;
      else
        nCovering++;

      if ( dPhi_dx ) {
        double* dPhi_A_dx = dPhi_dx + _dim * A;
        if ( covering )
          kernelFunctions[A]->computeKernelFunctionGradient( coord, dPhi_A_dx );
        else
          for ( int i = 0; i < _dim; i++ )
            dPhi_A_dx[i] = 0.0;
      }
    }
    return nCovering;
  }

  void ReproducingKernelApproximation::evaluateBasis( const double*         coord,
                                                      const KernelFunction& kernelFunction,
                                                      int                   completenessOrder,
                                                      BasisVector&          H,
                                                      BasisGradient*        dH_dx ) const
  {
    const Math::MonomialBasis& basis = _basisOfOrder[completenessOrder];
    Coordinates                x_minus_xA( _dim );
    for ( int i = 0; i < _dim; i++ )
      x_minus_xA( i ) = coord[i] - kernelFunction.getCenterCoordinates()[i];

    H.resize( basis.size() );
    basis.evaluate( x_minus_xA, H );
    if ( dH_dx ) {
      dH_dx->resize( basis.size(), _dim );
      basis.evaluateGradient( x_minus_xA, *dH_dx );
    }
  }

  void ReproducingKernelApproximation::assembleMomentMatrix(
    const double*                               coord,
    const std::vector< const KernelFunction* >& kernelFunctions,
    const double*                               phi,
    const double*                               dPhi_dx,
    int                                         completenessOrder,
    MomentMatrix&                               M,
    MomentMatrixGradient*                       dM_dx ) const
  {
    const int sizeH = _basisOfOrder[completenessOrder].size();
    M.setZero( sizeH, sizeH );
    if ( dM_dx )
      for ( int i = 0; i < _dim; i++ )
        ( *dM_dx )[i].setZero( sizeH, sizeH );

    BasisVector   H;
    BasisGradient dH_dx;
    for ( size_t B = 0; B < kernelFunctions.size(); B++ ) {
      if ( phi[B] == 0.0 )
        continue;
      evaluateBasis( coord, *kernelFunctions[B], completenessOrder, H, dM_dx ? &dH_dx : nullptr );

      // M += H H^T phi_B, and M_,i += ( H_,i H^T + H H_,i^T ) phi_B + H H^T phi_B,i, entry by entry: for the small
      // matrices of a basis, plain loops are faster than products of matrices of a runtime size
      const double phi_B = phi[B];
      for ( int b = 0; b < sizeH; b++ )
        for ( int a = 0; a < sizeH; a++ )
          M( a, b ) += H( a ) * H( b ) * phi_B;
      if ( dM_dx )
        for ( int i = 0; i < _dim; i++ ) {
          const double  dPhi_B_dxi = dPhi_dx[_dim * B + i];
          MomentMatrix& dM_dxi     = ( *dM_dx )[i];
          for ( int b = 0; b < sizeH; b++ )
            for ( int a = 0; a < sizeH; a++ )
              dM_dxi( a, b ) += ( dH_dx( a, i ) * H( b ) + H( a ) * dH_dx( b, i ) ) * phi_B +
                                H( a ) * H( b ) * dPhi_B_dxi;
        }
    }
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

    const int completenessOrder = getCorrectedCompletenessOrder( nCovering );
    const int sizeH             = _basisOfOrder[completenessOrder].size();

    MomentMatrix         M;
    MomentMatrixGradient dM_dx;
    assembleMomentMatrix( coord, kernelFunctions, phi, dPhi_dx, completenessOrder, M, dPsi_dx ? &dM_dx : nullptr );

    // b = M^-1 H0, with H0 = H(0) = [1, 0, ..., 0]^T
    const auto        MQr = factorizeMomentMatrix( M, coord, nCovering );
    const BasisVector b   = MQr.solve( BasisVector::Unit( sizeH, 0 ) );

    // b_,i = -M^-1 ( M_,i b ): the same for every kernel, one solve per dimension
    BasisGradient db_dx( sizeH, _dim );
    if ( dPsi_dx )
      for ( int i = 0; i < _dim; i++ )
        db_dx.col( i ) = -MQr.solve( dM_dx[i] * b );

    // Psi_A = b^T H phi_A, and Psi_A,i = ( b_,i^T H + b^T H_,i ) phi_A + b^T H phi_A,i
    BasisVector   H;
    BasisGradient dH_dx;
    for ( size_t A = 0; A < kernelFunctions.size(); A++ ) {
      const double phi_A = phi[A];
      if ( phi_A == 0.0 )
        continue; // not covering: zero value and gradient, already written
      evaluateBasis( coord, *kernelFunctions[A], completenessOrder, H, dPsi_dx ? &dH_dx : nullptr );

      const double bH = b.dot( H );
      Psi[A]          = bH * phi_A;
      if ( dPsi_dx )
        for ( int i = 0; i < _dim; i++ ) {
          const double dPhi_A_dxi = dPhi_dx[_dim * A + i];
          dPsi_dx[_dim * A + i]   = ( db_dx.col( i ).dot( H ) + b.dot( dH_dx.col( i ) ) ) * phi_A + bH * dPhi_A_dxi;
        }
    }
  }

} // namespace Marmot::Discretizations::Meshfree
