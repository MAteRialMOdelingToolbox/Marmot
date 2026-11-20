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
#include "Marmot/MarmotMeshfreeReproducingKernelApproximationImplicit.h"
#include "Marmot/MarmotMonomialBasisFunctions.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/src/Core/Matrix.h>
#include <cmath>

namespace Marmot::Discretizations::Meshfree {
  ReproducingKernelApproximationImplicit::ReproducingKernelApproximationImplicit( int dim, int completenessOrder )
    : ReproducingKernelApproximation( dim, completenessOrder )
  {
  }

  void ReproducingKernelApproximationImplicit::computeShapeFunctionsAndGradients(
    const double*                               coord,
    const std::vector< const KernelFunction* >& kernelFunctions,
    double*                                     shapeFunctionValues_,
    double*                                     shapeFunctionValueGradients_ ) const
  {
    // the kernel values are kept in the output array until the shape functions overwrite them; the implicit gradients
    // need no kernel gradients
    double*   phi       = shapeFunctionValues_;
    const int nCovering = evaluateKernelFunctions( coord, kernelFunctions, phi, nullptr );

    const int completenessOrder = getCorrectedCompletenessOrder( nCovering );
    const int sizeH             = _basisOfOrder[completenessOrder].size();

    MomentMatrix M;
    assembleMomentMatrix( coord, kernelFunctions, phi, nullptr, completenessOrder, M, nullptr );

    // first column for the values, the other columns for the gradients: they reproduce -d(x - x_I)^alpha/dx_i at
    // x_I, which is -1 for the linear monomial x_i and 0 otherwise; the linear monomials are the rows whose gradient
    // at the origin is nonzero
    using RightHandSides = Eigen::
      Matrix< double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, maxSizeOfBasis, 1 + maxDim >;
    RightHandSides H0 = RightHandSides::Zero( sizeH, 1 + _dim );
    H0( 0, 0 )        = 1.0;
    if ( completenessOrder > 0 ) {
      const Coordinates origin = Coordinates::Zero( _dim );
      BasisGradient     dHOrigin( sizeH, _dim );
      _basisOfOrder[completenessOrder].evaluateGradient( origin, dHOrigin );
      for ( int k = 0; k < sizeH; k++ )
        for ( int i = 0; i < _dim; i++ )
          if ( dHOrigin( k, i ) != 0.0 )
            H0( k, 1 + i ) = -1.0;
    }

    const RightHandSides bMat = factorizeMomentMatrix( M, coord, nCovering ).solve( H0 );

    Eigen::Map< Eigen::MatrixXd > shapeFunctionValueGradients( shapeFunctionValueGradients_,
                                                               _dim,
                                                               kernelFunctions.size() );
    shapeFunctionValueGradients.setZero();

    BasisVector H;
    for ( size_t A = 0; A < kernelFunctions.size(); A++ ) {
      const double phi_A = phi[A];
      if ( phi_A == 0.0 )
        continue; // not covering: zero value, already written, and zero gradient
      evaluateBasis( coord, *kernelFunctions[A], completenessOrder, H, nullptr );

      shapeFunctionValues_[A]              = bMat.col( 0 ).dot( H ) * phi_A;
      shapeFunctionValueGradients.col( A ) = bMat.rightCols( _dim ).transpose() * H * phi_A;
    }
  }

}; // namespace Marmot::Discretizations::Meshfree
