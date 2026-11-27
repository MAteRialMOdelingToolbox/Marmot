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
  ReproducingKernelApproximationImplicit::ReproducingKernelApproximationImplicit( int    dim,
                                                                                  int    completenessOrder,
                                                                                  double momentMatrixConditionLimit )
    : ReproducingKernelApproximation( dim, completenessOrder, momentMatrixConditionLimit )
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

    Eigen::MatrixXd H;
    MomentMatrix    M;
    const int       completenessOrder = chooseCompletenessOrderAndAssembleMomentMatrix( coord,
                                                                                  kernelFunctions,
                                                                                  phi,
                                                                                  nullptr,
                                                                                  nCovering,
                                                                                  H,
                                                                                  nullptr,
                                                                                  M,
                                                                                  nullptr );
    const int       sizeH             = _basisOfOrder[completenessOrder].size();

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

    if ( completenessOrder == 0 ) {
      // The implicit gradients of order 0 reproduce only the derivative of a constant: they vanish identically, and a
      // point evaluated with them would carry no strain. The point gets the Shepard gradients instead, the derivatives
      // of its values, as with explicit gradients: with b = 1 / sum_B phi_B, Psi_A,i = b phi_A,i + b_,i phi_A and
      // b_,i = -b^2 sum_B phi_B,i. Only here are the kernel gradients needed.
      const double b                    = bMat( 0, 0 );
      Coordinates  sumOfKernelGradients = Coordinates::Zero( _dim );
      for ( size_t A = 0; A < kernelFunctions.size(); A++ )
        if ( phi[A] != 0.0 ) {
          kernelFunctions[A]->computeKernelFunctionAndGradient( coord, shapeFunctionValueGradients.col( A ).data() );
          sumOfKernelGradients += shapeFunctionValueGradients.col( A );
        }
      const Coordinates db_dx = -b * b * sumOfKernelGradients;
      for ( size_t A = 0; A < kernelFunctions.size(); A++ )
        if ( phi[A] != 0.0 )
          shapeFunctionValueGradients.col( A ) = b * shapeFunctionValueGradients.col( A ) + phi[A] * db_dx;
    }

    for ( size_t A = 0; A < kernelFunctions.size(); A++ ) {
      const double phi_A = phi[A];
      if ( phi_A == 0.0 )
        continue; // not covering: zero value, already written, and zero gradient
      const auto H_A = H.col( A );

      shapeFunctionValues_[A] = bMat.col( 0 ).dot( H_A ) * phi_A;
      if ( completenessOrder > 0 )
        shapeFunctionValueGradients.col( A ) = bMat.rightCols( _dim ).transpose() * H_A * phi_A;
    }
  }

}; // namespace Marmot::Discretizations::Meshfree
