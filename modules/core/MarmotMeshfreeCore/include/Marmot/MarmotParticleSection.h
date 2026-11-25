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
#pragma once

#include "Marmot/MarmotMeshfreeApproximation.h"
#include "Marmot/MarmotMeshfreeKernelFunction.h"
#include "Marmot/MarmotMeshfreeQuadHexCell.h"
#include <Eigen/Core>
#include <array>
#include <vector>

/**
 * @file
 * @brief The integrals over the domain of a particle: the smoothed gradients of the shape functions, one evaluation
 * for the SQCNI and NSNI particles of both families.
 */

namespace Marmot::Discretizations::Meshfree {

  /**
   * @brief The integrals over the domain of a particle with a domain (a Lagrange cell): its volume and the smoothed
   * gradients of the shape functions over its smoothing domain.
   *
   * @details These are the integrals of the plane-strain and 3D particles: the smoothed gradient with one point per
   * face (the face center),
   * @f[
   *   \tilde G_{Bj} = \frac{1}{V_s} \sum_f N_B(\boldsymbol{Y}_f)\, n_{f,j}\, dA_f .
   * @f]
   * The second gradients of NSNI are smoothed alike, from the gradients at the face centers,
   * @f[
   *   \frac{\partial^2 N_B}{\partial Y_I \partial Y_J} \approx \frac{1}{V_s} \sum_f \frac{\partial
   *   N_B}{\partial Y_I}(\boldsymbol{Y}_f)\, n_{f,J}\, dA_f
   * @f]
   * (stored unsymmetrized, as `d2N_dYdY[J](I, B)`).
   *
   * @tparam nDim      The spatial dimension.
   * @tparam nVertices The number of vertices of the cell (4: Quad, 8: Hexa).
   */
  template < int nDim, int nVertices >
  class ParticleSection {

  public:
    using Cell = LagrangeCell< nDim, nVertices >;  ///< the domain of a particle
    using Vec  = Eigen::Matrix< double, nDim, 1 >; ///< a point or a vector

    /**
     * @brief The volume of a domain.
     * @param[in] cell The domain.
     * @return Its volume (area in plane strain).
     */
    static double volume( const Cell& cell ) { return cell.volume(); }

    /**
     * @brief The quadrature points of the faces of a smoothing domain: the face centers.
     * @param[in] cell The smoothing domain.
     * @return nDim x number of faces, face by face.
     */
    static Eigen::Matrix< double, nDim, Eigen::Dynamic > faceQuadraturePoints( const Cell& cell )
    {
      Eigen::Matrix< double, nDim, Eigen::Dynamic > points( nDim, cell.getNumberOfFaces() );
      for ( int f = 0; f < cell.getNumberOfFaces(); f++ )
        points.col( f ) = faceQuadraturePoint( cell, f + 1 );
      return points;
    }

    /**
     * @brief The smoothed gradients (and second gradients) of the shape functions over a smoothing domain, see the
     * class description.
     * @param[in] smoothingDomain The smoothing domain.
     * @param[in] approximation The meshfree approximation.
     * @param[in] kernelFunctions The kernel functions of the nodes.
     * @param[in] N The shape functions at the particle center, 1 x nNodes.
     * @param[out] dN_dY The smoothed gradients, nDim x nNodes.
     * @param[out] d2N_dYdY The smoothed second gradients, `d2N_dYdY[J](I, B)`; nullptr if not needed.
     */
    static void smoothShapeFunctionGradients( const Cell&                                 smoothingDomain,
                                              const Approximation&                        approximation,
                                              const std::vector< const KernelFunction* >& kernelFunctions,
                                              const Eigen::MatrixXd&                      N,
                                              Eigen::MatrixXd&                            dN_dY,
                                              std::array< Eigen::MatrixXd, nDim >*        d2N_dYdY );

  private:
    /**
     * @brief The quadrature point of a face.
     * @param[in] cell The smoothing domain.
     * @param[in] faceID The face (1-based).
     * @return The face center.
     */
    static Vec faceQuadraturePoint( const Cell& cell, int faceID ) { return cell.getFaceCenterCoordinates( faceID ); }
  };

  template < int nDim, int nVertices >
  void ParticleSection< nDim, nVertices >::smoothShapeFunctionGradients(
    const Cell&                                 smoothingDomain,
    const Approximation&                        approximation,
    const std::vector< const KernelFunction* >& kernelFunctions,
    const Eigen::MatrixXd&                      N,
    Eigen::MatrixXd&                            dN_dY,
    std::array< Eigen::MatrixXd, nDim >*        d2N_dYdY )
  {
    const int nNodes = static_cast< int >( N.cols() );
    dN_dY            = Eigen::MatrixXd::Zero( nDim, nNodes );
    if ( d2N_dYdY )
      for ( int J = 0; J < nDim; J++ )
        ( *d2N_dYdY )[J] = Eigen::MatrixXd::Zero( nDim, nNodes );

    Eigen::MatrixXd NBoundary( 1, nNodes );
    Eigen::MatrixXd dN_dY_Boundary = Eigen::MatrixXd::Zero( nDim, nNodes );

    for ( int f = 0; f < smoothingDomain.getNumberOfFaces(); f++ ) {
      const Vec Y = faceQuadraturePoint( smoothingDomain, f + 1 );

      if ( d2N_dYdY )
        approximation.computeShapeFunctionsAndGradients( Y.data(),
                                                         kernelFunctions,
                                                         NBoundary.data(),
                                                         dN_dY_Boundary.data() );
      else
        approximation.computeShapeFunctions( Y.data(), kernelFunctions, NBoundary.data() );

      const Vec n_dA = smoothingDomain.boundarySurfaceVector( f + 1 );
      dN_dY += n_dA * NBoundary;
      if ( d2N_dYdY )
        for ( int J = 0; J < nDim; J++ )
          ( *d2N_dYdY )[J] += n_dA( J ) * dN_dY_Boundary;
    }

    const double V = volume( smoothingDomain );
    dN_dY /= V;
    if ( d2N_dYdY )
      for ( int J = 0; J < nDim; J++ )
        ( *d2N_dYdY )[J] /= V;
  }

} // namespace Marmot::Discretizations::Meshfree
