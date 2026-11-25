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
 * @brief The integrals over the smoothing domain of a particle: its volume and the smoothed gradients of the shape
 * functions, one implementation for the SQCNI, SNNI, NSNI and SDI particles of both families.
 */

namespace Marmot::Discretizations::Meshfree {

  /**
   * @brief The integrals over the smoothing domain of a particle: its volume and the smoothed gradients of the shape
   * functions (stabilized conforming nodal integration, Chen et al., 2001).
   *
   * @details The smoothed gradient is the average of the gradient over the smoothing domain @f$ \Omega_s @f$, turned
   * into a boundary integral by the divergence theorem and integrated with one point per face @f$ f @f$, the face
   * center @f$ \boldsymbol{Y}_f @f$ (faceQuadraturePoints()),
   * @f[
   *   \tilde N_{A,i} = \frac{1}{|\Omega_s|}\int_{\Omega_s} \frac{\partial N_A}{\partial Y_i}\,dV
   *   = \frac{1}{|\Omega_s|}\oint_{\partial\Omega_s} N_A\,n_i\,dA
   *   \approx \frac{1}{|\Omega_s|} \sum_f N_A(\boldsymbol{Y}_f)\, n_{f,i}\, dA_f ,
   * @f]
   * with the coordinates @f$ \boldsymbol{Y} @f$ of the intermediate configuration, the volume @f$ |\Omega_s| @f$
   * (volume()) and the area vector @f$ \boldsymbol{n}_f\,dA_f @f$ of face @f$ f @f$
   * (LagrangeCell::boundarySurfaceVector()). The face rule is the only approximation, and only the values of the
   * shape functions at the face points are needed. NSNI smooths the second derivatives alike, from the (unsmoothed)
   * gradients at the same points,
   * @f[
   *   \hat N_{A,IJ} = \frac{1}{|\Omega_s|} \sum_f \frac{\partial N_A}{\partial Y_I}(\boldsymbol{Y}_f)\,
   *   n_{f,J}\, dA_f ,
   * @f]
   * stored unsymmetrized as `d2N_dYdY[J](I, A)`; NSNI::secondGradientOfNode() symmetrizes them.
   *
   * The particles take these integrals, and their evaluation points (faceQuadraturePoints()), from here and not from
   * the cell: the face rule of the smoothing is defined in this class only. *Section* is meant as for finite elements
   * (plane strain, 3D): the form of the integrals depends on it. Here the section is Cartesian, and the integrals carry
   * the weight 1. The class holds static functions only, so that a particle selects its integrals with its own template
   * parameters.
   *
   * @tparam nDim      The spatial dimension.
   * @tparam nVertices The number of vertices of the cell (4: Quad, 8: Hexa).
   */
  template < int nDim, int nVertices >
  class ParticleSectionIntegrals {

  public:
    using Cell = LagrangeCell< nDim, nVertices >;  ///< a smoothing domain (a Lagrange cell)
    using Vec  = Eigen::Matrix< double, nDim, 1 >; ///< a point or a vector

    /// @brief The number of quadrature points per face of a smoothing domain: one, the face center.
    static constexpr int nPointsPerFace = 1;

    /// @brief The quadrature points of one face, nDim x nPointsPerFace.
    using FacePoints = Eigen::Matrix< double, nDim, nPointsPerFace >;

    /**
     * @brief The volume of a smoothing domain.
     * @param[in] cell The smoothing domain.
     * @return Its volume @f$ |\Omega_s| @f$ (area in plane strain).
     */
    static double volume( const Cell& cell ) { return cell.volume(); }

    /**
     * @brief The number of face quadrature points of a smoothing domain.
     * @param[in] cell The smoothing domain.
     * @return The number of face quadrature points, `faceQuadraturePoints( cell ).cols()`.
     */
    static int numberOfFaceQuadraturePoints( const Cell& cell ) { return cell.getNumberOfFaces() * nPointsPerFace; }

    /**
     * @brief The quadrature points of one face of a smoothing domain: its center.
     * @param[in] cell The smoothing domain.
     * @param[in] faceID The face (1-based).
     * @return nDim x nPointsPerFace.
     */
    static FacePoints faceQuadraturePoints( const Cell& cell, int faceID )
    {
      return cell.getFaceCenterCoordinates( faceID );
    }

    /**
     * @brief The quadrature points of all faces of a smoothing domain.
     * @param[in] cell The smoothing domain.
     * @return nDim x numberOfFaceQuadraturePoints(), face by face.
     */
    static Eigen::Matrix< double, nDim, Eigen::Dynamic > faceQuadraturePoints( const Cell& cell )
    {
      Eigen::Matrix< double, nDim, Eigen::Dynamic > points( nDim, numberOfFaceQuadraturePoints( cell ) );
      for ( int f = 0; f < cell.getNumberOfFaces(); f++ )
        points.middleCols( f * nPointsPerFace, nPointsPerFace ) = faceQuadraturePoints( cell, f + 1 );
      return points;
    }

    /**
     * @brief The smoothed gradients (and second derivatives) of the shape functions over a smoothing domain, see the
     * class description. The points of a face share its area vector equally.
     * @param[in] smoothingDomain The smoothing domain.
     * @param[in] approximation The meshfree approximation.
     * @param[in] kernelFunctions The kernel functions of the nodes.
     * @param[in] N The shape functions at the particle center, 1 x nNodes; the Cartesian smoothing uses only its size
     * (the number of nodes).
     * @param[out] dN_dY The smoothed gradients, nDim x nNodes.
     * @param[out] d2N_dYdY The smoothed second derivatives, `d2N_dYdY[J](I, A)`; nullptr: first derivatives only
     * (SQCNI, SNNI, SDI); NSNI passes its array.
     */
    static void smoothShapeFunctionGradients( const Cell&                                 smoothingDomain,
                                              const Approximation&                        approximation,
                                              const std::vector< const KernelFunction* >& kernelFunctions,
                                              const Eigen::MatrixXd&                      N,
                                              Eigen::MatrixXd&                            dN_dY,
                                              std::array< Eigen::MatrixXd, nDim >*        d2N_dYdY );
  };

  template < int nDim, int nVertices >
  void ParticleSectionIntegrals< nDim, nVertices >::smoothShapeFunctionGradients(
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
      const FacePoints points = faceQuadraturePoints( smoothingDomain, f + 1 );
      const Vec        n_dA   = smoothingDomain.boundarySurfaceVector( f + 1 ) / double( nPointsPerFace );

      for ( int g = 0; g < nPointsPerFace; g++ ) {
        const Vec Y = points.col( g );

        if ( d2N_dYdY )
          approximation.computeShapeFunctionsAndGradients( Y.data(),
                                                           kernelFunctions,
                                                           NBoundary.data(),
                                                           dN_dY_Boundary.data() );
        else
          approximation.computeShapeFunctions( Y.data(), kernelFunctions, NBoundary.data() );

        dN_dY += n_dA * NBoundary;
        if ( d2N_dYdY )
          for ( int J = 0; J < nDim; J++ )
            ( *d2N_dYdY )[J] += n_dA( J ) * dN_dY_Boundary;
      }
    }

    const double V = volume( smoothingDomain );
    dN_dY /= V;
    if ( d2N_dYdY )
      for ( int J = 0; J < nDim; J++ )
        ( *d2N_dYdY )[J] /= V;
  }

} // namespace Marmot::Discretizations::Meshfree
