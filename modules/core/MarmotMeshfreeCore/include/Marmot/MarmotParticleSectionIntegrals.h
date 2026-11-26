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

#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMeshfreeApproximation.h"
#include "Marmot/MarmotMeshfreeKernelFunction.h"
#include "Marmot/MarmotMeshfreeQuadHexCell.h"
#include "Marmot/MarmotParticle.h"
#include "Marmot/MarmotParticleStiffness.h"
#include "Marmot/MarmotUtils.h"
#include <Eigen/Core>
#include <Fastor/Fastor.h>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

/**
 * @file
 * @brief The integrals over the smoothing domain of a particle in its section: its volume, the smoothed gradients of
 * the shape functions, its second moments and the face loads, one implementation for the SQCNI, SNNI, NSNI and SDI
 * particles of both families; and the hoop terms of the axisymmetric weak form. The axisymmetric derivation is on the
 * page "Axisymmetric Particles" of the documentation.
 */

namespace Marmot::Discretizations::Meshfree {

  /**
   * @brief The integrals over the smoothing domain of a particle in its section: its volume, the smoothed gradients of
   * the shape functions (stabilized conforming nodal integration, Chen et al., 2001), its second moments and the face
   * loads.
   *
   * @details The notation is that of the page "Meshfree Particles": @f$ N_{A,i} = \partial N_A / \partial Y_i @f$
   * with the coordinates @f$ \boldsymbol{Y} @f$ of the intermediate configuration, the volume @f$ |\Omega_s| @f$ of
   * the smoothing domain @f$ \Omega_s @f$ (volume()), and the area vector @f$ n_{f,i}\,dA_f @f$ of its face @f$ f @f$
   * (LagrangeCell::boundarySurfaceVector()). In 2D, @f$ dA_f @f$ is the length of the face (the area of unit
   * thickness in plane strain).
   *
   * The smoothed gradient is the average of the gradient over @f$ \Omega_s @f$, turned into a boundary integral by
   * the divergence theorem and integrated with a rule on each face (faceQuadrature()). In the Cartesian section (plane
   * strain, 3D) the rule is one point per face, the face center @f$ \boldsymbol{Y}_f @f$,
   * @f[
   *   \tilde N_{A,i} = \frac{1}{|\Omega_s|}\int_{\Omega_s} N_{A,i}\,dV
   *   = \frac{1}{|\Omega_s|}\oint_{\partial\Omega_s} N_A\,n_i\,dA
   *   \approx \frac{1}{|\Omega_s|} \sum_f N_A(\boldsymbol{Y}_f)\, n_{f,i}\, dA_f .
   * @f]
   * Here the face rule is the only approximation, and only the values of the shape functions at the face points are
   * needed. NSNI smooths the second derivatives alike, from the (unsmoothed) gradients at the same points,
   * @f[
   *   \hat N_{A,IJ} = \frac{1}{|\Omega_s|} \sum_f N_{A,I}(\boldsymbol{Y}_f)\, n_{f,J}\, dA_f ,
   * @f]
   * stored unsymmetrized as `d2N_dYdY[J](I, A)`; NSNI::secondGradientOfNode() symmetrizes them.
   *
   * In the axisymmetric section, @f$ \Omega_s @f$ is the ring swept by the smoothing domain of the meridian half plane
   * @f$ (r, z) = (Y_1, Y_2) @f$ (`Y( 0 )` is @f$ r @f$ in the code), and every integral carries the ring weight
   * @f$ 2\pi r @f$ (weight()): @f$ |\Omega_s| = 2\pi A_s \bar r @f$, with the area @f$ A_s @f$ and the radius
   * @f$ \bar r @f$ of the area centroid @f$ \bar{\boldsymbol{Y}} @f$ of the smoothing domain (the particle center).
   * The divergence theorem of the weighted gradient gains an area term, @f$ r N_{A,r} = (r N_A)_{,r} - N_A @f$, which
   * is integrated with the centroid; the faces are integrated with two Gauss points @f$ \boldsymbol{Y}_g @f$ each,
   * which share @f$ n_{f,i}\,dA_f @f$ equally:
   * @f[
   *   \tilde N_{A,i} = \frac{1}{|\Omega_s|} \sum_f \sum_{g=1}^{2} 2\pi r_g\, N_A(\boldsymbol{Y}_g)\,
   *   \frac{n_{f,i}\, dA_f}{2} - \delta_{ir}\, \frac{N_A(\bar{\boldsymbol{Y}})}{\bar r} .
   * @f]
   * The hoop operator @f$ \tilde H_A = N_A(\bar{\boldsymbol{Y}}) / \bar r @f$ (in Hoop::force() and
   * Hoop::addStiffnessOfNode()) is the ring average of @f$ N_A / r @f$ with the same area integral. Together they
   * satisfy the integration constraint of the ring, the discrete divergence theorem of the 3D field
   * @f$ N_A \boldsymbol{e}_r @f$, whose divergence is @f$ N_{A,r} + N_A / r @f$: for every smoothing domain,
   * @f[
   *   |\Omega_s| \bigl(\tilde N_{A,r} + \tilde H_A\bigr) = \sum_f \sum_{g=1}^{2} 2\pi r_g\, N_A(\boldsymbol{Y}_g)\,
   *   \frac{n_{f,r}\, dA_f}{2} , \qquad
   *   |\Omega_s|\, \tilde N_{A,z} = \sum_f \sum_{g=1}^{2} 2\pi r_g\, N_A(\boldsymbol{Y}_g)\, \frac{n_{f,z}\, dA_f}{2} ,
   * @f]
   * the area terms cancelling within the smoothing domain. Summed over smoothing domains that tile the body, the face
   * sums of the shared faces cancel, and what remains is the face sum over the boundary of the body. A face on the axis
   * has @f$ r_g = 0 @f$ and drops out: the axis is no boundary of the body of revolution. A face load is integrated
   * with the same rule (addPressure()), so that the load @f$ \boldsymbol{\sigma}\,\boldsymbol{n} @f$ of a uniform
   * stress balances the internal force of a linear field exactly. The second derivatives of NSNI are smoothed alike,
   * @f[
   *   \hat N_{A,IJ} = \frac{1}{|\Omega_s|} \sum_f \sum_{g=1}^{2} 2\pi \Big( r_g\, N_{A,I}(\boldsymbol{Y}_g)\, n_{f,J}
   *   - \delta_{Jr}\, N_A(\boldsymbol{Y}_g)\, n_{f,I} \Big) \frac{dA_f}{2} ,
   * @f]
   * and the second moments are those of the ring about the area centroid, @f$ 2\pi \int \boldsymbol{\xi} \otimes
   * \boldsymbol{\xi}\, r\, dA @f$.
   *
   * The particles take these integrals and their evaluation points (faceQuadraturePoints()) from here and not from the
   * cell: the face rule of the smoothing is defined in this class only. So do the SQCNI and NSNI particles whose
   * smoothing domain follows the deformation gradient (the update `DeformationGradient`, with which the smoothing
   * domain is the particle) take their face loads (addPressure()). The SDI particle
   * (DisplacementParticleSQCNIxSDI::computeDistributedLoad()) and the gradient-enhanced particle with another update of
   * its smoothing domain (GradientEnhancedFiniteStrainParticleSQCNI::faceQuadratureIntermediate()) integrate a
   * pressure at a face center of their own, which does not balance a uniform stress exactly. *Section* is meant as for
   * finite elements (plane strain, 3D, axisymmetric): the form of the integrals depends on it; in the Cartesian
   * section the integrals carry the weight 1. The class holds static functions only, so that a particle selects its
   * integrals with its own template parameters.
   *
   * @tparam section   The section.
   * @tparam nDim      The spatial dimension; 2 for the axisymmetric section.
   * @tparam nVertices The number of vertices of the cell (4: Quad, 8: Hexa).
   */
  template < Section section, int nDim, int nVertices >
  class ParticleSectionIntegrals {
    static_assert( section == Section::Cartesian || nDim == 2,
                   "an axisymmetric particle lives in the meridian half plane: nDim = 2" );

  public:
    using Cell = LagrangeCell< nDim, nVertices >;  ///< a smoothing domain (a Lagrange cell)
    using Vec  = Eigen::Matrix< double, nDim, 1 >; ///< a point or a vector

    /// @brief The section is axisymmetric.
    static constexpr bool isAxisymmetric = section == Section::Axisymmetric;

    /// @brief The number of quadrature points per face of a smoothing domain: one, the face center (Cartesian), or two
    /// Gauss points (axisymmetric).
    static constexpr int nPointsPerFace = isAxisymmetric ? 2 : 1;

    /// @brief The quadrature points of one face, nDim x nPointsPerFace.
    using FacePoints = Eigen::Matrix< double, nDim, nPointsPerFace >;

    /**
     * @brief The quadrature of one face: its points and the share of each point in the area vector of the face.
     * @details A face integral @f$ \int_f g\, w\, n_i\, dA @f$ is @f$ \sum_g g(\boldsymbol{Y}_g)\,
     * w(\boldsymbol{Y}_g)\, @f$ `n_dA( i )`, with the weight @f$ w @f$ of the section (weight()).
     */
    struct FaceQuadrature {
      FacePoints points; ///< the points @f$ \boldsymbol{Y}_g @f$, nDim x nPointsPerFace
      Vec        n_dA;   ///< the share of each point, @f$ n_{f,i}\,dA_f @f$ / nPointsPerFace (without the weight)
    };

    /**
     * @brief The volume of a smoothing domain.
     * @param[in] cell The smoothing domain.
     * @return Its volume @f$ |\Omega_s| @f$ (area in plane strain), or its ring volume @f$ 2\pi A_s \bar r @f$
     * (Pappus).
     */
    static double volume( const Cell& cell )
    {
      if constexpr ( isAxisymmetric )
        return 2.0 * Marmot::Math::Constants::Pi * cell.volume() * cell.centroid()( 0 );
      else
        return cell.volume();
    }

    /**
     * @brief The weight of the section at a point, the factor of the area element of a face.
     * @param[in] Y The point.
     * @return 1, or the ring @f$ 2\pi r @f$.
     */
    static double weight( const Vec& Y )
    {
      if constexpr ( isAxisymmetric )
        return 2.0 * Marmot::Math::Constants::Pi * Y( 0 );
      else
        return 1.0;
    }

    /**
     * @brief The number of face quadrature points of a smoothing domain.
     * @param[in] cell The smoothing domain.
     * @return The number of face quadrature points, `faceQuadraturePoints( cell ).cols()`.
     */
    static int numberOfFaceQuadraturePoints( const Cell& cell ) { return cell.getNumberOfFaces() * nPointsPerFace; }

    /**
     * @brief The quadrature points of one face of a smoothing domain: its center, or its two Gauss points.
     * @param[in] cell The smoothing domain.
     * @param[in] faceID The face (1-based).
     * @return nDim x nPointsPerFace: the face center (Cartesian), or the Gauss points @f$ \boldsymbol{Y}_f \mp
     * \boldsymbol{e}_f / (2\sqrt{3}) @f$ of the straight edge @f$ \boldsymbol{e}_f @f$ (axisymmetric).
     */
    static FacePoints faceQuadraturePoints( const Cell& cell, int faceID )
    {
      if constexpr ( isAxisymmetric ) {
        // the edge from the surface vector (e_1, -e_0) of a counter-clockwise quadrilateral
        const Vec  n_dA = cell.boundarySurfaceVector( faceID );
        const Vec  edge( -n_dA( 1 ), n_dA( 0 ) );
        const Vec  center = cell.getFaceCenterCoordinates( faceID );
        FacePoints points;
        for ( int g = 0; g < nPointsPerFace; g++ )
          points.col( g ) = center + ( g == 0 ? -1.0 : 1.0 ) / ( 2.0 * std::sqrt( 3.0 ) ) * edge;
        return points;
      }
      else
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
     * @brief The quadrature of one face of a cell: the rule of the smoothing, for the smoothing itself and for the
     * loads on the face.
     * @param[in] cell The cell.
     * @param[in] faceID The face (1-based).
     * @return The points (faceQuadraturePoints()) and the share of each in the area vector of the face.
     */
    static FaceQuadrature faceQuadrature( const Cell& cell, int faceID )
    {
      return { faceQuadraturePoints( cell, faceID ), cell.boundarySurfaceVector( faceID ) / double( nPointsPerFace ) };
    }

    /**
     * @brief The smoothed gradients (and second derivatives) of the shape functions over a smoothing domain, see the
     * class description.
     * @param[in] smoothingDomain The smoothing domain; its area centroid is the particle center.
     * @param[in] approximation The meshfree approximation.
     * @param[in] kernelFunctions The kernel functions of the nodes.
     * @param[in] N The shape functions at the particle center, 1 x nNodes: the area term of the ring (axisymmetric);
     * the Cartesian smoothing uses only its size (the number of nodes).
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

    /**
     * @brief The second moments of a domain about its area centroid in the section, @f$ \int \boldsymbol{\xi} \otimes
     * \boldsymbol{\xi}\, dV @f$.
     * @param[in] cell The domain.
     * @return The second moments of the cell, or of its ring (a 3 x 3 Gauss rule, exact for a quadrilateral).
     */
    static Eigen::Matrix< double, nDim, nDim > secondMoments( const Cell& cell );

    /**
     * @brief Refuse a particle the section cannot integrate; nothing to check for the Cartesian section.
     * @details An axisymmetric particle is refused with a vertex at @f$ r < 0 @f$ (beyond the round-off of a mesher,
     * relative to the diagonal of the bounding box), or with smoothing domains that do not tile the body (any update
     * of the smoothing domain but `DeformationGradient`): the integration constraint needs them. It may touch the axis
     * with a vertex or a face: the ring weight vanishes there, and its center has @f$ \bar r > 0 @f$.
     * @param[in] cell The undeformed domain of the particle.
     * @param[in] smoothingDomainsTileTheBody The smoothing domains follow the deformation of the particles
     * (`DeformationGradient`).
     * @param[in] particleLabel The label of the particle, for the message.
     * @throws std::invalid_argument for such a particle.
     */
    static void checkParticle( const Cell& cell, bool smoothingDomainsTileTheBody, int particleLabel );

    /**
     * @brief The distributed loads on the faces of an SQCNI particle in the section.
     * @tparam pressure The load type of `PRESSURE` of the particle family.
     * @tparam cwfCorrection The load type of `CWFCORRECTION` of the particle family.
     * @return `PRESSURE` and `CWFCORRECTION`; `PRESSURE` only for the axisymmetric section, for which the CWF
     * correction is not derived.
     */
    template < int pressure, int cwfCorrection >
    static const std::unordered_map< std::string, int >& supportedFaceLoads()
    {
      static const std::unordered_map< std::string, int >
        faceLoads = isAxisymmetric ? std::unordered_map< std::string, int >{ { "PRESSURE", pressure } }
                                   : std::unordered_map< std::string, int >{ { "PRESSURE", pressure },
                                                                             { "CWFCORRECTION", cwfCorrection } };
      return faceLoads;
    }

    /**
     * @brief Refuse a distributed load on a face that an SQCNI particle in the section does not support
     * (supportedFaceLoads()): the CWF correction of the axisymmetric section.
     * @param[in] isCWFCorrection The load is the CWF correction.
     * @param[in] caller The function that applies the load, for the message.
     * @throws std::invalid_argument for the CWF correction of an axisymmetric particle.
     */
    static void checkFaceLoadIsAvailable( bool isCWFCorrection, const char* caller )
    {
      if ( isAxisymmetric && isCWFCorrection )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << caller << ": CWFCORRECTION is not available for an axisymmetric particle" );
    }

    /**
     * @brief Add the follower pressure on one face of a particle to its external force, and its tangent.
     * @details With the face rule of the smoothing (faceQuadrature()), each point @f$ \boldsymbol{Y}_g @f$ of the face
     * carries the force of Nanson's formula in the plane of the section, from @f$ \boldsymbol{Y} @f$ to
     * @f$ \boldsymbol{x} @f$, times the weight of the section after the increment,
     * @f[
     *   \boldsymbol{f}_g = -\Delta J\, \Delta\boldsymbol{F}^{-\mathsf{T}}\, p\, \frac{\boldsymbol{n}_f\,dA_f}{n_g}\;
     *   w_g , \qquad f_{Ai} \mathrel{-}= \sum_g T_A(\boldsymbol{Y}_g)\, f_{g,i} ,
     * @f]
     * with @f$ w_g = 1 @f$ (Cartesian) or the ring of the point after the increment, @f$ w_g = 2\pi
     * x_r(\boldsymbol{Y}_g)
     * @f$ (Hoop::ringOfFacePoint()): the face moves with the motion of the particle, and its ring grows with the hoop
     * stretch of the point. As for the elements, a positive pressure pushes on the face (against the outward normal).
     * The tangent differentiates @f$ \Delta\boldsymbol{F} @f$ (through @f$ \tilde N_{B,L} @f$) and, axisymmetric,
     * @f$ w_g @f$ (Hoop::dRingOfFacePoint_dDeltaq_Br()). For a uniform stress, the pressure @f$ -\sigma_{nn} @f$ on
     * the boundary balances the internal force of a linear field exactly, because the face rule is that of the
     * integration constraint.
     * @tparam nodeBlockSize The number of dofs per node; the displacement comes first.
     * @tparam MaterialPoint The material point of the particle.
     * @param[in] face The face in the intermediate configuration.
     * @param[in] p The pressure.
     * @param[in] materialPoint The material point at the center, evaluated for the increment: its
     * @f$ \Delta\boldsymbol{F} @f$ in the plane of the section, and (axisymmetric) its @f$ \Delta F_{\theta\theta} @f$
     * and its position, the particle center @f$ \bar{\boldsymbol{Y}} @f$.
     * @param[in] N The trial functions at the center, 1 x nNodes.
     * @param[in] dN_dY Their smoothed gradients, nDim x nNodes.
     * @param[in] approximation The meshfree approximation, for the test functions at the points.
     * @param[in] kernelFunctions The kernel functions of the nodes.
     * @param[in,out] fExt The external force, node-wise.
     * @param[in,out] dFExt_ddQ Its tangent, column-major; nullptr for none (explicit).
     */
    template < int nodeBlockSize, class MaterialPoint >
    static void addPressure( const FaceQuadrature&                       face,
                             double                                      p,
                             const MaterialPoint&                        materialPoint,
                             const Eigen::MatrixXd&                      N,
                             const Eigen::MatrixXd&                      dN_dY,
                             const Approximation&                        approximation,
                             const std::vector< const KernelFunction* >& kernelFunctions,
                             double*                                     fExt,
                             double*                                     dFExt_ddQ );
  };

  /**
   * @brief Refuse VCI for a section it is not derived for.
   * @details VCI (Chen, Hillman and Rüter, 2013) is derived for the Cartesian section; for the ring it would have to be
   * derived again, with the ring measure and the hoop term, and the ring-weighted SQCNI satisfies the linear
   * integration constraint by itself.
   * @tparam section The section of the particle.
   * @param[in] caller The caller, for the message.
   * @throws std::invalid_argument for the axisymmetric section.
   */
  template < Section section >
  void checkVCIIsAvailable( [[maybe_unused]] const std::string& caller )
  {
    if constexpr ( section == Section::Axisymmetric )
      throw std::invalid_argument( caller +
                                   ": VCI is not available for an axisymmetric particle; its ring-weighted smoothing "
                                   "satisfies the linear integration constraint by itself" );
  }

  template < Section section, int nDim, int nVertices >
  void ParticleSectionIntegrals< section, nDim, nVertices >::smoothShapeFunctionGradients(
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
      const FaceQuadrature face = faceQuadrature( smoothingDomain, f + 1 );

      for ( int g = 0; g < nPointsPerFace; g++ ) {
        const Vec Y = face.points.col( g );

        if ( d2N_dYdY )
          approximation.computeShapeFunctionsAndGradients( Y.data(),
                                                           kernelFunctions,
                                                           NBoundary.data(),
                                                           dN_dY_Boundary.data() );
        else
          approximation.computeShapeFunctions( Y.data(), kernelFunctions, NBoundary.data() );

        // the point's share of the face, weighted by the section (the 2 pi r of the ring)
        const Vec ring_n_dA = weight( Y ) * face.n_dA;
        dN_dY += ring_n_dA * NBoundary;
        if ( d2N_dYdY ) {
          for ( int J = 0; J < nDim; J++ )
            ( *d2N_dYdY )[J] += ring_n_dA( J ) * dN_dY_Boundary;
          // the area term of the ring in the second derivatives, int dN/dY_I dA = oint N n_I ds
          if constexpr ( isAxisymmetric )
            ( *d2N_dYdY )[0] -= ( 2.0 * Marmot::Math::Constants::Pi * face.n_dA ) * NBoundary;
        }
      }
    }

    // the area term of the ring, -delta_ir 2 pi int N dA, with the centroid rule: -2 pi A_s N(center)
    if constexpr ( isAxisymmetric ) {
      const double A_s = smoothingDomain.volume(); // the area of the smoothing domain in the meridian half plane
      dN_dY.row( 0 ) -= ( 2.0 * Marmot::Math::Constants::Pi * A_s ) * N;
    }

    const double V = volume( smoothingDomain );
    dN_dY /= V;
    if ( d2N_dYdY )
      for ( int J = 0; J < nDim; J++ )
        ( *d2N_dYdY )[J] /= V;
  }

  template < Section section, int nDim, int nVertices >
  Eigen::Matrix< double, nDim, nDim > ParticleSectionIntegrals< section, nDim, nVertices >::secondMoments(
    const Cell& cell )
  {
    if constexpr ( isAxisymmetric ) {
      const Vec                           C       = cell.centroid();
      const std::array< double, 3 >       points  = { -std::sqrt( 0.6 ), 0.0, std::sqrt( 0.6 ) };
      const std::array< double, 3 >       weights = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
      Eigen::Matrix< double, nDim, nDim > moments = Eigen::Matrix< double, nDim, nDim >::Zero();
      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ ) {
          const Vec    s( points[i], points[j] );
          const Vec    Y  = cell.mapToPhysical( s );
          const Vec    xi = Y - C;
          const double dV = weight( Y ) * std::abs( cell.jacobian( s ).determinant() ) * weights[i] * weights[j];
          moments += xi * xi.transpose() * dV;
        }
      return moments;
    }
    else
      return cell.secondMoments();
  }

  template < Section section, int nDim, int nVertices >
  void ParticleSectionIntegrals< section, nDim, nVertices >::checkParticle(
    [[maybe_unused]] const Cell& cell,
    [[maybe_unused]] bool        smoothingDomainsTileTheBody,
    [[maybe_unused]] int         particleLabel )
  {
    if constexpr ( isAxisymmetric ) {
      // the round-off of a mesher, relative to the size of the particle, as for the axisymmetric elements
      constexpr double axisTolerance = 1e-12;
      const auto&      vertices      = cell.nodes();
      const double     h             = ( vertices.rowwise().maxCoeff() - vertices.rowwise().minCoeff() ).norm();
      Eigen::Index     i             = 0;
      const double     r             = vertices.row( 0 ).minCoeff( &i );
      if ( r < -axisTolerance * h )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << "axisymmetric particle " << particleLabel << ": vertex " << i + 1
                                     << " lies at r = " << r
                                     << " < 0; axisymmetric particles live in the meridian half plane r >= 0." );
      if ( !smoothingDomainsTileTheBody )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << "axisymmetric particle " << particleLabel
                                     << ": it needs the smoothing domain update DeformationGradient, whose smoothing "
                                        "domains tile the body, for the integration constraint of the ring." );
    }
  }

  /**
   * @brief The hoop terms of the weak form of an axisymmetric particle, shared by its formulations.
   *
   * @details The increment of the hoop stretch is @f$ \Delta F_{\theta\theta} = x_r / r = 1 + \Delta u_r / \bar r @f$
   * at the particle center, with the radial displacement increment @f$ \Delta u_r = N_B\, \Delta q_{Br} @f$ there and
   * its radius @f$ \bar r @f$ in the intermediate configuration. The virtual hoop strain is @f$ \delta d_{\theta\theta}
   * = \delta x_r / x_r = \delta u_r / (\bar r\, \Delta F_{\theta\theta}) @f$, so that the virtual work of the hoop
   * stress over the undeformed ring, @f$ \tau_{\theta\theta}\, \delta d_{\theta\theta}\, V_0 @f$, adds to the radial
   * residual of node @f$ A @f$ the term @f$ T_A\, \tau_{\theta\theta} / (\bar r\, \Delta F_{\theta\theta})\, V_0 @f$,
   * with the test function @f$ T_A @f$ of the node at the center: the term @f$ T_A / \bar r @f$ is the hoop operator
   * @f$ \tilde H_A = N_A(\bar{\boldsymbol{Y}}) / \bar r @f$ of ParticleSectionIntegrals with the test function (see
   * there for why it is the point value).
   * The faces of a particle carry the ring of their points (ringOfFacePoint()). The derivation is on the page
   * "Axisymmetric Particles" of the documentation.
   */
  namespace Hoop {

    /**
     * @brief The radius of the particle center in the intermediate configuration, where the shape functions are
     * evaluated.
     * @tparam MaterialPoint The material point of the particle.
     * @param[in] materialPoint The material point at the center.
     * @return @f$ \bar r = \bar Y_1 @f$.
     */
    template < class MaterialPoint >
    double radiusOfCenter( const MaterialPoint& materialPoint )
    {
      std::array< double, 2 > center;
      materialPoint.getCoordinatesAtCenter( center.data() );
      return center[0];
    }

    /**
     * @brief Refuse an increment that pulls the ring of a particle through the axis.
     * @param[in] deltaF33 The increment of the hoop stretch, @f$ \Delta F_{\theta\theta} = 1 + \Delta u_r / \bar r
     * @f$.
     * @param[in] r The radius @f$ \bar r @f$ of the center in the intermediate configuration.
     * @param[in] particleLabel The label of the particle, for the message.
     * @throws std::runtime_error if @p deltaF33 is not positive: the increment cannot be integrated, which a solver may
     * answer with a smaller increment.
     */
    inline void checkStretchIncrement( double deltaF33, double r, int particleLabel )
    {
      if ( !( deltaF33 > 0.0 ) )
        throw std::runtime_error( Marmot::Utilities::MakeString()
                                  << "axisymmetric particle " << particleLabel
                                  << ": the hoop stretch increment 1 + du_r / r = " << deltaF33 << " at r = " << r
                                  << " is not positive; the ring is pulled through the axis." );
    }

    /**
     * @brief Add the hoop stretch increment of the center, @f$ \Delta u_r / \bar r @f$, to the out-of-plane stretch
     * increment @f$ \Delta F_{33} @f$ of the material point (1 after its prepareYourself()), as the axisymmetric
     * elements do, and refuse it if it pulls the ring through the axis (checkStretchIncrement()).
     * @tparam MaterialPoint The 2D material point of the particle.
     * @param[in,out] materialPoint The material point at the center, with the in-plane increment added.
     * @param[in] du_r The radial displacement increment @f$ \Delta u_r @f$ at the center.
     * @throws std::runtime_error see checkStretchIncrement().
     */
    template < class MaterialPoint >
    void addStretchIncrement( MaterialPoint& materialPoint, double du_r )
    {
      const double r = radiusOfCenter( materialPoint );
      materialPoint.addToOutOfPlaneStretchIncrement( du_r / r );
      checkStretchIncrement( materialPoint.outOfPlaneStretchIncrement(), r, materialPoint.getMaterialPointNumber() );
    }

    /**
     * @brief The hoop part of the radial residual of node @f$ A @f$, @f$ T_A\, \tau_{\theta\theta} / (\bar r\, \Delta
     * F_{\theta\theta})\, V_0 @f$.
     * @param[in] T_A The test function of the node at the center.
     * @param[in] r The radius @f$ \bar r @f$ of the center.
     * @param[in] deltaF33 The increment of the hoop stretch.
     * @param[in] tau33 The hoop Kirchhoff stress.
     * @param[in] V0 The undeformed ring volume of the particle.
     * @return The force.
     */
    inline double force( double T_A, double r, double deltaF33, double tau33, double V0 )
    {
      return T_A / ( r * deltaF33 ) * tau33 * V0;
    }

    /**
     * @brief Add the hoop terms of the displacement rows to the stiffness of node @f$ A @f$ (see StiffnessOfNode; the
     * displacement comes first among the dofs of a node): the hoop row through the in-plane gradient and through the
     * radial displacement increment, and every displacement row through the hoop stretch in the stress,
     * @f[
     *   \frac{\partial r^U_{Ar}}{\partial (\partial \Delta u_k / \partial Y_l)} \mathrel{+}= c\,
     *   \frac{\partial \tau_{\theta\theta}}{\partial \Delta F_{kl}}, \qquad
     *   \frac{\partial r^U_{Ar}}{\partial \Delta u_r} \mathrel{+}= \frac{c}{\bar r} \Big( \frac{\partial
     *   \tau_{\theta\theta}}{\partial \Delta F_{\theta\theta}} - \frac{\tau_{\theta\theta}}{\Delta F_{\theta\theta}}
     *   \Big), \qquad
     *   \frac{\partial r^U_{Aj}}{\partial \Delta u_r} \mathrel{+}= \frac{V_0}{\bar r} \frac{\partial T_A}{\partial
     *   x_i} \frac{\partial \tau_{ij}}{\partial \Delta F_{\theta\theta}},
     * @f]
     * with @f$ c = T_A V_0 / (\bar r\, \Delta F_{\theta\theta}) @f$.
     * @tparam nDofPerNode The number of dofs per node.
     * @param[in] T_A The test function of the node at the center.
     * @param[in] dT_A_dx Its spatial gradient.
     * @param[in] r The radius @f$ \bar r @f$ of the center.
     * @param[in] deltaF33 The increment of the hoop stretch.
     * @param[in] V0 The undeformed ring volume of the particle.
     * @param[in] tau33 The hoop Kirchhoff stress.
     * @param[in] dTau_dDeltaF33 @f$ \partial \tau_{ij} / \partial \Delta F_{\theta\theta} @f$, in-plane @f$ ij @f$.
     * @param[in] dTau33_dDeltaF @f$ \partial \tau_{\theta\theta} / \partial \Delta F_{kl} @f$, in-plane @f$ kl @f$.
     * @param[in] dTau33_dDeltaF33 @f$ \partial \tau_{\theta\theta} / \partial \Delta F_{\theta\theta} @f$.
     * @param[in,out] ofNode The stiffness of the node; the hoop terms are added.
     */
    template < int nDofPerNode >
    void addStiffnessOfNode( double                                T_A,
                             const Fastor::Tensor< double, 2 >&    dT_A_dx,
                             double                                r,
                             double                                deltaF33,
                             double                                V0,
                             double                                tau33,
                             const Fastor::Tensor< double, 2, 2 >& dTau_dDeltaF33,
                             const Fastor::Tensor< double, 2, 2 >& dTau33_dDeltaF,
                             double                                dTau33_dDeltaF33,
                             StiffnessOfNode< 2, nDofPerNode >&    ofNode )
    {
      const double c = T_A * V0 / ( r * deltaF33 );
      for ( int l = 0; l < 2; l++ )
        for ( int k = 0; k < 2; k++ )
          ofNode.gradientBlock( l )( 0, k ) += c * dTau33_dDeltaF( k, l );
      ofNode.df_dIncrement( 0, 0 ) += c / r * ( dTau33_dDeltaF33 - tau33 / deltaF33 );
      for ( int j = 0; j < 2; j++ ) {
        double dfj_dDeltaF33 = 0.0;
        for ( int i = 0; i < 2; i++ )
          dfj_dDeltaF33 += dT_A_dx( i ) * dTau_dDeltaF33( i, j );
        ofNode.df_dIncrement( j, 0 ) += V0 / r * dfj_dDeltaF33;
      }
    }

    /**
     * @brief The ring of a point of a face of the particle after the increment, @f$ 2\pi x_r(\boldsymbol{Y}) @f$.
     * @details The face moves with the motion of the particle, @f$ \boldsymbol{x}(\boldsymbol{Y}) =
     * \bar{\boldsymbol{x}}
     * + \Delta\boldsymbol{F}\,(\boldsymbol{Y} - \bar{\boldsymbol{Y}}) @f$, the motion Nanson's formula of the face
     * assumes, with the center at @f$ \bar x_r = \bar r\,\Delta F_{\theta\theta} @f$. The ring @f$ 2\pi x_r @f$ is the
     * ring weight @f$ 2\pi r @f$ of the point times its hoop stretch increment @f$ x_r / r @f$, by which the area of
     * the face grows beyond that of the half plane. For a linear increment it is exact; at the center it is @f$ 2\pi
     * \bar r\,\Delta F_{\theta\theta} @f$. It needs no division by @f$ r @f$. A point on the axis has the ring 0 for an
     * increment that keeps the axis in place; the axis is no boundary and carries no load.
     * @param[in] Y The point, in the intermediate configuration.
     * @param[in] center The particle center @f$ \bar{\boldsymbol{Y}} @f$.
     * @param[in] deltaF33 @f$ \Delta F_{\theta\theta} @f$ of the particle.
     * @param[in] deltaF @f$ \Delta\boldsymbol{F} @f$ of the particle, in the half plane.
     * @return @f$ 2\pi x_r(\boldsymbol{Y}) @f$.
     */
    inline double ringOfFacePoint( const Eigen::Vector2d&                Y,
                                   const Eigen::Vector2d&                center,
                                   double                                deltaF33,
                                   const Fastor::Tensor< double, 2, 2 >& deltaF )
    {
      const Eigen::Vector2d xi  = Y - center;
      const double          x_r = center( 0 ) * deltaF33 + deltaF( 0, 0 ) * xi( 0 ) + deltaF( 0, 1 ) * xi( 1 );
      return 2.0 * Marmot::Math::Constants::Pi * x_r;
    }

    /**
     * @brief The derivative of ringOfFacePoint() with respect to the radial dof of node @f$ B @f$,
     * @f$ 2\pi \bigl(N_B(\bar{\boldsymbol{Y}}) + \tilde N_{B,j}\,(Y_j - \bar Y_j)\bigr) @f$: through the center
     * @f$ \bar x_r = \bar r + N_B\,\Delta q_{Br} @f$ and through @f$ \Delta F_{rj} = \delta_{rj} + \Delta q_{Br}\,
     * \tilde N_{B,j} @f$.
     * @param[in] Y The point, in the intermediate configuration.
     * @param[in] center The particle center @f$ \bar{\boldsymbol{Y}} @f$.
     * @param[in] N_B The trial function of node @f$ B @f$ at the center.
     * @param[in] dN_B_dY Its smoothed gradient @f$ \tilde N_{B,j} @f$ (2 values).
     * @return @f$ \partial (2\pi x_r) / \partial \Delta q_{Br} @f$.
     */
    inline double dRingOfFacePoint_dDeltaq_Br( const Eigen::Vector2d& Y,
                                               const Eigen::Vector2d& center,
                                               double                 N_B,
                                               const double*          dN_B_dY )
    {
      const Eigen::Vector2d xi = Y - center;
      return 2.0 * Marmot::Math::Constants::Pi * ( N_B + dN_B_dY[0] * xi( 0 ) + dN_B_dY[1] * xi( 1 ) );
    }

  } // namespace Hoop

  template < Section section, int nDim, int nVertices >
  template < int nodeBlockSize, class MaterialPoint >
  void ParticleSectionIntegrals< section, nDim, nVertices >::addPressure(
    const FaceQuadrature&                       face,
    double                                      p,
    const MaterialPoint&                        materialPoint,
    const Eigen::MatrixXd&                      N,
    const Eigen::MatrixXd&                      dN_dY,
    const Approximation&                        approximation,
    const std::vector< const KernelFunction* >& kernelFunctions,
    double*                                     fExt,
    double*                                     dFExt_ddQ )
  {
    using namespace Fastor;
    using namespace Marmot::Math::Tensors::FastorIndices;
    using TensorD = Tensor< double, nDim >;

    const int                     nNodes = static_cast< int >( N.cols() );
    Eigen::Map< Eigen::VectorXd > P( fExt, nNodes * nodeBlockSize );

    // Nanson's formula in the plane of the section, from Y to x
    const Tensor< double, nDim, nDim >             deltaF    = materialPoint.dx_dY();
    const Tensor< double, nDim, nDim >             deltaFInv = inverse( deltaF );
    const double                                   deltaJ    = determinant( deltaF );
    const Tensor< double, nDim, nDim, nDim, nDim > dFInv_dF  = -einsum< Ik, Ki, to_IikK >( deltaFInv, deltaFInv );

    // axisymmetric: the hoop stretch increment and the position of the particle center
    [[maybe_unused]] double deltaF33 = 1.0;
    [[maybe_unused]] Vec    center;
    if constexpr ( isAxisymmetric ) {
      deltaF33 = materialPoint.outOfPlaneStretchIncrement();
      materialPoint.getCoordinatesAtCenter( center.data() );
    }

    Eigen::MatrixXd T = Eigen::MatrixXd::Zero( 1, nNodes );

    for ( int g = 0; g < nPointsPerFace; g++ ) {
      const Vec Y = face.points.col( g );

      // a positive pressure pushes on the face: the traction acts against the outward normal
      const TensorD fY = -p * TensorD( face.n_dA.data() );

      TensorD                            f          = deltaJ * transpose( deltaFInv ) % fY;
      Tensor< double, nDim, nDim, nDim > df_dDeltaF = outer( f, transpose( deltaFInv ) ) +
                                                      deltaJ * einsum< IikK, Index< I_ > >( dFInv_dF, fY );

      // the ring of the point after the increment (a ring force, in the plane of the section times the ring)
      [[maybe_unused]] TensorD f_halfPlane;
      if constexpr ( isAxisymmetric ) {
        f_halfPlane       = f;
        const double ring = Hoop::ringOfFacePoint( Y, center, deltaF33, deltaF );
        f *= ring;
        df_dDeltaF *= ring;
      }

      T.setZero();
      approximation.computeShapeFunctions( Y.data(), kernelFunctions, T.data() );

      for ( int A = 0; A < nNodes; A++ ) {
        const int     idxA = nodeBlockSize * A;
        const TensorD r_U  = T( A ) * f;
        P.template segment< nDim >( idxA ) -= Eigen::Map< const Eigen::Matrix< double, nDim, 1 > >( r_U.data() );

        if ( !dFExt_ddQ )
          continue;

        Eigen::Map< Eigen::MatrixXd > K( dFExt_ddQ, nNodes * nodeBlockSize, nNodes * nodeBlockSize );
        for ( int B = 0; B < nNodes; B++ ) {
          const int  idxB    = nodeBlockSize * B;
          const auto dN_B_dY = TensorMap< const double, nDim >( dN_dY.col( B ).data() );

          const Tensor< double, nDim, nDim > df_ddQU_B = T( A ) * einsum< ijk, k >( df_dDeltaF, dN_B_dY );
          K.template block< nDim, nDim >( idxA, idxB ) -= Eigen::Map< Eigen::Matrix< double, nDim, nDim > >(
            torowmajor( df_ddQU_B ).data() );

          if constexpr ( isAxisymmetric ) {
            const double dRing_dDeltaq_Br = Hoop::dRingOfFacePoint_dDeltaq_Br( Y,
                                                                               center,
                                                                               N( B ),
                                                                               dN_dY.col( B ).data() );
            for ( int i = 0; i < nDim; i++ )
              K( idxA + i, idxB ) -= T( A ) * f_halfPlane( i ) * dRing_dDeltaq_Br;
          }
        }
      }
    }
  }

} // namespace Marmot::Discretizations::Meshfree
