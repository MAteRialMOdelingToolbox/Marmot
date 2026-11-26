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

#include "Marmot/DisplacementParticle.h"
#include "Marmot/MarmotMeshfreeApproximation.h"
#include "Marmot/MarmotParticleDomain.h" // Geometric component
#include "Marmot/MarmotParticleSectionIntegrals.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/src/Core/util/Constants.h>
/* Fastor/Fastor.h ends with `#pragma warning( default : ... )`; scope it so it neither overrides the build's warning
 * flags nor leaks into consumers. */
#ifdef _MSC_VER
#  pragma warning( push )
#endif
#include <Fastor/Fastor.h>
#ifdef _MSC_VER
#  pragma warning( pop )
#endif
#include <cmath>

namespace Marmot::Discretizations::Meshfree {

  /**
   * @class Marmot::Discretizations::Meshfree::DisplacementParticleSQCNI
   * @brief A displacement particle implementing the Stabilized Quasi-Conforming Nodal Integration (SQCNI) method.
   *
   * @details This class extends the basic DisplacementParticle to incorporate the SQCNI formulation,
   * which involves smoothing domains for integration and specific handling of deformation.
   * It uses a ParticleDomain object for its geometric representation via composition.
   *
   * The residual and the tangent are those of DisplacementParticle, evaluated with the shape functions
   * @f$ N_B @f$ at the particle center and the smoothed gradients over the smoothing domain @f$ \Omega_s @f$ of the
   * particle (in the intermediate configuration @f$ \boldsymbol{Y} @f$), computed by a one-point rule at the face
   * centers @f$ \boldsymbol{Y}_f @f$ of the smoothing domain,
   * @f[
   *   \frac{\partial N_B}{\partial Y_i} \approx \frac{1}{|\Omega_s|} \oint_{\partial\Omega_s} N_B\,n_i\,dA
   *   \approx \frac{1}{|\Omega_s|} \sum_f N_B(\boldsymbol{Y}_f)\, n_{f,i}\, dA_f
   * @f]
   * (ParticleSectionIntegrals::smoothShapeFunctionGradients()).
   * The test functions are the trial functions (unless VCI corrects the test gradients). How the smoothing domain
   * follows the deformation is set by the SmoothingDomainUpdateType: with the deformation gradient (SQCNI, the
   * smoothing domain conforms to the deformed particle), not at all (SNNI, it is only translated), with the rotation
   * only, or with the rotation and the principal stretches.
   *
   * The particle geometry also provides the faces for distributed loads: a follower pressure and the consistent weak
   * form (CWF) correction, see computeDistributedLoad().
   *
   * The volume, the smoothed gradients and the face loads are those of the section (ParticleSectionIntegrals): for
   * Section::Axisymmetric, the ring volume, the ring-weighted smoothed gradient with its area term, and two Gauss
   * points per face, both for the smoothing and for a load on the face. An axisymmetric particle needs the update type
   * `DeformationGradient` (smoothing domains that tile the body, for the integration constraint), refuses a vertex at
   * @f$ r < 0 @f$ (ParticleSectionIntegrals::checkParticle()), and supports the distributed load `PRESSURE` only. See
   * the page "Axisymmetric Particles" of the documentation.
   *
   * @tparam nDim The number of dimensions (e.g., 2 for 2D, 3 for 3D).
   * @tparam nVertices The number of vertices defining the particle's geometry.
   * @tparam section The section (Section::Cartesian: plane strain or 3D; Section::Axisymmetric: nDim = 2).
   */
  template < int nDim, int nVertices, Section section = Section::Cartesian >
  class DisplacementParticleSQCNI : public DisplacementParticle< nDim, section > {

    using TensorD  = Fastor::Tensor< double, nDim >;       ///< vector of size nDim
    using TensorDD = Fastor::Tensor< double, nDim, nDim >; ///< second-order tensor of size nDim

  protected:
    using ParentPointParticle = DisplacementParticle< nDim, section >; ///< the point particle base class
    using ParticleDomainType  = ParticleDomain< nDim, nVertices >;     ///< the particle geometry type

    ParticleDomainType _particleDomain; ///< geometry and smoothing domain of the particle

    /// the integrals over the smoothing domain
    using SectionIntegrals = ParticleSectionIntegrals< section, nDim, nVertices >;

    /**
     * @brief The smoothing domain of the particle.
     * @return The cell of the smoothing domain in the intermediate configuration.
     */
    const typename ParticleDomainType::LagrangeCellType& smoothingDomain() const
    {
      return _particleDomain.smoothingDomain();
    }

  public:
    /// how the smoothing domain follows the deformation
    using SmoothingDomainUpdateType = ParticleDomainType::SmoothingDomainUpdateType;

    /**
     * @brief Constructs a new DisplacementParticleSQCNI object.
     *
     * @param[in] elementID The unique identifier for the element.
     * @param[in] nodeCoordinates Pointer to an array of vertex coordinates (nDim * nVertices) in the undeformed
     * configuration; the center and the volume of the particle follow from them.
     * @param[in] nNodeCoordiantes The total number of coordinate values (nDim * nVertices).
     * @param[in] volume The volume of the particle. Must be 0, as the volume is computed from vertex coordinates.
     * @param[in] materialName The name of the material assigned to the particle.
     * @param[in] materialProperties Pointer to an array of material properties.
     * @param[in] sizeMaterialProperties The size of the material properties array.
     * @param[in] approximation The meshfree approximation object used for shape functions.
     * @param[in] smoothingVolumeUpdateType The strategy for updating the smoothing domain's volume.
     * @throws std::invalid_argument if the provided volume is not zero, and for an axisymmetric particle the section
     * cannot integrate (ParticleSectionIntegrals::checkParticle()).
     */
    DisplacementParticleSQCNI( int                                                          elementID,
                               const double*                                                nodeCoordinates,
                               int                                                          nNodeCoordiantes,
                               double                                                       volume,
                               const std::string&                                           materialName,
                               const double*                                                materialProperties,
                               int                                                          sizeMaterialProperties,
                               const Approximation&                                         approximation,
                               const typename ParticleDomainType::SmoothingDomainUpdateType smoothingVolumeUpdateType );

    // Override Particle interface methods to delegate to _particleDomain
    /**
     * @brief Vertices of the particle geometry in the intermediate configuration.
     * @param[out] coordinates The vertex coordinates (nDim x nVertices, column major).
     */
    void getVertexCoordinates( double* coordinates ) const override
    {
      Eigen::Map< Eigen::Matrix< double, nDim, nVertices > > coordinatesMap( coordinates );
      coordinatesMap = _particleDomain.getGeometryDeformedVertexCoordinates();
    }

    /**
     * @brief Vertices of the smoothing domain, used for visualization.
     * @param[out] coordinates The vertex coordinates (nDim x nVertices, column major).
     */
    void getVisualizationVertexCoordinates( double* coordinates ) const override
    {
      Eigen::Map< Eigen::Matrix< double, nDim, nVertices > > coordinatesMap( coordinates );
      coordinatesMap = _particleDomain.getSmoothingVertexCoordinates();
    }

    /**
     * @brief Number of vertices of the particle geometry.
     * @return nVertices.
     */
    int getNumberOfVertices() const override { return _particleDomain.getNumberOfVertices(); }

    /**
     * @brief Shape of the particle geometry.
     * @return The shape name, e.g. "quad4" or "hexa8".
     */
    std::string getParticleShape() const override { return _particleDomain.getParticleShape(); }

    /**
     * @brief Center of a face of the particle geometry (intermediate configuration).
     * @param[in] faceID The face ID (1-based).
     * @param[out] coordinates The coordinates (nDim values).
     */
    void getFaceCoordinates( int faceID, double* coordinates ) const override
    {
      Eigen::Map< Eigen::Matrix< double, nDim, 1 > > coordinatesMap( coordinates );
      coordinatesMap = _particleDomain.getFaceCenterCoordinates( faceID );
    }

    /**
     * @brief Evaluation points of the shape functions: the quadrature points of the faces of the smoothing domain (the
     * face centers, or two Gauss points per face of an axisymmetric particle,
     * ParticleSectionIntegrals::faceQuadraturePoints()).
     * @param[out] coordinates The coordinates (nDim x number of evaluation points, column major).
     */
    void getEvaluationCoordinates( double* coordinates ) const override
    {
      Eigen::Map< Eigen::Matrix< double, nDim, Eigen::Dynamic > > coordinatesMap( coordinates,
                                                                                  nDim,
                                                                                  getNumberOfEvaluationPoints() );
      coordinatesMap = SectionIntegrals::faceQuadraturePoints( smoothingDomain() );
    }

    /**
     * @brief Number of evaluation points.
     * @return The number of face quadrature points of the smoothing domain.
     */
    int getNumberOfEvaluationPoints() const override
    {
      return SectionIntegrals::numberOfFaceQuadraturePoints( smoothingDomain() );
    }

    /**
     * @brief Center of the particle in the intermediate configuration, from the material point.
     * @param[out] coordinates The coordinates (nDim values).
     */
    void getCenterCoordinates( double* coordinates ) const override
    {
      // Use the material point's center coordinates, which is physics-specific
      this->_mp->getCoordinatesAtCenter( coordinates );
    }

    /**
     * @brief Volume in the undeformed configuration.
     * @return @f$ V_0 @f$.
     */
    double getVolumeUndeformed() const override { return this->_mp->getVolumeUndeformed(); };

    /**
     * @brief Supported distributed loads, on the faces of the particle domain.
     * @return "PRESSURE" and "CWFCORRECTION"; "PRESSURE" only for an axisymmetric particle.
     */
    const std::unordered_map< std::string, int >& getSupportedDistributedLoadTypes() const override
    {
      return SectionIntegrals::template supportedFaceLoads< ParentPointParticle::Pressure,
                                                            ParentPointParticle::CWFCorrection >();
    };

    /**
     * @brief States: "vertex displacements" and "smoothing vertex displacements" (nDim x nVertices) of the geometry and
     * of the smoothing domain, otherwise the states of the material point.
     * @param[in] stateName Name of the state.
     * @param[in] qp Evaluation point (not used).
     * @return The view on the state.
     */
    Marmot::Utilities::StateView getStateView( const std::string& stateName, int qp ) const override
    {
      if ( stateName == "vertex displacements" )
        return Marmot::Utilities::StateView( const_cast< double* >(
                                               _particleDomain.getGeometryDeformedVertexDisplacements().data() ),
                                             nDim * nVertices );

      if ( stateName == "smoothing vertex displacements" )
        return Marmot::Utilities::StateView( const_cast< double* >(
                                               _particleDomain.getSmoothingDomainVertexDisplacements().data() ),
                                             nDim * nVertices );

      return ParentPointParticle::getStateView( stateName, qp );
    }

    /**
     * @brief Computes the deformed volume of the particle.
     * @return The deformed volume, calculated as undeformed volume multiplied by the determinant of the deformation
     * gradient @f$ \boldsymbol{F}_n @f$ of the last accepted state (i.e., the volume in the intermediate
     * configuration); for an axisymmetric particle the ring volume of its geometry.
     */
    virtual double getVolumeDeformed() const
    {
      if constexpr ( SectionIntegrals::isAxisymmetric )
        return SectionIntegrals::volume( _particleDomain.geometryDeformed() );
      else
        return this->getVolumeUndeformed() * Fastor::determinant( this->dY_dX() );
    }

    /**
     * @brief Accepts the increment (see DisplacementParticle::acceptStateAndPosition()) and moves the particle geometry
     * and the smoothing domain to the new intermediate configuration.
     */
    void acceptStateAndPosition() override
    {
      // First, call the DisplacementParticle's acceptStateAndPosition to update material point and GenericParticle's
      // center and volume.
      ParentPointParticle::acceptStateAndPosition();

      // Then, update the ParticleDomain's geometry based on the current deformation
      const auto F_physics          = Eigen::Matrix< double, nDim, nDim, Eigen::RowMajor >( this->dY_dX().data() );
      const auto centerDisplacement = Eigen::Matrix< double, nDim, 1 >( this->getDisplacementAtCenter().data() );

      _particleDomain.acceptStateAndPosition( F_physics, centerDisplacement );
    }

    /**
     * @brief Evaluate the shape functions with the assigned kernel functions, with smoothed gradients.
     *
     * N at the particle center, dN_dY as the boundary integral over the current smoothing domain, and the test
     * functions T = N, dT_dY = dN_dY.
     */
    void updateShapeFunctions() override
    {
      Eigen::Matrix< double, nDim, 1 > coords;
      this->getCenterCoordinates( coords.data() ); // Uses DisplacementParticleSQCNI's getCenterCoordinates

      this->_N = Eigen::MatrixXd::Zero( 1, this->_nNodes );

      // Compute N at the particle center (from GenericParticle)
      this->_meshfreeApproximation.computeShapeFunctions( coords.data(),
                                                          this->_assignedKernelFunctions,
                                                          this->_N.data() );

      // the smoothed gradients over the smoothing domain, in the section
      SectionIntegrals::smoothShapeFunctionGradients( smoothingDomain(),
                                                      this->_meshfreeApproximation,
                                                      this->_assignedKernelFunctions,
                                                      this->_N,
                                                      this->_dN_dY,
                                                      nullptr ); // no second derivatives (not NSNI)

      this->_T     = this->_N;
      this->_dT_dY = this->_dN_dY;
    }

    /**
     * @brief Computes the boundary integral part for the VCI (Variational Consistent Integration) test function.
     *
     * This method contributes to the R_AiC_RowMajor matrix, which is part of the VCI formulation.
     * It involves shape functions, monomial basis, and boundary surface vectors:
     * @f$ R_{AiC} \mathrel{+}= \sum_g T_A(\boldsymbol{Y}_g)\,P_C(\boldsymbol{Y}_g)\,(n_i\,dA)_g @f$, with the face
     * quadrature of faceQuadratureIntermediate().
     *
     * @param[in,out] R_AiC_RowMajor Pointer to the row-major matrix for VCI constraints.
     * @param[in] boundarySurfaceVector Pointer to the boundary surface vector (not used; the face of the particle
     * geometry is used instead).
     * @param[in] boundaryFaceID The ID of the boundary face.
     */
    void vci_compute_Test_P_BoundaryIntegral( double*       R_AiC_RowMajor,
                                              const double* boundarySurfaceVector,
                                              int           boundaryFaceID ) override
    {
      const auto face = faceQuadratureIntermediate( boundaryFaceID );
      for ( int g = 0; g < SectionIntegrals::nPointsPerFace; g++ ) {
        const Eigen::Matrix< double, nDim, 1 > Y_N       = face.points.col( g );
        const Eigen::Matrix< double, nDim, 1 > N_dAY     = SectionIntegrals::weight( Y_N ) * face.n_dA;
        Eigen::MatrixXd                        TBoundary = Eigen::MatrixXd::Zero( 1, this->_nNodes );

        this->_meshfreeApproximation.computeShapeFunctions( Y_N.data(),
                                                            this->_assignedKernelFunctions,
                                                            TBoundary.data() );

        // get P for the exact integration location at the boundary
        auto PBoundary = this->_P;
        Math::computeMonomialBasis( this->_vciOrder, Y_N, PBoundary );

        for ( int A = 0; A < this->_nNodes; A++ )
          for ( int i = 0; i < nDim; i++ )
            for ( int C = 0; C < this->_nVCIConstraints; C++ )
              R_AiC_RowMajor[A * ( nDim * this->_nVCIConstraints ) + i * this->_nVCIConstraints + C] += TBoundary( A ) *
                                                                                                        PBoundary( C ) *
                                                                                                        N_dAY[i];
      }
    };

    /**
     * @brief The quadrature of a face of the particle geometry in the intermediate configuration, for the loads on it.
     * @details The rule of the section (ParticleSectionIntegrals::faceQuadrature()) on the face of the particle
     * geometry; for the update type DeformationGradient (SQCNI) the geometry is the smoothing domain, so that a face
     * load is integrated with the rule of the smoothing.
     * @param[in] boundaryFaceID The ID of the boundary face.
     * @return The points and their share in the area vector of the face.
     */
    typename SectionIntegrals::FaceQuadrature faceQuadratureIntermediate( int boundaryFaceID ) const
    {
      return SectionIntegrals::faceQuadrature( _particleDomain.geometryDeformed(), boundaryFaceID );
    }

    /**
     * @brief Computes the distributed load and its derivative with respect to nodal displacements.
     *
     * This method handles different types of distributed loads, such as pressure,
     * and calculates the external force vector and its tangent matrix. On each point @f$ \boldsymbol{Y}_N @f$ of the
     * face quadrature (faceQuadratureIntermediate()), with its area vector @f$ \boldsymbol{N}\,dA_Y @f$, both types
     * assemble a face force @f$ \boldsymbol{f} @f$ as @f$ r_{Aj} \mathrel{-}= T_A(\boldsymbol{Y}_N)\,f_j @f$, together
     * with @f$ \partial f_j/\partial \Delta F_{kL}\,\partial N_B/\partial Y_L @f$ in the tangent:
     * - Pressure: @f$ \boldsymbol{f} = -\Delta J\,\Delta\boldsymbol{F}^{-\mathsf T}\,p\,\boldsymbol{N}\,dA_Y @f$
     *   (Nanson's formula from @f$ \boldsymbol{Y} @f$ to @f$ \boldsymbol{x} @f$), with @f$ p @f$ = load[0]; as for
     *   the elements, a positive pressure pushes on the face (against the outward normal); see
     *   ParticleSectionIntegrals::addPressure().
     * - CWFCorrection: @f$ \boldsymbol{f} = \boldsymbol{\tau}\,\Delta\boldsymbol{F}^{-\mathsf T}\,\boldsymbol{N}\,
     *   dA_Y / J_Y @f$, @f$ J_Y = \det\boldsymbol{F}_n @f$, i.e. the Cauchy traction
     *   @f$ \boldsymbol{\sigma}\,\boldsymbol{n}\,da @f$ of the particle's own stress on the face; its tangent
     *   contains the material part (from @f$ \partial\boldsymbol{\tau}/\partial\Delta\boldsymbol{F} @f$) and the
     *   geometric part (from @f$ \Delta\boldsymbol{F}^{-\mathsf T} @f$). load[0] (or no load) selects the
     *   corrected components as a bit mask: 1 = x, 2 = y, 4 = z, 0 = all.
     *
     * @param[in] type The type of distributed load (e.g., DisplacementParticle::Pressure).
     * @param[in] surfaceID The ID of the boundary face where the load is applied.
     * @param[in] load Pointer to the load value(s).
     * @param[in,out] fExt Pointer to the array where the external force vector will be accumulated.
     * @param[in,out] dExt_dQ Pointer to the array where the derivative of the external force with respect to nodal
     * displacements will be accumulated.
     * @param[in] timeNew The current time.
     * @param[in] dT The time increment.
     * For an axisymmetric particle, the face carries two Gauss points, the rule of the smoothing, and the pressure
     * force of each is that of the half plane times the ring of the point after the increment
     * (Hoop::ringOfFacePoint()), with its derivative in the tangent; CWFCorrection is not available.
     *
     * @throws std::invalid_argument if an invalid DistributedLoad type is specified, and for CWFCorrection of an
     * axisymmetric particle.
     */
    void computeDistributedLoad( int           type,
                                 int           surfaceID,
                                 const double* load,
                                 double*       fExt,
                                 double*       dExt_dQ,
                                 double        timeNew,
                                 double        dT ) const override;

    /**
     * @brief Distributed load for explicit time integration: the Pressure force of computeDistributedLoad(), without
     * the tangent.
     * @param[in] type The type of distributed load (only DisplacementParticle::Pressure).
     * @param[in] boundaryFaceID The ID of the boundary face where the load is applied.
     * @param[in] load Pointer to the load value (the pressure).
     * @param[in,out] fExt Pointer to the array where the external force vector will be accumulated.
     * @param[in] timeNew The current time.
     * @param[in] dT The time increment.
     * @throws std::invalid_argument if an invalid DistributedLoad type is specified.
     */
    void computeDistributedLoadExplicit( int           type,
                                         int           boundaryFaceID,
                                         const double* load,
                                         double*       fExt,
                                         double        timeNew,
                                         double        dT ) const override;
  };

  template < int nDim, int nVertices, Section section >
  DisplacementParticleSQCNI< nDim, nVertices, section >::DisplacementParticleSQCNI(
    int                  elementID,
    const double*        vertexCoordinates,
    int                  nVertexCoordinates,
    double               volume, // This volume must be 0, as it's computed from vertices
    const std::string&   materialName,
    const double*        materialProperties,
    int                  sizeMaterialProperties,
    const Approximation& approximation,
    const typename ParticleDomainType::SmoothingDomainUpdateType smoothingVolumeUpdateType )
    : DisplacementParticle< nDim, section >( elementID,
                                             ParticleDomainType::getCenterFromVertices(
                                               Eigen::Map< const Eigen::Matrix< double, nDim, nVertices > >(
                                                 vertexCoordinates ) )
                                               .data(),
                                             ParticleDomainType::CoordinatesSized::RowsAtCompileTime,
                                             SectionIntegrals::volume( typename ParticleDomainType::LagrangeCellType(
                                               Eigen::Map< const Eigen::Matrix< double, nDim, nVertices > >(
                                                 vertexCoordinates ) ) ),
                                             materialName,
                                             materialProperties,
                                             sizeMaterialProperties,
                                             approximation ),
      _particleDomain( vertexCoordinates, nVertexCoordinates, smoothingVolumeUpdateType )
  {
    SectionIntegrals::checkParticle( smoothingDomain(),
                                     smoothingVolumeUpdateType ==
                                       ParticleDomainType::SmoothingDomainUpdateType::DeformationGradient,
                                     elementID );
    if ( volume != 0 ) {
      throw std::
        invalid_argument( Marmot::Utilities::MakeString()
                          << __PRETTY_FUNCTION__
                          << ": volume argument must be zero for DisplacementParticleSQCNI, as volume is computed from "
                             "vertex coordinates." );
    }
  }

  template < int nDim, int nVertices, Section section >
  void DisplacementParticleSQCNI< nDim, nVertices, section >::computeDistributedLoad( int           type,
                                                                                      int           boundaryFaceID,
                                                                                      const double* load_,
                                                                                      double*       fExt,
                                                                                      double*       dFExt_ddQ,
                                                                                      double        timeNew,
                                                                                      double        dT ) const
  {

    SectionIntegrals::checkFaceLoadIsAvailable( type == DisplacementParticle< nDim, section >::CWFCorrection,
                                                __PRETTY_FUNCTION__ );

    switch ( type ) {

    case DisplacementParticle< nDim, section >::Pressure: {
      SectionIntegrals::template addPressure< nDim >( faceQuadratureIntermediate( boundaryFaceID ),
                                                      load_[0],
                                                      *this->_mp,
                                                      this->_N,
                                                      this->_dN_dY,
                                                      this->_meshfreeApproximation,
                                                      this->_assignedKernelFunctions,
                                                      fExt,
                                                      dFExt_ddQ );
      break;
    }
    case DisplacementParticle< nDim, section >::CWFCorrection: {
      const auto&   _mp           = this->_mp;
      const auto&   _nNodes       = this->_nNodes;
      constexpr int nodeBlockSize = nDim;

      Eigen::Map< Eigen::VectorXd > P( fExt, _nNodes * nodeBlockSize );
      Eigen::Map< Eigen::MatrixXd > K( dFExt_ddQ, _nNodes * nodeBlockSize, _nNodes * nodeBlockSize );

      using namespace Fastor;

      // one point with the weight 1, the Cartesian rule of the face: an axisymmetric particle has no CWF correction
      const auto                             face = faceQuadratureIntermediate( boundaryFaceID );
      const Eigen::Matrix< double, nDim, 1 > Y_N  = face.points.col( 0 );
      const TensorD                          N_dAY( face.n_dA.data() );

      const auto& S = _mp->response.S;
      const auto& t = _mp->tangents;

      const auto                         deltaF    = this->dx_dY();
      const Tensor< double, nDim, nDim > deltaFInv = inverse( deltaF );

      // the internal force integrates tau over the undeformed volume V0 = V_Y / J_Y, so the boundary term of the same
      // weak form is tau * deltaF^-T * N dA_Y / J_Y, with J_Y the Jacobian of the intermediate configuration
      const double  JY       = determinant( this->dY_dX() );
      const TensorD v        = transpose( deltaFInv ) % N_dAY / JY; // F^-T * N * dA / J_Y
      TensorD       traction = S % v;                               // tau * v

      Tensor< double, nDim, nDim, nDim > df_dDeltaF;
      df_dDeltaF.zeros();

      for ( int i = 0; i < nDim; ++i ) {
        for ( int m = 0; m < nDim; ++m ) {
          for ( int M = 0; M < nDim; ++M ) {

            double geo = 0.0;
            double mat = 0.0;

            for ( int j = 0; j < nDim; ++j ) {
              // Geometric stiffness part: derived from delta(F^-T)
              geo -= S( i, j ) * v( m ) * deltaFInv( M, j );

              // Material stiffness part: derived from delta(tau)
              mat += t.dS_dDeltaF( i, j, m, M ) * v( j );
            }

            df_dDeltaF( i, m, M ) = geo + mat;
          }
        }
      }

      // load[0] selects the corrected components as a bit mask (1 = x, 2 = y, 4 = z; 0 or no load = all): on a face
      // where only the normal displacement is prescribed (symmetry plane, frictionless platen) only that component of
      // the traction is a reaction, the tangential one is a natural (zero) condition
      const int mask = load_ ? static_cast< int >( std::lround( load_[0] ) ) : 0;
      TensorD   sel;
      for ( int i = 0; i < nDim; ++i )
        sel( i ) = ( mask == 0 || ( mask >> i ) & 1 ) ? 1.0 : 0.0;
      for ( int i = 0; i < nDim; ++i ) {
        traction( i ) *= sel( i );
        for ( int m = 0; m < nDim; ++m )
          for ( int M = 0; M < nDim; ++M )
            df_dDeltaF( i, m, M ) *= sel( i );
      }

      Eigen::RowVectorXd testBoundary = Eigen::RowVectorXd::Zero( this->_nNodes );

      this->_meshfreeApproximation.computeShapeFunctions( Y_N.data(),
                                                          this->_assignedKernelFunctions,
                                                          testBoundary.data() );

      for ( int A = 0; A < _nNodes; A++ ) {
        const int idxA_u = nodeBlockSize * A;

        TensorD r_U = testBoundary( A ) * traction;

        {
          using namespace Eigen;
          P.template segment< nDim >( idxA_u ) -= Map< Matrix< double, nDim, 1 > >( r_U.data() );
        }

        for ( int B = 0; B < _nNodes; B++ ) {
          const int  idxB_u  = nodeBlockSize * B;
          const auto dN_B_dY = TensorMap< const double, nDim >( this->_dN_dY.col( B ).data() );

          Tensor< double, nDim, nDim > df_ddQU_B;
          df_ddQU_B.zeros();

          // Contract df_dF with the reference shape function gradient
          for ( int i = 0; i < nDim; ++i ) {
            for ( int m = 0; m < nDim; ++m ) {
              for ( int M = 0; M < nDim; ++M ) {
                df_ddQU_B( i, m ) += df_dDeltaF( i, m, M ) * dN_B_dY( M );
              }
            }
          }

          df_ddQU_B *= testBoundary( A );

          {
            using namespace Eigen;
            K.template block< nDim, nDim >( idxA_u, idxB_u ) -= Map< Matrix< double, nDim, nDim > >(
              torowmajor( df_ddQU_B ).data() );
          }
        }
      }

      break;
    }
    default: {
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__ << ": invalid DistributedLoad type specified" );
    }
    }
  }

  template < int nDim, int nVertices, Section section >
  void DisplacementParticleSQCNI< nDim, nVertices, section >::computeDistributedLoadExplicit( int type,
                                                                                              int boundaryFaceID,
                                                                                              const double* load,
                                                                                              double*       fExt,
                                                                                              double        timeNew,
                                                                                              double        dT ) const
  {
    switch ( type ) {

    case DisplacementParticle< nDim, section >::Pressure: {
      SectionIntegrals::template addPressure< nDim >( faceQuadratureIntermediate( boundaryFaceID ),
                                                      load[0],
                                                      *this->_mp,
                                                      this->_N,
                                                      this->_dN_dY,
                                                      this->_meshfreeApproximation,
                                                      this->_assignedKernelFunctions,
                                                      fExt,
                                                      nullptr );
      break;
    }
    default: {
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__ << ": invalid DistributedLoad type specified" );
    }
    }
  }

} // namespace Marmot::Discretizations::Meshfree
