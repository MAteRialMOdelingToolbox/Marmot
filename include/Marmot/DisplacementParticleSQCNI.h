/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck,
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
 *
 * Matthias Neuner matthias.neuner@uibk.ac.at
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
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/src/Core/util/Constants.h>
#include <Fastor/Fastor.h>

namespace Marmot::Meshfree {

  /**
   * @brief A displacement particle implementing the Stabilized Quasi-Conforming Nodal Integration (SQCNI) method.
   *
   * This class extends the basic DisplacementParticle to incorporate the SQCNI formulation,
   * which involves smoothing domains for integration and specific handling of deformation.
   * It uses a ParticleDomain object for its geometric representation via composition.
   *
   * @tparam nDim The number of dimensions (e.g., 2 for 2D, 3 for 3D).
   * @tparam nVertices The number of vertices defining the particle's geometry.
   */
  template < int nDim, int nVertices >
  class DisplacementParticleSQCNI : public DisplacementParticle< nDim > {

    using TensorD  = Fastor::Tensor< double, nDim >;
    using TensorDD = Fastor::Tensor< double, nDim, nDim >;

  protected:
    using ParentPointParticle = DisplacementParticle< nDim >;
    using ParticleDomainType  = ParticleDomain< nDim, nVertices >;

    ParticleDomainType _particleDomain;

  public:
    using SmoothingDomainUpdateType = ParticleDomainType::SmoothingDomainUpdateType;

    /**
     * @brief Constructs a new DisplacementParticleSQCNI object.
     *
     * @param elementID The unique identifier for the element.
     * @param nodeCoordinates Pointer to an array of node coordinates (nDim * nVertices).
     * @param nNodeCoordiantes The total number of coordinate values (nDim * nVertices).
     * @param volume The volume of the particle. Must be 0, as the volume is computed from vertex coordinates.
     * @param materialName The name of the material assigned to the particle.
     * @param materialProperties Pointer to an array of material properties.
     * @param sizeMaterialProperties The size of the material properties array.
     * @param approximation The meshfree approximation object used for shape functions.
     * @param smoothingVolumeUpdateType The strategy for updating the smoothing domain's volume.
     * @throws std::invalid_argument if the provided volume is not zero.
     */
    DisplacementParticleSQCNI( int                                                          elementID,
                               const double*                                                nodeCoordinates,
                               int                                                          nNodeCoordiantes,
                               double                                                       volume,
                               const std::string&                                           materialName,
                               const double*                                                materialProperties,
                               int                                                          sizeMaterialProperties,
                               const MarmotMeshfreeApproximation&                           approximation,
                               const typename ParticleDomainType::SmoothingDomainUpdateType smoothingVolumeUpdateType );

    // Override MarmotParticle interface methods to delegate to _particleDomain
    virtual void getVertexCoordinates( double* coordinates ) const override
    {
      Eigen::Map< Eigen::Matrix< double, nDim, nVertices > > coordinatesMap( coordinates );
      coordinatesMap = _particleDomain.getGeometryDeformedVertexCoordinates();
    }

    virtual void getVisualizationVertexCoordinates( double* coordinates ) const override
    {
      Eigen::Map< Eigen::Matrix< double, nDim, nVertices > > coordinatesMap( coordinates );
      coordinatesMap = _particleDomain.getSmoothingVertexCoordinates();
    }

    virtual int getNumberOfVertices() const override { return _particleDomain.getNumberOfVertices(); }

    virtual std::string getParticleShape() const override { return _particleDomain.getParticleShape(); }

    virtual void getFaceCoordinates( int faceID, double* coordinates ) const override
    {
      Eigen::Map< Eigen::Matrix< double, nDim, 1 > > coordinatesMap( coordinates );
      coordinatesMap = _particleDomain.getFaceCenterCoordinates( faceID );
    }

    virtual void getEvaluationCoordinates( double* coordinates ) const override
    {
      Eigen::Map< Eigen::Matrix< double, nDim, Eigen::Dynamic > > coordinatesMap( coordinates,
                                                                                  nDim,
                                                                                  _particleDomain.getNumberOfFaces() );

      for ( int i = 0; i < _particleDomain.getNumberOfFaces(); i++ ) {
        coordinatesMap.col( i ) = _particleDomain.getSmoothingDomainFaceCenterCoordinates( i + 1 );
      }
    }

    virtual int getNumberOfEvaluationPoints() const override { return _particleDomain.getNumberOfFaces(); }

    virtual void getCenterCoordinates( double* coordinates ) const override
    {
      // Use the material point's center coordinates, which is physics-specific
      this->_mp->getCoordinatesAtCenter( coordinates );
    }

    virtual double getVolumeUndeformed() const override { return this->_mp->getVolumeUndeformed(); };

    virtual StateView getStateView( const std::string& stateName, int qp ) const override
    {
      if ( stateName == "vertex displacements" )
        return StateView( const_cast< double* >( _particleDomain.getGeometryDeformedVertexDisplacements().data() ),
                          nDim * nVertices );

      if ( stateName == "smoothing vertex displacements" )
        return StateView( const_cast< double* >( _particleDomain.getSmoothingDomainVertexDisplacements().data() ),
                          nDim * nVertices );

      return ParentPointParticle::getStateView( stateName, qp );
    }

    /**
     * @brief Computes the deformed volume of the particle.
     * @return The deformed volume, calculated as undeformed volume multiplied by the determinant of the deformation
     * gradient.
     */
    virtual double getVolumeDeformed() const
    {
      return this->getVolumeUndeformed() * Fastor::determinant( this->dY_dX() );
    }

    virtual void acceptStateAndPosition() override
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
     * @brief Assigns the meshfree kernel functions to the particle.
     *
     * This method also computes and updates the shape functions (N, dN_dY) and
     * their derivatives (T, dT_dY) based on the assigned kernel functions and
     * the current smoothing volume.
     *
     * @param kernelFunctions A vector of pointers to the meshfree kernel functions.
     */
    virtual void assignMeshfreeKernelFunctions(
      const std::vector< const MarmotMeshfreeKernelFunction* >& kernelFunctions ) override
    {
      this->_assignedKernelFunctions = kernelFunctions;
      this->_nNodes                  = this->_assignedKernelFunctions.size();

      Eigen::Matrix< double, nDim, 1 > coords;
      this->getCenterCoordinates( coords.data() ); // Uses DisplacementParticleSQCNI's getCenterCoordinates

      this->_N     = Eigen::MatrixXd::Zero( 1, this->_nNodes );
      this->_dN_dY = Eigen::MatrixXd::Zero( nDim, this->_nNodes );

      // Compute N at the particle center (from GenericParticle)
      this->_meshfreeApproximation.computeShapeFunctions( coords.data(),
                                                          this->_assignedKernelFunctions,
                                                          this->_N.data() );

      // Compute dN_dY using the SQCNI approach (boundary integral over smoothing domain)
      Eigen::MatrixXd smooth_NBoundary( 1, this->_nNodes );
      for ( int i = 0; i < _particleDomain.getNumberOfFaces(); i++ ) {

        auto smoothing_evaluation_point = _particleDomain.getSmoothingDomainFaceCenterCoordinates( i + 1 );
        auto smoothing_n_dA             = _particleDomain.getSmoothingBoundarySurfaceVector( i + 1 );

        this->_meshfreeApproximation.computeShapeFunctions( smoothing_evaluation_point.data(),
                                                            this->_assignedKernelFunctions,
                                                            smooth_NBoundary.data() );
        this->_dN_dY += smoothing_n_dA * smooth_NBoundary;
      }
      this->_dN_dY /= _particleDomain.getSmoothingVolume();

      this->_T     = this->_N;
      this->_dT_dY = this->_dN_dY;
    }

    /**
     * @brief Computes the boundary integral part for the VCI (Variational Consistent Integration) test function.
     *
     * This method contributes to the R_AiC_RowMajor matrix, which is part of the VCI formulation.
     * It involves shape functions, monomial basis, and boundary surface vectors.
     *
     * @param R_AiC_RowMajor Pointer to the row-major matrix for VCI constraints.
     * @param boundarySurfaceVector Pointer to the boundary surface vector.
     * @param boundaryFaceID The ID of the boundary face.
     */
    virtual void vci_compute_Test_P_BoundaryIntegral( double*       R_AiC_RowMajor,
                                                      const double* boundarySurfaceVector,
                                                      int           boundaryFaceID ) override
    {
      using namespace Fastor;

      const auto [N_dAY, Y_N]   = getBoundaryVectorIntermediate( boundaryFaceID );
      Eigen::MatrixXd TBoundary = Eigen::MatrixXd::Zero( 1, this->_nNodes );

      this->_meshfreeApproximation.computeShapeFunctions( Y_N.data(),
                                                          this->_assignedKernelFunctions,
                                                          TBoundary.data() );

      // get P for the exact integration location at the boundary
      auto PBoundary = this->_P;
      Math::computeMonomialBasis( this->_vciOrder, Eigen::Matrix< double, nDim, 1 >( Y_N.data() ), PBoundary );

      for ( int A = 0; A < this->_nNodes; A++ )
        for ( int i = 0; i < nDim; i++ )
          for ( int C = 0; C < this->_nVCIConstraints; C++ )
            R_AiC_RowMajor[A * ( nDim * this->_nVCIConstraints ) + i * this->_nVCIConstraints + C] += TBoundary( A ) *
                                                                                                      PBoundary( C ) *
                                                                                                      N_dAY[i];
    };

    /**
     * @brief Retrieves the boundary surface vector and the face center coordinates in the intermediate configuration.
     * @param boundaryFaceID The ID of the boundary face.
     * @return A tuple containing the boundary surface vector (N_dAY) and the face center coordinates (Y_N).
     */
    std::tuple< TensorD, TensorD > getBoundaryVectorIntermediate( int boundaryFaceID ) const;

    /**
     * @brief Computes the distributed load and its derivative with respect to nodal displacements.
     *
     * This method handles different types of distributed loads, such as pressure,
     * and calculates the external force vector and its tangent matrix.
     *
     * @param type The type of distributed load (e.g., DisplacementParticle::Pressure).
     * @param surfaceID The ID of the boundary face where the load is applied.
     * @param load Pointer to the load value(s).
     * @param fExt Pointer to the array where the external force vector will be accumulated.
     * @param dExt_dQ Pointer to the array where the derivative of the external force with respect to nodal
     * displacements will be accumulated.
     * @param timeNew The current time.
     * @param dT The time increment.
     * @throws std::invalid_argument if an invalid DistributedLoad type is specified.
     */
    virtual void computeDistributedLoad( int           type,
                                         int           surfaceID,
                                         const double* load,
                                         double*       fExt,
                                         double*       dExt_dQ,
                                         double        timeNew,
                                         double        dT ) const override;

    virtual void computeDistributedLoadExplicit( int           type,
                                                 int           boundaryFaceID,
                                                 const double* load,
                                                 double*       fExt,
                                                 double        timeNew,
                                                 double        dT ) const override;
  };

  template < int nDim, int nVertices >
  DisplacementParticleSQCNI< nDim, nVertices >::DisplacementParticleSQCNI(
    int                                elementID,
    const double*                      vertexCoordinates,
    int                                nVertexCoordinates,
    double                             volume, // This volume must be 0, as it's computed from vertices
    const std::string&                 materialName,
    const double*                      materialProperties,
    int                                sizeMaterialProperties,
    const MarmotMeshfreeApproximation& approximation,
    const typename ParticleDomainType::SmoothingDomainUpdateType smoothingVolumeUpdateType )
    : DisplacementParticle< nDim >( elementID,
                                    ParticleDomainType::getCenterFromVertices(
                                      Eigen::Map< const Eigen::Matrix< double, nDim, nVertices > >(
                                        vertexCoordinates ) )
                                      .data(),
                                    ParticleDomainType::CoordinatesSized::RowsAtCompileTime,
                                    ParticleDomainType::getVolumeFromVertices(
                                      Eigen::Map< const Eigen::Matrix< double, nDim, nVertices > >(
                                        vertexCoordinates ) ),
                                    materialName,
                                    materialProperties,
                                    sizeMaterialProperties,
                                    approximation ),
      _particleDomain( vertexCoordinates, nVertexCoordinates, smoothingVolumeUpdateType )
  {
    if ( volume != 0 ) {
      throw std::invalid_argument(
        MakeString() << __PRETTY_FUNCTION__
                     << ": volume argument must be zero for DisplacementParticleSQCNI, as volume is computed from "
                        "vertex coordinates." );
    }
  }

  template < int nDim, int nVertices >
  void DisplacementParticleSQCNI< nDim, nVertices >::computeDistributedLoad( int           type,
                                                                             int           boundaryFaceID,
                                                                             const double* load_,
                                                                             double*       fExt,
                                                                             double*       dFExt_ddQ,
                                                                             double        timeNew,
                                                                             double        dT ) const
  {

    switch ( type ) {

    case DisplacementParticle< nDim >::Pressure: {

      const auto&   _nNodes       = this->_nNodes; // From GenericParticle
      constexpr int nodeBlockSize = nDim;

      const auto [N_dAY, Y_N] = getBoundaryVectorIntermediate( boundaryFaceID );

      TensorD fY = N_dAY * load_[0];

      Eigen::Map< Eigen::VectorXd > P( fExt, _nNodes * nodeBlockSize );
      Eigen::Map< Eigen::MatrixXd > K( dFExt_ddQ, _nNodes * nodeBlockSize, _nNodes * nodeBlockSize );

      using namespace Fastor;
      using namespace FastorIndices;

      Tensor< double, nDim, nDim > Eye;
      Eye.eye();

      // apply Nanson's formula
      const auto deltaF = this->dx_dY(); // From DisplacementParticle

      const Tensor< double, nDim, nDim > deltaFInv = inverse( deltaF );
      const double                       deltaJ    = determinant( deltaF );

      const TensorD f = deltaJ * transpose( deltaFInv ) % fY;

      const Tensor< double, nDim, nDim, nDim, nDim > dFInv_dF = -einsum< Ik, Ki, to_IikK >( deltaFInv, deltaFInv );

      const Tensor< double, nDim, nDim, nDim > df_dDeltaF = outer( f, transpose( deltaFInv ) ) +
                                                            deltaJ * einsum< IikK, Index< I_ > >( dFInv_dF, fY );

      TensorD r_U( 0.0 );

      Eigen::MatrixXd testBoundary = Eigen::MatrixXd::Zero( 1, this->_nNodes ); // From GenericParticle

      this->_meshfreeApproximation.computeShapeFunctions( Y_N.data(),
                                                          this->_assignedKernelFunctions,
                                                          testBoundary.data() );

      for ( int A = 0; A < _nNodes; A++ ) {
        const int idxA_u = nodeBlockSize * A;

        r_U = testBoundary( A ) * f;

        {
          using namespace Eigen;
          P.template segment< nDim >( idxA_u ) -= Map< Matrix< double, nDim, 1 > >( r_U.data() );
        }

        for ( int B = 0; B < _nNodes; B++ ) {
          const int  idxB_u  = nodeBlockSize * B;
          const auto dN_B_dY = TensorMap< const double, nDim >( this->_dN_dY.col( B ).data() ); // From GenericParticle

          const Tensor< double, nDim, nDim > df_ddQU_B = testBoundary( A ) * einsum< ijk, k >( df_dDeltaF, dN_B_dY );

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
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid DistributedLoad type specified" );
    }
    }
  }

  template < int nDim, int nVertices >
  void DisplacementParticleSQCNI< nDim, nVertices >::computeDistributedLoadExplicit( int           type,
                                                                                     int           boundaryFaceID,
                                                                                     const double* load,
                                                                                     double*       fExt,
                                                                                     double        timeNew,
                                                                                     double        dT ) const
  {
    switch ( type ) {

    case DisplacementParticle< nDim >::Pressure: {

      const auto&   _nNodes       = this->_nNodes; // From GenericParticle
      constexpr int nodeBlockSize = nDim;

      const auto [N_dAY, Y_N] = getBoundaryVectorIntermediate( boundaryFaceID );

      TensorD fY = N_dAY * load[0];

      Eigen::Map< Eigen::VectorXd > P( fExt, _nNodes * nodeBlockSize );

      using namespace Fastor;
      using namespace FastorIndices;

      Tensor< double, nDim, nDim > Eye;
      Eye.eye();

      // apply Nanson's formula
      const auto deltaF = this->dx_dY(); // From DisplacementParticle

      const Tensor< double, nDim, nDim > deltaFInv = inverse( deltaF );
      const double                       deltaJ    = determinant( deltaF );

      const TensorD f = deltaJ * transpose( deltaFInv ) % fY;

      TensorD r_U( 0.0 );

      Eigen::MatrixXd testBoundary = Eigen::MatrixXd::Zero( 1, this->_nNodes ); // From GenericParticle

      this->_meshfreeApproximation.computeShapeFunctions( Y_N.data(),
                                                          this->_assignedKernelFunctions,
                                                          testBoundary.data() );

      for ( int A = 0; A < _nNodes; A++ ) {
        const int idxA_u = nodeBlockSize * A;

        r_U = testBoundary( A ) * f;

        {
          using namespace Eigen;
          P.template segment< nDim >( idxA_u ) -= Map< Matrix< double, nDim, 1 > >( r_U.data() );
        }
      }

      break;
    }
    default: {
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid DistributedLoad type specified" );
    }
    }
  }

  template < int nDim, int nVertices >
  std::tuple< typename DisplacementParticleSQCNI< nDim, nVertices >::TensorD,
              typename DisplacementParticleSQCNI< nDim, nVertices >::TensorD >
  DisplacementParticleSQCNI< nDim, nVertices >::getBoundaryVectorIntermediate( int boundaryFaceID ) const
  {

    TensorD N_dAY;
    TensorD Y;

    Eigen::Matrix< double, nDim, 1 > _Y_eigen;

    // the evaluation point depends: For real SQCNI, we do it on the smoothing domain boundary, for all others we do in
    // in the center of the deformed geometry

    if ( _particleDomain.smoothingVolumeUpdateType ==
         ParticleDomainType::SmoothingDomainUpdateType::DeformationGradient )
      _Y_eigen = _particleDomain.getSmoothingDomainFaceCenterCoordinates( boundaryFaceID );
    else
      _Y_eigen = _particleDomain.getFaceCenterCoordinates( boundaryFaceID );

    // N_dAY (boundary surface vector for distributed load) comes from the deformed geometry
    auto _N_dAY_eigen = _particleDomain.getFaceBoundaryVector( boundaryFaceID );

    for ( int i = 0; i < nDim; i++ ) {
      N_dAY[i] = _N_dAY_eigen[i];
      Y[i]     = _Y_eigen[i];
    }

    return { N_dAY, Y };
  }

} // namespace Marmot::Meshfree
