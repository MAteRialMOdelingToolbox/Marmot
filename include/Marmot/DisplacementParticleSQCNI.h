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
#include "Marmot/MarmotMeshfreeQuadHexCell.h"
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
   * It supports different strategies for updating the smoothing domain's volume.
   *
   * @tparam nDim The number of dimensions (e.g., 2 for 2D, 3 for 3D).
   * @tparam nVertices The number of vertices defining the particle's geometry.
   */
  template < int nDim, int nVertices >
  class DisplacementParticleSQCNI : public DisplacementParticle< nDim > {

    using TensorD          = Fastor::Tensor< double, nDim >;
    using TensorDD         = Fastor::Tensor< double, nDim, nDim >;
    using CoordinatesSized = Eigen::Matrix< double, nDim, 1 >;

  public:
    /**
     * @brief Defines how the smoothing domain's volume is updated.
     */
    enum SmoothingDomainUpdateType {
      None,                       ///< No update to the smoothing domain's deformation tensor (identity).
      DeformationGradient,        ///< Update using the total deformation gradient F.
      RotationOnly,               ///< Update using only the rotation part R from F = RU.
      RotationAndPrincipalStretch ///< Update using the rotation R and principal stretches from U.
    };

  protected:
    MarmotLagrangeCell< nDim, nVertices > _cellForGeometryUndeformed; ///< Cell representing the undeformed geometry.
    MarmotLagrangeCell< nDim, nVertices >
      _cellForGeometryIntermediate;                                   ///< Cell representing the intermediate geometry.
    MarmotLagrangeCell< nDim, nVertices > _cellForSmoothing;          ///< Cell representing the smoothing domain.

    const SmoothingDomainUpdateType _smoothingVolumeUpdateType;       ///< Type of update for the smoothing volume.

    const Eigen::Matrix< double, nDim, nVertices > _vertexCoordinates_Undeformed; ///< Undeformed vertex coordinates.
    Eigen::Matrix< double, nDim, nVertices >
      _vertex_displacements_smoothingDomain; ///< Displacements of vertices in the smoothing domain.

    using ParentPointParticle = DisplacementParticle< nDim >;

  public:
    /**
     * @brief Retrieves the current coordinates of the particle's vertices in the intermediate configuration.
     * @param coordinates Pointer to a double array where the vertex coordinates will be stored.
     */
    virtual void getVertexCoordinates( double* coordinates ) const override;

    /**
     * @brief Retrieves the current coordinates of the particle's vertices for visualization purposes.
     *        Currently, this is the same as getVertexCoordinates.
     * @param coordinates Pointer to a double array where the vertex coordinates will be stored.
     */
    virtual void getVisualizationVertexCoordinates( double* coordinates ) const override
    {
      getVertexCoordinates( coordinates );
    };

    /**
     * @brief Returns the number of vertices defining the particle.
     * @return The number of vertices.
     */
    virtual int getNumberOfVertices() const override { return nVertices; };

    /**
     * @brief Returns a string describing the shape of the particle's underlying cell.
     * @return A string representing the cell shape (e.g., "quad", "hex").
     */
    virtual std::string getParticleShape() const override { return _cellForGeometryUndeformed.getCellShape(); }

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
    DisplacementParticleSQCNI( int                                elementID,
                               const double*                      nodeCoordinates,
                               int                                nNodeCoordiantes,
                               double                             volume,
                               const std::string&                 materialName,
                               const double*                      materialProperties,
                               int                                sizeMaterialProperties,
                               const MarmotMeshfreeApproximation& approximation,
                               const SmoothingDomainUpdateType    smoothingVolumeUpdateType );

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
      const std::vector< const MarmotMeshfreeKernelFunction* >& kernelFunctions ) override;

    /**
     * @brief Get the smoothing volume of the particle.
     * @return The smoothing volume of the particle.
     * @details The smoothing volume is computed as the volume of the element in an updated configuration
     *          that is obtained by applying (parts of) the deformation gradient to the element in the undeformed
     *          configuration. In general, this is not consistent with the actual deformation of the particle.
     */
    virtual double getSmoothingVolume() const
    {
      const auto _smoothingVolume = _cellForSmoothing.volume();
      return _smoothingVolume;
    }

    /**
     * @brief Retrieves the coordinates of the particle's center in the current configuration.
     * @param coordinates Pointer to a double array where the center coordinates will be stored.
     */
    virtual void getCenterCoordinates( double* coordinates ) const override
    {
      this->_mp->getCoordinatesAtCenter( coordinates );
    }

    /**
     * @brief Computes the centroid coordinates from a given set of vertex coordinates.
     * @param vertexCoordinates An Eigen matrix containing the vertex coordinates.
     * @return An Eigen vector representing the centroid coordinates.
     */
    virtual CoordinatesSized getCenterFromVertices(
      const Eigen::Matrix< double, nDim, nVertices >& vertexCoordinates ) const
    {
      MarmotLagrangeCell< nDim, nVertices > _lagrangeCell( vertexCoordinates );
      return _lagrangeCell.centroid();
    }

    /**
     * @brief Computes the volume of a cell defined by a given set of vertex coordinates.
     * @param vertexCoordinates An Eigen matrix containing the vertex coordinates.
     * @return The computed volume.
     */
    virtual double getVolumeFromVertices( const Eigen::Matrix< double, nDim, nVertices >& vertexCoordinates ) const
    {
      MarmotLagrangeCell< nDim, nVertices > _lagrangeCell( vertexCoordinates );
      return _lagrangeCell.volume();
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

    /**
     * @brief Returns the number of required state variables for the particle.
     * @return The number of required state variables.
     */
    virtual int getNumberOfRequiredStateVars() const override
    {
      return DisplacementParticle< nDim >::getNumberOfRequiredStateVars();
    };

    /**
     * @brief Assigns state variables to the particle.
     * @param stateVars Pointer to an array of state variables.
     * @param nStateVars The number of state variables.
     */
    void assignStateVars( double* stateVars, int nStateVars ) override
    {
      DisplacementParticle< nDim >::assignStateVars( stateVars, nStateVars );
    }

    /**
     * @brief Provides a view into a specific state variable.
     * @param stateName The name of the state variable (e.g., "vertex displacements").
     * @param qp The quadrature point index (not used for "vertex displacements").
     * @return A StateView object providing access to the state variable data.
     */
    virtual StateView getStateView( const std::string& stateName, int qp ) const override
    {
      if ( stateName == "vertex displacements" ) {
        return StateView( const_cast< double* >( _vertex_displacements_smoothingDomain.data() ), nDim * nVertices );
      }
      return DisplacementParticle< nDim >::getStateView( stateName, qp );
    }

    /**
     * @brief Sets an initial condition for the particle.
     * @param conditionName The name of the initial condition.
     * @param value Pointer to the value(s) for the initial condition.
     */
    virtual void setInitialCondition( const std::string& conditionName, const double* value ) override
    {
      ParentPointParticle::setInitialCondition( conditionName, value );
    };

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

    /**
     * @brief Retrieves the boundary surface vector and the face center coordinates in the intermediate configuration.
     * @param boundaryFaceID The ID of the boundary face.
     * @return A tuple containing the boundary surface vector (N_dAY) and the face center coordinates (Y_N).
     */
    std::tuple< TensorD, TensorD > getBoundaryVectorIntermediate( int boundaryFaceID ) const;

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
                                                      int           boundaryFaceID )
    {
      using namespace Fastor;

      const auto [N_dAY, Y_N]   = getBoundaryVectorIntermediate( boundaryFaceID );
      Eigen::MatrixXd TBoundary = Eigen::MatrixXd::Zero( 1, ParentPointParticle::_nNodes );

      ParentPointParticle::_meshfreeApproximation.computeShapeFunctions( Y_N.data(),
                                                                         ParentPointParticle::_assignedKernelFunctions,
                                                                         TBoundary.data() );

      // get P for the exact integration location at the boundary
      auto PBoundary = ParentPointParticle::_P;
      Math::computeMonomialBasis( ParentPointParticle::_vciOrder, CoordinatesSized( Y_N.data() ), PBoundary );

      for ( int A = 0; A < ParentPointParticle::_nNodes; A++ )
        for ( int i = 0; i < nDim; i++ )
          for ( int C = 0; C < ParentPointParticle::_nVCIConstraints; C++ )
            R_AiC_RowMajor[A * ( nDim * ParentPointParticle::_nVCIConstraints ) +
                           i * ParentPointParticle::_nVCIConstraints + C] += TBoundary( A ) * PBoundary( C ) * N_dAY[i];
    };

    /**
     * @brief Retrieves the coordinates of the evaluation points (e.g., face centers) for the particle.
     * @param coordinates Pointer to a double array where the evaluation coordinates will be stored.
     */
    virtual void getEvaluationCoordinates( double* coordinates ) const
    {
      Eigen::Map< Eigen::Matrix< double, nDim, Eigen::Dynamic > > faceCenters( coordinates,
                                                                               nDim,
                                                                               getNumberOfEvaluationPoints() );

      for ( int i = 0; i < getNumberOfEvaluationPoints(); i++ ) {
        const auto faceCenterCoords = _cellForSmoothing.getFaceCenterCoordinates( i + 1 );
        faceCenters.col( i )        = faceCenterCoords;
      }
    }

    /**
     * @brief Retrieves the coordinates of the center of a specific face in the intermediate configuration.
     * @param faceID The ID of the face.
     * @param coordinates Pointer to a double array where the face center coordinates will be stored.
     */
    virtual void getFaceCoordinates( int faceID, double* coordinates ) const
    {
      const auto faceCenterCoords = _cellForGeometryIntermediate.getFaceCenterCoordinates( faceID );
      for ( int i = 0; i < nDim; i++ ) {
        coordinates[i] = faceCenterCoords( i );
      }
    }

    /**
     * @brief Returns the number of evaluation points for the particle.
     *        These are typically the face centers of the smoothing cell.
     * @return The number of evaluation points.
     */
    virtual int getNumberOfEvaluationPoints() const
    {
      // Technically, we also evaluate at the center for integratio, but
      // we assume that the center is captured also if at least one face is evaluated.
      return _cellForSmoothing.getNumberOfFaces();
    };

  private:
    /**
     * @brief Computes the total deformation tensor for the smoothing domain based on the configured update type.
     * @return An Eigen matrix representing the deformation tensor for the smoothing domain.
     */
    Eigen::Matrix< double, nDim, nDim > _computeSmoothingDomainDeformationTensorTotal();

    /**
     * @brief Updates the particle's position to the reference intermediate configuration.
     *
     * This method updates the vertex coordinates of the geometry and smoothing cells
     * based on the center displacement and the deformation gradient (or parts of it,
     * depending on the smoothing domain update type).
     */
    virtual void updateParticlePositionToReferenceIntermediate() override
    {
      const auto _centerDisplacement = Eigen::Matrix< double, nDim, 1 >( this->getDisplacementAtCenter().data() );

      const auto FIntermediate = Eigen::Matrix< double, nDim, nDim, Eigen::RowMajor >( this->dY_dX().data() );
      _cellForGeometryIntermediate.updateVertexCoordinates( _vertexCoordinates_Undeformed );
      _cellForGeometryIntermediate.applyDeformationGradient( FIntermediate );
      _cellForGeometryIntermediate.applyUniformDisplacement( _centerDisplacement );

      const auto FSmoothing = _computeSmoothingDomainDeformationTensorTotal();
      _cellForSmoothing.updateVertexCoordinates( _vertexCoordinates_Undeformed );
      _cellForSmoothing.applyUniformDisplacement( _centerDisplacement );
      _cellForSmoothing.applyDeformationGradient( FSmoothing );

      _vertex_displacements_smoothingDomain = _cellForSmoothing.nodes() - _vertexCoordinates_Undeformed;

      getCenterCoordinates( ParentPointParticle::_centerReferenceIntermediate.data() );
    };
  };

  template < int nDim, int nVertices >
  DisplacementParticleSQCNI< nDim, nVertices >::DisplacementParticleSQCNI(
    int                                                  elementID,
    const double*                                        vertexCoordinates,
    int                                                  nVertexCoordinates,
    double                                               volume,
    const std::string&                                   materialName,
    const double*                                        materialProperties,
    int                                                  sizeMaterialProperties,
    const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation,
    const SmoothingDomainUpdateType                      smoothingVolumeUpdateType )
    : DisplacementParticle< nDim >( elementID,
                                    getCenterFromVertices( Eigen::Map< const Eigen::Matrix< double, nDim, nVertices > >(
                                                             vertexCoordinates ) )
                                      .data(),
                                    CoordinatesSized::RowsAtCompileTime,
                                    getVolumeFromVertices( Eigen::Map< const Eigen::Matrix< double, nDim, nVertices > >(
                                      vertexCoordinates ) ),
                                    materialName,
                                    materialProperties,
                                    sizeMaterialProperties,
                                    approximation ),
      _cellForGeometryUndeformed( vertexCoordinates, nVertexCoordinates ),
      _cellForGeometryIntermediate( vertexCoordinates, nVertexCoordinates ),
      _cellForSmoothing( vertexCoordinates, nVertexCoordinates ),
      _smoothingVolumeUpdateType( smoothingVolumeUpdateType ),
      _vertexCoordinates_Undeformed( vertexCoordinates )
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

      const auto&   _nNodes       = DisplacementParticle< nDim >::_nNodes;
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
      const auto deltaF = this->dx_dY();

      const Tensor< double, nDim, nDim > deltaFInv = inverse( deltaF );
      const double                       deltaJ    = determinant( deltaF );

      const TensorD f = deltaJ * transpose( deltaFInv ) % fY;

      const Tensor< double, nDim, nDim, nDim, nDim > dFInv_dF = -einsum< Ik, Ki, to_IikK >( deltaFInv, deltaFInv );

      const Tensor< double, nDim, nDim, nDim > df_dDeltaF = outer( f, transpose( deltaFInv ) ) +
                                                            deltaJ * einsum< IikK, Index< I_ > >( dFInv_dF, fY );

      TensorD r_U( 0.0 );

      Eigen::MatrixXd testBoundary = Eigen::MatrixXd::Zero( 1, ParentPointParticle::_nNodes );

      ParentPointParticle::_meshfreeApproximation.computeShapeFunctions( Y_N.data(),
                                                                         ParentPointParticle::_assignedKernelFunctions,
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
          const auto dN_B_dY = TensorMap< const double, nDim >( DisplacementParticle< nDim >::_dN_dY.col( B ).data() );

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
  std::tuple< typename DisplacementParticleSQCNI< nDim, nVertices >::TensorD,
              typename DisplacementParticleSQCNI< nDim, nVertices >::TensorD >
  DisplacementParticleSQCNI< nDim, nVertices >::getBoundaryVectorIntermediate( int boundaryFaceID ) const
  {

    TensorD N_dAY;
    TensorD Y;

    auto _Y     = _cellForGeometryIntermediate.getFaceCenterCoordinates( boundaryFaceID );
    auto _N_dAY = _cellForGeometryIntermediate.boundarySurfaceVector( boundaryFaceID );

    for ( int i = 0; i < nDim; i++ ) {
      N_dAY[i] = _N_dAY[i];
      Y[i]     = _Y[i];
    }

    return { N_dAY, Y };
  }

  template < int nDim, int nVertices >
  void DisplacementParticleSQCNI< nDim, nVertices >::getVertexCoordinates( double* coordinates ) const
  {

    Eigen::Map< Eigen::Matrix< double, nDim, nVertices > > coordinatesDeformed( coordinates );

    coordinatesDeformed = _cellForGeometryIntermediate.nodes();
  }

  template < int nDim, int nVertices >
  Eigen::Matrix< double, nDim, nDim > DisplacementParticleSQCNI< nDim, nVertices >::
    _computeSmoothingDomainDeformationTensorTotal()
  {

    Eigen::Matrix< double, nDim, nDim > F;

    switch ( _smoothingVolumeUpdateType ) {

    case SmoothingDomainUpdateType::None: F.setIdentity(); break;

    case SmoothingDomainUpdateType::DeformationGradient: {

      const auto F_ = this->dY_dX();
      F             = Eigen::Map< Eigen::Matrix< double, nDim, nDim, Eigen::RowMajor > >( F_.data() );
      break;
    }

    case SmoothingDomainUpdateType::RotationOnly: {

      const auto F_ = this->dY_dX();
      F             = Eigen::Map< Eigen::Matrix< double, nDim, nDim, Eigen::RowMajor > >( F_.data() );

      Eigen::JacobiSVD< Eigen::MatrixXd > svd;
      svd.compute( F, Eigen::ComputeFullU | Eigen::ComputeFullV );

      Eigen::Matrix< double, nDim, nDim > R = svd.matrixU() * svd.matrixV().transpose();

      F = R;
      break;
    }

    case SmoothingDomainUpdateType::RotationAndPrincipalStretch: {
      const auto F_ = this->dY_dX();
      F             = Eigen::Map< Eigen::Matrix< double, nDim, nDim, Eigen::RowMajor > >( F_.data() );

      Eigen::JacobiSVD< Eigen::MatrixXd > svd;
      svd.compute( F, Eigen::ComputeFullU | Eigen::ComputeFullV );

      Eigen::Matrix< double, nDim, nDim > R = svd.matrixU() * svd.matrixV().transpose();

      Eigen::Matrix< double, nDim, nDim > U_ = Eigen::Matrix< double, nDim, nDim >::Identity();
      U_.diagonal()                          = ( R.transpose() * F ).diagonal();

      F = R * U_;

      break;
    }
    }

    return F;
  }

  template < int nDim, int nVertices >
  void DisplacementParticleSQCNI< nDim, nVertices >::assignMeshfreeKernelFunctions(
    const std::vector< const MarmotMeshfreeKernelFunction* >& kernelFunctions )
  {

    ParentPointParticle::_assignedKernelFunctions = kernelFunctions;

    ParentPointParticle::_nNodes = DisplacementParticle< nDim >::_assignedKernelFunctions.size();

    Eigen::Matrix< double, nDim, 1 > coords;
    this->getCenterCoordinates( coords.data() );

    ParentPointParticle::_N     = Eigen::MatrixXd::Zero( 1, ParentPointParticle::_nNodes );
    ParentPointParticle::_dN_dY = Eigen::MatrixXd::Zero( nDim, ParentPointParticle::_nNodes );

    ParentPointParticle::_meshfreeApproximation.computeShapeFunctions( coords.data(),
                                                                       ParentPointParticle::_assignedKernelFunctions,
                                                                       ParentPointParticle::_N.data() );

    Eigen::MatrixXd NBoundary( 1, ParentPointParticle::_nNodes );
    for ( int i = 0; i < _cellForSmoothing.getNumberOfFaces(); i++ ) {

      auto faceCenterCoords = _cellForSmoothing.getFaceCenterCoordinates( i + 1 );
      auto n                = _cellForSmoothing.boundarySurfaceVector( i + 1 );

      ParentPointParticle::_meshfreeApproximation.computeShapeFunctions( faceCenterCoords.data(),
                                                                         ParentPointParticle::_assignedKernelFunctions,
                                                                         NBoundary.data() );
      ParentPointParticle::_dN_dY += n * NBoundary;
    }

    ParentPointParticle::_dN_dY /= getSmoothingVolume();

    ParentPointParticle::_T     = ParentPointParticle::_N;
    ParentPointParticle::_dT_dY = ParentPointParticle::_dN_dY;
  }

} // namespace Marmot::Meshfree
