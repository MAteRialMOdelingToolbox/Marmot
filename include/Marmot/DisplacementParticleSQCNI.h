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

  template < int nDim, int nVertices >
  class DisplacementParticleSQCNI : public DisplacementParticle< nDim > {

    using TensorD          = Fastor::Tensor< double, nDim >;
    using TensorDD         = Fastor::Tensor< double, nDim, nDim >;
    using CoordinatesSized = Eigen::Matrix< double, nDim, 1 >;

  public:
    enum SmoothingDomainUpdateType { None, DeformationGradient, RotationOnly, RotationAndPrincipalStretch };

  protected:
    MarmotLagrangeCell< nDim, nVertices > _cellForGeometryUndeformed;
    MarmotLagrangeCell< nDim, nVertices > _cellForGeometryIntermediate;
    MarmotLagrangeCell< nDim, nVertices > _cellForSmoothing;

    const SmoothingDomainUpdateType _smoothingVolumeUpdateType;

    const Eigen::Matrix< double, nDim, nVertices > _vertexCoordinates_Undeformed;
    Eigen::Matrix< double, nDim, nVertices >       _vertex_displacements_smoothingDomain;

    using ParentPointParticle = DisplacementParticle< nDim >;

  public:
    virtual void getVertexCoordinates( double* coordinates ) const override;

    virtual void getVisualizationVertexCoordinates( double* coordinates ) const override
    {
      getVertexCoordinates( coordinates );
    };

    virtual int getNumberOfVertices() const override { return nVertices; };

    virtual std::string getParticleShape() const override { return _cellForGeometryUndeformed.getCellShape(); }

    DisplacementParticleSQCNI( int                                elementID,
                               const double*                      nodeCoordinates,
                               int                                nNodeCoordiantes,
                               double                             volume,
                               const std::string&                 materialName,
                               const double*                      materialProperties,
                               int                                sizeMaterialProperties,
                               const MarmotMeshfreeApproximation& approximation,
                               const SmoothingDomainUpdateType    smoothingVolumeUpdateType );

    virtual void assignMeshfreeKernelFunctions(
      const std::vector< const MarmotMeshfreeKernelFunction* >& kernelFunctions ) override;

    /// \brief Get the smoothing volume of the particle
    /// \return The smoothing volume of the particle
    /// \details The smoothing volume is computed as the volume of the element in an updated configuration
    ///         that is obtained by applying (parts of) the deformation gradient to the element in the undeformed
    ///         configuration. In general, this is not consistent with the actual deformation of the particle.
    ///
    virtual double getSmoothingVolume() const
    {
      const auto _smoothingVolume = _cellForSmoothing.volume();
      return _smoothingVolume;
    }

    virtual void getCenterCoordinates( double* coordinates ) const override
    {
      this->_mp->getCoordinatesAtCenter( coordinates );
    }

    virtual CoordinatesSized getCenterFromVertices(
      const Eigen::Matrix< double, nDim, nVertices >& vertexCoordinates ) const
    {
      MarmotLagrangeCell< nDim, nVertices > _lagrangeCell( vertexCoordinates );
      return _lagrangeCell.centroid();
    }

    virtual double getVolumeFromVertices( const Eigen::Matrix< double, nDim, nVertices >& vertexCoordinates ) const
    {
      MarmotLagrangeCell< nDim, nVertices > _lagrangeCell( vertexCoordinates );
      return _lagrangeCell.volume();
    }

    virtual double getVolumeDeformed() const
    {
      return this->getVolumeUndeformed() * Fastor::determinant( this->dY_dX() );
    }

    virtual int getNumberOfRequiredStateVars() const override
    {
      return DisplacementParticle< nDim >::getNumberOfRequiredStateVars();
    };

    void assignStateVars( double* stateVars, int nStateVars ) override
    {
      DisplacementParticle< nDim >::assignStateVars( stateVars, nStateVars );
    }

    virtual StateView getStateView( const std::string& stateName, int qp ) const override
    {
      if ( stateName == "vertex displacements" ) {
        return StateView( const_cast< double* >( _vertex_displacements_smoothingDomain.data() ), nDim * nVertices );
      }
      return DisplacementParticle< nDim >::getStateView( stateName, qp );
    }

    virtual void setInitialCondition( const std::string& conditionName, const double* value ) override
    {
      ParentPointParticle::setInitialCondition( conditionName, value );
    };

    virtual void computeDistributedLoad( int           type,
                                         int           surfaceID,
                                         const double* load,
                                         double*       fExt,
                                         double*       dExt_dQ,
                                         double        timeNew,
                                         double        dT ) const override;

    std::tuple< TensorD, TensorD > getBoundaryVectorIntermediate( int boundaryFaceID ) const;

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

    virtual void getFaceCoordinates( int faceID, double* coordinates ) const
    {
      const auto faceCenterCoords = _cellForGeometryIntermediate.getFaceCenterCoordinates( faceID );
      for ( int i = 0; i < nDim; i++ ) {
        coordinates[i] = faceCenterCoords( i );
      }
    }

    virtual int getNumberOfEvaluationPoints() const
    {
      // Technically, we also evaluate at the center for integratio, but
      // we assume that the center is captured also if at least one face is evaluated.
      return _cellForSmoothing.getNumberOfFaces();
    };

  private:
    Eigen::Matrix< double, nDim, nDim > _computeSmoothingDomainDeformationTensorTotal();

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
      _smoothingVolumeUpdateType( smoothingVolumeUpdateType ),
      _vertexCoordinates_Undeformed( vertexCoordinates ),
      _cellForGeometryUndeformed( vertexCoordinates, nVertexCoordinates ),
      _cellForGeometryIntermediate( vertexCoordinates, nVertexCoordinates ),
      _cellForSmoothing( vertexCoordinates, nVertexCoordinates )
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
    //
    for ( int i = 0; i < _cellForGeometryIntermediate.getNumberOfFaces(); i++ ) {

      auto faceCenterCoords = _cellForGeometryIntermediate.getFaceCenterCoordinates( i + 1 );
      auto n                = _cellForGeometryIntermediate.boundarySurfaceVector( i + 1 );

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
