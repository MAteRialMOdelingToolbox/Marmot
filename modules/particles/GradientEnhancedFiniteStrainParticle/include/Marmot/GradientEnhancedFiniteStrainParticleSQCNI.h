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
 * Thomas Mader    thomas.mader@boku.ac.at
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

#include "Marmot/GradientEnhancedFiniteStrainParticle.h"
#include "Marmot/MarmotMeshfreeApproximation.h"
#include "Marmot/MarmotMeshfreeQuadHexCell.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Fastor/Fastor.h>

namespace Marmot::Meshfree {

  /**
   * @brief Stabilized-conforming nodal integration (SQCNI) particle for gradient-enhanced
   *        finite-strain materials WITHOUT a micropolar continuum.
   *
   * The non-micropolar sibling of @ref GradientEnhancedMicropolarParticleSQCNI, and a near
   * verbatim port of it: apart from the node block losing the micro-rotation slot
   * (@f$n_\mathrm{dim}+1@f$ instead of @f$n_\mathrm{dim}+n_\mathrm{rot}+1@f$) nothing in the
   * smoothing-domain machinery is micropolar.
   *
   * What it adds over the plain @ref GradientEnhancedFiniteStrainParticle point particle:
   *  - the shape-function GRADIENTS come from a divergence (boundary) integral over the
   *    smoothing domain rather than from a point evaluation, which is what makes nodal
   *    integration stable;
   *  - the particle has VERTICES and FACES, so it can carry a distributed load.  This is what
   *    makes a constant-pressure confinement (a genuine triaxial test) possible at all --
   *    the point particle cannot take one, having no faces.
   *
   * `SmoothingDomainUpdateType` selects how the smoothing domain follows the deformation:
   * `DeformationGradient` (full SQCNI), `None` (SNNI, domain fixed), `RotationOnly` and
   * `RotationAndPrincipalStretch`.
   */
  template < int nDim, int nVertices >
  class GradientEnhancedFiniteStrainParticleSQCNI : public GradientEnhancedFiniteStrainParticle< nDim > {

    using TensorD          = Fastor::Tensor< double, nDim >;
    using CoordinatesSized = Eigen::Matrix< double, nDim, 1 >;

  public:
    enum SmoothingDomainUpdateType { None, DeformationGradient, RotationOnly, RotationAndPrincipalStretch };

  protected:
    using LagrangeCellType = MarmotLagrangeCell< nDim, nVertices >;

    const SmoothingDomainUpdateType _smoothingVolumeUpdateType;

    const Eigen::Matrix< double, nDim, nVertices > _vertexCoordinates_Undeformed;

    double* _vertexDisplacements_SmoothingDomain;

    using ParentPointParticle = GradientEnhancedFiniteStrainParticle< nDim >;

    /// \brief Build a Lagrange cell from the current (deformed) smoothing domain vertex coordinates
    LagrangeCellType _makeSmoothingDomainCell() const
    {
      Eigen::Matrix< double, nDim, nVertices > vertexCoordinates;
      getVertexCoordinates( vertexCoordinates.data() );
      return LagrangeCellType( vertexCoordinates.data(), nDim * nVertices );
    }

    /// \brief Build a Lagrange cell from the undeformed vertex coordinates
    LagrangeCellType _makeUndeformedCell() const
    {
      return LagrangeCellType( _vertexCoordinates_Undeformed.data(), nDim * nVertices );
    }

  public:
    virtual void getVertexCoordinates( double* coordinates ) const override;

    virtual void getVisualizationVertexCoordinates( double* coordinates ) const override
    {
      getVertexCoordinates( coordinates );
    };

    virtual int getNumberOfVertices() const override { return nVertices; };

    virtual std::string getParticleShape() const override { return _makeUndeformedCell().getCellShape(); }

    GradientEnhancedFiniteStrainParticleSQCNI( int                                elementID,
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
    virtual double getSmoothingVolume() const { return _makeSmoothingDomainCell().volume(); }

    virtual void getCenterCoordinates( double* coordinates ) const override
    {
      this->_mp.getCoordinatesAtCenter( coordinates );
    }

    virtual CoordinatesSized getCenterFromVertices(
      const Eigen::Matrix< double, nDim, nVertices >& vertexCoordinates ) const
    {
      return LagrangeCellType( vertexCoordinates.data(), nDim * nVertices ).centroid();
    }

    virtual double getVolumeFromVertices( const Eigen::Matrix< double, nDim, nVertices >& vertexCoordinates ) const
    {
      return LagrangeCellType( vertexCoordinates.data(), nDim * nVertices ).volume();
    }

    virtual double getVolumeDeformed() const
    {
      return this->getVolumeUndeformed() * Fastor::determinant( this->_mp.dY_dX() );
    }

    virtual int getNumberOfRequiredStateVars() const override
    {
      return GradientEnhancedFiniteStrainParticle< nDim >::getNumberOfRequiredStateVars() + nDim * nVertices;
    };

    void assignStateVars( double* stateVars, int nStateVars ) override
    {
      GradientEnhancedFiniteStrainParticle< nDim >::assignStateVars( stateVars, nStateVars - nDim * nVertices );
      _vertexDisplacements_SmoothingDomain = stateVars + nStateVars - nDim * nVertices;
    }

    virtual StateView getStateView( const std::string& stateName, int qp ) const override
    {
      if ( stateName == "vertex displacements" ) {
        return { _vertexDisplacements_SmoothingDomain, nDim * nVertices };
      }
      return GradientEnhancedFiniteStrainParticle< nDim >::getStateView( stateName, qp );
    }

    virtual void computeDistributedLoad( int           type,
                                         int           surfaceID,
                                         const double* load,
                                         double*       fExt,
                                         double*       dExt_dQ,
                                         double        timeNew,
                                         double        dT ) const override;

    std::tuple< TensorD, TensorD > getIntermediateConfBoundaryVector( int boundaryFaceID ) const;

    virtual void vci_compute_Test_P_BoundaryIntegral( double*       R_AiC_RowMajor,
                                                      const double* boundarySurfaceVector,
                                                      int           boundaryFaceID )
    {
      using namespace Fastor;

      const auto [N_dAY, Y_N]   = getIntermediateConfBoundaryVector( boundaryFaceID );
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
      const auto cell = _makeSmoothingDomainCell();

      Eigen::Map< Eigen::Matrix< double, nDim, Eigen::Dynamic > > faceCenters( coordinates,
                                                                               nDim,
                                                                               cell.getNumberOfFaces() );

      for ( int i = 0; i < cell.getNumberOfFaces(); i++ )
        faceCenters.col( i ) = cell.getFaceCenterCoordinates( i + 1 );
    }

    virtual void getFaceCoordinates( int faceID, double* coordinates ) const
    {
      Eigen::Map< Eigen::Matrix< double, nDim, 1 > > faceCenter( coordinates );
      faceCenter = _makeSmoothingDomainCell().getFaceCenterCoordinates( faceID );
    }

    virtual int getNumberOfEvaluationPoints() const { return _makeUndeformedCell().getNumberOfFaces(); };

  private:
    void _updateVertexDisplacementsFromMaterialPointDeformation();

    virtual void updateParticlePositionToReferenceIntermediate() override
    {
      _updateVertexDisplacementsFromMaterialPointDeformation();

      getCenterCoordinates( ParentPointParticle::_centerReferenceIntermediate.data() );
    };
  };

  template < int nDim, int nVertices >
  GradientEnhancedFiniteStrainParticleSQCNI< nDim, nVertices >::GradientEnhancedFiniteStrainParticleSQCNI(
    int                                                  elementID,
    const double*                                        vertexCoordinates,
    int                                                  nVertexCoordinates,
    double                                               volume,
    const std::string&                                   materialName,
    const double*                                        materialProperties,
    int                                                  sizeMaterialProperties,
    const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation,
    const SmoothingDomainUpdateType                      smoothingVolumeUpdateType )
    : GradientEnhancedFiniteStrainParticle< nDim >( elementID,
                                                  getCenterFromVertices(
                                                    Eigen::Map< const Eigen::Matrix< double, nDim, nVertices > >(
                                                      vertexCoordinates ) )
                                                    .data(),
                                                  CoordinatesSized::RowsAtCompileTime,
                                                  getVolumeFromVertices(
                                                    Eigen::Map< const Eigen::Matrix< double, nDim, nVertices > >(
                                                      vertexCoordinates ) ),
                                                  materialName,
                                                  materialProperties,
                                                  sizeMaterialProperties,
                                                  approximation ),
      _smoothingVolumeUpdateType( smoothingVolumeUpdateType ),
      _vertexCoordinates_Undeformed( vertexCoordinates )

  {
  }

  template < int nDim, int nVertices >
  void GradientEnhancedFiniteStrainParticleSQCNI< nDim, nVertices >::computeDistributedLoad( int           type,
                                                                                           int           boundaryFaceID,
                                                                                           const double* load_,
                                                                                           double*       fExt,
                                                                                           double*       dFExt_ddQ,
                                                                                           double        timeNew,
                                                                                           double        dT ) const
  {

    switch ( type ) {

    case GradientEnhancedFiniteStrainParticle< nDim >::Pressure: {

      const auto&   _mp           = GradientEnhancedFiniteStrainParticle< nDim >::_mp;
      const auto&   _nNodes       = GradientEnhancedFiniteStrainParticle< nDim >::_nNodes;
      constexpr int nodeBlockSize = nDim + 1;

      const auto [N_dAY, Y_N] = getIntermediateConfBoundaryVector( boundaryFaceID );

      TensorD fY = N_dAY * load_[0];

      Eigen::Map< Eigen::VectorXd > P( fExt, _nNodes * nodeBlockSize );
      Eigen::Map< Eigen::MatrixXd > K( dFExt_ddQ, _nNodes * nodeBlockSize, _nNodes * nodeBlockSize );

      using namespace Fastor;
      using namespace FastorIndices;

      Tensor< double, nDim, nDim > Eye;
      Eye.eye();

      // apply Nanson's formula
      const auto deltaF = _mp.dx_dY();

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

        /* const double T_A = GradientEnhancedFiniteStrainParticle< nDim >::_T( A ); */
        r_U = testBoundary( A ) * f;

        {
          using namespace Eigen;
          P.template segment< nDim >( idxA_u ) -= Map< Matrix< double, nDim, 1 > >( r_U.data() );
        }

        for ( int B = 0; B < _nNodes; B++ ) {
          const int  idxB_u  = nodeBlockSize * B;
          const auto dN_B_dY = TensorMap< const double, nDim >(
            GradientEnhancedFiniteStrainParticle< nDim >::_dN_dY.col( B ).data() );

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
    case GradientEnhancedFiniteStrainParticle< nDim >::CWFCorrection: {
      const auto&   _mp           = GradientEnhancedFiniteStrainParticle< nDim >::_mp;
      const auto&   _nNodes       = GradientEnhancedFiniteStrainParticle< nDim >::_nNodes;
      constexpr int nodeBlockSize = nDim + 1;

      const auto [N_dAY, Y_N] = getIntermediateConfBoundaryVector( boundaryFaceID );

      Eigen::Map< Eigen::VectorXd > P( fExt, _nNodes * nodeBlockSize );
      Eigen::Map< Eigen::MatrixXd > K( dFExt_ddQ, _nNodes * nodeBlockSize, _nNodes * nodeBlockSize );

      using namespace Fastor;
      using namespace FastorIndices;

      Tensor< double, nDim, nDim > Eye;
      Eye.eye();

      const auto& S = _mp.response.S;
      // const auto& t = _mp.tangents;
      // apply Nanson's formula
      const auto                         deltaF    = _mp.dx_dY();
      const Tensor< double, nDim, nDim > deltaFInv = inverse( deltaF );
      const double                       deltaJ    = determinant( deltaF );

      const TensorD n_dA = deltaJ * transpose( deltaFInv ) % N_dAY;

      const Tensor< double, nDim, nDim, nDim, nDim > dFInv_dF = -einsum< Ik, Ki, to_IikK >( deltaFInv, deltaFInv );

      const Tensor< double, nDim, nDim, nDim > dndA_dDeltaF = outer( n_dA, transpose( deltaFInv ) ) +
                                                              deltaJ * einsum< IikK, Index< I_ > >( dFInv_dF, N_dAY );

      TensorD r_U( 0.0 );

      Eigen::MatrixXd testBoundary = Eigen::MatrixXd::Zero( 1, ParentPointParticle::_nNodes );

      ParentPointParticle::_meshfreeApproximation.computeShapeFunctions( Y_N.data(),
                                                                         ParentPointParticle::_assignedKernelFunctions,
                                                                         testBoundary.data() );

      for ( int A = 0; A < _nNodes; A++ ) {
        const int idxA_u = nodeBlockSize * A;

        r_U = testBoundary( A ) * S % n_dA;

        {
          using namespace Eigen;
          P.template segment< nDim >( idxA_u ) -= Map< Matrix< double, nDim, 1 > >( r_U.data() );
        }
        // for ( int B = 0; B < _nNodes; B++ ) {
        //   const int  idxB_u  = nodeBlockSize * B;
        //   const auto dN_B_dY = TensorMap< const double, nDim >(
        //     GradientEnhancedFiniteStrainParticle< nDim >::_dN_dY.col( B ).data() );

        //  const auto dS_dqU_B = evaluate ( + einsum < ijkl, l > ( t.dS_dDeltaF, dN_B_dY ));
        //  const auto dS_dqU_B_ndA = testBoundary(A) * dS_dqU_B % n_dA;

        //  const Tensor< double, nDim, nDim > dndA_ddQU_B = testBoundary( A ) * einsum< ijk, k >( dndA_dDeltaF, dN_B_dY
        //  ); const Tensor< double, nDim, nDim > S_dndA_ddQU_B = S % dndA_ddQU_B;

        //  const Tensor< double, nDim, nDim > cwfTangent = dS_dqU_B_ndA + S_dndA_ddQU_B;

        //  {
        //    using namespace Eigen;
        //    K.template block< nDim, nDim >( idxA_u, idxB_u ) -= Map< Matrix< double, nDim, nDim > >(
        //      torowmajor( cwfTangent ).data() );
        //  }
        //}
      }
      break;
    }
    default: {
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid DistributedLoad type specified" );
    }
    }
  }

  template < int nDim, int nVertices >
  std::tuple< typename GradientEnhancedFiniteStrainParticleSQCNI< nDim, nVertices >::TensorD,
              typename GradientEnhancedFiniteStrainParticleSQCNI< nDim, nVertices >::TensorD >
  GradientEnhancedFiniteStrainParticleSQCNI< nDim, nVertices >::getIntermediateConfBoundaryVector(
    int boundaryFaceID ) const
  {

    TensorD N_dAY;
    TensorD Y;

    if ( _smoothingVolumeUpdateType == DeformationGradient ) {
      // For the full SQCNI case, we actually operate on the deformed element.
      // This means, that the deformed smoothing domain is consistent with the physical domain of the particle.
      const auto cell = _makeSmoothingDomainCell();

      const Eigen::Matrix< double, nDim, 1 > n_dA = cell.boundarySurfaceVector( boundaryFaceID );
      const Eigen::Matrix< double, nDim, 1 > y    = cell.getFaceCenterCoordinates( boundaryFaceID );

      for ( int i = 0; i < nDim; i++ ) {
        N_dAY[i] = n_dA[i];
        Y[i]     = y[i];
      }
    }
    else {
      // In all other cases, the deformed smoothing domain is not consistent with the physical domain of the particle.
      // Accordingly, it is not reasonable to let the boundary vector be acting on the deformed smoothing domain.
      // Instead, we compute the boundary vector in the undeformed configuration and apply the deformation gradient
      // using Nanson's formula. Also, the let the origin of the boundary segment be the center of the particle in the
      // deformed configuration.

      const auto cellUndeformed = _makeUndeformedCell();

      const Eigen::Matrix< double, nDim, 1 > N_dA0_eigen = cellUndeformed.boundarySurfaceVector( boundaryFaceID );
      const Eigen::Matrix< double, nDim, 1 > Y0          = cellUndeformed.getFaceCenterCoordinates( boundaryFaceID );

      const auto FY = this->_mp.dY_dX();

      // apply Nanson's formula to shift to Y configuration:
      using namespace Fastor;
      N_dAY = determinant( FY ) * transpose( inverse( FY ) ) % TensorD( N_dA0_eigen.data() );

      // and the center of the boundary segment in the deformed configuration
      this->_mp.getCenterDisplacement( Y.data() );
      Y += TensorD( Y0.data() );
    }

    return { N_dAY, Y };
  }

  template < int nDim, int nVertices >
  void GradientEnhancedFiniteStrainParticleSQCNI< nDim, nVertices >::getVertexCoordinates( double* coordinates ) const
  {

    Eigen::Map< Eigen::Matrix< double, nDim, nVertices > > coordinatesDeformed( coordinates );
    Eigen::Map< Eigen::Matrix< double, nDim, nVertices > > vertexDisplacementsSmoothingDomain(
      _vertexDisplacements_SmoothingDomain );

    coordinatesDeformed = _vertexCoordinates_Undeformed + vertexDisplacementsSmoothingDomain;
  }

  template < int nDim, int nVertices >
  void GradientEnhancedFiniteStrainParticleSQCNI< nDim,
                                                nVertices >::_updateVertexDisplacementsFromMaterialPointDeformation()
  {

    Eigen::Matrix< double, nDim, nDim > F;

    switch ( _smoothingVolumeUpdateType ) {

    case SmoothingDomainUpdateType::None: F.setIdentity(); break;

    case SmoothingDomainUpdateType::DeformationGradient: {

      Fastor::Tensor< double, nDim, nDim > F_ = GradientEnhancedFiniteStrainParticle< nDim >::_mp.dY_dX();
      F = Eigen::Map< Eigen::Matrix< double, nDim, nDim, Eigen::RowMajor > >( F_.data() );
      break;
    }

    case SmoothingDomainUpdateType::RotationOnly: {

      Fastor::Tensor< double, nDim, nDim > F_ = GradientEnhancedFiniteStrainParticle< nDim >::_mp.dY_dX();
      F = Eigen::Map< Eigen::Matrix< double, nDim, nDim, Eigen::RowMajor > >( F_.data() );

      Eigen::JacobiSVD< Eigen::MatrixXd > svd;
      svd.compute( F, Eigen::ComputeFullU | Eigen::ComputeFullV );

      Eigen::Matrix< double, nDim, nDim > R = svd.matrixU() * svd.matrixV().transpose();

      F = R;
      break;
    }

    case SmoothingDomainUpdateType::RotationAndPrincipalStretch: {
      Fastor::Tensor< double, nDim, nDim > F_ = GradientEnhancedFiniteStrainParticle< nDim >::_mp.dY_dX();
      F = Eigen::Map< Eigen::Matrix< double, nDim, nDim, Eigen::RowMajor > >( F_.data() );

      Eigen::JacobiSVD< Eigen::MatrixXd > svd;
      svd.compute( F, Eigen::ComputeFullU | Eigen::ComputeFullV );

      Eigen::Matrix< double, nDim, nDim > R = svd.matrixU() * svd.matrixV().transpose();

      Eigen::Matrix< double, nDim, nDim > U_ = Eigen::Matrix< double, nDim, nDim >::Identity();
      U_.diagonal()                          = ( R.transpose() * F ).diagonal();

      F = R * U_;

      break;
    }
    }

    Eigen::Matrix< double, nDim, nVertices > coordinatesRelToCenter0 = _vertexCoordinates_Undeformed;
    for ( int i = 0; i < nDim; i++ ) {
      coordinatesRelToCenter0.row( i ).array() -= this->_centerCoordinatesUndeformed( i );
    }

    Eigen::Matrix< double, nDim, nVertices > vertexCoordinatesDeformed = F * coordinatesRelToCenter0;

    Eigen::Matrix< double, nDim, 1 > _mp_displacement;
    GradientEnhancedFiniteStrainParticle< nDim >::_mp.getCenterDisplacement( _mp_displacement.data() );

    for ( int i = 0; i < nDim; i++ ) {
      vertexCoordinatesDeformed.row( i ).array() += this->_centerCoordinatesUndeformed( i ) + _mp_displacement( i );
    }

    Eigen::Map< Eigen::Matrix< double, nDim, nVertices > > vertexDisplacements( _vertexDisplacements_SmoothingDomain );
    vertexDisplacements = vertexCoordinatesDeformed - _vertexCoordinates_Undeformed;
  }

  template < int nDim, int nVertices >
  void GradientEnhancedFiniteStrainParticleSQCNI< nDim, nVertices >::assignMeshfreeKernelFunctions(
    const std::vector< const MarmotMeshfreeKernelFunction* >& kernelFunctions )
  {

    ParentPointParticle::_assignedKernelFunctions = kernelFunctions;

    ParentPointParticle::_nNodes = GradientEnhancedFiniteStrainParticle< nDim >::_assignedKernelFunctions.size();

    Eigen::Matrix< double, nDim, 1 > coords;
    ParentPointParticle::_mp.getCoordinatesAtCenter( coords.data() );

    ParentPointParticle::_N     = Eigen::MatrixXd::Zero( 1, ParentPointParticle::_nNodes );
    ParentPointParticle::_dN_dY = Eigen::MatrixXd::Zero( nDim, ParentPointParticle::_nNodes );

    ParentPointParticle::_meshfreeApproximation.computeShapeFunctions( coords.data(),
                                                                       ParentPointParticle::_assignedKernelFunctions,
                                                                       ParentPointParticle::_N.data() );

    Eigen::MatrixXd NBoundary( 1, ParentPointParticle::_nNodes );

    const auto cell = _makeSmoothingDomainCell();

    for ( int i = 0; i < cell.getNumberOfFaces(); i++ ) {

      const Eigen::Matrix< double, nDim, 1 > n_dA       = cell.boundarySurfaceVector( i + 1 );
      const Eigen::Matrix< double, nDim, 1 > faceCenter = cell.getFaceCenterCoordinates( i + 1 );

      ParentPointParticle::_meshfreeApproximation.computeShapeFunctions( faceCenter.data(),
                                                                         ParentPointParticle::_assignedKernelFunctions,
                                                                         NBoundary.data() );

      ParentPointParticle::_dN_dY += n_dA * NBoundary;
    }

    ParentPointParticle::_dN_dY /= cell.volume();

    ParentPointParticle::_T     = ParentPointParticle::_N;
    ParentPointParticle::_dT_dY = ParentPointParticle::_dN_dY;
  }

} // namespace Marmot::Meshfree
