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
#include "Marmot/GenericSDIParticle.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Fastor/Fastor.h>
#include <stdexcept>

namespace Marmot::Meshfree {

  /**
   * @class Marmot::Meshfree::DisplacementParticleSQCNIxSDI
   * @brief Displacement particle with subdomain integration (SDI): smoothed (SQCNI/SNNI) gradients on each subdomain
   * of the particle.
   *
   * @details The particle geometry is uniformly subdivided by GenericSDIParticle (into 4 quadrilaterals in 2D, 8
   * hexahedra in 3D). Each subdomain @f$ s @f$ owns a DisplacementMaterialPoint (see MaterialPointType) at its center,
   * with its undeformed volume @f$ V_{0,s} @f$, and its own shape functions: @f$ N_B @f$ at the subdomain center and
   * the gradients smoothed over the subdomain's smoothing domain (whose update follows the SmoothingDomainUpdateType,
   * as in DisplacementParticleSQCNI). Each subdomain contributes the residual and the tangent of DisplacementParticle
   * (Kirchhoff stress, geometric stiffness, Newmark-beta inertia),
   * @f[
   *   r_{Aj} = \sum_s \left( \frac{\partial T^s_A}{\partial x_i}\,\tau^s_{ij}
   *     + \rho_0\,a^s_j\,T^s_A \right) V_{0,s},
   * @f]
   * so that the tangent is exact. The central displacement and deformation gradient of the whole particle, computed by
   * GenericSDIParticle from the smoothed gradient over the whole particle, move the particle geometry and the
   * subdomains at acceptStateAndPosition() and enter the pressure load.
   *
   * Properties (in this order): "VCI order" (from GenericSDIParticle), "newmark-beta beta", "newmark-beta gamma".
   *
   * @tparam nDim The number of dimensions (2 or 3).
   * @tparam nVertices The number of vertices of the particle geometry (4 or 8).
   */
  template < int nDim, int nVertices >
  class DisplacementParticleSQCNIxSDI : public GenericSDIParticle< nDim, nVertices > {

    using TensorD  = GenericSDIParticle< nDim, nVertices >::TensorD;  ///< vector of size nDim
    using TensorDD = GenericSDIParticle< nDim, nVertices >::TensorDD; ///< second-order tensor of size nDim
    /// the material points, one per subdomain (in the order of the subdomains)
    std::vector< std::unique_ptr< MaterialPointType< nDim > > > _subdomainMaterialPoints;

    /// number of vertex and center displacement values (not used by this class)
    constexpr int static nStateVarsParticle = nDim * nVertices + nDim; // vertex displacements + center displacement

    double _newmark_beta;  ///< Newmark parameter @f$ \beta @f$ (property "newmark-beta beta", default 0)
    double _newmark_gamma; ///< Newmark parameter @f$ \gamma @f$ (property "newmark-beta gamma", default 0)

    /// static vector of valid properties
    inline static const std::vector< std::string > _validProperties = {
      "newmark-beta beta",
      "newmark-beta gamma",
    };

  public:
    /// Body load types.
    enum BodyLoadTypes {
      BodyForce, ///< body force per unit undeformed volume ("BODYFORCE")
    };

    /// Distributed load types.
    enum DistributedLoadTypes {
      Pressure ///< follower pressure on a particle face ("PRESSURE")
    };

    /**
     * @brief Supported body loads.
     * @return "BODYFORCE".
     */
    const std::unordered_map< std::string, int >& getSupportedBodyLoadTypes() const override
    {
      static const std::unordered_map< std::string, int > _supportedBodyLoadTypes = { { "BODYFORCE", BodyForce } };
      return _supportedBodyLoadTypes;
    };

    /**
     * @brief Supported distributed loads.
     * @return "PRESSURE".
     */
    const std::unordered_map< std::string, int >& getSupportedDistributedLoadTypes() const override
    {
      static const std::unordered_map< std::string, int > _supportedDistributedLoadTypes = { { "PRESSURE", Pressure } };
      return _supportedDistributedLoadTypes;
    };

    static constexpr int nDofPerNodeU = nDim; ///< dofs per node of the displacement field

    /**
     * @brief Sets a property of the subdomains ("newmark-beta beta" or "newmark-beta gamma").
     * @param[in] propertyName Name of the property.
     * @param[in] property Its value.
     * @throws std::runtime_error for an unknown name.
     */
    virtual void setPropertyOnSubdomains( const std::string& propertyName, const double* property ) override
    {
      if ( propertyName == "newmark-beta beta" ) {
        _newmark_beta = property[0];
      }
      else if ( propertyName == "newmark-beta gamma" ) {
        _newmark_gamma = property[0];
      }
      else {
        throw std::runtime_error( "Property " + propertyName + " not supported!" );
      }
    };

    /// \brief Get the names of the properties
    /// \return The names of the properties
    virtual std::vector< std::string > getSubdomainPropertyNames() const override { return _validProperties; };

    /**
     * @brief Number of dofs per node.
     * @return nDim.
     */
    virtual int getNBaseDof() const override { return nDofPerNodeU; };

    /**
     * @brief Initializes the material points of all subdomains.
     */
    void initializeYourselfOnSubdomains() override
    {
      for ( auto& mp : _subdomainMaterialPoints ) {
        mp->initializeYourself();
      }
    };

    /**
     * @brief Body load: a body force @f$ \boldsymbol{b} @f$ per unit undeformed volume (dead load), integrated over
     * the subdomains, @f$ P_{Ai} \mathrel{-}= \sum_s T^s_A\,b_i\,V^s_0 @f$; the tangent is zero.
     * @param[in] type The body load type (BodyForce).
     * @param[in] load The body force vector (nDim values).
     * @param[in,out] fExt Load vector, the contribution is added.
     * @param[in,out] dExt_dQ Tangent (not modified).
     * @param[in] timeNew Time at the end of the increment.
     * @param[in] dT Time increment.
     * @throws std::invalid_argument for another load type.
     */
    virtual void computeBodyLoad( int           type,
                                  const double* load,
                                  double*       fExt,
                                  double*       dExt_dQ,
                                  double        timeNew,
                                  double        dT ) const override
    {
      switch ( type ) {
      case BodyForce: {
        // integrated over the subdomains, with their test functions
        for ( size_t s = 0; s < this->_subDomainShapeFunctions.size(); s++ ) {
          const double V0 = _subdomainMaterialPoints[s]->getVolumeUndeformed();
          const auto&  sd = this->_subDomainShapeFunctions[s];
          for ( int A = 0; A < this->_nNodes; A++ )
            for ( int i = 0; i < nDofPerNodeU; i++ )
              fExt[nDofPerNodeU * A + i] -= sd.T( A ) * load[i] * V0;
        }
        break;
      }
      default: throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid body load type" );
      }
    }

    /**
     * @brief Node fields.
     * @return "displacement".
     */
    virtual const std::vector< std::string >& getFields() const override
    {
      static const std::vector< std::string > nodeFields = { "displacement" };
      return nodeFields;
    };

    /**
     * @brief Constructs the particle, its subdomains and one material point per subdomain.
     * @param[in] elementID Label of the particle (also the label of the material points).
     * @param[in] nodeCoordinates Vertex coordinates of the particle geometry in the undeformed configuration
     * (nDim * nVertices values).
     * @param[in] nNodeCoordiantes Number of coordinates (nDim * nVertices).
     * @param[in] volume Volume argument, passed to GenericSDIParticle.
     * @param[in] materialName Name of the finite-strain material.
     * @param[in] materialProperties Material properties.
     * @param[in] sizeMaterialProperties Number of material properties.
     * @param[in] approximation The meshfree approximation for the shape functions.
     * @param[in] smoothingVolumeUpdateType How the smoothing domains follow the deformation.
     * @throws std::invalid_argument for an unknown material.
     */
    DisplacementParticleSQCNIxSDI(
      int                                                                    elementID,
      const double*                                                          nodeCoordinates,
      int                                                                    nNodeCoordiantes,
      double                                                                 volume,
      const std::string&                                                     materialName,
      const double*                                                          materialProperties,
      int                                                                    sizeMaterialProperties,
      const MarmotMeshfreeApproximation&                                     approximation,
      const GenericSDIParticle< nDim, nVertices >::SmoothingDomainUpdateType smoothingVolumeUpdateType );

    /**
     * @brief Volume in the undeformed configuration, the sum over the subdomains.
     * @return @f$ V_0 = \sum_s V_{0,s} @f$.
     */
    virtual double getVolumeUndeformed() const override
    {

      double V0 = 0.0;
      for ( const auto& mp : _subdomainMaterialPoints ) {
        V0 += mp->getVolumeUndeformed();
      }
      return V0;
    }

    /**
     * @brief Volume in the intermediate configuration (last accepted state), the sum over the subdomains.
     * @return @f$ \sum_s V_{0,s}\,\det\boldsymbol{F}_{n,s} @f$.
     */
    virtual double getVolumeDeformed() const
    {
      double volDeformed = 0.0;
      for ( const auto& mp : _subdomainMaterialPoints ) {
        volDeformed += mp->getVolumeUndeformed() * determinant( mp->dY_dX() );
      }
      return volDeformed;
    }

    /**
     * @brief Volume of a subdomain in the intermediate configuration, from its material point.
     * @param[in] subdomain The subdomain (one of the subdomains of this particle).
     * @return @f$ V_{0,s}\,\det\boldsymbol{F}_{n,s} @f$.
     */
    virtual double getSubdomainVolume( const ParticleDomain< nDim, nVertices >& subdomain ) const override
    {
      int subdomainIndex = -1;
      for ( size_t i = 0; i < this->_subDomains.size(); i++ ) {
        if ( &subdomain == &( this->_subDomains[i] ) ) {
          subdomainIndex = static_cast< int >( i );
          break;
        }
      }

      assert( subdomainIndex != -1 && "Subdomain not found!" );

      const auto& mp = _subdomainMaterialPoints[subdomainIndex];

      return mp->getVolumeUndeformed() * determinant( mp->dY_dX() );
    };

    /**
     * @brief Accepts the increment of the material points of all subdomains.
     */
    virtual void acceptStateAndPositionOnSubdomains() override
    {
      for ( auto& mp : _subdomainMaterialPoints ) {
        mp->acceptStateAndPosition();
      }
    };

    /**
     * @brief Number of state variables of the subdomains: the state of each material point, padded to a multiple of 8.
     * @return The number of state variables.
     */
    virtual int getNumberOfRequiredStateVarsOnSubdomains() const override
    {
      int nStateVars = 0;

      for ( const auto& mp : _subdomainMaterialPoints ) {
        nStateVars += this->paddedStateVarSize( mp->getNumberOfRequiredStateVars() );
      }

      return nStateVars;
    };

    /**
     * @brief Assigns consecutive (padded) blocks of the state vector to the material points of the subdomains.
     * @param[in] stateVars State vector of the subdomains.
     * @param[in] nStateVars Its length.
     * @throws std::runtime_error if the length does not match getNumberOfRequiredStateVarsOnSubdomains().
     */
    virtual void assignStateVarsOnSubdomains( double* stateVars, int nStateVars ) override
    {
      int offset = 0;

      for ( auto& mp : _subdomainMaterialPoints ) {

        int nStateVarsSubParticle = mp->getNumberOfRequiredStateVars();
        mp->assignStateVars( stateVars + offset, nStateVarsSubParticle );
        offset += this->paddedStateVarSize( nStateVarsSubParticle ); // every subdomain block starts aligned
      }

      if ( offset != nStateVars ) {
        throw std::runtime_error( "Error: Number of state variables does not match!" );
      }
    }

    /**
     * @brief State of the material point of a subdomain.
     * @param[in] stateName Name of the state (see DisplacementMaterialPoint::getStateView()).
     * @param[in] subdomainIndex Index of the subdomain.
     * @return The view on the state.
     */
    virtual StateView getStateViewOnSubdomains( const std::string& stateName, int subdomainIndex ) const override
    {
      return _subdomainMaterialPoints[subdomainIndex]->getStateView( stateName );
    }

    /**
     * @brief Updates the material points of all subdomains with the increment dQ and assembles their residuals and
     * tangents (see the class description), including the Newmark-beta update of their velocities and accelerations.
     * @param[in] dQ Nodal displacement increments of the current step (nDim values per node).
     * @param[in,out] fInt Residual, the contribution is added.
     * @param[in,out] dFInt_ddQ Tangent @f$ \partial r/\partial\Delta q @f$, the contribution is added.
     * @param[in] timeNew Time at the end of the increment.
     * @param[in] dT Time increment.
     */
    virtual void computePhysicsKernelsOnSubdomains( const double* dQ,
                                                    double*       fInt,
                                                    double*       dFInt_ddQ,
                                                    double        timeNew,
                                                    double        dT ) override;

    /**
     * @brief Follower pressure on a face of the particle geometry.
     * @details @f$ \boldsymbol{f} = \Delta J_c\,\Delta\boldsymbol{F}_c^{-\mathsf T}\,p\,\boldsymbol{N}\,dA_Y @f$
     * with the central incremental deformation gradient @f$ \Delta\boldsymbol{F}_c @f$ and the boundary vector of the
     * face in the intermediate configuration, assembled as @f$ r_{Aj} \mathrel{-}= T_A\,f_j @f$ with @f$ T_A @f$ at
     * the face center of the smoothing domain; the tangent uses the smoothed gradient over the whole particle.
     * @param[in] type The distributed load type (Pressure).
     * @param[in] surfaceID The face ID (1-based).
     * @param[in] load The pressure @f$ p @f$ (load[0]).
     * @param[in,out] fExt Load vector, the contribution is added.
     * @param[in,out] dExt_dQ Tangent, the contribution is added.
     * @param[in] timeNew Time at the end of the increment.
     * @param[in] dT Time increment.
     * @throws std::invalid_argument for an unsupported type.
     */
    virtual void computeDistributedLoad( int           type,
                                         int           surfaceID,
                                         const double* load,
                                         double*       fExt,
                                         double*       dExt_dQ,
                                         double        timeNew,
                                         double        dT ) const override;
    /**
     * @brief VCI boundary term @f$ R_{AiC} \mathrel{+}= T_A(\boldsymbol{Y}_N)\,P_C(\boldsymbol{Y}_N)\,(N\,dA_Y)_i @f$
     * at the face center @f$ \boldsymbol{Y}_N @f$ of the particle geometry, with its boundary vector in the
     * intermediate configuration.
     * @param[in,out] R_AiC_RowMajor VCI matrix (nNodes x nDim x nVCIConstraints, row major), the contribution is
     * added.
     * @param[in] boundarySurfaceVector Boundary surface vector (not used; the face of the particle geometry is used).
     * @param[in] boundaryFaceID The face ID (1-based).
     */
    virtual void vci_compute_Test_P_BoundaryIntegral( double*       R_AiC_RowMajor,
                                                      const double* boundarySurfaceVector,
                                                      int           boundaryFaceID ) override
    {
      const auto [N_dAY, Y_N]   = this->getIntermediateConfigurationBoundaryVector( boundaryFaceID,
                                                                                  this->_particleDomainMain );
      Eigen::MatrixXd TBoundary = Eigen::MatrixXd::Zero( 1, this->_nNodes );

      this->_meshfreeApproximation.computeShapeFunctions( Y_N.data(),
                                                          this->_assignedKernelFunctions,
                                                          TBoundary.data() );

      Eigen::VectorXd                  PBoundary( this->_nVCIConstraints );
      Eigen::Matrix< double, nDim, 1 > Y_N_coords( Y_N.data() );
      Math::computeMonomialBasis( this->_vciOrder, Y_N_coords, PBoundary );

      for ( int A = 0; A < this->_nNodes; A++ ) {
        for ( int i = 0; i < nDim; i++ ) {
          for ( int C = 0; C < this->_nVCIConstraints; C++ ) {

            R_AiC_RowMajor[A * ( nDim * this->_nVCIConstraints ) + i * this->_nVCIConstraints + C] += TBoundary( 0,
                                                                                                                 A ) *
                                                                                                      PBoundary( C ) *
                                                                                                      N_dAY[i];
            // TBoundary( 0, A ) * PBoundary( C ) * boundarySurfaceVector[i];// N_dAY[i];
          }
        }
      }
    };

    /**
     * @brief Evaluation points: for each face of the particle, the corresponding face centers of the smoothing
     * domains of the subdomains adjacent to that face.
     * @param[out] coordinates The coordinates (nDim x getNumberOfEvaluationPoints(), column major).
     */
    virtual void getEvaluationCoordinates( double* coordinates ) const override
    {
      const int nEvalPoints = this->getNumberOfEvaluationPoints();

      Eigen::Map< Eigen::Matrix< double, nDim, Eigen::Dynamic > > coordinatesMap( coordinates, nDim, nEvalPoints );

      int i = 0;
      for ( int f = 0; f < this->_particleDomainMain.getNumberOfFaces(); f++ ) {
        const auto subcellsAttachedToFace = this->_particleDomainMain.getSubCellIndicesOnParentFace( f + 1 );
        for ( size_t sc = 0; sc < subcellsAttachedToFace.size(); sc++ ) {
          const int subcellIndex  = subcellsAttachedToFace[sc];
          coordinatesMap.col( i ) = this->_subDomains[subcellIndex].getSmoothingDomainFaceCenterCoordinates( f + 1 );
          i++;
        }
      }
    }

    /**
     * @brief Number of evaluation points, see getEvaluationCoordinates().
     * @return The number of evaluation points.
     */
    virtual int getNumberOfEvaluationPoints() const override
    {

      int nEvalPoints = 0;
      for ( int i = 0; i < this->_particleDomainMain.getNumberOfFaces(); i++ ) {
        const auto subcellsAttachedToFace = this->_particleDomainMain.getSubCellIndicesOnParentFace( i + 1 );
        nEvalPoints += static_cast< int >( subcellsAttachedToFace.size() );
      }
      return nEvalPoints;
    };

    /**
     * @brief Initial conditions are not supported.
     * @param[in] conditionName Name of the initial condition.
     * @param[in] value Its values.
     * @throws std::invalid_argument always.
     */
    virtual void setInitialCondition( const std::string& conditionName, const double* value ) override
    {
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid initial condition" );
    };
  };

  template < int nDim, int nVertices >
  DisplacementParticleSQCNIxSDI< nDim, nVertices >::DisplacementParticleSQCNIxSDI(
    int                                                                    elementID,
    const double*                                                          vertexCoordinates,
    int                                                                    nVertexCoordinates,
    double                                                                 volume,
    const std::string&                                                     materialName,
    const double*                                                          materialProperties,
    int                                                                    sizeMaterialProperties,
    const Marmot::Meshfree::MarmotMeshfreeApproximation&                   approximation,
    const GenericSDIParticle< nDim, nVertices >::SmoothingDomainUpdateType smoothingVolumeUpdateType )
    : GenericSDIParticle< nDim, nVertices >( elementID,
                                             vertexCoordinates,
                                             nVertexCoordinates,
                                             volume,
                                             approximation,
                                             smoothingVolumeUpdateType ),
      _newmark_beta( 0. ),
      _newmark_gamma( 0. )
  {
    MarmotMaterialSection section( materialName, materialProperties, sizeMaterialProperties );

    for ( size_t i = 0; i < this->_subDomains.size(); i++ ) {

      const auto& initialParticleDomain = this->_subDomains[i];

      const double subV0 = initialParticleDomain.getVolumeUndeformed();
      const auto   X0    = initialParticleDomain.getCenterCoordinates();

      _subdomainMaterialPoints.push_back(
        std::make_unique< MaterialPointType< nDim > >( elementID, X0.data(), 1, subV0 ) );
    }

    for ( auto& mp : _subdomainMaterialPoints ) {
      mp->assignMaterial( section );
    }
  }

  template < int nDim, int nVertices >
  void DisplacementParticleSQCNIxSDI< nDim, nVertices >::computeDistributedLoad( int           type,
                                                                                 int           boundaryFaceID,
                                                                                 const double* load_,
                                                                                 double*       fExt,
                                                                                 double*       dFExt_ddQ,
                                                                                 double        timeNew,
                                                                                 double        dT ) const
  {

    switch ( type ) {

    case DisplacementParticle< nDim >::Pressure: {

      constexpr int nodeBlockSize = nDim;

      const auto [N_dAY, Y_N] = this->getIntermediateConfigurationBoundaryVector( boundaryFaceID,
                                                                                  this->_particleDomainMain );

      const auto T_Boundary = this->evaluateShapeFunctionsOnFace( this->_particleDomainMain, boundaryFaceID );

      // const auto dT_dY = dN_dY;
      const auto [_, dN_dY] = this->evaluateShapeFunctionsAndDerivativesForParticleDomain( this->_particleDomainMain );

      TensorD fY = N_dAY * load_[0];

      Eigen::Map< Eigen::VectorXd > P( fExt, this->_nNodes * nodeBlockSize );
      Eigen::Map< Eigen::MatrixXd > K( dFExt_ddQ, this->_nNodes * nodeBlockSize, this->_nNodes * nodeBlockSize );

      using namespace Fastor;
      using namespace FastorIndices;

      Tensor< double, nDim, nDim > Eye;
      Eye.eye();

      // apply Nanson's formula
      const TensorDD deltaF    = transpose( TensorDD( this->_centralDeformationGradientDelta.data() ) );
      const TensorDD deltaFInv = inverse( deltaF );
      const double   deltaJ    = determinant( deltaF );

      const TensorD f = deltaJ * transpose( deltaFInv ) % fY;

      const Tensor< double, nDim, nDim, nDim, nDim > dFInv_dF = -einsum< Ik, Ki, to_IikK >( deltaFInv, deltaFInv );

      const Tensor< double, nDim, nDim, nDim > df_dDeltaF = outer( f, transpose( deltaFInv ) ) +
                                                            deltaJ * einsum< IikK, Index< I_ > >( dFInv_dF, fY );

      TensorD r_U( 0.0 );

      for ( int A = 0; A < this->_nNodes; A++ ) {
        const int idxA_u = nodeBlockSize * A;

        r_U = T_Boundary( A ) * f;

        {
          using namespace Eigen;
          P.template segment< nDim >( idxA_u ) -= Map< Matrix< double, nDim, 1 > >( r_U.data() );
        }

        for ( int B = 0; B < this->_nNodes; B++ ) {
          const int  idxB_u  = nodeBlockSize * B;
          const auto dN_B_dY = Tensor< double, nDim >( dN_dY.col( B ).data() );

          const Tensor< double, nDim, nDim > df_ddQU_B = T_Boundary( A ) * einsum< ijk, k >( df_dDeltaF, dN_B_dY );

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
  void DisplacementParticleSQCNIxSDI< nDim, nVertices >::computePhysicsKernelsOnSubdomains( const double* dQ,
                                                                                            double*       fInt,
                                                                                            double*       dFInt_ddQ,
                                                                                            double        timeNew,
                                                                                            double        dT )
  {
    using namespace Marmot::FastorIndices;
    using namespace Fastor;
    using to_jk = Fastor::OIndex< j_, k_ >;

    const static Tensor< double, nDim, nDim > I(
      ( Eigen::Matrix< double, nDim, nDim >() << Eigen::Matrix< double, nDim, nDim >::Identity() ).finished().data() );

    constexpr int nodeBlockSize = nDim;

    for ( size_t mpNumber = 0; mpNumber < this->_subDomainShapeFunctions.size(); mpNumber++ ) {
      auto& mp = this->_subdomainMaterialPoints[mpNumber];
      auto& sd = this->_subDomainShapeFunctions[mpNumber];

      Tensor< double, nDim > du( 0.0 );

      Tensor< double, nDim, nDim > du_dY( 0.0 );

      for ( int B = 0; B < this->_nNodes; B++ ) {

        const int idxB_u = nodeBlockSize * B;

        const double N_B     = sd.N( B );
        const auto   dN_B_dY = Tensor< double, nDim >( sd.dN_dY.col( B ).data() ); // works because ColumnMajor of Eigen

        const auto dQU = Tensor< double, nDim >( dQ + idxB_u );

        du += N_B * dQU;

        du_dY += einsum< i, j >( dQU, dN_B_dY );
      }

      mp->prepareYourself( timeNew, dT );
      mp->incrementDeformation( du, du_dY );
      mp->computeYourself( timeNew, dT );

      const double density0 = mp->getDensityUndeformed();

      auto v = mp->getVelocity();
      auto a = mp->getAcceleration();

      Tensor< double, nDim, nDim > da_ddu( 0.0 );
      Marmot::TimeIntegration::newmarkBetaIntegration< nDim >( du.data(),
                                                               v.data(),
                                                               a.data(),
                                                               dT,
                                                               this->_newmark_beta,
                                                               this->_newmark_gamma,
                                                               da_ddu.data() );
      mp->setVelocity( v );
      mp->setAcceleration( a );

      Tensor< double, nDim > r_U( 0.0 );

      Tensor< double, nDim, nDim > k_UU( 0.0 );

      const auto& S = mp->response.S;

      const double V0 = mp->getVolumeUndeformed();

      const auto& t = mp->tangents;

      Eigen::Map< Eigen::VectorXd > P( fInt, this->_nNodes * nodeBlockSize );
      Eigen::Map< Eigen::MatrixXd > K( dFInt_ddQ, this->_nNodes * nodeBlockSize, this->_nNodes * nodeBlockSize );

      // clang-format off
      for ( int A = 0; A < this->_nNodes; A++ ) {

        const double T_A = sd.T( A );
        const auto                   dT_A_dY = TensorMap< const double, nDim >( sd.dT_dY.col( A ).data() );
        const Tensor< double, nDim > dT_A_dx = einsum< ji, j >( inv( mp->dx_dY() ), dT_A_dY );

        const int idxA_u = nodeBlockSize * A;
        r_U = ( +einsum< i, ij >( dT_A_dx, S ) ) * V0;

        // add inertia
        r_U += density0 * a * T_A * V0;

        {
          using namespace Eigen;
          P.template segment< nDim >( idxA_u ) += Map< Matrix< double, nDim, 1 > >( r_U.data() );
        }

        for ( int B = 0; B < this->_nNodes; B++ ) {

          const int idxB_u = nodeBlockSize * B;

          const double                 N_B     = sd.N( B );
          const auto dN_B_dY = TensorMap< const double, nDim >( sd.dN_dY.col(B).data() );
          const auto dN_B_dx = evaluate( einsum< ji, j >( inv( mp->dx_dY() ), dN_B_dY ) );

          // aux stiffness tensors
          const auto dS_dqU_B = evaluate ( + einsum < ijkl, l > ( t.dS_dDeltaF, dN_B_dY ) );
          k_UU  = ( + einsum< i, ijk        > ( dT_A_dx, dS_dqU_B )                       ) * V0;

          k_UU += ( - einsum< k, ij, i, to_jk >( dT_A_dx, S, dN_B_dx ) ) * V0;

          k_UU += density0 * da_ddu * T_A * N_B * V0;

          {
              using namespace Eigen;
              // TODO: check if we can use transpose instead of torowmajor:
              K.template block< nDim, nDim >( idxA_u, idxB_u ) += Map< Matrix< double, nDim, nDim > >( torowmajor( k_UU ).data() );
          }
        }
      }

      // clang-format on
    }
  }

} // namespace Marmot::Meshfree
