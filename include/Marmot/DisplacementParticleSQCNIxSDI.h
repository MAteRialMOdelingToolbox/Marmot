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
#include "Marmot/GenericSDIParticle.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Fastor/Fastor.h>
#include <stdexcept>

namespace Marmot::Meshfree {

  template < int nDim, int nVertices >
  class DisplacementParticleSQCNIxSDI : public GenericSDIParticle< nDim, nVertices > {

    using TensorD  = GenericSDIParticle< nDim, nVertices >::TensorD;
    using TensorDD = GenericSDIParticle< nDim, nVertices >::TensorDD;
    std::vector< std::unique_ptr< MaterialPointType< nDim > > > _subdomainMaterialPoints;

    constexpr int static nStateVarsParticle = nDim * nVertices + nDim; // vertex displacements + center displacement

    double _newmark_beta;
    double _newmark_gamma;

    /// static vector of valid properties
    inline static const std::vector< std::string > _validProperties = {
      "newmark-beta beta",
      "newmark-beta gamma",
    };

  public:
    enum BodyLoadTypes {
      BodyForce,
    };

    enum DistributedLoadTypes { Pressure };

    const std::unordered_map< std::string, int >& getSupportedBodyLoadTypes() const override
    {
      static const std::unordered_map< std::string, int > _supportedBodyLoadTypes = { { "BODYFORCE", BodyForce } };
      return _supportedBodyLoadTypes;
    };

    const std::unordered_map< std::string, int >& getSupportedDistributedLoadTypes() const override
    {
      static const std::unordered_map< std::string, int > _supportedDistributedLoadTypes = { { "PRESSURE", Pressure } };
      return _supportedDistributedLoadTypes;
    };

    static constexpr int nDofPerNodeU = nDim; // Displacement   field U

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

    virtual int getNBaseDof() const override { return nDofPerNodeU; };

    void initializeYourselfOnSubdomains() override
    {
      for ( auto& mp : _subdomainMaterialPoints ) {
        mp->initializeYourself();
      }
    };

    virtual void computeBodyLoad( int           type,
                                  const double* load,
                                  double*       fExt,
                                  double*       dExt_dQ,
                                  double        timeNew,
                                  double        dT ) const override
    {
      throw std::runtime_error( "Not implemented yet!" );
    }

    virtual const std::vector< std::string >& getFields() const override
    {
      static const std::vector< std::string > nodeFields = { "displacement" };
      return nodeFields;
    };

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

    virtual double getVolumeUndeformed() const
    {

      double V0 = 0.0;
      for ( const auto& mp : _subdomainMaterialPoints ) {
        V0 += mp->getVolumeUndeformed();
      }
      return V0;
    }

    virtual double getVolumeDeformed() const
    {
      double volDeformed = 0.0;
      for ( const auto& mp : _subdomainMaterialPoints ) {
        volDeformed += mp->getVolumeUndeformed() * determinant( mp->dY_dX() );
      }
      return volDeformed;
    }

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

    virtual void acceptStateAndPositionOnSubdomains() override
    {
      for ( auto& mp : _subdomainMaterialPoints ) {
        mp->acceptStateAndPosition();
      }
    };

    virtual int getNumberOfRequiredStateVarsOnSubdomains() const override
    {
      int nStateVars = 0;

      for ( const auto& mp : _subdomainMaterialPoints ) {
        nStateVars += mp->getNumberOfRequiredStateVars();
      }

      return nStateVars;
    };

    virtual void assignStateVarsOnSubdomains( double* stateVars, int nStateVars ) override
    {
      int offset = 0;

      for ( auto& mp : _subdomainMaterialPoints ) {

        int nStateVarsSubParticle = mp->getNumberOfRequiredStateVars();
        mp->assignStateVars( stateVars + offset, nStateVarsSubParticle );
        offset += nStateVarsSubParticle;
      }

      if ( offset != nStateVars ) {
        throw std::runtime_error( "Error: Number of state variables does not match!" );
      }
    }

    virtual StateView getStateViewOnSubdomains( const std::string& stateName, int subdomainIndex ) const
    {
      return _subdomainMaterialPoints[subdomainIndex]->getStateView( stateName );
    }

    virtual void computePhysicsKernelsOnSubdomains( const double* dQ,
                                                    double*       fInt,
                                                    double*       dFInt_ddQ,
                                                    double        timeNew,
                                                    double        dT );

    virtual void computeDistributedLoad( int           type,
                                         int           surfaceID,
                                         const double* load,
                                         double*       fExt,
                                         double*       dExt_dQ,
                                         double        timeNew,
                                         double        dT ) const override;

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

    virtual int getNumberOfEvaluationPoints() const
    {

      int nEvalPoints = 0;
      for ( int i = 0; i < this->_particleDomainMain.getNumberOfFaces(); i++ ) {
        const auto subcellsAttachedToFace = this->_particleDomainMain.getSubCellIndicesOnParentFace( i + 1 );
        nEvalPoints += static_cast< int >( subcellsAttachedToFace.size() );
      }
      return nEvalPoints;
    };

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
