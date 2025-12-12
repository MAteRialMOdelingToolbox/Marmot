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
#include "Marmot/MarmotParticle.h"
#include "Marmot/MarmotParticleDomain.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Fastor/Fastor.h>
#include <stdexcept>

namespace Marmot::Meshfree {

  template < int nDim, int nVertices >
  class DisplacementParticleSQCNIxSDI : public MarmotParticle {

    using TensorD                = Fastor::Tensor< double, nDim >;
    using TensorDD               = Fastor::Tensor< double, nDim, nDim >;
    using VertexCoordinatesSized = Eigen::Matrix< double, nDim, nVertices >;
    using CoordinatesSized       = Eigen::Matrix< double, nDim, 1 >;

    struct SubDomain {

      ParticleDomain< nDim, nVertices > particleDomain;

      std::unique_ptr< MaterialPointType< nDim > > materialPoint;

      Eigen::MatrixXd N;
      Eigen::MatrixXd dN_dY;

      Eigen::MatrixXd T;
      Eigen::MatrixXd dT_dY;

      Eigen::VectorXd P;
      Eigen::MatrixXd P_Gradient;
    };

    ParticleDomain< nDim, nVertices > _particleDomainMain;
    std::vector< SubDomain >          _subDomains;

  public:
    using SmoothingDomainUpdateType = ParticleDomain< nDim, nVertices >::SmoothingDomainUpdateType;
    // enum SmoothingVolumeUpdateType { None, DeformationGradient, RotationOnly, RotationAndPrincipalStretch };

  private:
    constexpr int static nStateVarsParticle = nDim * nVertices + nDim; // vertex displacements + center displacement

    int    _elementID;
    double _newmark_beta;
    double _newmark_gamma;
    int    _nNodes;
    int    _vciOrder;
    int    _nVCIConstraints;

    const MarmotMeshfreeApproximation& _meshfreeApproximation;
    using JacobianSized = Eigen::Matrix< double, nDim, nDim >;

    Eigen::Map< CoordinatesSized > _centerDisplacement;
    Eigen::Map< JacobianSized >    _centralDeformationGradient;

    TensorDD _dx_dY_center;
    TensorD  _du_center;

    std::vector< const MarmotMeshfreeKernelFunction* > _assignedKernelFunctions;

    /// static vector of valid properties
    inline static const std::vector< std::string > _validProperties = {
      "newmark-beta beta",
      "newmark-beta gamma",
      "VCI order",
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

    using Material = MarmotMaterialFiniteStrain;

    using ForceSized = Eigen::Matrix< double, nDim, 1 >;

    virtual void setProperties( const double* properties, int nProperties ) override
    {
      if ( nProperties != static_cast< int >( _validProperties.size() ) ) {
        throw std::runtime_error( "Number of properties does not match!" );
      }

      for ( int i = 0; i < nProperties; i++ ) {
        setProperty( _validProperties[i], &properties[i] );
      }
    };

    virtual void setProperty( const std::string& propertyName, const double* property )
    {
      if ( propertyName == "newmark-beta beta" ) {
        _newmark_beta = property[0];
      }
      else if ( propertyName == "newmark-beta gamma" ) {
        _newmark_gamma = property[0];
      }
      else if ( propertyName == "VCI order" ) {
        _vciOrder = static_cast< int >( property[0] );
        this->setVCIOrder( _vciOrder );
      }
      else {
        throw std::runtime_error( "Property " + propertyName + " not supported!" );
      }
    };

    /// \brief Get the names of the properties
    /// \return The names of the properties
    virtual std::vector< std::string > getPropertyNames() const { return _validProperties; };

    virtual void getVertexCoordinates( double* coordinates ) const override;

    virtual void getFaceCoordinates( int faceID, double* coordinates ) const
    {
      Eigen::Map< Eigen::Matrix< double, nDim, 1 > > coordinatesMap( coordinates );
      coordinatesMap = _particleDomainMain.getFaceCenterCoordinates( faceID );
    }

    virtual void getCenterCoordinates( double* coordinates ) const override
    {
      Eigen::Map< Eigen::Matrix< double, nDim, 1 > > centerCoordinatesMap( coordinates );
      centerCoordinatesMap = _particleDomainMain.getCenterCoordinates();
    }

    virtual void getVisualizationVertexCoordinates( double* coordinates ) const override
    {
      getVertexCoordinates( coordinates );
    };

    virtual int getNumberOfVertices() const override { return nVertices; };

    virtual int getNBaseDof() const override { return nDofPerNodeU; }

    void initializeYourself() override
    {
      for ( auto& sd : _subDomains ) {
        sd.materialPoint->initializeYourself();
      }

      _centerDisplacement.setZero();
      _centralDeformationGradient.setIdentity();

      _dx_dY_center.eye();
      _du_center.zeros();
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
    virtual std::string getParticleShape() const override { return _particleDomainMain.getParticleShape(); };

    DisplacementParticleSQCNIxSDI( int                                elementID,
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

    virtual double getVolumeUndeformed() const
    {

      double V0 = 0.0;
      for ( const auto& sd : _subDomains ) {
        V0 += sd.materialPoint->getVolumeUndeformed();
      }
      return V0;
    }

    virtual double getVolumeDeformed() const
    {
      double volDeformed = 0.0;
      for ( const auto& sd : _subDomains ) {
        volDeformed += sd.materialPoint->getVolumeUndeformed() * determinant( sd.materialPoint->dY_dX() );
      }
      return volDeformed;
    }

    // That goes to the generic sdi particle
    virtual void acceptStateAndPosition() override
    {
      _centerDisplacement += Eigen::Matrix< double, nDim, 1 >( _du_center.data() );
      _du_center.zeros();

      Eigen::Map< Eigen::Matrix< double, nDim, nDim, Eigen::RowMajor > > dx_dY_map( _dx_dY_center.data() );
      _centralDeformationGradient = dx_dY_map * _centralDeformationGradient;
      _dx_dY_center.eye();

      _particleDomainMain.acceptStateAndPosition( _centralDeformationGradient, _centerDisplacement );

      // make that cleaner
      const auto newParticleDomains = _particleDomainMain.uniformSubdivided( 1 );
      for ( size_t i = 0; i < _subDomains.size(); i++ ) {
        _subDomains[i].particleDomain = newParticleDomains[i];
      }

      for ( auto& sd : _subDomains ) {
        sd.materialPoint->acceptStateAndPosition();
      }
    };

    virtual int getNumberOfRequiredStateVars() const override
    {
      int nStateVars = 0;

      nStateVars += nDim;        // center displacement
      nStateVars += nDim * nDim; // central deformation gradient

      for ( const auto& sd : _subDomains ) {
        nStateVars += sd.materialPoint->getNumberOfRequiredStateVars();
      }

      return nStateVars;
    };

    void assignStateVars( double* stateVars, int nStateVars ) override
    {

      int offset = 0;

      new ( &_centerDisplacement ) Eigen::Map< CoordinatesSized >( stateVars + offset );
      offset += nDim;

      new ( &_centralDeformationGradient ) Eigen::Map< JacobianSized >( stateVars + offset );
      offset += nDim * nDim;

      for ( auto& sd : _subDomains ) {
        int nStateVarsSubParticle = sd.materialPoint->getNumberOfRequiredStateVars();
        sd.materialPoint->assignStateVars( stateVars + offset, nStateVarsSubParticle );
        offset += nStateVarsSubParticle;
      }

      if ( offset != nStateVars ) {
        throw std::runtime_error( "Error: Number of state variables does not match!" );
      }
    }

    virtual StateView getStateView( const std::string& stateName, int qp ) const override
    {
      if ( stateName == "vertex displacements" )
        return StateView( const_cast< double* >( _particleDomainMain.getSmoothingDomainVertexDisplacements().data() ),
                          nDim * nVertices );

      return _subDomains[qp].materialPoint->getStateView( stateName );
    }

    virtual void computePhysicsKernels( const double* dQ,
                                        double*       fInt,
                                        double*       dFInt_ddQ,
                                        double        timeNew,
                                        double        dT ) override;

    virtual void computeDistributedLoad( int           type,
                                         int           surfaceID,
                                         const double* load,
                                         double*       fExt,
                                         double*       dExt_dQ,
                                         double        timeNew,
                                         double        dT ) const override;

    std::tuple< TensorD, TensorD > getIntermediateConfBoundaryVector( int boundaryFaceID ) const;

    virtual int getDimension() const override { return nDim; };

    virtual void getInterpolationVector( double* vec, const double* coordinates ) const override
    {
      _meshfreeApproximation.computeShapeFunctions( coordinates, _assignedKernelFunctions, vec );
    };

    virtual int vci_getNumberOfConstraints() override { return _nVCIConstraints; }

    virtual void vci_compute_Test_P_BoundaryIntegral( double*       R_AiC_RowMajor,
                                                      const double* boundarySurfaceVector,
                                                      int           boundaryFaceID ){

      // using namespace Fastor;

      // const auto [N_dAY, Y_N]   = getIntermediateConfBoundaryVector( boundaryFaceID );
      // Eigen::MatrixXd TBoundary = Eigen::MatrixXd::Zero( 1, _nNodes );

      // _meshfreeApproximation.computeShapeFunctions( Y_N.data(), _assignedKernelFunctions, TBoundary.data() );

      // // get P for the exact integration location at the boundary
      // auto PBoundary = Eigen::VectorXd( _nVCIConstraints );
      // Math::computeMonomialBasis( _vciOrder, CoordinatesSized( Y_N.data() ), PBoundary );

      // for ( int A = 0; A < _nNodes; A++ )
      //   for ( int i = 0; i < nDim; i++ )
      //     for ( int C = 0; C < _nVCIConstraints; C++ )
      //       R_AiC_RowMajor[A * ( nDim * _nVCIConstraints ) + i * _nVCIConstraints + C] += TBoundary( A ) *
      //                                                                                     PBoundary( C ) * N_dAY[i];
    };

    virtual void vci_compute_TestGradient_P_Integral( double* R_AiC_RowMajor ) override{
      // dimensions of R_AiC_RowMajor: _nNodes x nDim x _nVCIConstraints
      // for ( const auto& sd : _subDomains)
      //   for ( int A = 0; A < _nNodes; A++ )
      //     for ( int i = 0; i < nDim; i++ )
      //       for ( int C = 0; C < _nVCIConstraints; C++ )
      //         R_AiC_RowMajor[A * ( nDim * _nVCIConstraints ) + i * _nVCIConstraints + C] += sd.dT_dY( i, A ) *
      //                                                                                       sd.P( C ) *
      //                                                                                       sd.V_IntermediateReference;
    };

    virtual void vci_compute_Test_PGradient_Integral( double* R_AiC_RowMajor ) override{
      // dimensions of R_AiC_RowMajor: _nNodes x nDim x _nVCIConstraints
      // for ( const auto& sd : _subDomains)
      //   for ( int A = 0; A < _nNodes; A++ )
      //     for ( int i = 0; i < nDim; i++ )
      //       for ( int C = 0; C < _nVCIConstraints; C++ )
      //         R_AiC_RowMajor[A * ( nDim * _nVCIConstraints ) + i * _nVCIConstraints + C] += sd.T( A ) *
      //                                                                                       sd.P_Gradient( C, i ) *
      //                                                                                       sd.V_IntermediateReference;
    };

    virtual void vci_compute_MMatrix( double* mMatrix_ACD_RowMajor ) override{
      // dimensions of R_AiC_RowMajor: _nNodes x nDim x _nVCIConstraints

      // for ( const auto& sd : _subDomains)
      //   for ( int A = 0; A < _nNodes; A++ ) {
      //     const double R_A = _assignedKernelFunctions[A]->isInSupport( sd.center_IntermediateReference.data() ) ? 1.0
      //                                                                                                           :
      //                                                                                                           0.0;
      //     // const double R_A = 1.0;

      //     for ( int C = 0; C < _nVCIConstraints; C++ )
      //       for ( int D = 0; D < _nVCIConstraints; D++ )
      //         mMatrix_ACD_RowMajor[A * ( _nVCIConstraints * _nVCIConstraints ) + C * _nVCIConstraints +
      //                              D] += R_A * sd.P( C ) * sd.P( D ) * sd.V_IntermediateReference;
      //   }
    };

    virtual void vci_assignTestFunctionCorrectionTerms( const double* eta_AiC_RowMajor ) override{

      // for ( auto& sd : _subDomains)
      //   for ( int A = 0; A < _nNodes; A++ ) {
      //     const double R_A = _assignedKernelFunctions[A]->isInSupport( sd.center_IntermediateReference.data() ) ? 1.0
      //                                                                                                           :
      //                                                                                                           0.0;
      //     // const double R_A = 1.0;
      //     for ( int i = 0; i < nDim; i++ ) {
      //       for ( int C = 0; C < _nVCIConstraints; C++ ) {
      //         sd.dT_dY( i, A ) += eta_AiC_RowMajor[A * ( nDim * _nVCIConstraints ) + i * _nVCIConstraints + C] * R_A
      //         *
      //                             sd.P( C );
      //       }
      //     }
      //   }
    };

    virtual void getEvaluationCoordinates( double* coordinates ) const override
    {
      Eigen::Map< Eigen::Matrix< double, nDim, Eigen::Dynamic > >
        coordinatesMap( coordinates, nDim, _particleDomainMain.getNumberOfFaces() );

      for ( int i = 0; i < _particleDomainMain.getNumberOfFaces(); i++ ) {
        coordinatesMap.col( i ) = _particleDomainMain.getSmoothingDomainFaceCenterCoordinates( i + 1 );
      }
    }

    virtual int getNumberOfEvaluationPoints() const { return _particleDomainMain.getNumberOfFaces(); };

    virtual void setInitialCondition( const std::string& conditionName, const double* value ) override
    {
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid initial condition" );
    };

  private:
    /// \brief Evaluate the shape functions for a vertex-shaped domain
    /// \details This function evaluates the shape functions for a vertex-shaped
    ///         domain using the vertex coordinates of the particle.
    ///         \param particleDomain The particle domain for which to evaluate shape functions.
    ///         \return The shape functions (at center) and their gradients computed from smoothing around the domain.
    std::tuple< Eigen::MatrixXd, Eigen::MatrixXd > evaluateShapeFunctionsForParticleDomain(
      const ParticleDomain< nDim, nVertices >& particleDomain ) const;

    void setVCIOrder( int order )
    {
      _vciOrder        = order;
      _nVCIConstraints = ( order + 1 ) * ( order + 2 ) / 2; // number of VCI constraints for polynomial basis of order
    };
  };

  template < int nDim, int nVertices >
  DisplacementParticleSQCNIxSDI< nDim, nVertices >::DisplacementParticleSQCNIxSDI(
    int                                                  elementID,
    const double*                                        vertexCoordinates,
    int                                                  nVertexCoordinates,
    double                                               volume,
    const std::string&                                   materialName,
    const double*                                        materialProperties,
    int                                                  sizeMaterialProperties,
    const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation,
    const SmoothingDomainUpdateType                      smoothingVolumeUpdateType )
    : _elementID( elementID ),
      _newmark_beta( 0. ),
      _newmark_gamma( 0. ),
      _meshfreeApproximation( approximation ),
      _vciOrder( 0 ),
      _particleDomainMain( vertexCoordinates, nVertexCoordinates, smoothingVolumeUpdateType ),
      _centerDisplacement( nullptr ),
      _centralDeformationGradient( nullptr )
  {
    _dx_dY_center.eye();
    _du_center.zeros();

    const auto initialParticleDomains = _particleDomainMain.uniformSubdivided( 1 );
    for ( size_t i = 0; i < initialParticleDomains.size(); i++ ) {

      const auto& initialParticleDomain = initialParticleDomains[i];

      const double subV0 = initialParticleDomain.getVolumeUndeformed();
      const auto   X0    = initialParticleDomain.getCenterCoordinates();

      _subDomains.push_back(
        SubDomain{ .particleDomain = initialParticleDomain,
                   .materialPoint = std::make_unique< MaterialPointType< nDim > >( elementID, X0.data(), 1, subV0 ) } );
    }

    int                   materialCode = MarmotLibrary::MarmotMaterialFactory::getMaterialCodeFromName( materialName );
    MarmotMaterialSection section( materialCode, materialProperties, sizeMaterialProperties );
    for ( auto& sd : _subDomains ) {
      sd.materialPoint->assignMaterial( section );
    }

    this->setVCIOrder( _vciOrder );
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

      const auto [testBoundary, dN_dY] = evaluateShapeFunctionsForParticleDomain( _particleDomainMain );
      const auto dT_dY                 = dN_dY;

      const auto [N_dAY, Y_N] = getIntermediateConfBoundaryVector( boundaryFaceID );

      TensorD fY = N_dAY * load_[0];

      Eigen::Map< Eigen::VectorXd > P( fExt, _nNodes * nodeBlockSize );
      Eigen::Map< Eigen::MatrixXd > K( dFExt_ddQ, _nNodes * nodeBlockSize, _nNodes * nodeBlockSize );

      using namespace Fastor;
      using namespace FastorIndices;

      Tensor< double, nDim, nDim > Eye;
      Eye.eye();

      // apply Nanson's formula
      // const auto deltaF = _mp->dx_dY(); // TODO replace with

      const Tensor< double, nDim, nDim > deltaFInv = inverse( _dx_dY_center );
      const double                       deltaJ    = determinant( _dx_dY_center );

      const TensorD f = deltaJ * transpose( deltaFInv ) % fY;

      const Tensor< double, nDim, nDim, nDim, nDim > dFInv_dF = -einsum< Ik, Ki, to_IikK >( deltaFInv, deltaFInv );

      const Tensor< double, nDim, nDim, nDim > df_dDeltaF = outer( f, transpose( deltaFInv ) ) +
                                                            deltaJ * einsum< IikK, Index< I_ > >( dFInv_dF, fY );

      TensorD r_U( 0.0 );

      for ( int A = 0; A < _nNodes; A++ ) {
        const int idxA_u = nodeBlockSize * A;

        r_U = testBoundary( A ) * f;

        {
          using namespace Eigen;
          P.template segment< nDim >( idxA_u ) -= Map< Matrix< double, nDim, 1 > >( r_U.data() );
        }

        for ( int B = 0; B < _nNodes; B++ ) {
          const int  idxB_u  = nodeBlockSize * B;
          const auto dN_B_dY = TensorMap< const double, nDim >( dN_dY.col( B ).data() );

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

  // that will go to the generic sdi particle
  template < int nDim, int nVertices >
  std::tuple< typename DisplacementParticleSQCNIxSDI< nDim, nVertices >::TensorD,
              typename DisplacementParticleSQCNIxSDI< nDim, nVertices >::TensorD >
  DisplacementParticleSQCNIxSDI< nDim, nVertices >::getIntermediateConfBoundaryVector( int boundaryFaceID ) const
  {

    TensorD N_dAY;
    TensorD Y;

    Eigen::Matrix< double, nDim, 1 > _Y_eigen;

    // the evaluation point depends: For real SQCNI, we do it on the smoothing domain boundary, for all others we do in
    // in the center of the deformed geometry

    if ( _particleDomainMain.smoothingVolumeUpdateType == SmoothingDomainUpdateType::DeformationGradient )
      _Y_eigen = _particleDomainMain.getSmoothingDomainFaceCenterCoordinates( boundaryFaceID );
    else
      throw std::invalid_argument( "not implemented" );

    // N_dAY (boundary surface vector for distributed load) comes from the deformed geometry
    auto _N_dAY_eigen = _particleDomainMain.getFaceBoundaryVector( boundaryFaceID );

    for ( int i = 0; i < nDim; i++ ) {
      N_dAY[i] = _N_dAY_eigen[i];
      Y[i]     = _Y_eigen[i];
    }

    return { N_dAY, Y };
  }

  // that will later go to general sdi particle
  template < int nDim, int nVertices >
  void DisplacementParticleSQCNIxSDI< nDim, nVertices >::getVertexCoordinates( double* coordinates ) const
  {
    Eigen::Map< Eigen::Matrix< double, nDim, nVertices > > coordinatesMap( coordinates );
    coordinatesMap = _particleDomainMain.getGeometryDeformedVertexCoordinates();
  }

  /// That will later go to general subdomain particle.
  template < int nDim, int nVertices >
  void DisplacementParticleSQCNIxSDI< nDim, nVertices >::assignMeshfreeKernelFunctions(
    const std::vector< const MarmotMeshfreeKernelFunction* >& kernelFunctions )
  {

    _assignedKernelFunctions = kernelFunctions;
    _nNodes                  = kernelFunctions.size();

    for ( size_t mpNumber = 0; mpNumber < _subDomains.size(); mpNumber++ ) {

      auto& sd = _subDomains[mpNumber];

      const auto [N, dN_dY] = evaluateShapeFunctionsForParticleDomain( sd.particleDomain );

      CoordinatesSized mpCenter;
      sd.materialPoint->getCoordinatesAtCenter( mpCenter.data() );

      sd.N     = N;
      sd.dN_dY = dN_dY;

      sd.T     = N;
      sd.dT_dY = dN_dY;

      sd.P.resize( _nVCIConstraints );
      sd.P_Gradient.resize( _nVCIConstraints, nDim );

      Math::computeMonomialBasis( _vciOrder, mpCenter, sd.P );
      Math::computeMonomialBasisGradient( _vciOrder, mpCenter, sd.P_Gradient );
    }
  }

  template < int nDim, int nVertices >
  std::tuple< Eigen::MatrixXd, Eigen::MatrixXd > DisplacementParticleSQCNIxSDI< nDim, nVertices >::
    evaluateShapeFunctionsForParticleDomain( const ParticleDomain< nDim, nVertices >& particleDomain ) const

  {

    Eigen::Matrix< double, nDim, 1 > coords;
    // Use the center of the provided particleDomain
    coords = particleDomain.getCenterCoordinates();

    Eigen::MatrixXd N     = Eigen::MatrixXd::Zero( 1, this->_nNodes );
    Eigen::MatrixXd dN_dY = Eigen::MatrixXd::Zero( nDim, this->_nNodes );

    // Compute N at the particle center (from GenericParticle)
    this->_meshfreeApproximation.computeShapeFunctions( coords.data(), this->_assignedKernelFunctions, N.data() );

    // Compute dN_dY using the SQCNI approach (boundary integral over smoothing domain)
    Eigen::MatrixXd smooth_NBoundary( 1, this->_nNodes );
    for ( int i = 0; i < particleDomain.getNumberOfFaces(); i++ ) {

      auto smoothing_evaluation_point = particleDomain.getSmoothingDomainFaceCenterCoordinates( i + 1 );
      auto smoothing_n_dA             = particleDomain.getSmoothingBoundarySurfaceVector( i + 1 );

      this->_meshfreeApproximation.computeShapeFunctions( smoothing_evaluation_point.data(),
                                                          this->_assignedKernelFunctions,
                                                          smooth_NBoundary.data() );
      dN_dY += smoothing_n_dA * smooth_NBoundary;
    }
    dN_dY /= particleDomain.getSmoothingVolume();

    return std::make_tuple( N, dN_dY );
  }

  template < int nDim, int nVertices >
  void DisplacementParticleSQCNIxSDI< nDim, nVertices >::computePhysicsKernels( const double* dQ,
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

    // update central deformation and displacement.
    _dx_dY_center.eye();
    _du_center.zeros();
    {
      const auto [N, dN_dY] = evaluateShapeFunctionsForParticleDomain( _particleDomainMain );

      Tensor< double, nDim, nDim > du_dY( 0.0 );
      Tensor< double, nDim >       du( 0.0 );

      for ( int B = 0; B < _nNodes; B++ ) {

        const int idxB_u = nodeBlockSize * B;

        const auto dN_B_dY = Tensor< double, nDim >( dN_dY.col( B ).data() ); // works because ColumnMajor of Eigen

        const auto dQU = Tensor< double, nDim >( dQ + idxB_u );

        du_dY += einsum< i, j >( dQU, dN_B_dY );
        du += N( B ) * dQU;
      }
      _dx_dY_center += du_dY;
      _du_center += du;
    }

    for ( auto& sd : _subDomains ) {

      Tensor< double, nDim > du( 0.0 );

      Tensor< double, nDim, nDim > du_dY( 0.0 );

      for ( int B = 0; B < _nNodes; B++ ) {

        const int idxB_u = nodeBlockSize * B;

        const double N_B     = sd.N( B );
        const auto   dN_B_dY = Tensor< double, nDim >( sd.dN_dY.col( B ).data() ); // works because ColumnMajor of Eigen

        const auto dQU = Tensor< double, nDim >( dQ + idxB_u );

        du += N_B * dQU;

        du_dY += einsum< i, j >( dQU, dN_B_dY );
      }

      sd.materialPoint->prepareYourself( timeNew, dT );
      sd.materialPoint->incrementDeformation( du, du_dY );
      sd.materialPoint->computeYourself( timeNew, dT );

      const double density0 = sd.materialPoint->getDensityUndeformed();

      auto v = sd.materialPoint->getVelocity();
      auto a = sd.materialPoint->getAcceleration();

      Tensor< double, nDim, nDim > da_ddu( 0.0 );
      Marmot::TimeIntegration::newmarkBetaIntegration< nDim >( du.data(),
                                                               v.data(),
                                                               a.data(),
                                                               dT,
                                                               this->_newmark_beta,
                                                               this->_newmark_gamma,
                                                               da_ddu.data() );
      sd.materialPoint->setVelocity( v );
      sd.materialPoint->setAcceleration( a );

      Tensor< double, nDim > r_U( 0.0 );

      Tensor< double, nDim, nDim > k_UU( 0.0 );

      const auto& S = sd.materialPoint->response.S;

      const double V0 = sd.materialPoint->getVolumeUndeformed();

      const auto& t = sd.materialPoint->tangents;

      Eigen::Map< Eigen::VectorXd > P( fInt, _nNodes * nodeBlockSize );
      Eigen::Map< Eigen::MatrixXd > K( dFInt_ddQ, _nNodes * nodeBlockSize, _nNodes * nodeBlockSize );

      // clang-format off
      for ( int A = 0; A < _nNodes; A++ ) {

        const double T_A = sd.T( A );
        const auto                   dT_A_dY = TensorMap< const double, nDim >( sd.dT_dY.col( A ).data() );
        const Tensor< double, nDim > dT_A_dx = einsum< ji, j >( inv( sd.materialPoint->dx_dY() ), dT_A_dY );

        const int idxA_u = nodeBlockSize * A;
        r_U = ( +einsum< i, ij >( dT_A_dx, S ) ) * V0;

        // add inertia
        r_U += density0 * a * T_A * V0;

        {
          using namespace Eigen;
          P.template segment< nDim >( idxA_u ) += Map< Matrix< double, nDim, 1 > >( r_U.data() );
        }

        for ( int B = 0; B < _nNodes; B++ ) {

          const int idxB_u = nodeBlockSize * B;

          const double                 N_B     = sd.N( B );
          const auto dN_B_dY = TensorMap< const double, nDim >( sd.dN_dY.col(B).data() );
          const auto dN_B_dx = evaluate( einsum< ji, j >( inv( sd.materialPoint->dx_dY() ), dN_B_dY ) );

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
