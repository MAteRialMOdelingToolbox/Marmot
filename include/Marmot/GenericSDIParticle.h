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

#include "Marmot/MarmotMeshfreeApproximation.h"
#include "Marmot/MarmotMonomialBasisFunctions.h"
#include "Marmot/MarmotParticle.h"
#include "Marmot/MarmotParticleDomain.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Fastor/Fastor.h>
#include <stdexcept>

namespace Marmot::Meshfree {

  template < int nDim, int nVertices >
  class GenericSDIParticle : public MarmotParticle {

  protected:
    using TensorD                = Fastor::Tensor< double, nDim >;
    using TensorDD               = Fastor::Tensor< double, nDim, nDim >;
    using VertexCoordinatesSized = Eigen::Matrix< double, nDim, nVertices >;
    using CoordinatesSized       = Eigen::Matrix< double, nDim, 1 >;
    using JacobianSized          = Eigen::Matrix< double, nDim, nDim >;

    struct SubDomainShapeFunctions {
      Eigen::MatrixXd N;
      Eigen::MatrixXd dN_dY;

      Eigen::MatrixXd T;
      Eigen::MatrixXd dT_dY;

      Eigen::VectorXd P;
      Eigen::MatrixXd P_Gradient;
    };

    constexpr int static nStateVarsParticle = nDim * nVertices + nDim; // vertex displacements + center displacement

    int _elementID;
    int _nNodes;
    int _vciOrder;
    int _nVCIConstraints;

    const MarmotMeshfreeApproximation& _meshfreeApproximation;

    // goest to the generic sdi particle
    ParticleDomain< nDim, nVertices > _particleDomainMain;
    // goes to the generic sdi particle
    std::vector< SubDomainShapeFunctions > _subDomainShapeFunctions;
    // gos to the generic sdi particle
    std::vector< ParticleDomain< nDim, nVertices > > _subDomains;

    Eigen::Map< CoordinatesSized > _centerDisplacement;
    Eigen::Map< JacobianSized >    _centralDeformationGradient;

    TensorDD _dx_dY_center;
    TensorD  _du_center;

    std::vector< const MarmotMeshfreeKernelFunction* > _assignedKernelFunctions;

    /// static vector of valid properties
    inline static const std::vector< std::string > _validProperties = {
      "VCI order",
    };

  public:
    using SmoothingDomainUpdateType = ParticleDomain< nDim, nVertices >::SmoothingDomainUpdateType;

    virtual void setProperties( const double* properties, int nProperties ) override {};

    virtual void setProperty( const std::string& propertyName, const double* property )
    {
      if ( propertyName == "VCI order" ) {
        _vciOrder = static_cast< int >( property[0] );
        this->setVCIOrder( _vciOrder );
      }

      setPropertyOnSubdomains( propertyName, property );
    };

    virtual void setPropertyOnSubdomains( const std::string& propertyName, const double* property ) = 0;

    /// \brief Get the names of the properties
    /// \return The names of the properties
    virtual std::vector< std::string > getPropertyNames() const
    {
      auto validProperties = _validProperties;
      validProperties.insert( validProperties.end(),
                              getSubdomainPropertyNames().begin(),
                              getSubdomainPropertyNames().end() );
      return validProperties;
    };

    virtual std::vector< std::string > getSubdomainPropertyNames() const = 0;

    virtual void getVertexCoordinates( double* coordinates ) const override;

    virtual void getFaceCoordinates( int faceID, double* coordinates ) const override final
    {
      Eigen::Map< Eigen::Matrix< double, nDim, 1 > > coordinatesMap( coordinates );
      coordinatesMap = _particleDomainMain.getFaceCenterCoordinates( faceID );
    }

    virtual void getCenterCoordinates( double* coordinates ) const override final
    {
      Eigen::Map< Eigen::Matrix< double, nDim, 1 > > centerCoordinatesMap( coordinates );
      centerCoordinatesMap = _particleDomainMain.getCenterCoordinates();
    }

    virtual void getVisualizationVertexCoordinates( double* coordinates ) const override
    {
      getVertexCoordinates( coordinates );
    };

    virtual int getNumberOfVertices() const override final { return nVertices; };

    virtual double getSubdomainVolume( const ParticleDomain< nDim, nVertices >& subdomain ) const = 0;

    void initializeYourself() override
    {

      _centerDisplacement.setZero();
      _centralDeformationGradient.setIdentity();

      _dx_dY_center.eye();
      _du_center.zeros();

      initializeYourselfOnSubdomains();
    };

    virtual void initializeYourselfOnSubdomains() = 0;

    virtual std::string getParticleShape() const override final { return _particleDomainMain.getParticleShape(); };

    GenericSDIParticle( int                                elementID,
                        const double*                      nodeCoordinates,
                        int                                nNodeCoordiantes,
                        double                             volume,
                        const MarmotMeshfreeApproximation& approximation,
                        const SmoothingDomainUpdateType    smoothingVolumeUpdateType );

    virtual void assignMeshfreeKernelFunctions(
      const std::vector< const MarmotMeshfreeKernelFunction* >& kernelFunctions ) override;

    virtual void acceptStateAndPosition() override
    {
      _centerDisplacement += Eigen::Matrix< double, nDim, 1 >( _du_center.data() );
      _du_center.zeros();

      Eigen::Map< Eigen::Matrix< double, nDim, nDim, Eigen::RowMajor > > dx_dY_map( _dx_dY_center.data() );
      _centralDeformationGradient = dx_dY_map * _centralDeformationGradient;
      _dx_dY_center.eye();

      _particleDomainMain.acceptStateAndPosition( _centralDeformationGradient, _centerDisplacement );

      _subDomains = _particleDomainMain.uniformSubdivided();

      acceptStateAndPositionOnSubdomains();
    };

    virtual void acceptStateAndPositionOnSubdomains() = 0;

    virtual int getNumberOfRequiredStateVars() const override
    {
      int nStateVars = 0;

      nStateVars += nDim;        // center displacement
      nStateVars += nDim * nDim; // central deformation gradient
                                 //

      nStateVars += getNumberOfRequiredStateVarsOnSubdomains();

      return nStateVars;
    };

    virtual int getNumberOfRequiredStateVarsOnSubdomains() const = 0;

    void assignStateVars( double* stateVars, int nStateVars ) override final
    {

      int offset = 0;

      new ( &_centerDisplacement ) Eigen::Map< CoordinatesSized >( stateVars + offset );
      offset += nDim;

      new ( &_centralDeformationGradient ) Eigen::Map< JacobianSized >( stateVars + offset );
      offset += nDim * nDim;

      if ( offset > nStateVars ) {
        throw std::runtime_error( "Error: Number of state variables does not match!" );
      }

      assignStateVarsOnSubdomains( stateVars + offset, nStateVars - offset );
    }

    virtual void assignStateVarsOnSubdomains( double* stateVars, int nStateVars ) = 0;

    virtual StateView getStateView( const std::string& stateName, int qp ) const override final
    {
      if ( stateName == "vertex displacements" )
        return StateView( const_cast< double* >( _particleDomainMain.getSmoothingDomainVertexDisplacements().data() ),
                          nDim * nVertices );

      return getStateViewOnSubdomains( stateName, qp );
    }

    virtual StateView getStateViewOnSubdomains( const std::string& stateName, int subdomainIndex ) const = 0;

    virtual void computePhysicsKernels( const double* dQ,
                                        double*       fInt,
                                        double*       dFInt_ddQ,
                                        double        timeNew,
                                        double        dT ) override final;

    virtual void computePhysicsKernelsOnSubdomains( const double* dQ,
                                                    double*       fInt,
                                                    double*       dFInt_ddQ,
                                                    double        timeNew,
                                                    double        dT ) = 0;

    std::tuple< TensorD, TensorD > getIntermediateConfigurationBoundaryVector(
      int                                      boundaryFaceID,
      const ParticleDomain< nDim, nVertices >& particleDomain ) const;

    virtual int getDimension() const override final { return nDim; };

    virtual void getInterpolationVector( double* vec, const double* coordinates ) const override final
    {
      _meshfreeApproximation.computeShapeFunctions( coordinates, _assignedKernelFunctions, vec );
    };

    virtual int vci_getNumberOfConstraints() override final { return _nVCIConstraints; }

    virtual void vci_compute_Test_P_BoundaryIntegral( double*       R_AiC_RowMajor,
                                                      const double* boundarySurfaceVector,
                                                      int           boundaryFaceID ) {

      // using namespace Fastor;

      // const auto [N_dAY, Y_N]   = getIntermediateConfigurationBoundaryVector( boundaryFaceID );
      // Eigen::MatrixXd TBoundary = Eigen::MatrixXd::Zero( 1, _nNodes );

      // _meshfreeApproximation.computeShapeFunctions( Y_N.data(), _assignedKernelFunctions, TBoundary.data() );

      // // get P for the exact integration location at the boundary
      // auto PBoundary = Eigen::VectorXd( _nVCIConstraints );
      // Math::computeMonomialBasis( _vciOrder, CoordinatesSized( Y_N.data() ), PBoundary );

      // for ( int A = 0; A < _nNodes; A++ )
      //   for ( int i = 0; i < nDim; i++ )
      //     for ( int C = 0; C < _nVCIConstraints; C++ )
      //       R_AiC_RowMajor[A * ( nDim * _nVCIConstraints ) + i * _nVCIConstraints + C] += TBoundary( A ) * PBoundary(
      //       C ) * N_dAY[i];

    };

    virtual void vci_compute_TestGradient_P_Integral( double* R_AiC_RowMajor ) override {
      // // dimensions of R_AiC_RowMajor: _nNodes x nDim x _nVCIConstraints
      // for ( const auto& sd : _subDomains)
      //   for ( int A = 0; A < _nNodes; A++ )
      //     for ( int i = 0; i < nDim; i++ )
      //       for ( int C = 0; C < _nVCIConstraints; C++ )
      //         R_AiC_RowMajor[A * ( nDim * _nVCIConstraints ) + i * _nVCIConstraints + C] += sd.dT_dY( i, A ) *
      //                                                                                       sd.P( C ) *
      //                                                                                       getSubdomainVolume( sd );
    };

    virtual void vci_compute_Test_PGradient_Integral( double* R_AiC_RowMajor ) override {
      // dimensions of R_AiC_RowMajor: _nNodes x nDim x _nVCIConstraints
      // for ( const auto& sd : _subDomains)
      //   for ( int A = 0; A < _nNodes; A++ )
      //     for ( int i = 0; i < nDim; i++ )
      //       for ( int C = 0; C < _nVCIConstraints; C++ )
      //         R_AiC_RowMajor[A * ( nDim * _nVCIConstraints ) + i * _nVCIConstraints + C] += sd.T( A ) *
      //                                                                                       sd.P_Gradient( C, i ) *
      //                                                                                       getSubdomainVolume( sd );
    };

    virtual void vci_compute_MMatrix( double* mMatrix_ACD_RowMajor ) override {
      // // dimensions of R_AiC_RowMajor: _nNodes x nDim x _nVCIConstraints

      // for ( const auto& sd : _subDomains)
      //   for ( int A = 0; A < _nNodes; A++ ) {
      //     const double R_A = _assignedKernelFunctions[A]->isInSupport( sd.center_IntermediateReference.data() ) ? 1.0
      //                                                                                                           :
      //                                                                                                           0.0;
      //     // const double R_A = 1.0;

      //     for ( int C = 0; C < _nVCIConstraints; C++ )
      //       for ( int D = 0; D < _nVCIConstraints; D++ )
      //         mMatrix_ACD_RowMajor[A * ( _nVCIConstraints * _nVCIConstraints ) + C * _nVCIConstraints +
      //                              D] += R_A * sd.P( C ) * sd.P( D ) * getSubdomainVolume( sd );
      //   }
    };

    virtual void vci_assignTestFunctionCorrectionTerms( const double* eta_AiC_RowMajor ) override {

      // for ( auto& sd : _subDomains)
      //   for ( int A = 0; A < _nNodes; A++ ) {
      //     const double R_A = _assignedKernelFunctions[A]->isInSupport( sd.center_IntermediateReference.data() ) ? 1.0
      //                                                                                                           :
      //                                                                                                           0.0;
      //     // const double R_A = 1.0;
      //     for ( int i = 0; i < nDim; i++ ) {
      //       for ( int C = 0; C < _nVCIConstraints; C++ ) {
      //         sd.dT_dY( i, A ) += eta_AiC_RowMajor[A * ( nDim * _nVCIConstraints ) + i * _nVCIConstraints + C] * R_A
      //         * sd.P( C );
      //       }
      //     }
      //   }
    };

  protected:
    /// \brief Evaluate the shape functions for a vertex-shaped domain
    /// \details This function evaluates the shape functions for a vertex-shaped
    ///         domain using the vertex coordinates of the particle.
    ///         \param particleDomain The particle domain for which to evaluate shape functions.
    ///         \return The shape functions (at center) and their gradients computed from smoothing around the domain.
    std::tuple< Eigen::MatrixXd, Eigen::MatrixXd > evaluateShapeFunctionsAndDerivativesForParticleDomain(
      const ParticleDomain< nDim, nVertices >& particleDomain ) const;

    std::tuple< Eigen::MatrixXd, Eigen::MatrixXd > evaluateShapeFunctionsAndDerivativesOnFace(
      const ParticleDomain< nDim, nVertices >& particleDomain,
      int                                      faceID ) const;

    Eigen::MatrixXd _smoothDerivativeShapeFunctionsForParticleDomain(
      const ParticleDomain< nDim, nVertices >& particleDomain ) const;

    void setVCIOrder( int order )
    {
      _vciOrder        = order;
      _nVCIConstraints = ( order + 1 ) * ( order + 2 ) / 2; // number of VCI constraints for polynomial basis of order
    };
  };

  template < int nDim, int nVertices >
  GenericSDIParticle< nDim, nVertices >::GenericSDIParticle(
    int                                                  elementID,
    const double*                                        vertexCoordinates,
    int                                                  nVertexCoordinates,
    double                                               volume,
    const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation,
    const SmoothingDomainUpdateType                      smoothingVolumeUpdateType )
    : _elementID( elementID ),
      _meshfreeApproximation( approximation ),
      _vciOrder( 0 ),
      _particleDomainMain( vertexCoordinates, nVertexCoordinates, smoothingVolumeUpdateType ),
      _centerDisplacement( nullptr ),
      _centralDeformationGradient( nullptr )
  {
    _dx_dY_center.eye();
    _du_center.zeros();

    _subDomains = _particleDomainMain.uniformSubdivided();
    for ( size_t i = 0; i < _subDomains.size(); i++ ) {
      _subDomainShapeFunctions.push_back( SubDomainShapeFunctions() );
    }

    this->setVCIOrder( _vciOrder );
  }

  // that will go to the generic sdi particle
  template < int nDim, int nVertices >
  std::tuple< typename GenericSDIParticle< nDim, nVertices >::TensorD, typename GenericSDIParticle< nDim, nVertices >::TensorD > GenericSDIParticle<
    nDim,
    nVertices >::getIntermediateConfigurationBoundaryVector( int                                      boundaryFaceID,
                                                             const ParticleDomain< nDim, nVertices >& particleDomain )
    const
  {

    TensorD N_dAY;
    TensorD Y;

    Eigen::Matrix< double, nDim, 1 > _Y_eigen;

    // the evaluation point depends: For real SQCNI, we do it on the smoothing domain boundary, for all others we do in
    // in the center of the deformed geometry

    if ( particleDomain.smoothingVolumeUpdateType == SmoothingDomainUpdateType::DeformationGradient )
      _Y_eigen = particleDomain.getSmoothingDomainFaceCenterCoordinates( boundaryFaceID );
    else
      throw std::invalid_argument( "not implemented" );

    // N_dAY (boundary surface vector for distributed load) comes from the deformed geometry
    auto _N_dAY_eigen = particleDomain.getFaceBoundaryVector( boundaryFaceID );

    for ( int i = 0; i < nDim; i++ ) {
      N_dAY[i] = _N_dAY_eigen[i];
      Y[i]     = _Y_eigen[i];
    }

    return { N_dAY, Y };
  }

  template < int nDim, int nVertices >
  void GenericSDIParticle< nDim, nVertices >::getVertexCoordinates( double* coordinates ) const
  {
    Eigen::Map< Eigen::Matrix< double, nDim, nVertices > > coordinatesMap( coordinates );
    coordinatesMap = _particleDomainMain.getGeometryDeformedVertexCoordinates();
  }

  template < int nDim, int nVertices >
  void GenericSDIParticle< nDim, nVertices >::assignMeshfreeKernelFunctions(
    const std::vector< const MarmotMeshfreeKernelFunction* >& kernelFunctions )
  {

    _assignedKernelFunctions = kernelFunctions;
    _nNodes                  = kernelFunctions.size();

    for ( size_t mpNumber = 0; mpNumber < _subDomainShapeFunctions.size(); mpNumber++ ) {

      auto&       mpl = _subDomainShapeFunctions[mpNumber];
      const auto& sd  = _subDomains[mpNumber];

      const auto [N, dN_dY] = evaluateShapeFunctionsAndDerivativesForParticleDomain( sd );

      const auto mpCenter = sd.getCenterCoordinates();

      mpl.N     = N;
      mpl.dN_dY = dN_dY;

      mpl.T     = N;
      mpl.dT_dY = dN_dY;

      mpl.P.resize( _nVCIConstraints );
      mpl.P_Gradient.resize( _nVCIConstraints, nDim );

      Math::computeMonomialBasis( _vciOrder, mpCenter, mpl.P );
      Math::computeMonomialBasisGradient( _vciOrder, mpCenter, mpl.P_Gradient );
    }
  }

  template < int nDim, int nVertices >
  std::tuple< Eigen::MatrixXd, Eigen::MatrixXd > GenericSDIParticle< nDim, nVertices >::
    evaluateShapeFunctionsAndDerivativesForParticleDomain(
      const ParticleDomain< nDim, nVertices >& particleDomain ) const

  {

    Eigen::Matrix< double, nDim, 1 > coords;
    // Use the center of the provided particleDomain
    coords = particleDomain.getCenterCoordinates();

    Eigen::MatrixXd N = Eigen::MatrixXd::Zero( 1, this->_nNodes );

    // Compute N at the particle center (from GenericParticle)
    this->_meshfreeApproximation.computeShapeFunctions( coords.data(), this->_assignedKernelFunctions, N.data() );

    Eigen::MatrixXd dN_dY = _smoothDerivativeShapeFunctionsForParticleDomain( particleDomain );

    return std::make_tuple( N, dN_dY );
  }

  // this will go to generic sdi particle
  template < int nDim, int nVertices >
  std::tuple< Eigen::MatrixXd, Eigen::MatrixXd > GenericSDIParticle< nDim, nVertices >::
    evaluateShapeFunctionsAndDerivativesOnFace( const ParticleDomain< nDim, nVertices >& particleDomain,
                                                int                                      faceID ) const
  {

    Eigen::Matrix< double, nDim, 1 > coordsFaceCenter;
    // Use the center of the provided particleDomain
    coordsFaceCenter = particleDomain.getSmoothingDomainFaceCenterCoordinates( faceID );

    Eigen::MatrixXd N     = Eigen::MatrixXd::Zero( 1, this->_nNodes );
    Eigen::MatrixXd dN_dY = Eigen::MatrixXd::Zero( nDim, this->_nNodes );

    // Compute N at the particle center (from GenericParticle)
    this->_meshfreeApproximation.computeShapeFunctions( coordsFaceCenter.data(),
                                                        this->_assignedKernelFunctions,
                                                        N.data() );

    dN_dY = _smoothDerivativeShapeFunctionsForParticleDomain( particleDomain );

    return std::make_tuple( N, dN_dY );
  }

  // this will go to generic sdi particle
  template < int nDim, int nVertices >
  Eigen::MatrixXd GenericSDIParticle< nDim, nVertices >::_smoothDerivativeShapeFunctionsForParticleDomain(
    const ParticleDomain< nDim, nVertices >& particleDomain ) const
  {
    Eigen::MatrixXd dN_dY = Eigen::MatrixXd::Zero( nDim, this->_nNodes );

    // Compute dN_dY using the SNNI/SCNI approach (boundary integral over smoothing domain)
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

    return dN_dY;
  }

  template < int nDim, int nVertices >
  void GenericSDIParticle< nDim, nVertices >::computePhysicsKernels( const double* dQ,
                                                                     double*       fInt,
                                                                     double*       dFInt_ddQ,
                                                                     double        timeNew,
                                                                     double        dT )
  {
    using namespace Fastor;
    using namespace Marmot::FastorIndices;
    constexpr int nodeBlockSize = nDim;
    // update central deformation and displacement.
    _dx_dY_center.eye();
    _du_center.zeros();
    {
      const auto [N, dN_dY] = evaluateShapeFunctionsAndDerivativesForParticleDomain( _particleDomainMain );

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

    computePhysicsKernelsOnSubdomains( dQ, fInt, dFInt_ddQ, timeNew, dT );
  }

} // namespace Marmot::Meshfree
