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

#include "Marmot/DisplacementMaterialPoint.h"
#include "Marmot/DisplacementParticle.h"
#include "Marmot/MarmotGeometryElement.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotMeshfreeApproximation.h"
#include "Marmot/MarmotParticle.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/src/Core/Matrix.h>
#include <Fastor/Fastor.h>

namespace Marmot::Meshfree {

  template < int nDim, int nVertices >
  class DisplacementParticleSQCNIxSDI : public MarmotParticle {

    using TensorD                = Fastor::Tensor< double, nDim >;
    using TensorDD               = Fastor::Tensor< double, nDim, nDim >;
    using VertexCoordinatesSized = Eigen::Matrix< double, nDim, nVertices >;
    using CoordinatesSized       = Eigen::Matrix< double, nDim, 1 >;

    std::vector< std::unique_ptr< MaterialPoints::DisplacementMaterialPoint< nDim > > > _mps;

    struct SubIntegrationDomain {

      Eigen::MatrixXd N;
      Eigen::MatrixXd dN_dY;

      Eigen::MatrixXd T;
      Eigen::MatrixXd dT_dY;

      double           V_IntermediateReference;
      CoordinatesSized center_IntermediateReference;

      Eigen::VectorXd P;
      Eigen::MatrixXd P_Gradient;
    };

    std::vector< SubIntegrationDomain > _subIntegrationDomains;

  public:
    enum SmoothingVolumeUpdateType { None, DeformationGradient, RotationOnly, RotationAndPrincipalStretch };

  private:
    constexpr int static nStateVarsParticle = nDim * nVertices + nDim; // vertex displacements + center displacement

    int    _elementID;
    double _newmark_beta;
    double _newmark_gamma;
    int    _nNodes;
    int    _vciOrder;
    int    _nVCIConstraints;

    // const SmoothingVolumeUpdateType                _smoothingVolumeUpdateType;
    const Eigen::Matrix< double, nDim, nVertices > _vertexCoordinates_Undeformed;
    CoordinatesSized                               _centerCoordinates_Undeformed;

    const MarmotMeshfreeApproximation&   _meshfreeApproximation;
    Eigen::Map< VertexCoordinatesSized > _vertexDisplacements_Intermediate;
    Eigen::Map< CoordinatesSized >       _centerDisplacement;

    MarmotGeometryElement< nDim, nVertices > _geoElement;
    double                                   _V0_from_geoElement;
    TensorDD                                 _dx_dY_center;

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

      Eigen::Map< Eigen::Matrix< double, nDim, 1 > > segmentCenter( coordinates );

      Eigen::Matrix< double, nDim, nVertices > vertexCoordinates;
      getVertexCoordinates( vertexCoordinates.data() );

      segmentCenter = 0.5 * ( vertexCoordinates.col( faceID % nVertices ) +
                              vertexCoordinates.col( ( faceID - 1 ) % nVertices ) );
    }

    virtual void getCenterCoordinates( double* coordinates ) const override
    {
      const VertexCoordinatesSized vertexCoordinates_Intermediate = _vertexCoordinates_Undeformed +
                                                                    _vertexDisplacements_Intermediate;
      const CoordinatesSized centerCoordinates_Intermediate = _geoElement.N( { 0, 0 } ) *
                                                              vertexCoordinates_Intermediate.transpose();

      for ( int i = 0; i < nDim; i++ ) {
        coordinates[i] = centerCoordinates_Intermediate[i];
      }
    }

    /// \brief Get the visualization vertex coordinates
    /// \param coordinates The coordinates of the vertices
    /// \details This function is used for visualization purposes. It returns the coordinates of the vertices in the
    /// deformed configuration.
    virtual void getVisualizationVertexCoordinates( double* coordinates ) const override
    {
      getVertexCoordinates( coordinates );
    };

    /// \brief Get the number of vertices
    /// \return The number of vertices
    /// \details This function returns the number of vertices of the particle.
    virtual int getNumberOfVertices() const override { return nVertices; };

    /// \brief Get the number of base degrees of freedom per attached node (meshfree kernel function)
    /// \return The number of base degrees of freedom
    /// \details This function returns the number of base degrees of freedom per attached node (meshfree kernel
    /// function).
    virtual int getNBaseDof() const override { return nDofPerNodeU; }

    void initializeYourself() override
    {
      for ( auto& mp : _mps ) {
        mp->initializeYourself();
      }

      _dx_dY_center.eye();
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

    /// \brief Get the names of the fields
    /// \return The names of the fields
    /// \details This function returns the names of the fields of the particle. The fields are:
    /// - "displacement"
    virtual const std::vector< std::string >& getFields() const override
    {
      static const std::vector< std::string > nodeFields = { "displacement" };
      return nodeFields;
    };

    /// \brief Get the shape of the particle
    /// \return The shape of the particle
    /// \details This function returns the shape of the particle. The shape is defined by the geometry element.
    virtual std::string getParticleShape() const override { return _geoElement.getElementShape(); }

    DisplacementParticleSQCNIxSDI( int                                elementID,
                                   const double*                      nodeCoordinates,
                                   int                                nNodeCoordiantes,
                                   double                             volume,
                                   const std::string&                 materialName,
                                   const double*                      materialProperties,
                                   int                                sizeMaterialProperties,
                                   const MarmotMeshfreeApproximation& approximation,
                                   const SmoothingVolumeUpdateType    smoothingVolumeUpdateType );

    /// \brief Assign the meshfree kernel functions
    /// \param kernelFunctions The meshfree kernel functions
    /// \details This function assigns the meshfree kernel functions to the particle. The kernel functions are used to
    /// compute the shape functions and their gradients.
    virtual void assignMeshfreeKernelFunctions(
      const std::vector< const MarmotMeshfreeKernelFunction* >& kernelFunctions ) override;

    /// \brief Get the smoothing volume of the particle
    /// \return The smoothing volume of the particle
    /// \details The smoothing volume is computed as the volume of the element in an updated configuration
    ///         that is obtained by applying (parts of) the deformation gradient to the element in the undeformed
    ///         configuration. In general, this is not consistent with the actual deformation of the particle.
    ///
    virtual double getSmoothingVolume( const VertexCoordinatesSized& vertices ) const
    {

      MarmotGeometryElement< nDim, nVertices > _geometryElementForSmoothing;

      _geometryElementForSmoothing.assignNodeCoordinates( vertices.data() );

      const auto   dNd_dXi_center  = _geometryElementForSmoothing.dNdXi( Eigen::Matrix< double, nDim, 1 >::Zero() );
      const double smoothingVolume = _geometryElementForSmoothing.Jacobian( dNd_dXi_center ).determinant() *
                                     std::pow( 2, nDim );

      return smoothingVolume;
    }

    CoordinatesSized getCenterFromVertices( const VertexCoordinatesSized& vertices ) const
    {
      const CoordinatesSized centerValues = _geoElement.N( { 0, 0 } ) * vertices.transpose();
      return centerValues;
    }

    /// \brief Get the volume of the particle in the undeformed configuration
    /// \return The volume of the particle in the undeformed configuration
    /// \details This function returns the volume of the particle in the undeformed configuration. The volume is
    /// computed using the geometry element. The volume is computed as the determinant of the Jacobian of the geometry
    /// element.
    virtual double getVolumeUndeformed() const { return this->_V0_from_geoElement; };

    /// \brief Get the volume of the particle in the deformed configuration
    /// \return The volume of the particle in the deformed configuration
    /// \details This function returns the volume of the particle in the deformed configuration. The volume is
    /// computed using the geometry element. The volume is
    /// computed as the determinant of the Jacobian of the geometry
    /// element.
    virtual double getVolumeDeformed() const
    {
      double volDeformed = 0.0;
      for ( const auto& subDomain : _subIntegrationDomains ) {
        // Compute the volume of the subdomain in the deformed configuration
        volDeformed += subDomain.V_IntermediateReference;
      }
      return volDeformed;
    }

    ///
    virtual void acceptStateAndPosition() override
    {

      for ( auto& mp : _mps ) {
        mp->acceptStateAndPosition();
      }

      _updateVertexDisplacementsFromMaterialPointDeformation();
      _dx_dY_center.eye();
    };

    virtual int getNumberOfRequiredStateVars() const override
    {
      // return DisplacementParticle< nDim >::getNumberOfRequiredStateVars() + 8;
      int nStateVars = 0;
      nStateVars += nStateVarsParticle;

      for ( const auto& mp : _mps ) {
        nStateVars += mp->getNumberOfRequiredStateVars();
      }

      return nStateVars;
    };

    /// \brief Assign the state variables
    /// \param stateVars The state variables
    /// \param nStateVars The number of state variables
    /// \details This function assigns the state variables to the particle.
    void assignStateVars( double* stateVars, int nStateVars ) override
    {

      int offset = 0;
      new ( &_vertexDisplacements_Intermediate ) Eigen::Map< VertexCoordinatesSized >( stateVars + offset );
      offset += nDim * nVertices;
      new ( &_centerDisplacement ) Eigen::Map< CoordinatesSized >( stateVars + offset );
      offset += nDim;

      for ( auto& mp : _mps ) {
        int nStateVarsSubParticle = mp->getNumberOfRequiredStateVars();
        mp->assignStateVars( stateVars + offset, nStateVarsSubParticle );
        offset += nStateVarsSubParticle;
      }

      if ( offset != nStateVars ) {
        throw std::runtime_error( "Error: Number of state variables does not match!" );
      }
    }

    /// \brief Get the state view of the particle
    /// \param stateName The name of the state variable
    /// \param qp The quadrature point
    /// \return The state view of the particle
    virtual StateView getStateView( const std::string& stateName, int qp ) const override
    {
      if ( stateName == "vertex displacements" ) {
        return StateView( (double*)_vertexDisplacements_Intermediate.data(), nDim * nVertices );
      }
      return _mps[qp]->getStateView( stateName );
    }

    /// \brief Compute the physics kernels
    /// \param dQ The field variable increments
    /// \param fInt The internal forces
    /// \param dFInt_ddQ The internal force gradients wrt. to the field variables
    /// \param timeNew The new time
    /// \param dT The time increment
    /// \details This function computes the physics kernels for the particle. The physics kernels are used to compute
    /// the internal forces and their gradients. The physics kernels are computed using the meshfree approximation.
    /// The physics kernels are computed for each material point in the particle.
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
                                                      int           boundaryFaceID )
    {

      using namespace Fastor;

      const auto [N_dAY, Y_N]   = getIntermediateConfBoundaryVector( boundaryFaceID );
      Eigen::MatrixXd TBoundary = Eigen::MatrixXd::Zero( 1, _nNodes );

      _meshfreeApproximation.computeShapeFunctions( Y_N.data(), _assignedKernelFunctions, TBoundary.data() );

      // get P for the exact integration location at the boundary
      auto PBoundary = Eigen::VectorXd( _nVCIConstraints );
      Math::computeMonomialBasis( _vciOrder, CoordinatesSized( Y_N.data() ), PBoundary );

      for ( int A = 0; A < _nNodes; A++ )
        for ( int i = 0; i < nDim; i++ )
          for ( int C = 0; C < _nVCIConstraints; C++ )
            R_AiC_RowMajor[A * ( nDim * _nVCIConstraints ) + i * _nVCIConstraints + C] += TBoundary( A ) *
                                                                                          PBoundary( C ) * N_dAY[i];
    };

    virtual void vci_compute_TestGradient_P_Integral( double* R_AiC_RowMajor ) override
    {
      // dimensions of R_AiC_RowMajor: _nNodes x nDim x _nVCIConstraints
      for ( const auto& sd : _subIntegrationDomains )
        for ( int A = 0; A < _nNodes; A++ )
          for ( int i = 0; i < nDim; i++ )
            for ( int C = 0; C < _nVCIConstraints; C++ )
              R_AiC_RowMajor[A * ( nDim * _nVCIConstraints ) + i * _nVCIConstraints + C] += sd.dT_dY( i, A ) *
                                                                                            sd.P( C ) *
                                                                                            sd.V_IntermediateReference;
    };

    virtual void vci_compute_Test_PGradient_Integral( double* R_AiC_RowMajor ) override
    {
      // dimensions of R_AiC_RowMajor: _nNodes x nDim x _nVCIConstraints
      for ( const auto& sd : _subIntegrationDomains )
        for ( int A = 0; A < _nNodes; A++ )
          for ( int i = 0; i < nDim; i++ )
            for ( int C = 0; C < _nVCIConstraints; C++ )
              R_AiC_RowMajor[A * ( nDim * _nVCIConstraints ) + i * _nVCIConstraints + C] += sd.T( A ) *
                                                                                            sd.P_Gradient( C, i ) *
                                                                                            sd.V_IntermediateReference;
    };

    virtual void vci_compute_MMatrix( double* mMatrix_ACD_RowMajor ) override
    {
      // dimensions of R_AiC_RowMajor: _nNodes x nDim x _nVCIConstraints

      for ( const auto& sd : _subIntegrationDomains )
        for ( int A = 0; A < _nNodes; A++ ) {
          const double R_A = _assignedKernelFunctions[A]->isInSupport( sd.center_IntermediateReference.data() ) ? 1.0
                                                                                                                : 0.0;
          // const double R_A = 1.0;

          for ( int C = 0; C < _nVCIConstraints; C++ )
            for ( int D = 0; D < _nVCIConstraints; D++ )
              mMatrix_ACD_RowMajor[A * ( _nVCIConstraints * _nVCIConstraints ) + C * _nVCIConstraints +
                                   D] += R_A * sd.P( C ) * sd.P( D ) * sd.V_IntermediateReference;
        }
    };

    virtual void vci_assignTestFunctionCorrectionTerms( const double* eta_AiC_RowMajor ) override
    {

      for ( auto& sd : _subIntegrationDomains )
        for ( int A = 0; A < _nNodes; A++ ) {
          const double R_A = _assignedKernelFunctions[A]->isInSupport( sd.center_IntermediateReference.data() ) ? 1.0
                                                                                                                : 0.0;
          // const double R_A = 1.0;
          for ( int i = 0; i < nDim; i++ ) {
            for ( int C = 0; C < _nVCIConstraints; C++ ) {
              sd.dT_dY( i, A ) += eta_AiC_RowMajor[A * ( nDim * _nVCIConstraints ) + i * _nVCIConstraints + C] * R_A *
                                  sd.P( C );
            }
          }
        }
    };

    virtual void getEvaluationCoordinates( double* coordinates ) const
    {

      Eigen::Map< Eigen::Matrix< double, nDim, 8 > > segmentCenters( coordinates );

      Eigen::Matrix< double, nDim, nVertices > vertexCoordinates;
      getVertexCoordinates( vertexCoordinates.data() );

      for ( int i = 0; i < nVertices; i++ ) {
        segmentCenters.col( i * 2 )     = ( 0.25 * vertexCoordinates.col( ( i + 1 ) % nVertices ) +
                                        0.75 * vertexCoordinates.col( i ) );
        segmentCenters.col( i * 2 + 1 ) = ( 0.75 * vertexCoordinates.col( ( i + 1 ) % nVertices ) +
                                            0.25 * vertexCoordinates.col( i ) );
      }
    }

    virtual int getNumberOfEvaluationPoints() const
    {
      return nVertices * 2; // 2 evaluation points per segment
    };

    virtual void setInitialCondition( const std::string& conditionName, const double* value ) override
    {
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid initial condition" );
    };

  private:
    /// \brief Update the vertex displacements from the material point deformation
    /// \details This function updates the vertex displacements of the particle
    ///         by applying the deformation gradient of the material points
    ///         to the vertex coordinates in the undeformed configuration.
    void _updateVertexDisplacementsFromMaterialPointDeformation();

    /// \brief Evaluate the shape functions for a vertex-shaped domain
    /// \details This function evaluates the shape functions for a vertex-shaped
    ///         domain using the vertex coordinates of the particle.
    ///         \param vertexCoordinates The vertex coordinates of the particle
    ///         \return The shape functions (at center) and their gradients computed from smoothing around the domain.
    std::tuple< Eigen::MatrixXd, Eigen::MatrixXd > evaluateShapeFunctionsForVertexShapedDomain(
      const VertexCoordinatesSized& vertexCoordinates ) const;

    /// \brief Compute the 3x3 coordinates from the 2x2 coordinates
    /// \details This function computes the 3x3 coordinates from the 2x2 coordinates
    ///        using the geometry element. The 3x3 coordinates are used to compute the
    ///        shape functions and their gradients.
    Eigen::Matrix< double, nDim, 9 > _compute3x3From2x2(
      const Eigen::Matrix< double, nDim, nVertices >& coordinates2x2 ) const;

    /// \brief Split the 3x3 coordinates into sub-particles
    /// \details This function splits the 3x3 coordinates into sub-particles
    ///         \param coordinates3x3 The 3x3 coordinates of the particle
    ///         \return The sub-particles
    ///         \details This function splits the 3x3 coordinates into sub-particles
    ///         using the geometry element. The sub-particles are used to compute
    ///         the shape functions and their gradients.
    std::vector< Eigen::Matrix< double, nDim, nVertices > > _split3x3ToSubParticles(
      const Eigen::Matrix< double, nDim, 9 >& coordinates3x3 ) const;

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
    const SmoothingVolumeUpdateType                      smoothingVolumeUpdateType )
    : _elementID( elementID ),
      _newmark_beta( 0. ),
      _newmark_gamma( 0. ),
      _vertexCoordinates_Undeformed( vertexCoordinates ),
      _centerCoordinates_Undeformed( getCenterFromVertices( _vertexCoordinates_Undeformed ) ),
      _meshfreeApproximation( approximation ),
      _vertexDisplacements_Intermediate( nullptr ),
      _centerDisplacement( nullptr ),
      _vciOrder( 0 )
  {

    _geoElement.assignNodeCoordinates( _vertexCoordinates_Undeformed.data() );
    const auto dNd_dXi_center = _geoElement.dNdXi( Eigen::Matrix< double, nDim, 1 >::Zero() );
    _V0_from_geoElement       = _geoElement.Jacobian( dNd_dXi_center ).determinant() * std::pow( 2, nDim );

    const auto subvertices       = _compute3x3From2x2( _vertexCoordinates_Undeformed );
    const auto subDomainVertices = _split3x3ToSubParticles( subvertices );

    for ( const auto& subDomainVertices : subDomainVertices ) {

      const auto   mpCenter = getCenterFromVertices( subDomainVertices );
      const double mpVolume = getSmoothingVolume( subDomainVertices );

      auto theprt = std::make_unique< MaterialPoints::DisplacementMaterialPoint2D >( elementID,
                                                                                     mpCenter.data(),
                                                                                     1,
                                                                                     mpVolume );
      _mps.push_back( std::move( theprt ) );
    }

    int                   materialCode = MarmotLibrary::MarmotMaterialFactory::getMaterialCodeFromName( materialName );
    MarmotMaterialSection section( materialCode, materialProperties, sizeMaterialProperties );

    for ( auto& mp : _mps ) {
      mp->assignMaterial( section );
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

      Eigen::Matrix< double, nDim, nVertices > vertexCoordinates;
      getVertexCoordinates( vertexCoordinates.data() );

      const auto [testBoundary, dN_dY] = evaluateShapeFunctionsForVertexShapedDomain( vertexCoordinates );
      // const auto dT_dY                 = _applyVCICorrectionTermsToShapeFunctionGradients( dN_dY );
      const auto dT_dY = dN_dY;

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

  template < int nDim, int nVertices >
  std::tuple< typename DisplacementParticleSQCNIxSDI< nDim, nVertices >::TensorD,
              typename DisplacementParticleSQCNIxSDI< nDim, nVertices >::TensorD >
  DisplacementParticleSQCNIxSDI< nDim, nVertices >::getIntermediateConfBoundaryVector( int boundaryFaceID ) const
  {

    TensorD N_dAY;
    TensorD Y;

    // if ( _smoothingVolumeUpdateType == DeformationGradient ) {
    // For the full SQCNI case, we actually operate on the deformed element.
    // This means, that the deformed smoothing domain is consistent with the physical domain of the particle.
    Eigen::Matrix< double, nDim, nVertices > vertexCoordinates;
    getVertexCoordinates( vertexCoordinates.data() );

    Eigen::Vector2d t; // tangent vector
    Eigen::Vector2d y; // origin of the boundary segment

    if ( boundaryFaceID > nVertices || boundaryFaceID < 1 )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid boundaryFaceID specified" );

    t = ( vertexCoordinates.col( boundaryFaceID % nVertices ) -
          vertexCoordinates.col( ( boundaryFaceID - 1 ) % nVertices ) );
    y = 0.5 * ( vertexCoordinates.col( boundaryFaceID % nVertices ) +
                vertexCoordinates.col( ( boundaryFaceID - 1 ) % nVertices ) );

    N_dAY = { t( 1 ), -t( 0 ) };
    Y     = { y( 0 ), y( 1 ) };

    return { N_dAY, TensorD( Y.data() ) };
  }

  template < int nDim, int nVertices >
  void DisplacementParticleSQCNIxSDI< nDim, nVertices >::getVertexCoordinates( double* coordinates ) const
  {

    Eigen::Map< Eigen::Matrix< double, nDim, nVertices > > coordinatesDeformed( coordinates );
    Eigen::Map< Eigen::Matrix< double, nDim, nVertices > > vertexDisplacements( _vertexDisplacements_Intermediate );

    coordinatesDeformed = _vertexCoordinates_Undeformed + vertexDisplacements;
  }

  /// \brief Update the vertex displacements from the central deformation gradient and the applied displacement
  /// \details This function updates the vertex displacements of the particle
  template < int nDim, int nVertices >
  void DisplacementParticleSQCNIxSDI< nDim, nVertices >::_updateVertexDisplacementsFromMaterialPointDeformation()
  {

    // CoordinatesSized centerCoordinates_Undeformed = _geoElement.N( { 0, 0 } ) *
    //                                                 _vertexCoordinates_Undeformed.transpose();

    const auto       dx_dY = Eigen::Map< Eigen::Matrix< double, nDim, nDim, Eigen::RowMajor > >( _dx_dY_center.data() );
    CoordinatesSized centerCoordinates_Intermediate;
    getCenterCoordinates( centerCoordinates_Intermediate.data() );

    VertexCoordinatesSized relative_VertexCoordinates_Intermediate = _vertexCoordinates_Undeformed +
                                                                     _vertexDisplacements_Intermediate;
    for ( int i = 0; i < nDim; i++ ) {
      relative_VertexCoordinates_Intermediate.row( i ).array() -= centerCoordinates_Intermediate( i );
    }

    VertexCoordinatesSized relative_vertexCoordinates_Deformed = dx_dY * relative_VertexCoordinates_Intermediate;

    VertexCoordinatesSized vertexCoordinates_Deformed;
    for ( int i = 0; i < nDim; i++ ) {
      vertexCoordinates_Deformed.row( i ).array() = relative_vertexCoordinates_Deformed.row( i ).array() +
                                                    _centerCoordinates_Undeformed( i ) + _centerDisplacement( i );
    }

    _vertexDisplacements_Intermediate = vertexCoordinates_Deformed - _vertexCoordinates_Undeformed;
  }

  template < int nDim, int nVertices >

  void DisplacementParticleSQCNIxSDI< nDim, nVertices >::assignMeshfreeKernelFunctions(
    const std::vector< const MarmotMeshfreeKernelFunction* >& kernelFunctions )
  {

    _assignedKernelFunctions = kernelFunctions;
    _nNodes                  = kernelFunctions.size();

    // This should go to prepare yourself!
    const VertexCoordinatesSized vertexCoordinates_Intermediate = _vertexCoordinates_Undeformed +
                                                                  _vertexDisplacements_Intermediate;

    const auto vertexCoordinatesDeformed3x3 = _compute3x3From2x2( vertexCoordinates_Intermediate );
    const auto subSmoothingDomains          = _split3x3ToSubParticles( vertexCoordinatesDeformed3x3 );

    _subIntegrationDomains.clear();

    for ( size_t mpNumber = 0; mpNumber < _mps.size(); mpNumber++ ) {

      auto& mp = _mps[mpNumber];

      const auto& smoothingDomain = subSmoothingDomains[mpNumber];

      const auto [N, dN_dY] = evaluateShapeFunctionsForVertexShapedDomain( smoothingDomain );

      CoordinatesSized center;
      mp->getCoordinatesAtCenter( center.data() );

      auto sd = SubIntegrationDomain{ .N                       = N,
                                      .dN_dY                   = dN_dY,
                                      .T                       = N,
                                      .dT_dY                   = dN_dY,
                                      .V_IntermediateReference = mp->getVolumeUndeformed() *
                                                                 Fastor::determinant( mp->dY_dX() ),
                                      .center_IntermediateReference = center };

      //
      sd.P.resize( _nVCIConstraints );
      sd.P_Gradient.resize( _nVCIConstraints, nDim );

      Math::computeMonomialBasis( _vciOrder, sd.center_IntermediateReference, sd.P );
      Math::computeMonomialBasisGradient( _vciOrder, sd.center_IntermediateReference, sd.P_Gradient );

      _subIntegrationDomains.push_back( sd );
    }
  }

  template < int nDim, int nVertices >

  std::tuple< Eigen::MatrixXd, Eigen::MatrixXd > DisplacementParticleSQCNIxSDI< nDim, nVertices >::
    evaluateShapeFunctionsForVertexShapedDomain( const VertexCoordinatesSized& vertexCoordinates ) const

  {

    int nNodes = _assignedKernelFunctions.size();

    const auto center = getCenterFromVertices( vertexCoordinates );

    Eigen::MatrixXd N( 1, nNodes );
    N.setZero();
    Eigen::MatrixXd dN_dY( nDim, nNodes );
    dN_dY.setZero();

    _meshfreeApproximation.computeShapeFunctions( center.data(), _assignedKernelFunctions, N.data() );

    Eigen::MatrixXd NBoundary( 1, nNodes );

    for ( int i = 0; i < nVertices; i++ ) {

      Eigen::Vector2d t;
      Eigen::Vector2d n;
      Eigen::Vector2d segmentCenter;
      t             = ( vertexCoordinates.col( ( i + 1 ) % nVertices ) - vertexCoordinates.col( i ) );
      segmentCenter = ( vertexCoordinates.col( ( i + 1 ) % nVertices ) + vertexCoordinates.col( i ) ) / 2;
      n << t( 1 ), -t( 0 );

      _meshfreeApproximation.computeShapeFunctions( segmentCenter.data(), _assignedKernelFunctions, NBoundary.data() );

      dN_dY += n * NBoundary;
    }

    dN_dY /= getSmoothingVolume( vertexCoordinates );

    return std::make_tuple( N, dN_dY );
  }

  template < int nDim, int nVertices >
  Eigen::Matrix< double, nDim, 9 > DisplacementParticleSQCNIxSDI< nDim, nVertices >::_compute3x3From2x2(
    const Eigen::Matrix< double, nDim, nVertices >& coordinates2x2 ) const
  {

    MarmotGeometryElement< nDim, 4 > _geometryElementForInterpolating;

    Eigen::Matrix< double, nDim, 9 > coordinates3x3;
    coordinates3x3.setZero();

    for ( int i = 0; i < 4; i++ ) {
      coordinates3x3.col( i ) = coordinates2x2.col( i );
    }

    coordinates3x3.col( 4 ) = _geometryElementForInterpolating.N( { 0, -1 } ) * coordinates2x2.transpose();
    coordinates3x3.col( 5 ) = _geometryElementForInterpolating.N( { 1, 0 } ) * coordinates2x2.transpose();
    coordinates3x3.col( 6 ) = _geometryElementForInterpolating.N( { 0, 1 } ) * coordinates2x2.transpose();
    coordinates3x3.col( 7 ) = _geometryElementForInterpolating.N( { -1, 0 } ) * coordinates2x2.transpose();
    coordinates3x3.col( 8 ) = _geometryElementForInterpolating.N( { 0, 0 } ) * coordinates2x2.transpose();

    return coordinates3x3;
  }

  template < int nDim, int nVertices >
  std::vector< Eigen::Matrix< double, nDim, nVertices > >

  DisplacementParticleSQCNIxSDI< nDim, nVertices >::_split3x3ToSubParticles(
    const Eigen::Matrix< double, nDim, 9 >& coordinates3x3 ) const
  {

    // (3)--(6)--(2)
    //  | 4  | 3  |
    // (7)--(8)--(5)
    //  | 1  | 2  |
    // (0)--(4)--(1)

    Eigen::Matrix< double, nDim, nVertices > resParticle1;
    Eigen::Matrix< double, nDim, nVertices > resParticle2;
    Eigen::Matrix< double, nDim, nVertices > resParticle3;
    Eigen::Matrix< double, nDim, nVertices > resParticle4;

    resParticle1.col( 0 ) = coordinates3x3.col( 0 );
    resParticle1.col( 1 ) = coordinates3x3.col( 4 );
    resParticle1.col( 2 ) = coordinates3x3.col( 8 );
    resParticle1.col( 3 ) = coordinates3x3.col( 7 );

    resParticle2.col( 0 ) = coordinates3x3.col( 4 );
    resParticle2.col( 1 ) = coordinates3x3.col( 1 );
    resParticle2.col( 2 ) = coordinates3x3.col( 5 );
    resParticle2.col( 3 ) = coordinates3x3.col( 8 );

    resParticle3.col( 0 ) = coordinates3x3.col( 8 );
    resParticle3.col( 1 ) = coordinates3x3.col( 5 );
    resParticle3.col( 2 ) = coordinates3x3.col( 2 );

    resParticle3.col( 3 ) = coordinates3x3.col( 6 );

    resParticle4.col( 0 ) = coordinates3x3.col( 7 );
    resParticle4.col( 1 ) = coordinates3x3.col( 8 );
    resParticle4.col( 2 ) = coordinates3x3.col( 6 );
    resParticle4.col( 3 ) = coordinates3x3.col( 3 );

    return { resParticle1, resParticle2, resParticle3, resParticle4 };
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

    const VertexCoordinatesSized vertexCoordinates_Intermediate = _vertexCoordinates_Undeformed +
                                                                  _vertexDisplacements_Intermediate;

    _dx_dY_center.eye();
    {
      const auto [N, dN_dY] = evaluateShapeFunctionsForVertexShapedDomain( vertexCoordinates_Intermediate );

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
      _centerDisplacement += Eigen::Map< CoordinatesSized >( du.data() );
    }

    for ( size_t mpNumber = 0; mpNumber < _mps.size(); mpNumber++ ) {

      auto& mp = _mps[mpNumber];

      const auto& subDomain = _subIntegrationDomains[mpNumber];

      Tensor< double, nDim > du( 0.0 );

      Tensor< double, nDim, nDim > du_dY( 0.0 );

      for ( int B = 0; B < _nNodes; B++ ) {

        const int idxB_u = nodeBlockSize * B;

        const double N_B     = subDomain.N( B );
        const auto   dN_B_dY = Tensor< double, nDim >(
          subDomain.dN_dY.col( B ).data() ); // works because ColumnMajor of Eigen

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

      Eigen::Map< Eigen::VectorXd > P( fInt, _nNodes * nodeBlockSize );
      Eigen::Map< Eigen::MatrixXd > K( dFInt_ddQ, _nNodes * nodeBlockSize, _nNodes * nodeBlockSize );

      // clang-format off
      for ( int A = 0; A < _nNodes; A++ ) {

        const double T_A = subDomain.T( A );
        const auto                   dT_A_dY = TensorMap< const double, nDim >( subDomain.dT_dY.col( A ).data() );
        const Tensor< double, nDim > dT_A_dx = einsum< ji, j >( inv( mp->dx_dY() ), dT_A_dY );

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

          const double                 N_B     = subDomain.N( B );
          const auto dN_B_dY = TensorMap< const double, nDim >( subDomain.dN_dY.col(B).data() );
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
