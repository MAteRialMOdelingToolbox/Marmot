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

#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMeshfreeApproximation.h"
#include "Marmot/MarmotMonomialBasisFunctions.h"
#include "Marmot/MarmotParticle.h"
#include "Marmot/MarmotParticleDomain.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Fastor/Fastor.h>
#include <stdexcept>
#include <vector> // Explicitly include vector for clarity

namespace Marmot::Meshfree {

  /**
   * @brief A generic particle class for the subdomain integration (SDI) method.
   *
   * This class extends MarmotParticle and provides common functionalities for
   * particles used in meshfree methods, particularly those involving
   * smoothing domains and a meshfree approximation.
   * It handles particle domain management, shape function evaluation,
   * and state variable management for central displacement and deformation gradient.
   *
   * @tparam nDim The number of dimensions (e.g., 2 for 2D, 3 for 3D).
   * @tparam nVertices The number of vertices defining the particle's geometry.
   */
  template < int nDim, int nVertices >
  class GenericSDIParticle : public MarmotParticle {

  protected:
    // Type aliases for improved readability
    using TensorD  = Fastor::Tensor< double, nDim >;       ///< Alias for a Fastor tensor of dimension nDim.
    using TensorDD = Fastor::Tensor< double, nDim, nDim >; ///< Alias for a Fastor tensor of dimension nDim x nDim.
    using VertexCoordinatesSized = Eigen::Matrix< double, nDim, nVertices >; ///< Alias for Eigen matrix storing vertex
                                                                             ///< coordinates.
    using CoordinatesSized = Eigen::Matrix< double, nDim, 1 >; ///< Alias for Eigen vector storing coordinates.
    using JacobianSized = Eigen::Matrix< double, nDim, nDim >; ///< Alias for Eigen matrix storing Jacobian/deformation
                                                               ///< gradient.
    using KernelFunctionVector = std::vector< const MarmotMeshfreeKernelFunction* >; ///< Alias for a vector of kernel
                                                                                     ///< function pointers.
    using ParticleDomainType = ParticleDomain< nDim, nVertices >; ///< Alias for the particle domain type.

    /**
     * @brief Structure to hold shape functions and their gradients for a subdomain.
     */
    struct SubDomainShapeFunctions {
      Eigen::MatrixXd N;          ///< Shape function values.
      Eigen::MatrixXd dN_dY;      ///< Shape function gradients with respect to reference coordinates.

      Eigen::MatrixXd T;          ///< Test function values (initially same as N).
      Eigen::MatrixXd dT_dY;      ///< Test function gradients (initially same as dN_dY).

      Eigen::VectorXd P;          ///< Monomial basis function values.
      Eigen::MatrixXd P_Gradient; ///< Monomial basis function gradients.
    };

    /// @brief Static constant for the number of state variables per particle.
    /// This includes vertex displacements and center displacement.
    constexpr static int nStateVarsParticle = nDim * nVertices + nDim; // vertex displacements + center displacement

    int _elementID;       ///< The ID of the element this particle belongs to.
    int _nNodes;          ///< The number of nodes (kernel functions) influencing this particle.
    int _vciOrder;        ///< The order of the VCI (Variational Consistent Integration) polynomial basis.
    int _nVCIConstraints; ///< The number of VCI constraints, derived from _vciOrder.

    const MarmotMeshfreeApproximation& _meshfreeApproximation;  ///< Reference to the meshfree approximation object.

    ParticleDomainType                     _particleDomainMain; ///< The main smoothing domain of the particle.
    std::vector< SubDomainShapeFunctions > _subDomainShapeFunctions; ///< Shape functions for each subdomain.
    std::vector< ParticleDomainType >      _subDomains;              ///< Subdomains for integration.

    Eigen::Map< CoordinatesSized > _centerDisplacement;              ///< Mapped central displacement vector.
    Eigen::Map< JacobianSized >    _centralDeformationGradient;      ///< Mapped central deformation gradient tensor.
    Eigen::Map< JacobianSized >    _centralDeformationGradientDelta; ///< Mapped central deformation gradient tensor.

    KernelFunctionVector _assignedKernelFunctions; ///< Pointers to the kernel functions assigned to this particle.

    /// @brief Static vector of valid properties for this particle type.
    inline static const std::vector< std::string > _validProperties = {
      "VCI order",
    };

  public:
    using SmoothingDomainUpdateType = ParticleDomainType::SmoothingDomainUpdateType; ///< Alias for smoothing domain
                                                                                     ///< update type.

    /**
     * @brief Sets multiple properties of the particle.
     * @param properties Pointer to an array of property values.
     * @param nProperties The number of properties in the array.
     */
    virtual void setProperties( [[maybe_unused]] const double* properties,
                                [[maybe_unused]] int           nProperties ) override{};

    /**
     * @brief Sets a single property of the particle by name.
     * @param propertyName The name of the property to set.
     * @param property Pointer to the value of the property.
     */
    virtual void setProperty( const std::string& propertyName, const double* property ) override
    {
      if ( propertyName == "VCI order" ) {
        _vciOrder = static_cast< int >( property[0] );
        this->setVCIOrder( _vciOrder );
      }

      setPropertyOnSubdomains( propertyName, property );
    };

    /**
     * @brief Sets a property on all subdomains.
     * @param propertyName The name of the property.
     * @param property Pointer to the property value.
     */
    virtual void setPropertyOnSubdomains( const std::string& propertyName, const double* property ) = 0;

    /**
     * @brief Get the names of the properties supported by this particle.
     * @return A vector of strings containing the property names.
     */
    virtual std::vector< std::string > getPropertyNames() const override
    {
      auto validProperties = _validProperties;
      validProperties.insert( validProperties.end(),
                              getSubdomainPropertyNames().begin(),
                              getSubdomainPropertyNames().end() );
      return validProperties;
    };

    /**
     * @brief Get the names of the properties supported by the subdomains.
     * @return A vector of strings containing the subdomain property names.
     */
    virtual std::vector< std::string > getSubdomainPropertyNames() const = 0;

    /**
     * @brief Retrieves the current coordinates of the particle's vertices.
     * @param coordinates Pointer to an array where the vertex coordinates will be stored.
     */
    virtual void getVertexCoordinates( double* coordinates ) const override;

    /**
     * @brief Retrieves the coordinates of the center of a specific face.
     * @param faceID The ID of the face.
     * @param coordinates Pointer to an array where the face center coordinates will be stored.
     */
    virtual void getFaceCoordinates( int faceID, double* coordinates ) const override final
    {
      Eigen::Map< CoordinatesSized > coordinatesMap( coordinates );
      coordinatesMap = _particleDomainMain.getFaceCenterCoordinates( faceID );
    }

    /**
     * @brief Retrieves the coordinates of the particle's center.
     * @param coordinates Pointer to an array where the center coordinates will be stored.
     */
    virtual void getCenterCoordinates( double* coordinates ) const override final
    {
      Eigen::Map< CoordinatesSized > centerCoordinatesMap( coordinates );
      centerCoordinatesMap = _particleDomainMain.getCenterCoordinates();
    }

    /**
     * @brief Retrieves the coordinates of the vertices for visualization purposes.
     * @param coordinates Pointer to an array where the visualization vertex coordinates will be stored.
     */
    virtual void getVisualizationVertexCoordinates( double* coordinates ) const override
    {
      getVertexCoordinates( coordinates );
    };

    /**
     * @brief Get the total number of vertices defining the particle's geometry.
     * @return The number of vertices.
     */
    virtual int getNumberOfVertices() const override final { return nVertices; };

    /**
     * @brief Calculates the volume of a given subdomain.
     * @param subdomain The particle domain for which to calculate the volume.
     * @return The volume of the subdomain.
     */
    virtual double getSubdomainVolume( const ParticleDomainType& subdomain ) const = 0;

    /**
     * @brief Initializes the particle's state, setting displacements to zero and
     *        deformation gradient to identity.
     */
    void initializeYourself() override
    {
      _centerDisplacement.setZero();
      _centralDeformationGradient.setIdentity();
      _centralDeformationGradientDelta.setIdentity();

      initializeYourselfOnSubdomains();
    };

    /**
     * @brief Initializes the state of all subdomains.
     */
    virtual void initializeYourselfOnSubdomains() = 0;

    /**
     * @brief Get the shape string of the particle domain.
     * @return A string representing the particle's shape.
     */
    virtual std::string getParticleShape() const override final { return _particleDomainMain.getParticleShape(); };

    /**
     * @brief Constructor for GenericSDIParticle.
     * @param elementID The ID of the element this particle belongs to.
     * @param vertexCoordinates Pointer to an array of initial vertex coordinates.
     * @param nVertexCoordinates The number of vertex coordinates (nDim * nVertices).
     * @param volume The initial volume of the particle.
     * @param approximation Reference to the meshfree approximation object.
     * @param smoothingVolumeUpdateType The type of smoothing domain update to use.
     */
    GenericSDIParticle( int                                elementID,
                        const double*                      vertexCoordinates,
                        [[maybe_unused]] int               nVertexCoordinates,
                        [[maybe_unused]] double            volume,
                        const MarmotMeshfreeApproximation& approximation,
                        const SmoothingDomainUpdateType    smoothingVolumeUpdateType );

    /**
     * @brief Assigns a vector of meshfree kernel functions to the particle.
     * @param kernelFunctions A vector of pointers to the kernel functions.
     */
    virtual void assignMeshfreeKernelFunctions( const KernelFunctionVector& kernelFunctions ) override;

    /**
     * @brief Accepts the current incremental state and updates the particle's position and deformation.
     *        This typically involves adding incremental displacements and updating the deformation gradient.
     */
    virtual void acceptStateAndPosition() override
    {
      _centralDeformationGradient = _centralDeformationGradientDelta * _centralDeformationGradient;
      _centralDeformationGradientDelta.setIdentity();

      _particleDomainMain.acceptStateAndPosition( _centralDeformationGradient, _centerDisplacement );
      for ( auto& sd : _subDomains )
        sd.acceptStateAndPosition( _centralDeformationGradient, _centerDisplacement );

      acceptStateAndPositionOnSubdomains();
    };

    /**
     * @brief Accepts the current incremental state and updates the subdomains' positions and deformations.
     */
    virtual void acceptStateAndPositionOnSubdomains() = 0;

    /**
     * @brief Get the total number of required state variables for this particle.
     * @return The number of required state variables.
     */
    virtual int getNumberOfRequiredStateVars() const override
    {
      int nStateVars = 0;

      nStateVars += nDim;        // center displacement
      nStateVars += nDim * nDim; // central deformation gradient
      nStateVars += nDim * nDim; // central deformation gradient delta

      nStateVars += getNumberOfRequiredStateVarsOnSubdomains();

      return nStateVars;
    };

    /**
     * @brief Get the number of required state variables for the subdomains.
     * @return The number of required state variables for subdomains.
     */
    virtual int getNumberOfRequiredStateVarsOnSubdomains() const = 0;

    /**
     * @brief Assigns a block of memory to store the particle's state variables.
     *        This maps the internal Eigen::Map members to the provided memory block.
     * @param stateVars Pointer to the memory block.
     * @param nStateVars The total number of state variables available in the block.
     */
    void assignStateVars( double* stateVars, int nStateVars ) override final
    {
      int offset = 0;

      new ( &_centerDisplacement ) Eigen::Map< CoordinatesSized >( stateVars + offset );
      offset += nDim;

      new ( &_centralDeformationGradient ) Eigen::Map< JacobianSized >( stateVars + offset );
      offset += nDim * nDim;

      new ( &_centralDeformationGradientDelta ) Eigen::Map< JacobianSized >( stateVars + offset );
      offset += nDim * nDim;

      if ( offset > nStateVars ) {
        throw std::runtime_error( "Error: Number of state variables does not match!" );
      }

      assignStateVarsOnSubdomains( stateVars + offset, nStateVars - offset );
    }

    /**
     * @brief Assigns a block of memory to store the subdomains' state variables.
     * @param stateVars Pointer to the memory block for subdomains.
     * @param nStateVars The total number of state variables available for subdomains.
     */
    virtual void assignStateVarsOnSubdomains( double* stateVars, int nStateVars ) = 0;

    /**
     * @brief Provides a view into a specific state variable.
     * @param stateName The name of the state variable (e.g., "vertex displacements").
     * @param qp The index of the quadrature point or subdomain (unused for "vertex displacements").
     * @return A StateView object providing access to the state variable data.
     */
    virtual StateView getStateView( const std::string& stateName, int qp ) const override final
    {
      if ( stateName == "vertex displacements" )
        return StateView( const_cast< double* >( _particleDomainMain.getGeometryDeformedVertexDisplacements().data() ),
                          nDim * nVertices );

      if ( stateName == "smoothing vertex displacements" )
        return StateView( const_cast< double* >( _particleDomainMain.getSmoothingDomainVertexDisplacements().data() ),
                          nDim * nVertices );

      return getStateViewOnSubdomains( stateName, qp );
    }

    /**
     * @brief Provides a view into a specific state variable for a given subdomain.
     * @param stateName The name of the state variable.
     * @param subdomainIndex The index of the subdomain.
     * @return A StateView object providing access to the state variable data.
     */
    virtual StateView getStateViewOnSubdomains( const std::string& stateName, int subdomainIndex ) const = 0;

    /**
     * @brief Computes the physics kernels (internal forces and their derivatives) for the particle.
     * @param dQ Incremental nodal displacements.
     * @param fInt Internal force vector.
     * @param dFInt_ddQ Stiffness matrix (derivative of internal forces with respect to incremental displacements).
     * @param timeNew The current simulation time.
     * @param dT The time step size.
     */
    virtual void computePhysicsKernels( const double*           dQ,
                                        double*                 fInt,
                                        double*                 dFInt_ddQ,
                                        [[maybe_unused]] double timeNew,
                                        [[maybe_unused]] double dT ) override final;

    /**
     * @brief Computes the physics kernels for all subdomains.
     * @param dQ Incremental nodal displacements.
     * @param fInt Internal force vector.
     * @param dFInt_ddQ Stiffness matrix (derivative of internal forces with respect to incremental displacements).
     * @param timeNew The current simulation time.
     * @param dT The time step size.
     */
    virtual void computePhysicsKernelsOnSubdomains( const double* dQ,
                                                    double*       fInt,
                                                    double*       dFInt_ddQ,
                                                    double        timeNew,
                                                    double        dT ) = 0;

    /**
     * @brief Retrieves the intermediate configuration boundary vector and evaluation point for a given face.
     * @param boundaryFaceID The ID of the boundary face.
     * @param particleDomain The particle domain to consider.
     * @return A tuple containing the boundary surface vector (TensorD) and the evaluation point (TensorD).
     */
    std::tuple< TensorD, TensorD > getIntermediateConfigurationBoundaryVector(
      int                       boundaryFaceID,
      const ParticleDomainType& particleDomain ) const;

    /**
     * @brief Get the dimension of the problem (e.g., 2 for 2D, 3 for 3D).
     * @return The dimension.
     */
    virtual int getDimension() const override final { return nDim; };

    /**
     * @brief Computes the interpolation vector (shape functions) at a given coordinate.
     * @param vec Pointer to an array where the interpolation vector will be stored.
     * @param coordinates Pointer to the coordinates at which to evaluate.
     */
    virtual void getInterpolationVector( double* vec, const double* coordinates ) const override final
    {
      _meshfreeApproximation.computeShapeFunctions( coordinates, _assignedKernelFunctions, vec );
    };

    /**
     * @brief Get the number of VCI constraints.
     * @return The number of VCI constraints.
     */
    virtual int vci_getNumberOfConstraints() override final { return _nVCIConstraints; }

    /**
     * @brief Computes the boundary integral part of the VCI test function P.
     * @param R_AiC_RowMajor Pointer to the result matrix (row-major).
     * @param boundarySurfaceVector Pointer to the boundary surface vector.
     * @param boundaryFaceID The ID of the boundary face.
     */
    virtual void vci_compute_Test_P_BoundaryIntegral( [[maybe_unused]] double*       R_AiC_RowMajor,
                                                      [[maybe_unused]] const double* boundarySurfaceVector,
                                                      [[maybe_unused]] int           boundaryFaceID ) override{

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

    /**
     * @brief Computes the integral of the VCI test function gradient P.
     * @param R_AiC_RowMajor Pointer to the result matrix (row-major).
     */
    virtual void vci_compute_TestGradient_P_Integral( [[maybe_unused]] double* R_AiC_RowMajor ) override{
      // // dimensions of R_AiC_RowMajor: _nNodes x nDim x _nVCIConstraints
      // for ( const auto& sd : _subDomains)
      //   for ( int A = 0; A < _nNodes; A++ )
      //     for ( int i = 0; i < nDim; i++ )
      //       for ( int C = 0; C < _nVCIConstraints; C++ )
      //         R_AiC_RowMajor[A * ( nDim * _nVCIConstraints ) + i * _nVCIConstraints + C] += sd.dT_dY( i, A ) *
      //                                                                                       sd.P( C ) *
      //                                                                                       getSubdomainVolume( sd );
    };

    /**
     * @brief Computes the integral of the VCI test function P gradient.
     * @param R_AiC_RowMajor Pointer to the result matrix (row-major).
     */
    virtual void vci_compute_Test_PGradient_Integral( [[maybe_unused]] double* R_AiC_RowMajor ) override{
      // dimensions of R_AiC_RowMajor: _nNodes x nDim x _nVCIConstraints
      // for ( const auto& sd : _subDomains)
      //   for ( int A = 0; A < _nNodes; A++ )
      //     for ( int i = 0; i < nDim; i++ )
      //       for ( int C = 0; C < _nVCIConstraints; C++ )
      //         R_AiC_RowMajor[A * ( nDim * _nVCIConstraints ) + i * _nVCIConstraints + C] += sd.T( A ) *
      //                                                                                       sd.P_Gradient( C, i ) *
      //                                                                                       getSubdomainVolume( sd );
    };

    /**
     * @brief Computes the M-matrix for VCI.
     * @param mMatrix_ACD_RowMajor Pointer to the result matrix (row-major).
     */
    virtual void vci_compute_MMatrix( [[maybe_unused]] double* mMatrix_ACD_RowMajor ) override{
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

    /**
     * @brief Assigns test function correction terms for VCI.
     * @param eta_AiC_RowMajor Pointer to the correction terms matrix (row-major).
     */
    virtual void vci_assignTestFunctionCorrectionTerms( [[maybe_unused]] const double* eta_AiC_RowMajor ) override{

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
    /**
     * @brief Evaluates the shape functions and their derivatives for a given particle domain.
     * @details This function evaluates the shape functions and their gradients at the center
     *          of the provided particle domain, using a smoothing approach.
     * @param particleDomain The particle domain for which to evaluate shape functions.
     * @return A tuple containing the shape functions (N) and their gradients (dN_dY) as Eigen matrices.
     */
    std::tuple< Eigen::MatrixXd, Eigen::MatrixXd > evaluateShapeFunctionsAndDerivativesForParticleDomain(
      const ParticleDomainType& particleDomain ) const;

    /**
     * @brief Evaluates the shape functions on a specific face of a particle domain.
     * @param particleDomain The particle domain.
     * @param faceID The ID of the face.
     * @return An Eigen::MatrixXd containing the shape functions (N) on the specified
     * face.
     */
    Eigen::MatrixXd evaluateShapeFunctionsOnFace( const ParticleDomainType& particleDomain, int faceID ) const;

    /**
     * @brief Evaluates the shape functions and their derivatives on a specific face of a particle domain.
     * @param particleDomain The particle domain.
     * @param faceID The ID of the face.
     * @return A tuple containing the shape functions (N) and their gradients (dN_dY) as Eigen matrices.
     */
    std::tuple< Eigen::MatrixXd, Eigen::MatrixXd > evaluateShapeFunctionsAndDerivativesOnFace(
      const ParticleDomainType& particleDomain,
      int                       faceID ) const;

    /**
     * @brief Computes the smoothed derivatives of shape functions for a particle domain.
     * @details This method uses a boundary integral approach (SNNI/SCNI) over the smoothing domain
     *          to compute the derivatives of the shape functions.
     * @param particleDomain The particle domain for which to compute smoothed derivatives.
     * @return An Eigen::MatrixXd containing the smoothed derivatives of shape functions (dN_dY).
     */
    Eigen::MatrixXd _smoothDerivativeShapeFunctionsForParticleDomain( const ParticleDomainType& particleDomain ) const;

    /**
     * @brief Sets the VCI order and updates the number of VCI constraints.
     * @param order The desired VCI order.
     */
    void setVCIOrder( int order )
    {
      _vciOrder = order;
      // Calculate the number of VCI constraints for a complete polynomial basis of order 'order'
      // in 'nDim' dimensions. This is given by the binomial coefficient C(order + nDim, nDim).
      long long num_constraints = 1;
      for ( int i = 0; i < nDim; ++i ) {
        num_constraints = num_constraints * ( order + 1 + i ) / ( i + 1 );
      }
      _nVCIConstraints = static_cast< int >( num_constraints );
    };
  };

  /**
   * @brief Constructor for GenericSDIParticle.
   * @tparam nDim The number of dimensions.
   * @tparam nVertices The number of vertices.
   * @param elementID The ID of the element this particle belongs to.
   * @param vertexCoordinates Pointer to an array of initial vertex coordinates.
   * @param nVertexCoordinates The number of vertex coordinates (nDim * nVertices).
   * @param volume The initial volume of the particle.
   * @param approximation Reference to the meshfree approximation object.
   * @param smoothingVolumeUpdateType The type of smoothing domain update to use.
   */
  template < int nDim, int nVertices >
  GenericSDIParticle< nDim, nVertices >::GenericSDIParticle(
    int                                                  elementID,
    const double*                                        vertexCoordinates,
    [[maybe_unused]] int                                 nVertexCoordinates, // Marked as unused
    [[maybe_unused]] double                              volume,             // Marked as unused
    const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation,
    const SmoothingDomainUpdateType                      smoothingVolumeUpdateType )
    : _elementID( elementID ),
      _nNodes( 0 ),          // Initialized to 0, will be set in assignMeshfreeKernelFunctions
      _vciOrder( 0 ),        // Initialized here, then set by setVCIOrder
      _nVCIConstraints( 0 ), // Initialized here, then set by setVCIOrder
      _meshfreeApproximation( approximation ),
      _particleDomainMain( vertexCoordinates, nVertexCoordinates, smoothingVolumeUpdateType ),
      _centerDisplacement( nullptr ),             // Initialized to nullptr, will be re-mapped
      _centralDeformationGradient( nullptr ),     // Initialized to nullptr, will be re-mapped
      _centralDeformationGradientDelta( nullptr ) // Initialized to nullptr, will be re-mapped
  {
    _subDomains = _particleDomainMain.uniformSubdivided();
    _subDomainShapeFunctions.reserve( _subDomains.size() ); // Pre-allocate memory
    for ( size_t i = 0; i < _subDomains.size(); i++ ) {
      _subDomainShapeFunctions.push_back( SubDomainShapeFunctions() );
    }

    this->setVCIOrder( _vciOrder ); // Call setVCIOrder to correctly initialize _nVCIConstraints
  }

  /**
   * @brief Retrieves the intermediate configuration boundary vector and evaluation point for a given face.
   * @tparam nDim The number of dimensions.
   * @tparam nVertices The number of vertices.
   * @param boundaryFaceID The ID of the boundary face.
   * @param particleDomain The particle domain to consider.
   * @return A tuple containing the boundary surface vector (TensorD) and the evaluation point (TensorD).
   */
  template < int nDim, int nVertices >
  std::tuple< typename GenericSDIParticle< nDim, nVertices >::TensorD, typename GenericSDIParticle< nDim, nVertices >::TensorD > GenericSDIParticle<
    nDim,
    nVertices >::getIntermediateConfigurationBoundaryVector( int                       boundaryFaceID,
                                                             const ParticleDomainType& particleDomain ) const
  {

    TensorD N_dAY;
    TensorD Y;

    CoordinatesSized _Y_eigen; // Use alias

    _Y_eigen = particleDomain.getFaceCenterCoordinates( boundaryFaceID );

    // N_dAY (boundary surface vector for distributed load) comes from the deformed geometry
    auto _N_dAY_eigen = particleDomain.getFaceBoundaryVector( boundaryFaceID );

    for ( int i = 0; i < nDim; i++ ) {
      N_dAY[i] = _N_dAY_eigen[i];
      Y[i]     = _Y_eigen[i];
    }

    return { N_dAY, Y };
  }

  /**
   * @brief Retrieves the current coordinates of the particle's vertices.
   * @tparam nDim The number of dimensions.
   * @tparam nVertices The number of vertices.
   * @param coordinates Pointer to an array where the vertex coordinates will be stored.
   */
  template < int nDim, int nVertices >
  void GenericSDIParticle< nDim, nVertices >::getVertexCoordinates( double* coordinates ) const
  {
    Eigen::Map< VertexCoordinatesSized > coordinatesMap( coordinates ); // Use alias
    coordinatesMap = _particleDomainMain.getGeometryDeformedVertexCoordinates();
  }

  /**
   * @brief Assigns a vector of meshfree kernel functions to the particle.
   * @tparam nDim The number of dimensions.
   * @tparam nVertices The number of vertices.
   * @param kernelFunctions A vector of pointers to the kernel functions.
   */
  template < int nDim, int nVertices >
  void GenericSDIParticle< nDim, nVertices >::assignMeshfreeKernelFunctions(
    const KernelFunctionVector& kernelFunctions ) // Use alias
  {

    _assignedKernelFunctions = kernelFunctions;
    _nNodes                  = static_cast< int >( kernelFunctions.size() ); // Cast to int

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

  /**
   * @brief Evaluates the shape functions and their derivatives for a given particle domain.
   * @tparam nDim The number of dimensions.
   * @tparam nVertices The number of vertices.
   * @param particleDomain The particle domain for which to evaluate shape functions.
   * @return A tuple containing the shape functions (N) and their gradients (dN_dY) as Eigen matrices.
   */
  template < int nDim, int nVertices >
  std::tuple< Eigen::MatrixXd, Eigen::MatrixXd > GenericSDIParticle< nDim, nVertices >::
    evaluateShapeFunctionsAndDerivativesForParticleDomain( const ParticleDomainType& particleDomain ) const // Use alias

  {

    CoordinatesSized coords; // Use alias
    // Use the center of the provided particleDomain
    coords = particleDomain.getCenterCoordinates();

    Eigen::MatrixXd N = Eigen::MatrixXd::Zero( 1, this->_nNodes );

    // Compute N at the particle center (from GenericParticle)
    this->_meshfreeApproximation.computeShapeFunctions( coords.data(), this->_assignedKernelFunctions, N.data() );

    Eigen::MatrixXd dN_dY = _smoothDerivativeShapeFunctionsForParticleDomain( particleDomain );

    return std::make_tuple( N, dN_dY );
  }

  /**
   * @brief Evaluates the shape functions on a specific face of a particle domain.
   * @tparam nDim The number of dimensions.
   * @tparam nVertices The number of vertices.
   * @param particleDomain The particle domain.
   * @param faceID The ID of the face.
   * @return The shape functions (N)  as Eigen matrices.
   */
  template < int nDim, int nVertices >
  Eigen::MatrixXd GenericSDIParticle< nDim, nVertices >::evaluateShapeFunctionsOnFace(
    const ParticleDomainType& particleDomain, // Use alias
    int                       faceID ) const
  {

    CoordinatesSized coordsFaceCenter; // Use alias
    // Use the center of the provided particleDomain
    coordsFaceCenter = particleDomain.getSmoothingDomainFaceCenterCoordinates( faceID );

    Eigen::MatrixXd N = Eigen::MatrixXd::Zero( 1, this->_nNodes );

    // Compute N at the particle center (from GenericParticle)
    this->_meshfreeApproximation.computeShapeFunctions( coordsFaceCenter.data(),
                                                        this->_assignedKernelFunctions,
                                                        N.data() );
    return N;
  }

  /**
   * @brief Evaluates the shape functions and their derivatives on a specific face of a particle domain.
   * @tparam nDim The number of dimensions.
   * @tparam nVertices The number of vertices.
   * @param particleDomain The particle domain.
   * @param faceID The ID of the face.
   * @return A tuple containing the shape functions (N) and their gradients (dN_dY) as Eigen matrices.
   */
  template < int nDim, int nVertices >
  std::tuple< Eigen::MatrixXd, Eigen::MatrixXd > GenericSDIParticle< nDim, nVertices >::
    evaluateShapeFunctionsAndDerivativesOnFace( const ParticleDomainType& particleDomain, // Use alias
                                                int                       faceID ) const
  {

    CoordinatesSized coordsFaceCenter; // Use alias
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

  /**
   * @brief Computes the smoothed derivatives of shape functions for a particle domain.
   * @tparam nDim The number of dimensions.
   * @tparam nVertices The number of vertices.
   * @param particleDomain The particle domain for which to compute smoothed derivatives.
   * @return An Eigen::MatrixXd containing the smoothed derivatives of shape functions (dN_dY).
   */
  template < int nDim, int nVertices >
  Eigen::MatrixXd GenericSDIParticle< nDim, nVertices >::_smoothDerivativeShapeFunctionsForParticleDomain(
    const ParticleDomainType& particleDomain ) const // Use alias
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

  /**
   * @brief Computes the physics kernels (internal forces and their derivatives) for the particle.
   * @tparam nDim The number of dimensions.
   * @tparam nVertices The number of vertices.
   * @param dQ Incremental nodal displacements.
   * @param fInt Internal force vector.
   * @param dFInt_ddQ Stiffness matrix (derivative of internal forces with respect to incremental displacements).
   * @param timeNew The current simulation time.
   * @param dT The time step size.
   */
  template < int nDim, int nVertices >
  void GenericSDIParticle< nDim, nVertices >::computePhysicsKernels(
    const double*           dQ,
    double*                 fInt,
    double*                 dFInt_ddQ,
    [[maybe_unused]] double timeNew, // Marked as unused
    [[maybe_unused]] double dT )     // Marked as unused
  {
    using namespace Fastor;
    using namespace Marmot::FastorIndices;
    constexpr int nodeBlockSize = nDim;
    // update central deformation and displacement.
    TensorD  _du_center( 0.0 );
    TensorDD _dx_dY_center;
    _dx_dY_center.eye();
    {
      const auto [N, dN_dY] = evaluateShapeFunctionsAndDerivativesForParticleDomain( _particleDomainMain );

      TensorDD du_dY( 0.0 ); // Use alias
      TensorD  du( 0.0 );    // Use alias

      for ( int B = 0; B < _nNodes; B++ ) {

        const int idxB_u = nodeBlockSize * B;

        const TensorD dN_B_dY = TensorD(
          dN_dY.col( B ).data() ); // works because Eigen is ColumnMajor, Fastor is RowMajor.
                                   // This implicitly transposes dN_dY.col(B) into a row vector for Fastor.

        const TensorD dQU = TensorD( dQ + idxB_u ); // Use alias

        du_dY += einsum< i, j >( dQU, dN_B_dY );
        du += N( B ) * dQU;
      }
      _dx_dY_center += du_dY;
      _du_center += du;
    }

    _centerDisplacement += CoordinatesSized( _du_center.data() );

    Eigen::Map< Eigen::Matrix< double, nDim, nDim, Eigen::RowMajor > > dx_dY_map(
      _dx_dY_center.data() ); // Use RowMajor for Fastor compatibility

    _centralDeformationGradientDelta = dx_dY_map;

    computePhysicsKernelsOnSubdomains( dQ, fInt, dFInt_ddQ, timeNew, dT );
  }

} // namespace Marmot::Meshfree
