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
#include <tuple>  // Explicitly include tuple for std::tuple usage
#include <vector> // Explicitly include vector for clarity

namespace Marmot::Meshfree {

  /**
   * @class Marmot::Meshfree::GenericSDIParticle
   * @brief A generic, physics-independent base class for particles with subdomain integration (SDI).
   *
   * The particle domain (a ParticleDomain, e.g., a quadrilateral or a hexahedron) is uniformly subdivided once into
   * @f$ 2^{n_{\text{dim}}} @f$ subdomains (ParticleDomain::uniformSubdivided()). Each subdomain @f$ s @f$ is an
   * integration point of its own: the derived class attaches a material point to it and integrates the weak form as
   * a sum over the subdomains. The shape functions of a subdomain are evaluated at its center, and their gradients
   * are smoothed over its smoothing domain @f$ \Omega_s @f$ by the divergence theorem, with one point (the face
   * center @f$ \boldsymbol{Y}_f @f$) per face @f$ f @f$ and the boundary surface vector
   * @f$ \boldsymbol{n}_f\,dA_f @f$,
   * @f[
   *   \frac{\partial N_A}{\partial Y_i}\bigg|_s \approx \frac{1}{|\Omega_s|} \sum_f N_A(\boldsymbol{Y}_f)\,
   *   n_{f,i}\, dA_f
   * @f]
   * (see _smoothDerivativeShapeFunctionsForParticleDomain()). The test functions start as copies of the trial
   * functions and are changed only by the VCI correction.
   *
   * The particle as a whole moves with its **center**: in computePhysicsKernels(), the displacement increment and
   * the deformation gradient increment of the center are computed from the shape functions of the main domain,
   * @f[
   *   \Delta\boldsymbol{u}_c = \sum_B N_B\, \Delta\boldsymbol{q}_B, \qquad
   *   \Delta\boldsymbol{F}_c = \boldsymbol{I} + \sum_B \Delta\boldsymbol{q}_B \otimes
   *   \frac{\partial N_B}{\partial \boldsymbol{Y}},
   * @f]
   * and in acceptStateAndPosition() the total central deformation gradient
   * @f$ \boldsymbol{F}_c \leftarrow \Delta\boldsymbol{F}_c\, \boldsymbol{F}_c @f$ and the total center displacement
   * @f$ \boldsymbol{u}_c @f$ move the main domain and all subdomains (ParticleDomain::acceptStateAndPosition()).
   *
   * **State variable layout.** The block given to assignStateVars() is
   * | offset | size | content |
   * |---|---|---|
   * | 0 | nDim | center displacement @f$ \boldsymbol{u}_c @f$ |
   * | nDim | nDim * nDim | central deformation gradient @f$ \boldsymbol{F}_c @f$ (column-major) |
   * | nDim + nDim * nDim | nDim * nDim | its increment @f$ \Delta\boldsymbol{F}_c @f$ (column-major) |
   * | paddedStateVarSize(nStateVarsCenter) | ... | the states of the subdomains (derived class) |
   *
   * The center block of nStateVarsCenter doubles is padded to a multiple of 8 doubles (64 bytes), and derived
   * classes pad the block of every subdomain in the same way (paddedStateVarSize()), so that each subdomain block
   * starts at a multiple of 8 doubles from the start of the particle block.
   *
   * Derived classes implement the physics on the subdomains (the ...OnSubdomains() methods, getSubdomainVolume(),
   * the loads and the VCI boundary integral).
   *
   * @tparam nDim The number of dimensions (2 or 3).
   * @tparam nVertices The number of vertices defining the particle's geometry (4 or 8).
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
     * @struct SubDomainShapeFunctions
     * @brief Shape functions, test functions and VCI basis of one subdomain.
     */
    struct SubDomainShapeFunctions {
      Eigen::MatrixXd N;          ///< Shape function values at the subdomain center (1 x nNodes).
      Eigen::MatrixXd dN_dY;      ///< Smoothed shape function gradients with respect to the reference intermediate
                                  ///< coordinates (nDim x nNodes).

      Eigen::MatrixXd T;          ///< Test function values (initially same as N).
      Eigen::MatrixXd dT_dY;      ///< Test function gradients (dN_dY plus the VCI correction).

      Eigen::VectorXd P;          ///< VCI monomial basis at the subdomain center.
      Eigen::MatrixXd P_Gradient; ///< Gradient of the VCI monomial basis at the subdomain center
                                  ///< (nConstraints x nDim).
    };

    /// @brief Size of vertex displacements plus center displacement (not used by GenericSDIParticle itself).
    constexpr static int nStateVarsParticle = nDim * nVertices + nDim; // vertex displacements + center displacement
    /// @brief Number of state variables of the center: center displacement, central deformation gradient and its
    /// increment (unpadded).
    constexpr static int nStateVarsCenter = nDim + 2 * nDim * nDim;

    int _elementID;       ///< The number (ID) of the particle.
    int _nNodes;          ///< The number of nodes (kernel functions) influencing this particle.
    int _vciOrder;        ///< The order of the VCI (Variationally Consistent Integration) polynomial basis.
    int _nVCIConstraints; ///< The number of VCI constraints, derived from _vciOrder.

    const MarmotMeshfreeApproximation& _meshfreeApproximation;  ///< Reference to the meshfree approximation object.

    ParticleDomainType                     _particleDomainMain; ///< Main domain (geometry, smoothing domain).
    std::vector< SubDomainShapeFunctions > _subDomainShapeFunctions; ///< Shape functions for each subdomain.
    std::vector< ParticleDomainType >      _subDomains;              ///< Subdomains (uniform subdivision).

    /// Total center displacement @f$ \boldsymbol{u}_c @f$ (mapped state variables).
    Eigen::Map< CoordinatesSized > _centerDisplacement;
    /// Total central deformation gradient @f$ \boldsymbol{F}_c @f$ of the last accepted state (mapped state
    /// variables).
    Eigen::Map< JacobianSized > _centralDeformationGradient;
    /// Central deformation gradient increment @f$ \Delta\boldsymbol{F}_c @f$ with respect to the last accepted
    /// state (mapped state variables).
    Eigen::Map< JacobianSized > _centralDeformationGradientDelta;

    KernelFunctionVector _assignedKernelFunctions; ///< Pointers to the kernel functions assigned to this particle.

    /// @brief Static vector of valid properties for this particle type.
    inline static const std::vector< std::string > _validProperties = {
      "VCI order",
    };

  public:
    using SmoothingDomainUpdateType = ParticleDomainType::SmoothingDomainUpdateType; ///< Alias for smoothing domain
                                                                                     ///< update type.

    /**
     * @brief Sets multiple properties of the particle, in the order of getPropertyNames().
     * @param[in] properties Pointer to an array of property values.
     * @param[in] nProperties The number of properties in the array.
     * @throws std::invalid_argument if @p nProperties does not match the number of property names.
     */
    virtual void setProperties( const double* properties, int nProperties ) override
    {
      const std::vector< std::string > propertyNames = getPropertyNames();

      if ( nProperties != static_cast< int >( propertyNames.size() ) )
        throw std::invalid_argument(
          "GenericSDIParticle::setProperties: number of properties (" + std::to_string( nProperties ) +
          ") does not match number of supported properties (" + std::to_string( propertyNames.size() ) + ")." );

      for ( int i = 0; i < nProperties; ++i )
        setProperty( propertyNames[i], &properties[i] );
    };
    /**
     * @brief Sets a single property of the particle by name.
     * @details "VCI order" is handled here (setVCIOrder()), all other names are passed to setPropertyOnSubdomains().
     * @param[in] propertyName The name of the property to set.
     * @param[in] property Pointer to the value of the property.
     */
    virtual void setProperty( const std::string& propertyName, const double* property ) override
    {
      if ( propertyName == "VCI order" ) {
        _vciOrder = static_cast< int >( property[0] );
        this->setVCIOrder( _vciOrder );
        return;
      }

      setPropertyOnSubdomains( propertyName, property );
    };

    /**
     * @brief Sets a property on all subdomains (physics-specific properties).
     * @param[in] propertyName The name of the property.
     * @param[in] property Pointer to the property value.
     */
    virtual void setPropertyOnSubdomains( const std::string& propertyName, const double* property ) = 0;

    /**
     * @brief Get the names of the properties supported by this particle.
     * @return "VCI order" followed by getSubdomainPropertyNames().
     */
    virtual std::vector< std::string > getPropertyNames() const override
    {
      auto       validProperties     = _validProperties;
      const auto subdomainProperties = getSubdomainPropertyNames(); // returned by value: iterate over one copy
      validProperties.insert( validProperties.end(), subdomainProperties.begin(), subdomainProperties.end() );
      return validProperties;
    };

    /**
     * @brief Get the names of the properties supported by the subdomains.
     * @return A vector of strings containing the subdomain property names.
     */
    virtual std::vector< std::string > getSubdomainPropertyNames() const = 0;

    /**
     * @brief Retrieves the coordinates of the particle's vertices (deformed geometry of the main domain).
     * @param[out] coordinates Pointer to an array of nDim * nVertices values.
     */
    virtual void getVertexCoordinates( double* coordinates ) const override;

    /**
     * @brief Retrieves the coordinates of the center of a specific face (deformed geometry of the main domain).
     * @param[in] faceID The ID of the face (1-based).
     * @param[out] coordinates Pointer to an array where the face center coordinates will be stored.
     */
    virtual void getFaceCoordinates( int faceID, double* coordinates ) const override final
    {
      Eigen::Map< CoordinatesSized > coordinatesMap( coordinates );
      coordinatesMap = _particleDomainMain.getFaceCenterCoordinates( faceID );
    }

    /**
     * @brief Retrieves the coordinates of the particle's center (centroid of the deformed geometry).
     * @param[out] coordinates Pointer to an array where the center coordinates will be stored.
     */
    virtual void getCenterCoordinates( double* coordinates ) const override final
    {
      Eigen::Map< CoordinatesSized > centerCoordinatesMap( coordinates );
      centerCoordinatesMap = _particleDomainMain.getCenterCoordinates();
    }

    /**
     * @brief Retrieves the coordinates of the vertices for visualization purposes (same as
     *        getVertexCoordinates()).
     * @param[out] coordinates Pointer to an array where the visualization vertex coordinates will be stored.
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
     * @brief Calculates the volume of a given subdomain, the integration weight of the VCI integrals.
     * @param[in] subdomain The particle domain (an element of _subDomains) for which to calculate the volume.
     * @return The volume of the subdomain.
     */
    virtual double getSubdomainVolume( const ParticleDomainType& subdomain ) const = 0;

    /**
     * @brief Initializes the particle's state, setting the center displacement to zero and the central
     *        deformation gradient and its increment to identity, then initializes the subdomains.
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
     * @brief Constructor for GenericSDIParticle; builds the main domain and its uniform subdivision.
     * @param[in] elementID The number (ID) of the particle.
     * @param[in] vertexCoordinates Pointer to an array of initial vertex coordinates.
     * @param[in] nVertexCoordinates The number of vertex coordinates (nDim * nVertices).
     * @param[in] volume Unused; the volume follows from the vertex coordinates.
     * @param[in] approximation Reference to the meshfree approximation object; it must outlive the particle.
     * @param[in] smoothingVolumeUpdateType The type of smoothing domain update to use.
     */
    GenericSDIParticle( int                                elementID,
                        const double*                      vertexCoordinates,
                        [[maybe_unused]] int               nVertexCoordinates,
                        [[maybe_unused]] double            volume,
                        const MarmotMeshfreeApproximation& approximation,
                        const SmoothingDomainUpdateType    smoothingVolumeUpdateType );

    /**
     * @brief Assigns the meshfree kernel functions and evaluates, for every subdomain, N at its center, the
     *        smoothed gradients dN_dY, the test functions (T = N, dT_dY = dN_dY) and the VCI basis at its center.
     * @param[in] kernelFunctions A vector of pointers to the kernel functions.
     */
    virtual void assignMeshfreeKernelFunctions( const KernelFunctionVector& kernelFunctions ) override;

    /**
     * @brief Accepts the current incremental state and updates the particle's position and deformation.
     * @details @f$ \boldsymbol{F}_c \leftarrow \Delta\boldsymbol{F}_c \boldsymbol{F}_c @f$, the increment is reset
     *          to identity, the main domain and all subdomains are moved with @f$ \boldsymbol{F}_c @f$ and the
     *          center displacement, and finally acceptStateAndPositionOnSubdomains() is called.
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
     * @brief Accepts the state of the subdomains (e.g., of their material points).
     */
    virtual void acceptStateAndPositionOnSubdomains() = 0;

    /**
     * @brief Get the total number of required state variables for this particle.
     * @return The padded size of the center block plus getNumberOfRequiredStateVarsOnSubdomains().
     */
    virtual int getNumberOfRequiredStateVars() const override
    {
      return paddedStateVarSize( nStateVarsCenter ) + getNumberOfRequiredStateVarsOnSubdomains();
    };

    /**
     * @brief The size of a block of state variables, rounded up to a multiple of 8 doubles (64 bytes).
     *
     * The blocks of the particle itself and of every subdomain start at such a multiple, so that the state of each
     * subdomain is aligned as the state of a stand-alone material point: Fastor may use aligned SIMD stores on the
     * tensor maps into it (seen with GCC 14.2 -O3, a segmentation fault for a block starting at an odd double).
     * The offsets are relative to the start of the particle's block, so the absolute alignment also requires the
     * host framework to pass an aligned block.
     * @param[in] n The unpadded number of state variables.
     * @return @p n rounded up to a multiple of 8.
     */
    static constexpr int paddedStateVarSize( int n ) { return ( n + 7 ) / 8 * 8; }

    /**
     * @brief Get the number of required state variables for the subdomains.
     * @return The number of required state variables for subdomains, each subdomain block padded with
     *         paddedStateVarSize().
     */
    virtual int getNumberOfRequiredStateVarsOnSubdomains() const = 0;

    /**
     * @brief Assigns a block of memory to store the particle's state variables.
     *        This maps the internal Eigen::Map members to the provided memory block (layout: see the class
     *        documentation) and passes the rest, starting at paddedStateVarSize(nStateVarsCenter), to
     *        assignStateVarsOnSubdomains().
     * @param[in] stateVars Pointer to the memory block.
     * @param[in] nStateVars The total number of state variables available in the block.
     * @throws std::runtime_error if the block is smaller than the padded center block.
     */
    void assignStateVars( double* stateVars, int nStateVars ) override final
    {
      int offset = 0;

      new ( &_centerDisplacement ) Eigen::Map< CoordinatesSized >( stateVars + offset );
      offset += nDim;

      new ( &_centralDeformationGradient ) Eigen::Map< JacobianSized >( stateVars + offset );
      offset += nDim * nDim;

      new ( &_centralDeformationGradientDelta ) Eigen::Map< JacobianSized >( stateVars + offset );
      offset = paddedStateVarSize( nStateVarsCenter );

      if ( offset > nStateVars ) {
        throw std::runtime_error( "Error: Number of state variables does not match!" );
      }

      assignStateVarsOnSubdomains( stateVars + offset, nStateVars - offset );
    }

    /**
     * @brief Assigns a block of memory to store the subdomains' state variables.
     * @param[in] stateVars Pointer to the memory block for subdomains.
     * @param[in] nStateVars The total number of state variables available for subdomains.
     */
    virtual void assignStateVarsOnSubdomains( double* stateVars, int nStateVars ) = 0;

    /**
     * @brief Provides a view into a specific state variable.
     * @details "vertex displacements" and "smoothing vertex displacements" (nDim * nVertices values each) refer to
     *          the main domain; all other names are passed to getStateViewOnSubdomains().
     * @param[in] stateName The name of the state variable (e.g., "vertex displacements").
     * @param[in] qp The index of the subdomain (unused for the vertex displacements).
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
     * @param[in] stateName The name of the state variable.
     * @param[in] subdomainIndex The index of the subdomain.
     * @return A StateView object providing access to the state variable data.
     */
    virtual StateView getStateViewOnSubdomains( const std::string& stateName, int subdomainIndex ) const = 0;

    /**
     * @brief Computes the physics kernels (internal forces and their derivatives) for the particle.
     * @details Updates the center kinematics from the main domain (the center displacement is incremented by
     *          @f$ \sum_B N_B \Delta\boldsymbol{q}_B @f$ and @f$ \Delta\boldsymbol{F}_c @f$ is set, see the class
     *          documentation), then calls computePhysicsKernelsOnSubdomains(). As for all Marmot entities, the
     *          state variables are updated in place, and the host restores them to the last accepted state before
     *          each evaluation (EdelweissMeshfree does so in every computePhysicsKernels call): the center
     *          displacement thus stays the total one, which acceptStateAndPosition() applies to the undeformed
     *          domain.
     * @param[in] dQ Incremental nodal displacements (the first nDim dofs of each node are used).
     * @param[in,out] fInt Internal force vector.
     * @param[in,out] dFInt_ddQ Stiffness matrix (derivative of internal forces with respect to incremental
     * displacements).
     * @param[in] timeNew The current simulation time.
     * @param[in] dT The time step size.
     */
    virtual void computePhysicsKernels( const double*           dQ,
                                        double*                 fInt,
                                        double*                 dFInt_ddQ,
                                        [[maybe_unused]] double timeNew,
                                        [[maybe_unused]] double dT ) override final;

    /**
     * @brief Computes the physics kernels for all subdomains.
     * @param[in] dQ Incremental nodal displacements.
     * @param[in,out] fInt Internal force vector.
     * @param[in,out] dFInt_ddQ Stiffness matrix (derivative of internal forces with respect to incremental
     * displacements).
     * @param[in] timeNew The current simulation time.
     * @param[in] dT The time step size.
     */
    virtual void computePhysicsKernelsOnSubdomains( const double* dQ,
                                                    double*       fInt,
                                                    double*       dFInt_ddQ,
                                                    double        timeNew,
                                                    double        dT ) = 0;

    /**
     * @brief Retrieves the intermediate configuration boundary vector and evaluation point for a given face.
     * @details Both are taken from the deformed geometry of @p particleDomain: the boundary surface vector
     *          @f$ \boldsymbol{N}\,dA_Y @f$ and the face center.
     * @param[in] boundaryFaceID The ID of the boundary face (1-based).
     * @param[in] particleDomain The particle domain to consider.
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
     * @param[out] vec Pointer to an array where the interpolation vector will be stored.
     * @param[in] coordinates Pointer to the coordinates at which to evaluate.
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
     * @brief Not available here: the VCI boundary integral @f$ \int T_A P_C n_i\,dA @f$ is physics-specific.
     * @param[in,out] R_AiC_RowMajor Pointer to the result matrix (row-major).
     * @param[in] boundarySurfaceVector Pointer to the boundary surface vector.
     * @param[in] boundaryFaceID The ID of the boundary face.
     * @throws std::runtime_error always; derived classes must override it.
     */
    virtual void vci_compute_Test_P_BoundaryIntegral( [[maybe_unused]] double*       R_AiC_RowMajor,
                                                      [[maybe_unused]] const double* boundarySurfaceVector,
                                                      [[maybe_unused]] int           boundaryFaceID ) override
    {

      // This method depends on dY_dX, which is physics-specific.
      // It must be implemented in derived classes.
      throw std::runtime_error( "Error: GenericSDIParticle::vci_compute_Test_P_BoundaryIntegral not implemented. "
                                "Must be implemented in a physics-specific derived class." );
    };

    /**
     * @brief Accumulates the VCI integral of the test function gradient times the basis over the subdomains,
     *        @f$ R_{AiC} \mathrel{+}= \sum_s \partial_i T_A|_s\, P_C|_s\, V_s @f$ (@f$ V_s @f$: getSubdomainVolume()).
     * @param[in,out] R_AiC_RowMajor Pointer to the result matrix (row-major, nNodes x nDim x nConstraints).
     */
    virtual void vci_compute_TestGradient_P_Integral( [[maybe_unused]] double* R_AiC_RowMajor ) override
    {
      // // dimensions of R_AiC_RowMajor: _nNodes x nDim x _nVCIConstraints
      for ( size_t s = 0; s < _subDomains.size(); ++s ) { //( const auto& sd : _subDomains)
        const auto&  sd  = _subDomains[s];
        const auto&  sdf = _subDomainShapeFunctions[s];
        const double vol = getSubdomainVolume( sd );

        for ( int A = 0; A < _nNodes; A++ )
          for ( int i = 0; i < nDim; i++ )
            for ( int C = 0; C < _nVCIConstraints; C++ )
              R_AiC_RowMajor[A * ( nDim * _nVCIConstraints ) + i * _nVCIConstraints + C] += sdf.dT_dY( i, A ) *
                                                                                            sdf.P( C ) * vol;
        // getSubdomainVolume( sd );
      }
    };

    /**
     * @brief Accumulates the VCI integral of the test function times the basis gradient over the subdomains,
     *        @f$ R_{AiC} \mathrel{+}= \sum_s T_A|_s\, \partial_i P_C|_s\, V_s @f$.
     * @param[in,out] R_AiC_RowMajor Pointer to the result matrix (row-major, nNodes x nDim x nConstraints).
     */
    virtual void vci_compute_Test_PGradient_Integral( [[maybe_unused]] double* R_AiC_RowMajor ) override
    {
      // dimensions of R_AiC_RowMajor: _nNodes x nDim x _nVCIConstraints
      for ( size_t s = 0; s < _subDomains.size(); ++s ) { // ( const auto& sd : _subDomains)

        const auto&  sd  = _subDomains[s];
        const auto&  sdf = _subDomainShapeFunctions[s];
        const double vol = getSubdomainVolume( sd );

        for ( int A = 0; A < _nNodes; A++ )
          for ( int i = 0; i < nDim; i++ )
            for ( int C = 0; C < _nVCIConstraints; C++ )
              R_AiC_RowMajor[A * ( nDim * _nVCIConstraints ) + i * _nVCIConstraints + C] += sdf.T( A ) *
                                                                                            sdf.P_Gradient( C, i ) *
                                                                                            vol;
        // getSubdomainVolume( sd );
      }
    };

    /**
     * @brief Accumulates the VCI moment matrix over the subdomains,
     *        @f$ M_{ACD} \mathrel{+}= \sum_s \chi_A\, P_C|_s\, P_D|_s\, V_s @f$.
     * @details @f$ \chi_A @f$ is 1 if the center of the main domain (not of the subdomain) is in the support of
     *          kernel function @f$ A @f$, else 0.
     * @param[in,out] mMatrix_ACD_RowMajor Pointer to the result matrix (row-major, nNodes x nConstraints x
     * nConstraints).
     */
    virtual void vci_compute_MMatrix( [[maybe_unused]] double* mMatrix_ACD_RowMajor ) override
    {
      // // dimensions of R_AiC_RowMajor: _nNodes x nDim x _nVCIConstraints
      auto particleCenter = _particleDomainMain.getCenterCoordinates();

      for ( size_t s = 0; s < _subDomains.size(); ++s ) { // ( const auto& sd : _subDomains)
        const auto&  sd     = _subDomains[s];
        const auto&  sdf    = _subDomainShapeFunctions[s];
        const double vol    = getSubdomainVolume( sd );
        auto         center = sd.getCenterCoordinates();

        for ( int A = 0; A < _nNodes; A++ ) {
          const double R_A = _assignedKernelFunctions[A]->isInSupport( particleCenter.data() ) ? 1.0 : 0.0;
          // const double R_A = _assignedKernelFunctions[A]->isInSupport( center.data() ) ? 1.0 : 0.0;
          //  const double R_A = 1.0;

          for ( int C = 0; C < _nVCIConstraints; C++ )
            for ( int D = 0; D < _nVCIConstraints; D++ )
              mMatrix_ACD_RowMajor[A * ( _nVCIConstraints * _nVCIConstraints ) + C * _nVCIConstraints +
                                   D] += R_A * sdf.P( C ) * sdf.P( D ) * vol; // getSubdomainVolume( sd );
        }
      }
    };

    /**
     * @brief Assigns test function correction terms for VCI,
     *        @f$ \partial_i T_A|_s = \partial_i N_A|_s + \chi_A \sum_C \eta_{AiC} P_C|_s @f$ for every subdomain.
     * @details The correction is applied to the uncorrected gradients, so repeated calls do not accumulate;
     *          @f$ \chi_A @f$ as in vci_compute_MMatrix().
     * @param[in] eta_AiC_RowMajor Pointer to the correction terms matrix (row-major, nNodes x nDim x
     * nConstraints).
     */
    virtual void vci_assignTestFunctionCorrectionTerms( [[maybe_unused]] const double* eta_AiC_RowMajor ) override
    {
      auto particleCenter = _particleDomainMain.getCenterCoordinates();

      for ( size_t s = 0; s < _subDomains.size(); ++s ) { // ( auto& sd : _subDomains)
        const auto& sd     = _subDomains[s];
        auto&       sdf    = _subDomainShapeFunctions[s];
        auto        center = sd.getCenterCoordinates();

        for ( int A = 0; A < _nNodes; A++ ) {
          const double R_A = _assignedKernelFunctions[A]->isInSupport( particleCenter.data() ) ? 1.0 : 0.0;
          // const double R_A = _assignedKernelFunctions[A]->isInSupport( center.data() ) ? 1.0 : 0.0;
          //  const double R_A = 1.0;
          for ( int i = 0; i < nDim; i++ ) {
            double correction = 0.0;

            for ( int C = 0; C < _nVCIConstraints; C++ ) {
              correction += eta_AiC_RowMajor[A * ( nDim * _nVCIConstraints ) + i * _nVCIConstraints + C] * R_A *
                            sdf.P( C );
              // sdf.dT_dY( i, A ) += eta_AiC_RowMajor[A * ( nDim * _nVCIConstraints ) + i * _nVCIConstraints + C] * R_A
              // * sdf.P( C );
            }
            // sdf.dT_dY(i, A ) += correction;
            sdf.dT_dY( i, A ) = sdf.dN_dY( i, A ) + correction;
          }
        }
      }
    };

  protected:
    /**
     * @brief Evaluates the shape functions and their derivatives for a given particle domain.
     * @details N is evaluated at the center of the provided particle domain (the centroid of its deformed
     *          geometry), dN_dY is the smoothed gradient of _smoothDerivativeShapeFunctionsForParticleDomain().
     * @param[in] particleDomain The particle domain for which to evaluate shape functions.
     * @return A tuple containing the shape functions (N) and their gradients (dN_dY) as Eigen matrices.
     */
    std::tuple< Eigen::MatrixXd, Eigen::MatrixXd > evaluateShapeFunctionsAndDerivativesForParticleDomain(
      const ParticleDomainType& particleDomain ) const;

    /**
     * @brief Evaluates the shape functions at the face center of the (deformed) geometry of a particle domain.
     * @details This is the point at which getIntermediateConfigurationBoundaryVector() returns the face center and
     *          at which the VCI boundary term is evaluated; it differs from the face center of the smoothing domain
     *          unless the smoothing domain follows the full deformation gradient.
     * @param[in] particleDomain The particle domain.
     * @param[in] faceID The ID of the face (1-based).
     * @return An Eigen::MatrixXd containing the shape functions (N) on the specified
     * face.
     */
    Eigen::MatrixXd evaluateShapeFunctionsOnFace( const ParticleDomainType& particleDomain, int faceID ) const;

    /**
     * @brief Evaluates the shape functions at the face center of the smoothing domain, together with the smoothed
     *        gradients of the whole particle domain.
     * @param[in] particleDomain The particle domain.
     * @param[in] faceID The ID of the face (1-based).
     * @return A tuple containing the shape functions (N) and their gradients (dN_dY) as Eigen matrices.
     */
    std::tuple< Eigen::MatrixXd, Eigen::MatrixXd > evaluateShapeFunctionsAndDerivativesOnFace(
      const ParticleDomainType& particleDomain,
      int                       faceID ) const;

    /**
     * @brief Computes the smoothed derivatives of shape functions for a particle domain.
     * @details This method uses a boundary integral approach (SNNI/SCNI) over the smoothing domain
     *          to compute the derivatives of the shape functions, with one point per face,
     *          @f$ \partial N_A / \partial Y_i = \frac{1}{|\Omega_s|}\sum_f N_A(\boldsymbol{Y}_f)\, n_{f,i}\,dA_f @f$,
     *          where @f$ \boldsymbol{Y}_f @f$ is the face center and @f$ \boldsymbol{n}_f\,dA_f @f$ the boundary
     *          surface vector of the smoothing domain and @f$ |\Omega_s| @f$ its volume.
     * @param[in] particleDomain The particle domain for which to compute smoothed derivatives.
     * @return An Eigen::MatrixXd containing the smoothed derivatives of shape functions (dN_dY).
     */
    Eigen::MatrixXd _smoothDerivativeShapeFunctionsForParticleDomain( const ParticleDomainType& particleDomain ) const;

    /**
     * @brief Sets the VCI order and updates the number of VCI constraints.
     * @details The number of constraints is @f$ \binom{k + n_{\text{dim}}}{n_{\text{dim}}} @f$ for the order
     *          @f$ k @f$. Unlike GenericParticle::setVCIOrder(), the basis is not evaluated here but in
     *          assignMeshfreeKernelFunctions() (per subdomain), which must therefore be called after the order has
     *          been set.
     * @param[in] order The desired VCI order.
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

  /*
   * (Definition; documented at the declaration.)
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

  /*
   * (Definition; documented at the declaration.)
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

  /*
   * (Definition; documented at the declaration.)
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

  /*
   * (Definition; documented at the declaration.)
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

  /*
   * (Definition; documented at the declaration.)
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

  /*
   * (Definition; documented at the declaration.)
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

    // the face center of the geometry, consistent with getIntermediateConfigurationBoundaryVector()
    const CoordinatesSized coordsFaceCenter = particleDomain.getFaceCenterCoordinates( faceID );

    Eigen::MatrixXd N = Eigen::MatrixXd::Zero( 1, this->_nNodes );

    this->_meshfreeApproximation.computeShapeFunctions( coordsFaceCenter.data(),
                                                        this->_assignedKernelFunctions,
                                                        N.data() );
    return N;
  }

  /*
   * (Definition; documented at the declaration.)
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

  /*
   * (Definition; documented at the declaration.)
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

  /*
   * (Definition; documented at the declaration.)
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

    // the state is restored to the accepted one before every call (see the state contract of MarmotParticle), so
    // this adds the increment of the step to the accepted center displacement once per iteration, not cumulatively
    Eigen::Map< CoordinatesSized > du_center_eigen( _du_center.data() ); // Use alias
    _centerDisplacement += du_center_eigen;

    Eigen::Map< Eigen::Matrix< double, nDim, nDim, Eigen::RowMajor > > dx_dY_map(
      _dx_dY_center.data() ); // Use RowMajor for Fastor compatibility

    _centralDeformationGradientDelta = dx_dY_map;

    computePhysicsKernelsOnSubdomains( dQ, fInt, dFInt_ddQ, timeNew, dT );
  }

} // namespace Marmot::Meshfree
