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

#include "Marmot/MarmotMeshfreeKernelFunction.h"
#include "Marmot/MarmotUtils.h"
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace Marmot::Discretizations::Meshfree {

  class StiffnessVisitor;

  /**
   * @brief The section of a particle: the measure of the integrals over its domain (see Particle::getSection() and
   * ParticleSectionIntegrals).
   */
  enum class Section {
    Cartesian,    ///< the domain itself: plane strain (unit thickness) in 2D, the volume in 3D
    Axisymmetric, ///< the ring the domain sweeps about the axis, @f$ dV = 2\pi r\,dA @f$ in the meridian half plane
                  ///< @f$ (r, z) = (Y_1, Y_2) @f$
  };

  /**
   * @class Marmot::Discretizations::Meshfree::Particle
   * @brief Abstract interface of a particle of the reproducing kernel particle method (RKPM).
   *
   * A particle is the integration entity of the meshfree layer: it represents a piece of the body (its domain,
   * commonly a point, a quadrilateral or a hexahedron) and it evaluates the weak form on that piece. Unlike a finite
   * element, a particle has no fixed nodes. It interacts with the degrees of freedom, which live on the nodes of the
   * meshfree kernel functions (KernelFunction), via the shape functions of a meshfree approximation
   * (Approximation). The host framework (e.g., EdelweissMeshfree) finds the kernel functions whose
   * support covers the particle and hands them over with assignMeshfreeKernelFunctions(); the particle then
   * evaluates the shape functions @f$ N_A @f$ and (possibly smoothed) gradients @f$ \partial N_A / \partial Y @f$
   * of these @f$ A = 1 \dots n_{\text{nodes}} @f$ nodes.
   *
   * Similar to finite elements, a particle computes the internal forces and their derivatives
   * (computePhysicsKernels()), body loads (computeBodyLoad()) and distributed loads (computeDistributedLoad()).
   * Concrete particles work in an updated Lagrangian setting: gradients are taken with respect to the coordinates
   * @f$ Y @f$ of the last accepted (reference intermediate) configuration, and acceptStateAndPosition() moves the
   * particle to the new configuration at the end of a converged increment.
   *
   * A typical life cycle driven by the host framework is
   * -# creation via Marmot::Registry::ParticleFactory::createParticle(),
   * -# setProperties() (e.g., the order of the variationally consistent integration and time integration
   *    parameters),
   * -# assignStateVars() with a block of getNumberOfRequiredStateVars() doubles, initializeYourself(),
   * -# per increment: assignMeshfreeKernelFunctions(), or updateShapeFunctions() if the kernel functions that cover
   *    the particle did not change, optionally the VCI methods (vci_...), then repeatedly computePhysicsKernels() and
   *    the load methods, and after convergence acceptStateAndPosition().
   *
   * That is an implicit increment. The calls of an explicit increment, in their order, and the three contracts that
   * follow from it (a zero-length increment leaves the rates alone, an accept is idempotent, the material's tangent is
   * required) are written down once, on the page "Meshfree Particles" of the documentation, section "The explicit life
   * cycle of a particle". describeStiffness() is not one of them: it leaves the particle as it was found (and is of an
   * accepted state after acceptStateAndPosition()).
   *
   * **State contract:** the state variable block is a *trial* copy. Before every call of computePhysicsKernels(),
   * the host restores it to the values committed by the last acceptStateAndPosition(), and @c dQ is the total
   * increment since that state (not the Newton correction). Implementations may therefore update their state in
   * place during computePhysicsKernels() (e.g. add the increment of the center displacement or of a nonlocal field
   * to its committed value, or integrate velocities and accelerations), since every call starts again from the
   * committed state; a cutback is simply a restore. EdelweissMeshfree implements this by copying the committed
   * block into the trial block before each computePhysicsKernels() and back in acceptStateAndPosition().
   *
   * For the residual vector and the stiffness matrix, unlike elements and cells, no (field-)blocked layout is
   * possible. This is due to the fact that the number of nodes per particle is not fixed and may vary, and
   * accordingly, the dofIndicesPermutationPattern used to designate the structure of a blocked storage may vary.
   * Accordingly, for performance reasons, no dofIndicesPermutationPattern is provided, and the residual vector and
   * the stiffness matrix are computed in a non-blocked, node-wise layout. Example: [node_1_displacement,
   * node_1_temperature, node_2_displacement, node_2_temperature, ...], in which nodes designate the nodes of the
   * kernel functions attached to the particle, in the order in which they were assigned. The size of both is thus
   * @f$ n_{\text{nodes}} \cdot @f$ getNBaseDof().
   *
   * For the stiffness matrix, a column-major layout is used.
   *
   * The particle also offers the local contributions to the variationally consistent integration (VCI) of Chen,
   * Hillman and Rüter (2013), see the vci_... methods. The global assembly and solution of the VCI correction is
   * done by the host framework.
   *
   * Methods added to this interface are appended at the end of the class, so that the vtable slots of the existing
   * ones do not move: a library and a host built against different versions then fail on the one new call instead of
   * silently calling the wrong method. This is why methods that belong together by subject may stand apart.
   */

  class Particle {

  public:
    /// @brief Default constructor.
    Particle() = default;

    /// @brief Virtual destructor.
    virtual ~Particle() = default;

    /**
     * @brief Assign all the properties of the particle.
     * @details The properties need to follow the order of the property names returned by getPropertyNames().
     * @param[in] properties Array of property values.
     * @param[in] nProperties Number of values in @p properties; must equal the number of property names.
     */
    virtual void setProperties( const double* properties, int nProperties ) = 0;

    /**
     * @brief Assign a single property of the particle by its name.
     * @param[in] propertyName Name of the property, one of getPropertyNames().
     * @param[in] property Pointer to the value of the property.
     */
    virtual void setProperty( const std::string& propertyName, const double* property ) = 0;

    /**
     * @brief Get the names of all the valid properties of the particle.
     * @details The order of the names is the order expected by setProperties().
     * @return The property names.
     */
    virtual std::vector< std::string > getPropertyNames() const = 0;

    /**
     * @brief Assign the meshfree kernel functions (the nodes) that interact with this particle.
     * @details Implementations store the kernel functions and evaluate the shape functions and their gradients
     * (trial and test functions) for the current configuration of the particle (updateShapeFunctions()). The order of
     * @p kernelFunctions defines the node order of the residual vector and the stiffness matrix.
     * @param[in] kernelFunctions Kernel functions whose support covers the particle.
     */
    virtual void assignMeshfreeKernelFunctions( const std::vector< const KernelFunction* >& kernelFunctions ) = 0;

    /**
     * @brief Assign the memory block that holds the state variables of the particle.
     * @param[in] stateVars Pointer to a block of (at least) getNumberOfRequiredStateVars() doubles, owned by the
     * host framework.
     * @param[in] nStateVars Size of the block.
     */
    virtual void assignStateVars( double* stateVars, int nStateVars ) = 0;

    /// @brief Initialize the state variables (e.g., the deformation gradient to identity) and the material.
    virtual void initializeYourself() = 0;

    /**
     * @brief Accept the current state after a converged increment and move the particle to the new configuration.
     * @details The configuration reached becomes the reference intermediate configuration @f$ Y @f$ of the next
     * increment (updated position, volume, geometry and, for VCI, the monomial basis). It must be idempotent: a host
     * may accept twice without an increment in between (calls 2 and 3 of an explicit increment, see the page
     * "Meshfree Particles", section "The explicit life cycle of a particle").
     */
    virtual void acceptStateAndPosition() = 0;

    /**
     * @brief Get a view on a named state (e.g., a state variable of the material point).
     * @param[in] stateName Name of the state.
     * @param[in] qp Index of the evaluation point or subdomain; ignored by particles with a single material point.
     * @return A Marmot::Utilities::StateView (pointer and size) into the state storage.
     */
    virtual Marmot::Utilities::StateView getStateView( const std::string& stateName, int qp ) const = 0;

    /// @brief Get the number of state variables (doubles) required by the particle.
    /// @return The size of the block to be passed to assignStateVars().
    virtual int getNumberOfRequiredStateVars() const = 0;

    /// @brief Get the spatial dimension of the particle.
    /// @return 2 or 3.
    virtual int getDimension() const = 0;

    /// @brief Get the number of vertices of the particle domain (1 for a point particle).
    /// @return The number of vertices.
    virtual int getNumberOfVertices() const = 0;

    /// @brief Get the volume of the particle in the undeformed (initial) configuration.
    /// @return The undeformed volume.
    virtual double getVolumeUndeformed() const = 0;

    /// @brief Get the shape of the particle (e.g., point, quad, hexa) in Ensight Gold notation.
    /// @return The shape name.
    virtual std::string getParticleShape() const = 0;

    /// @brief Get the coordinates of the vertices of the particle (e.g., the nodes of a quadrilateral) for
    /// visualization.
    /// @param[out] coordinates Array of getDimension() * getNumberOfVertices() values, vertex by vertex.
    virtual void getVisualizationVertexCoordinates( double* coordinates ) const = 0;

    /// @brief Get the coordinates of the vertices of the particle (e.g., the nodes of a quadrilateral), e.g., for
    /// testing the coverage of the particle by shape functions.
    /// @param[out] coordinates Array of getDimension() * getNumberOfVertices() values, vertex by vertex.
    virtual void getVertexCoordinates( double* coordinates ) const = 0;

    /// @brief Get the coordinates of the center of a face of the particle, e.g., for boundary loads.
    /// @param[in] faceID ID of the face (1-based).
    /// @param[out] coordinates Array of getDimension() values.
    virtual void getFaceCoordinates( int faceID, double* coordinates ) const = 0;

    /// @brief Get the coordinates of the center of the particle.
    /// @param[out] coordinates Array of getDimension() values.
    virtual void getCenterCoordinates( double* coordinates ) const = 0;

    /// @brief Get the number of dofs per attached node (e.g., 3 for displacement, 1 for temperature, = 4 in total).
    /// @details The actual number of dofs results from the number of attached nodes and the number of dofs per node.
    /// @return The number of dofs per node.
    virtual int getNBaseDof() const = 0;

    /// @brief Get the fields on which the particle lives (e.g., displacement, temperature, ...).
    /// @details For a particle, the fields are the same for each attached node.
    /// @return The field names.
    virtual const std::vector< std::string >& getFields() const = 0;

    /**
     * @brief Compute the internal force vector and its derivative with respect to the dof increment.
     * @details Both are accumulated (+=) in the node-wise layout described in the class documentation.
     * @param[in] dQ Increment of the nodal dofs since the last accepted state, node-wise (the state variables are
     * restored to the accepted state before each call, see the state contract in the class documentation).
     * @param[in,out] fInt Internal force vector of size @f$ n_{\text{nodes}} \cdot @f$ getNBaseDof().
     * @param[in,out] dFInt_ddQ Stiffness matrix, column-major, square of the same size.
     * @param[in] timeNew Time at the end of the increment.
     * @param[in] dT Time increment.
     */
    virtual void computePhysicsKernels( const double* dQ,
                                        double*       fInt,
                                        double*       dFInt_ddQ,
                                        double        timeNew,
                                        double        dT ) = 0;

    /**
     * @brief Compute a body load and its derivative with respect to the dof increment.
     * @param[in] type Load type, a value of getSupportedBodyLoadTypes().
     * @param[in] load Load values (e.g., the body force vector).
     * @param[in,out] fExt External force vector (node-wise).
     * @param[in,out] dExt_dQ Its derivative (column-major).
     * @param[in] timeNew Time at the end of the increment.
     * @param[in] dT Time increment.
     */
    virtual void computeBodyLoad( int           type,
                                  const double* load,
                                  double*       fExt,
                                  double*       dExt_dQ,
                                  double        timeNew,
                                  double        dT ) const = 0;

    /**
     * @brief Compute a distributed (surface) load on a face of the particle and its derivative.
     * @param[in] type Load type, a value of getSupportedDistributedLoadTypes().
     * @param[in] boundaryFaceID ID of the loaded face (1-based).
     * @param[in] load Load values (e.g., the pressure).
     * @param[in,out] fExt External force vector (node-wise).
     * @param[in,out] dExt_dQ Its derivative (column-major).
     * @param[in] timeNew Time at the end of the increment.
     * @param[in] dT Time increment.
     */
    virtual void computeDistributedLoad( int           type,
                                         int           boundaryFaceID,
                                         const double* load,
                                         double*       fExt,
                                         double*       dExt_dQ,
                                         double        timeNew,
                                         double        dT ) const = 0;

    /**
     * @brief Explicit dynamics: update the kinematics and the material state for the dof increment.
     * @details Call 1 of an explicit increment; a zero-length increment leaves the rates alone. The call order of an
     * explicit increment and its contracts: the page "Meshfree Particles", section "The explicit life cycle of a
     * particle".
     * @param[in] dQ Increment of the nodal dofs, node-wise.
     * @param[in] timeNew Time at the end of the increment.
     * @param[in] dT Time increment.
     * @throws std::runtime_error in the default implementation: a particle without it would be integrated with its
     * mass and without its force, silently.
     */
    virtual void updatePhysicsExplicit( [[maybe_unused]] const double* dQ,
                                        [[maybe_unused]] double        timeNew,
                                        [[maybe_unused]] double        dT )
    {
      throw std::runtime_error( "Particle::updatePhysicsExplicit: not implemented for this particle." );
    };

    /**
     * @brief Explicit dynamics: compute the internal force vector for the state set by updatePhysicsExplicit().
     * @details Call 5 of an explicit increment, on the kernel functions of the new configuration. The call order of an
     * explicit increment and its contracts: the page "Meshfree Particles", section "The explicit life cycle of a
     * particle".
     * @param[in,out] fInt Internal force vector (node-wise).
     * @throws std::runtime_error in the default implementation: a particle without it would be integrated with its
     * mass and without its force, silently.
     */
    virtual void computePhysicsKernelsExplicit( [[maybe_unused]] double* fInt )
    {
      throw std::runtime_error( "Particle::computePhysicsKernelsExplicit: not implemented for this particle." );
    };

    /**
     * @brief Explicit dynamics: compute a body load (no tangent).
     * @param[in] type Load type, a value of getSupportedBodyLoadTypes().
     * @param[in] load Load values.
     * @param[in,out] fExt External force vector (node-wise).
     * @param[in] timeNew Time at the end of the increment.
     * @param[in] dT Time increment.
     * @throws std::runtime_error in the default implementation: a particle without it would leave the load out
     * silently. EdelweissMeshfree does not call it: it loads the particles of an explicit analysis through
     * computeBodyLoad() and discards the tangent.
     */
    virtual void computeBodyLoadExplicit( [[maybe_unused]] int           type,
                                          [[maybe_unused]] const double* load,
                                          [[maybe_unused]] double*       fExt,
                                          [[maybe_unused]] double        timeNew,
                                          [[maybe_unused]] double        dT ) const
    {
      throw std::runtime_error( "Particle::computeBodyLoadExplicit: not implemented for this particle." );
    };

    /**
     * @brief Explicit dynamics: compute a distributed load (no tangent).
     * @param[in] type Load type, a value of getSupportedDistributedLoadTypes().
     * @param[in] boundaryFaceID ID of the loaded face (1-based).
     * @param[in] load Load values.
     * @param[in,out] fExt External force vector (node-wise).
     * @param[in] timeNew Time at the end of the increment.
     * @param[in] dT Time increment.
     * @throws std::runtime_error in the default implementation: a particle without it would leave the load out
     * silently. EdelweissMeshfree does not call it: it loads the particles of an explicit analysis through
     * computeDistributedLoad() and discards the tangent.
     */
    virtual void computeDistributedLoadExplicit( [[maybe_unused]] int           type,
                                                 [[maybe_unused]] int           boundaryFaceID,
                                                 [[maybe_unused]] const double* load,
                                                 [[maybe_unused]] double*       fExt,
                                                 [[maybe_unused]] double        timeNew,
                                                 [[maybe_unused]] double        dT ) const
    {
      throw std::runtime_error( "Particle::computeDistributedLoadExplicit: not implemented for this particle." );
    };

    /**
     * @brief Compute the lumped (row-sum) inertia contribution of the particle.
     * @details Call 6 of an explicit increment. The call order of an explicit increment and its contracts: the page
     * "Meshfree Particles", section "The explicit life cycle of a particle".
     * @param[in,out] mLumped Lumped mass vector (node-wise).
     * @throws std::runtime_error in the default implementation.
     */
    virtual void computeLumpedInertia( [[maybe_unused]] double* mLumped ) const
    {
      throw std::runtime_error( "Particle::computeLumpedInertia: not implemented for this particle." );
    };

    /**
     * @brief Compute the lumped momentum contribution of the particle.
     * @details Call 5 of an explicit increment, with the internal force. The call order of an explicit increment and
     * its contracts: the page "Meshfree Particles", section "The explicit life cycle of a particle".
     * @param[in,out] pLumped Lumped momentum vector (node-wise).
     * @throws std::runtime_error in the default implementation.
     */
    virtual void computeLumpedMomentum( [[maybe_unused]] double* pLumped ) const
    {
      throw std::runtime_error( "Particle::computeLumpedMomentum: not implemented for this particle." );
    };

    /**
     * @brief The shape functions of the assigned kernel functions, evaluated exactly at arbitrary coordinates,
     * @f$ N_B(\boldsymbol{x}) @f$.
     * @details For a location that is not one of the particle's own points, e.g. where a host interpolates a field of
     * the kernel functions somewhere in the body. At its own points (its center, its vertices, its face centers) a
     * particle knows a field only by the value and the gradient at its center: its own interpolation there is
     * getInterpolationVectorLinearizedAboutCenter(), which differs from this one away from the center.
     * @param[out] vec Array of @f$ n_{\text{nodes}} @f$ shape function values.
     * @param[in] coordinates Evaluation point (getDimension() values).
     */
    virtual void getInterpolationVector( double* vec, const double* coordinates ) const = 0;

    /// @brief Get the supported body load types.
    /// @return A map from the (upper case) load type name to the @c type integer of computeBodyLoad().
    virtual const std::unordered_map< std::string, int >& getSupportedBodyLoadTypes() const = 0;

    /// @brief Get the supported distributed load types.
    /// @return A map from the (upper case) load type name to the @c type integer of computeDistributedLoad().
    virtual const std::unordered_map< std::string, int >& getSupportedDistributedLoadTypes() const = 0;

    /// @brief Get the coordinates of the evaluation points of the particle (where the shape functions are
    /// evaluated), e.g., for testing the coverage of the particle by shape functions.
    /// @param[out] coordinates Array of getDimension() * getNumberOfEvaluationPoints() values, point by point.
    virtual void getEvaluationCoordinates( double* coordinates ) const = 0;

    /// @brief Get the number of evaluation points of the particle.
    /// @return The number of evaluation points.
    virtual int getNumberOfEvaluationPoints() const = 0;

    /**
     * @name Variationally consistent integration (VCI)
     * The VCI method of Chen, Hillman and Rüter (2013) corrects the gradients of the test functions,
     * @f$ \partial T_A / \partial Y_i = \partial N_A / \partial Y_i + \chi_A \sum_C \eta_{AiC}\, P_C @f$, such that
     * the integration by parts is exact under the numerical integration for a polynomial basis @f$ P_C(Y) @f$
     * (the complete monomials up to the VCI order); @f$ \chi_A @f$ is the indicator of the support of kernel
     * function @f$ A @f$,
     * @f[
     *   \int_\Omega \frac{\partial T_A}{\partial Y_i} P_C \, dV + \int_\Omega T_A \frac{\partial P_C}{\partial Y_i} \,
     * dV = \int_{\partial\Omega} T_A P_C n_i \, dA .
     * @f]
     * The particles provide their local contributions to these integrals, all arrays are row-major, sized
     * @f$ n_{\text{nodes}} \times n_{\text{dim}} \times n_C @f$ (index @f$ AiC @f$) or
     * @f$ n_{\text{nodes}} \times n_C \times n_C @f$ (index @f$ ACD @f$), and contributions are accumulated (+=).
     * The host framework assembles them over all particles per node @f$ A @f$ and solves
     * @f$ \sum_D M_{ACD}\, \eta_{AiD} = R_{AiC} @f$ with the residual
     * @f$ R_{AiC} = \int_{\partial\Omega} T_A P_C n_i\,dA - \int_\Omega \partial_i T_A\, P_C\,dV
     * - \int_\Omega T_A\, \partial_i P_C\,dV @f$.
     * @{
     */

    /// @brief Get the number of VCI constraints @f$ n_C @f$ (the size of the polynomial basis).
    /// @return The number of VCI constraints.
    virtual int vci_getNumberOfConstraints() = 0;

    /**
     * @brief Accumulate the boundary integral @f$ \int_{\partial\Omega} T_A P_C n_i \, dA @f$ of a boundary face.
     * @param[in,out] R_AiC_RowMajor Row-major array @f$ n_{\text{nodes}} \times n_{\text{dim}} \times n_C @f$.
     * @param[in] boundarySurfaceVector Outward boundary surface vector @f$ n\,dA @f$ given by the host.
     * @param[in] boundaryFaceID ID of the particle face on the boundary (1-based).
     */
    virtual void vci_compute_Test_P_BoundaryIntegral( double*       R_AiC_RowMajor,
                                                      const double* boundarySurfaceVector,
                                                      int           boundaryFaceID ) = 0;

    /**
     * @brief Accumulate the volume integral @f$ \int_\Omega \partial T_A / \partial Y_i \, P_C \, dV @f$.
     * @param[in,out] R_AiC_RowMajow Row-major array @f$ n_{\text{nodes}} \times n_{\text{dim}} \times n_C @f$.
     */
    virtual void vci_compute_TestGradient_P_Integral( double* R_AiC_RowMajow ) = 0;

    /**
     * @brief Accumulate the volume integral @f$ \int_\Omega T_A \, \partial P_C / \partial Y_i \, dV @f$.
     * @param[in,out] R_AiC_RowMajor Row-major array @f$ n_{\text{nodes}} \times n_{\text{dim}} \times n_C @f$.
     */
    virtual void vci_compute_Test_PGradient_Integral( double* R_AiC_RowMajor ) = 0;

    /**
     * @brief Accumulate the moment matrix @f$ M_{ACD} = \int_\Omega \chi_A P_C P_D \, dV @f$ of the correction.
     * @param[in,out] M_ACD_RowMajor Row-major array @f$ n_{\text{nodes}} \times n_C \times n_C @f$.
     */
    virtual void vci_compute_MMatrix( double* M_ACD_RowMajor ) = 0;

    /**
     * @brief Assign the correction terms to the gradients of the test functions.
     * @param[in] eta_AiC_RowMajor Row-major array @f$ n_{\text{nodes}} \times n_{\text{dim}} \times n_C @f$ of the
     * coefficients @f$ \eta_{AiC} @f$ for the nodes of this particle.
     */
    virtual void vci_assignTestFunctionCorrectionTerms( const double* eta_AiC_RowMajor ) = 0;

    /** @} */

    /**
     * @brief Set an initial condition (e.g., "geostaticstress").
     * @param[in] conditionName Name of the initial condition.
     * @param[in] value Values of the initial condition.
     */
    virtual void setInitialCondition( const std::string& conditionName, const double* value ) = 0;

    /**
     * @brief Explicit dynamics: compute the lumped (row-sum) damping contribution of the particle.
     * @details The coefficient of each field's FIRST time derivative, node-wise and in the layout of
     * computeLumpedInertia(), which holds the coefficient of the SECOND: the nonlocal viscosity of a gradient-enhanced
     * field, which integrates a first-order field on its own and damps a second-order one. The default adds nothing,
     * i.e. no damping, which is the case of a displacement field. Call 7 of an explicit increment; the call order and
     * its contracts: the page "Meshfree Particles", section "The explicit life cycle of a particle".
     * @param[in,out] cLumped Lumped damping vector (node-wise); entries left untouched mean no damping.
     */
    virtual void computeLumpedDamping( [[maybe_unused]] double* cLumped ) const {};

    /**
     * @brief Describe the stiffness of the particle to a visitor, node by node: the shape functions at its evaluation
     * point, then for every node the derivatives of its force with respect to the increment the particle interpolates
     * (see StiffnessVisitor).
     * @details The stiffness is the tangent of computePhysicsKernels() at a zero increment, without the inertia: the
     * material point is evaluated with the configuration of the last acceptStateAndPosition() and the working copy of
     * the state. It is the stiffness of an accepted state after acceptStateAndPosition(), where an explicit host asks
     * for it (in the middle of an increment, the working copy holds the advanced state, and the stiffness describes
     * neither). Wherever it is called, the particle is left as it was found (state variables, and the response and the
     * tangents the material point keeps). The default throws: a particle that does not describe its stiffness refuses
     * a host that needs it, e.g. for the bound of the critical time step of an explicit scheme
     * (addBoundOfAbsoluteStiffnessTimes() in MarmotParticleStiffnessVisitors.h).
     * @param[in,out] visitor The visitor.
     * @param[in] timeNew The time of the zero increment the material point is evaluated with.
     * @param[in] dT The time increment of the zero increment.
     */
    virtual void describeStiffness( [[maybe_unused]] StiffnessVisitor& visitor,
                                    [[maybe_unused]] double            timeNew,
                                    [[maybe_unused]] double            dT )
    {
      throw std::runtime_error( "Particle::describeStiffness: not implemented for this particle." );
    };

    /**
     * @brief The shape functions linearized about the center of the particle,
     * @f$ G_B(\boldsymbol{Y}) = N_B(\boldsymbol{Y}_c) + \partial N_B / \partial \boldsymbol{Y} \cdot (\boldsymbol{Y} -
     * \boldsymbol{Y}_c) @f$: the particle's own interpolation at its own points.
     * @details A particle knows every field it carries only by its value at its center, @f$ \sum_B N_B
     * \boldsymbol{Q}_B @f$, and its (possibly smoothed) gradient there, @f$ \sum_B \boldsymbol{Q}_B \otimes \partial
     * N_B / \partial \boldsymbol{Y} @f$, with the shape functions and gradients that drive its kinematics. At a point
     * @f$ \boldsymbol{Y} @f$ of its own (its center, a vertex, a face center) its interpolation of a field is therefore
     * the linearization about the center,
     * @f[
     *   G_B = N_B(\boldsymbol{Y}_c) + \frac{\partial N_B}{\partial \boldsymbol{Y}} \cdot (\boldsymbol{Y} -
     *   \boldsymbol{Y}_c) ,
     * @f]
     * for every field: for the displacement, the motion of the point in the next increment,
     * @f$ \Delta\boldsymbol{x} = \sum_B G_B\, \Delta\boldsymbol{Q}_B @f$ (the particle moves as a whole, its center
     * with @f$ N_B @f$ and its domain with the homogeneous deformation gradient increment about the center); for a
     * nonlocal field or a temperature, the value the particle carries at the point. It is not getInterpolationVector()
     * at the point, @f$ N_B(\boldsymbol{Y}) @f$, except at the center: that one is exact at an arbitrary location, not
     * the particle's view of its own points. A constraint that measures a field at a point of the particle and
     * distributes its force with the same vector, @f$ \boldsymbol{f}_B = G_B\,\boldsymbol{f} @f$, has a symmetric
     * stiffness and is conservative (the work of the force on the nodes is the work on the point). The vector sums to
     * one (a rigid translation moves the point as the nodes). It has no default, unlike the methods next to it: the
     * only candidate, the shape functions at the point, is not the particle's own interpolation.
     * @param[out] vec Array of @f$ n_{\text{nodes}} @f$ values.
     * @param[in] coordinates The point in the reference intermediate configuration (getDimension() values).
     */
    virtual void getInterpolationVectorLinearizedAboutCenter( double* vec, const double* coordinates ) const = 0;

    /**
     * @brief Evaluate the shape functions anew where the particle is now, with the kernel functions it holds.
     * @details The shape functions of a particle are evaluated in its configuration (at its center, on its smoothing
     * domain or on its subdomains): the trial functions @f$ N_A @f$, their (possibly smoothed) gradients
     * @f$ \partial N_A / \partial Y_i @f$, and the test functions @f$ T_A = N_A @f$, which drops any VCI correction.
     * acceptStateAndPosition() moves the particle, so they have to be evaluated anew in every increment. If the kernel
     * functions that cover the particle changed, the host hands over the new ones with assignMeshfreeKernelFunctions(),
     * which stores them and then calls this method. If they did not change, the host calls this method alone: the
     * particle keeps its kernel functions, and with them its nodes and the layout of its residual and stiffness.
     */
    virtual void updateShapeFunctions() = 0;

    /**
     * @brief The section of the particle: what its integrals measure. A host that integrates other entities with the
     * particles (finite elements, contact points on the faces of the particles) integrates them over the same
     * measure: over the rings about the axis for Section::Axisymmetric, whose particles live in the meridian half
     * plane @f$ (r, z) @f$.
     * @return Section::Cartesian, unless a particle says otherwise.
     */
    virtual Section getSection() const { return Section::Cartesian; }

  protected:
    /**
     * @brief The shape functions linearized about the center (see getInterpolationVectorLinearizedAboutCenter()),
     * @f$ G_B = N_B + \partial N_B / \partial \boldsymbol{Y} \cdot (\boldsymbol{Y} - \boldsymbol{Y}_c) @f$, from the
     * shape functions and gradients that drive the particle's kinematics.
     * @param[out] vec Array of @p nNodes values.
     * @param[in] N Shape functions at the center, @p nNodes values.
     * @param[in] dN_dY Their gradients, column-major @p nDim x @p nNodes.
     * @param[in] centerCoordinates The center @f$ \boldsymbol{Y}_c @f$ (@p nDim values).
     * @param[in] coordinates The point @f$ \boldsymbol{Y} @f$ (@p nDim values).
     * @param[in] nDim The number of dimensions.
     * @param[in] nNodes The number of nodes.
     */
    static void interpolationVectorLinearizedAboutCenter( double*       vec,
                                                          const double* N,
                                                          const double* dN_dY,
                                                          const double* centerCoordinates,
                                                          const double* coordinates,
                                                          int           nDim,
                                                          int           nNodes )
    {
      for ( int B = 0; B < nNodes; B++ ) {
        vec[B] = N[B];
        for ( int i = 0; i < nDim; i++ )
          vec[B] += dN_dY[B * nDim + i] * ( coordinates[i] - centerCoordinates[i] );
      }
    }
  };

} // namespace Marmot::Discretizations::Meshfree
