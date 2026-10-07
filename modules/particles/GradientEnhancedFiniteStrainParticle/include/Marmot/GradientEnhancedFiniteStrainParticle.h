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

#include "Marmot/GradientEnhancedFiniteStrainMaterialPoint.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrain.h"
#include "Marmot/MarmotMeshfreeApproximation.h"
#include "Marmot/MarmotMonomialBasisFunctions.h"
#include "Marmot/MarmotParticle.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotUtils.h"
#include "Marmot/NewmarkBetaIntegrator.h"
#include <vector>

namespace Marmot::Meshfree {

  /**
   * @class Marmot::Meshfree::GradientEnhancedFiniteStrainParticle
   * @brief Nodally integrated point particle for gradient-enhanced (implicit-gradient) finite-strain materials
   *        without a micropolar continuum.
   *
   * The particle couples the displacement field @f$ \boldsymbol{u} @f$ to a single scalar nonlocal field
   * @f$ \bar{N} @f$ (nodal fields `displacement` and `nonlocal damage`, node block
   * @f$ n_\mathrm{dim}+1 @f$, dofs ordered node by node as @f$ \{u_1,\dots,u_{n_\mathrm{dim}},\bar{N}\} @f$). It owns
   * a single Marmot::MaterialPoints::GradientEnhancedFiniteStrainMaterialPoint, which drives a
   * MarmotMaterialGradientEnhancedFiniteStrain, and is integrated with that single point (direct nodal integration).
   *
   * **Kinematics.** The dofs are the increments of the current step. The trial shape functions @f$ N_B @f$ and
   * gradients @f$ \partial N_B/\partial\boldsymbol{Y} @f$ are those of the meshfree approximation at the particle
   * center in the intermediate reference configuration @f$ \boldsymbol{Y} @f$ (the last accepted configuration),
   * evaluated in assignMeshfreeKernelFunctions(). They give
   * @f[
   *   \Delta\boldsymbol{u} = N_B\,\Delta\boldsymbol{q}^U_B,\quad
   *   \Delta F_{ij} = \delta_{ij} + \Delta q^U_{Bi}\,\frac{\partial N_B}{\partial Y_j},\quad
   *   \Delta\bar{N} = N_B\,\Delta q^N_B,\quad
   *   \frac{\partial\Delta\bar{N}}{\partial Y_j} = \Delta q^N_B\,\frac{\partial N_B}{\partial Y_j},
   * @f]
   * and the material point is evaluated at @f$ \boldsymbol{F} = \Delta\boldsymbol{F}\,\boldsymbol{F}_n @f$.
   *
   * **Weak forms.** With the test functions @f$ T_A @f$ (identical to @f$ N_A @f$, except that VCI corrects their
   * gradients), the undeformed particle volume @f$ V_0 @f$, @f$ c = R^2 @f$ from the material and
   * @f$ \partial(\cdot)/\partial x_i = \partial(\cdot)/\partial Y_j\,\Delta F^{-1}_{ji} @f$,
   * @f$ \partial(\cdot)/\partial X_i = \partial(\cdot)/\partial Y_j\,F_{n,ji} @f$, the residuals per node A are
   * @f[
   *   r^U_{Aj} = \frac{\partial T_A}{\partial x_i}\,\tau_{ij}\,V_0 + \rho_0\,a_j\,T_A\,V_0 ,\qquad
   *   r^N_A = \Bigl( T_A\,\Delta\bar{N} + c\,\frac{\partial T_A}{\partial X_i}\frac{\partial\Delta\bar{N}}{\partial
   * X_i}
   *           - T_A\,\Delta L \Bigr) V_0 ,
   * @f]
   * i.e. the Helmholtz equation @f$ \bar{N} - c\,\nabla_X^2\bar{N} = L @f$ is solved in the undeformed configuration
   * for the INCREMENT of the nonlocal field, with the change @f$ \Delta L @f$ of the local driving force
   * (GradientEnhancedFiniteStrainMaterialPoint::response) as its source: the nodal values of the reproducing kernels
   * carry only the increment of the step, and the kernels are rebuilt on the moved nodes every increment, so the total
   * form @f$ T_A(\bar{N} - L) + c\,\nabla_X T_A\cdot\nabla_X\bar{N} @f$ would need @f$ \nabla_X\bar{N} @f$ as an
   * accumulated state of the particle (tested: it changed the nonlocal field by less than the step dependence of the
   * kinematics, without bringing it closer to the result of a single increment). The acceleration @f$ \boldsymbol{a}
   * @f$ follows from the Newmark-beta update of @f$ \Delta\boldsymbol{u} @f$ (@f$ \beta = 0 @f$ switches the inertia
   * off).
   *
   * **Tangent.** With @f$ \partial\boldsymbol{\tau}/\partial\Delta\boldsymbol{F} @f$ etc. from the material point,
   * @f[
   *   K^{UU}_{AjBk} = \Bigl( \frac{\partial T_A}{\partial x_i}\frac{\partial\tau_{ij}}{\partial\Delta F_{kL}}
   *     \frac{\partial N_B}{\partial Y_L} - \frac{\partial T_A}{\partial x_k}\,\tau_{ij}\,
   *     \frac{\partial N_B}{\partial x_i} + \rho_0\,\frac{\partial a_j}{\partial\Delta u_k}\,T_A N_B \Bigr) V_0 ,\quad
   *   K^{UN}_{AjB} = \frac{\partial T_A}{\partial x_i}\frac{\partial\tau_{ij}}{\partial\bar{N}}\,N_B\,V_0 ,
   * @f]
   * @f[
   *   K^{NU}_{ABk} = -T_A\,\frac{\partial L}{\partial\Delta F_{kL}}\frac{\partial N_B}{\partial Y_L}\,V_0 ,\qquad
   *   K^{NN}_{AB} = \Bigl( T_A N_B \bigl( 1 - \frac{\partial L}{\partial\bar{N}} \bigr)
   *     + c\,\frac{\partial T_A}{\partial X_i}\frac{\partial N_B}{\partial X_i} \Bigr) V_0 ,
   * @f]
   * as in the MPM cell and the finite element.
   *
   * **VCI.** The particle implements the variationally consistent integration hooks of MarmotParticle: the test
   * function gradients are corrected by @f$ \partial T_A/\partial Y_i \mathrel{+}= \eta_{AiC}\,P_C(\boldsymbol{Y}) @f$
   * with a monomial basis @f$ P_C @f$ of order `VCI order` (Chen, Hillman, Rüter, 2013). The correction applies to
   * both weak forms.
   *
   * The point particle has no faces and supports no distributed loads; use GradientEnhancedFiniteStrainParticleSQCNI
   * for boundary loads. computeBodyLoad() applies a body force per unit undeformed volume, also for the derived
   * particles.
   *
   * @tparam nDim Spatial dimension (2: plane strain, 3: 3D).
   */
  template < int nDim >
  class GradientEnhancedFiniteStrainParticle : public Marmot::Meshfree::MarmotParticle {

  protected:
    Eigen::Matrix< double, nDim, 1 > _centerCoordinatesUndeformed; ///< center @f$ \boldsymbol{X}_c @f$, undeformed
    Eigen::Matrix< double, nDim, 1 > _centerReferenceIntermediate; ///< center in the intermediate reference
                                                                   ///< configuration @f$ \boldsymbol{Y} @f$
    double _volReferenceIntermediate; ///< volume @f$ V_Y = V_0 \det\boldsymbol{F}_n @f$ (updated at each accept)

    MaterialPoints::GradientEnhancedFiniteStrainMaterialPoint< nDim >* __mp;     ///< the owned material point
    MaterialPoints::GradientEnhancedFiniteStrainMaterialPoint< nDim >& _mp;      ///< reference to the material point
    const MarmotMeshfreeApproximation&                 _meshfreeApproximation;   ///< the meshfree approximation
    std::vector< const MarmotMeshfreeKernelFunction* > _assignedKernelFunctions; ///< kernel functions of the nodes

    double _newmark_beta;  ///< Newmark-beta parameter @f$ \beta @f$ (property `newmark-beta beta`)
    double _newmark_gamma; ///< Newmark-beta parameter @f$ \gamma @f$ (property `newmark-beta gamma`)

    /// The number of currently assigned nodes (= meshfree kernel functions)
    int _nNodes;

    int             _vciOrder;        ///< order of the VCI polynomial basis (property `VCI order`)
    int             _nVCIConstraints; ///< number of monomials of degree <= _vciOrder
    Eigen::VectorXd _P;               ///< VCI monomial basis @f$ P_C @f$ at the particle center
    Eigen::MatrixXd _P_Gradient;      ///< its gradient @f$ \partial P_C/\partial Y_i @f$ (C x i)

    /// The vector of trial shape functions
    Eigen::MatrixXd _N;
    /// The matrix of trial shape function gradients
    Eigen::MatrixXd _dN_dY;

    /// The vector of test shape functions
    Eigen::MatrixXd _T;
    /// The matrix of test shape function gradients
    Eigen::MatrixXd _dT_dY;

    /// static vector of valid properties
    inline static const std::vector< std::string > _validProperties = {
      "newmark-beta beta",
      "newmark-beta gamma",
      "VCI order",
    };

  public:
    /// @brief Body load types.
    enum BodyLoadTypes {
      BodyForce, ///< body force per unit undeformed volume (`BODYFORCE`)
    };

    /// @brief Distributed load types.
    enum DistributedLoadTypes {
      Pressure,     ///< follower pressure on a face (`PRESSURE`); implemented by the SQCNI particles
      CWFCorrection ///< consistent weak form boundary correction (`CWFCORRECTION`); implemented by the SQCNI particles
    };

    /**
     * @brief Supported body loads.
     * @return `BODYFORCE`.
     */
    const std::unordered_map< std::string, int >& getSupportedBodyLoadTypes() const override
    {
      static const std::unordered_map< std::string, int > _supportedBodyLoadTypes = { { "BODYFORCE", BodyForce } };
      return _supportedBodyLoadTypes;
    };

    /**
     * @brief Supported distributed loads: none, a point particle has no faces. The particles with a domain
     * (GradientEnhancedFiniteStrainParticleSQCNI and its derivatives) support `PRESSURE` and `CWFCORRECTION`.
     * @return An empty map.
     */
    const std::unordered_map< std::string, int >& getSupportedDistributedLoadTypes() const override
    {
      static const std::unordered_map< std::string, int > _supportedDistributedLoadTypes = {};
      return _supportedDistributedLoadTypes;
    };

    static constexpr int nDofPerNodeU = nDim;                    ///< dofs per node of the displacement field U
    static constexpr int nDofPerNodeN = 1;                       ///< dofs per node of the nonlocal field N

    using Material = MarmotMaterialGradientEnhancedFiniteStrain; ///< the consumed material interface

    using ForceSized = Eigen::Matrix< double, nDim, 1 >;         ///< force vector of size nDim

    /**
     * @brief Set all properties, in the order of getPropertyNames().
     * @param[in] properties  `newmark-beta beta`, `newmark-beta gamma`, `VCI order`.
     * @param[in] nProperties Number of properties; must be 3.
     * @throws std::runtime_error for a wrong number of properties.
     */
    virtual void setProperties( const double* properties, int nProperties ) override
    {
      if ( nProperties != static_cast< int >( _validProperties.size() ) ) {
        std::ostringstream oss;
        oss << "Error in " << __PRETTY_FUNCTION__ << ": ";
        oss << "Expected " << _validProperties.size() << " properties, but got " << nProperties << ". ";
        oss << "Valid properties are: ";
        for ( const auto& prop : _validProperties ) {
          oss << prop << ", ";
        }
        throw std::runtime_error( oss.str() );
      }

      for ( int i = 0; i < nProperties; i++ ) {
        setProperty( _validProperties[i], &properties[i] );
      }
    };

    /**
     * @brief Set a single property.
     *
     * Setting `VCI order` also resizes and evaluates the VCI monomial basis.
     *
     * @param[in] propertyName `newmark-beta beta`, `newmark-beta gamma` or `VCI order`.
     * @param[in] property     Pointer to the value.
     * @throws std::runtime_error for an unknown property.
     */
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
        std::ostringstream oss;
        oss << "Property " << propertyName << " not supported! Valid properties are: ";
        for ( const auto& prop : _validProperties ) {
          oss << prop << ", ";
        }
        throw std::runtime_error( oss.str() );
      }
    };

    /**
     * @brief Names of the properties.
     * @return `newmark-beta beta`, `newmark-beta gamma`, `VCI order`.
     */
    virtual std::vector< std::string > getPropertyNames() const { return _validProperties; };

    /**
     * @brief Number of state variables: those of the material point (including the material).
     * @return Required size of the state variable vector.
     */
    virtual int getNumberOfRequiredStateVars() const override { return _mp.getNumberOfRequiredStateVars(); };

    /**
     * @brief Assign the state variable vector to the material point.
     * @param[in,out] stateVars  State variable vector.
     * @param[in]     nStateVars Its size.
     */
    void assignStateVars( double* stateVars, int nStateVars ) override { _mp.assignStateVars( stateVars, nStateVars ); }

    /**
     * @brief Assign the kernel functions of the nodes and evaluate the shape functions.
     *
     * @f$ N_B @f$ and @f$ \partial N_B/\partial\boldsymbol{Y} @f$ are evaluated by the meshfree approximation at the
     * particle center of the last accepted state; the test functions are set equal to them (before any VCI
     * correction).
     *
     * @param[in] kernelFunctions Kernel functions of the nodes that support the particle.
     */
    void assignMeshfreeKernelFunctions(
      const std::vector< const MarmotMeshfreeKernelFunction* >& kernelFunctions ) override
    {
      _assignedKernelFunctions = kernelFunctions;

      _nNodes = _assignedKernelFunctions.size();

      Eigen::Matrix< double, nDim, 1 > coords;
      _mp.getCoordinatesAtCenter( coords.data() );

      _N     = Eigen::MatrixXd::Zero( 1, _nNodes );
      _dN_dY = Eigen::MatrixXd::Zero( nDim, _nNodes );

      _meshfreeApproximation.computeShapeFunctionsAndGradients( coords.data(),
                                                                _assignedKernelFunctions,
                                                                _N.data(),
                                                                _dN_dY.data() );

      _T     = _N;
      _dT_dY = _dN_dY;
    }

    /**
     * @brief Number of dofs per node.
     * @return @f$ n_\mathrm{dim} + 1 @f$.
     */
    virtual int getNBaseDof() const { return nDofPerNodeU + nDofPerNodeN; }

    /**
     * @brief Fields per node.
     * @return `displacement`, `nonlocal damage`.
     */
    virtual const std::vector< std::string >& getFields() const override
    {
      static const std::vector< std::string > nodeFields = { "displacement", "nonlocal damage" };
      return nodeFields;
    };

    /**
     * @brief Coordinates of the single vertex: the material point position of the last accepted state.
     * @param[out] coordinates Coordinates (nDim values).
     */
    virtual void getVertexCoordinates( double* coordinates ) const override { _mp.getVertexCoordinates( coordinates ); }

    /**
     * @brief Not available: a point particle has no faces.
     * @param[in]  faceID      Face id.
     * @param[out] coordinates Unused.
     * @throws std::runtime_error always.
     */
    virtual void getFaceCoordinates( int faceID, double* coordinates ) const override
    {
      throw std::runtime_error( "Error: GradientEnhancedFiniteStrainParticle::getFaceCoordinates not implemented." );
    }

    /**
     * @brief Center coordinates, identical to getVertexCoordinates().
     * @param[out] coordinates Coordinates (nDim values).
     */
    virtual void getCenterCoordinates( double* coordinates ) const override { getVertexCoordinates( coordinates ); }

    /**
     * @brief Vertex coordinates for visualization, identical to getVertexCoordinates().
     * @param[out] coordinates Coordinates (nDim values).
     */
    virtual void getVisualizationVertexCoordinates( double* coordinates ) const override
    {
      getVertexCoordinates( coordinates );
    };

    /**
     * @brief Number of vertices.
     * @return 1.
     */
    virtual int getNumberOfVertices() const override { return 1; };

    /**
     * @brief Shape of the particle.
     * @return "point".
     */
    virtual std::string getParticleShape() const override { return "point"; }

    /**
     * @brief Spatial dimension.
     * @return nDim.
     */
    virtual int getDimension() const override { return nDim; };

    /**
     * @brief Construct a point particle and its material point, and assign the material.
     * @param[in] elementID           Label of the particle (also of its material point).
     * @param[in] nodeCoordinates     Undeformed center coordinates (nDim values).
     * @param[in] nNodeCoordiantes    Number of coordinates; must equal nDim.
     * @param[in] volume              Undeformed volume @f$ V_0 @f$.
     * @param[in] materialName        Name of a MarmotMaterialGradientEnhancedFiniteStrain material.
     * @param[in] materialProperties  Material properties.
     * @param[in] nMaterialProperties Number of material properties.
     * @param[in] approximation       Meshfree approximation used for the shape functions.
     * @throws std::invalid_argument for a wrong number of coordinates or an unsuitable material.
     */
    GradientEnhancedFiniteStrainParticle( int                                elementID,
                                          const double*                      nodeCoordinates,
                                          int                                nNodeCoordiantes,
                                          double                             volume,
                                          const std::string&                 materialName,
                                          const double*                      materialProperties,
                                          int                                nMaterialProperties,
                                          const MarmotMeshfreeApproximation& approximation );

    /**
     * @brief Initialize the material point (unit deformation gradient, material state, density).
     */
    void initializeYourself() override { _mp.initializeYourself(); };

    /**
     * @brief Accept the increment and move the intermediate reference configuration.
     *
     * Accepts the state of the material point (@f$ \boldsymbol{F}_n \leftarrow \Delta\boldsymbol{F}\,\boldsymbol{F}_n
     * @f$, @f$ \boldsymbol{u} \leftarrow \boldsymbol{u} + \Delta\boldsymbol{u} @f$), resets its increment, updates
     * @f$ V_Y = V_0\det\boldsymbol{F}_n @f$ and the center @f$ \boldsymbol{Y} @f$, and re-evaluates the VCI basis
     * there. The shape functions are re-evaluated only when the host calls assignMeshfreeKernelFunctions() again.
     */
    virtual void acceptStateAndPosition() override
    {
      _mp.acceptStateAndPosition();
      _mp.prepareYourself( 0, 0 );

      updateVolumeToReferenceIntermediate();
      updateParticlePositionToReferenceIntermediate();

      Math::computeMonomialBasis( _vciOrder, _centerReferenceIntermediate, _P );
      Math::computeMonomialBasisGradient( _vciOrder, _centerReferenceIntermediate, _P_Gradient );
    };

    /**
     * @brief Residual and tangent of the particle for the increment @p dQ.
     *
     * Interpolates @f$ \Delta\boldsymbol{u} @f$, @f$ \Delta\boldsymbol{F} @f$ and @f$ \Delta\bar{N} @f$ from @p dQ,
     * resets and increments the material point, evaluates it, updates velocity and acceleration by Newmark-beta, and
     * adds the residuals @f$ r^U, r^N @f$ and the tangent blocks given in the class description. The results are ADDED
     * to @p fInt and @p dFInt_ddQ (dofs node by node, @f$ \{u_1,\dots,u_{n_\mathrm{dim}},\bar{N}\} @f$).
     *
     * @param[in]     dQ        Increment of the nodal dofs since the last accepted state.
     * @param[in,out] fInt      Internal force vector.
     * @param[in,out] dFInt_ddQ Tangent, column-major.
     * @param[in]     timeNew   Time at the end of the increment.
     * @param[in]     dT        Time increment.
     */
    virtual void computePhysicsKernels( const double* dQ,
                                        double*       fInt,
                                        double*       dFInt_ddQ,
                                        double        timeNew,
                                        double        dT ) override;

    /**
     * @brief Body load: a body force @f$ \boldsymbol{b} @f$ per unit undeformed volume (dead load) on the
     * displacement field, @f$ P_{Ai} \mathrel{-}= T_A\,b_i\,V_0 @f$ (the host's sign convention for external loads,
     * as in the cells); the tangent is zero and the nonlocal rows are not loaded. Also used by the derived SQCNI /
     * NSNI particles.
     * @param[in]     type    Load type (BodyForce).
     * @param[in]     load    Body force vector (nDim values).
     * @param[in,out] fExt    Load vector, the contribution is added.
     * @param[in,out] dExt_dQ Load tangent (untouched).
     * @throws std::invalid_argument for another load type.
     * @param[in]     timeNew Time at the end of the increment.
     * @param[in]     dT      Time increment.
     */
    virtual void computeBodyLoad( int           type,
                                  const double* load,
                                  double*       fExt,
                                  double*       dExt_dQ,
                                  double        timeNew,
                                  double        dT ) const override;

    /**
     * @brief Distributed load: a point particle has no faces, so there is none (see
     * getSupportedDistributedLoadTypes() and GradientEnhancedFiniteStrainParticleSQCNI::computeDistributedLoad).
     * @throws std::invalid_argument always.
     * @param[in]     type      Load type.
     * @param[in]     surfaceID Face id.
     * @param[in]     load      Load values.
     * @param[in,out] fExt      Load vector (untouched).
     * @param[in,out] dExt_dQ   Load tangent (untouched).
     * @param[in]     timeNew   Time at the end of the increment.
     * @param[in]     dT        Time increment.
     */
    virtual void computeDistributedLoad( int           type,
                                         int           surfaceID,
                                         const double* load,
                                         double*       fExt,
                                         double*       dExt_dQ,
                                         double        timeNew,
                                         double        dT ) const override;

    /**
     * @brief Access a state variable.
     *
     * `vertex displacements` is the displacement of the material point (the single vertex); every other name is
     * forwarded to GradientEnhancedFiniteStrainMaterialPoint::getStateView.
     *
     * @param[in] stateName Name of the state variable.
     * @param[in] qp        Evaluation point (unused, there is only one).
     * @return View on the state variable.
     */
    virtual StateView getStateView( const std::string& stateName, int qp ) const override;

    /**
     * @brief Meshfree shape functions of the assigned nodes at a point.
     * @param[out] vec         Shape functions (one per assigned node).
     * @param[in]  coordinates Coordinates of the point.
     */
    virtual void getInterpolationVector( double* vec, const double* coordinates ) const override
    {
      _meshfreeApproximation.computeShapeFunctions( coordinates, _assignedKernelFunctions, vec );
    };

    // VCI:

    /**
     * @brief Number of VCI constraints (monomials) per node and direction.
     * @return _nVCIConstraints.
     */
    virtual int vci_getNumberOfConstraints() override { return _nVCIConstraints; }

    /**
     * @brief Add the boundary term @f$ T_A\,P_C\,n_i\,dA_Y @f$ of the VCI integration constraint.
     *
     * The undeformed boundary surface vector @f$ \boldsymbol{n}\,dA_0 @f$ is mapped to the intermediate reference
     * configuration by Nanson's formula with @f$ \boldsymbol{F}_n @f$; @f$ T_A @f$ and @f$ P_C @f$ are taken at the
     * particle center.
     *
     * @param[in,out] R_AiC_RowMajor        Constraint residual @f$ R_{AiC} @f$, row-major.
     * @param[in]     boundarySurfaceVector Undeformed boundary surface vector @f$ \boldsymbol{n}\,dA_0 @f$.
     * @param[in]     boundaryFaceID        Face id (unused).
     */
    virtual void vci_compute_Test_P_BoundaryIntegral( double*       R_AiC_RowMajor,
                                                      const double* boundarySurfaceVector,
                                                      int           boundaryFaceID ) override
    {
      using namespace Fastor;

      Tensor< double, nDim > n_dA0( boundarySurfaceVector ); // undeformed boundary surface vector N dA_0

      // apply Nanson's formula
      const Tensor< double, nDim, nDim > FInv = inverse( _mp.dY_dX() );
      const double                       J    = determinant( _mp.dY_dX() );

      const Tensor< double, nDim > n_dAY = J * transpose( FInv ) % n_dA0;

      for ( int A = 0; A < _nNodes; A++ )
        for ( int i = 0; i < nDim; i++ )
          for ( int C = 0; C < _nVCIConstraints; C++ )
            R_AiC_RowMajor[A * ( nDim * _nVCIConstraints ) + i * _nVCIConstraints + C] += _T( A ) * _P( C ) * n_dAY[i];
    };

    /**
     * @brief Add the domain term @f$ \partial T_A/\partial Y_i\,P_C\,V_Y @f$ of the VCI integration constraint.
     * @param[in,out] R_AiC_RowMajor Constraint residual @f$ R_{AiC} @f$, row-major.
     */
    virtual void vci_compute_TestGradient_P_Integral( double* R_AiC_RowMajor ) override
    {
      for ( int A = 0; A < _nNodes; A++ )
        for ( int i = 0; i < nDim; i++ )
          for ( int C = 0; C < _nVCIConstraints; C++ )
            R_AiC_RowMajor[A * ( nDim * _nVCIConstraints ) + i * _nVCIConstraints + C] += _dT_dY( i, A ) * _P( C ) *
                                                                                          _volReferenceIntermediate;
    };

    /**
     * @brief Add the domain term @f$ T_A\,\partial P_C/\partial Y_i\,V_Y @f$ of the VCI integration constraint.
     * @param[in,out] R_AiC_RowMajor Constraint residual @f$ R_{AiC} @f$, row-major.
     */
    virtual void vci_compute_Test_PGradient_Integral( double* R_AiC_RowMajor ) override
    {
      for ( int A = 0; A < _nNodes; A++ )
        for ( int i = 0; i < nDim; i++ )
          for ( int C = 0; C < _nVCIConstraints; C++ )
            R_AiC_RowMajor[A * ( nDim * _nVCIConstraints ) + i * _nVCIConstraints + C] += _T( A ) *
                                                                                          _P_Gradient( C, i ) *
                                                                                          _volReferenceIntermediate;
    };

    /**
     * @brief Add the contribution @f$ \chi_A\,P_C P_D\,V_Y @f$ to the VCI moment matrix of each node.
     *
     * @f$ \chi_A = 1 @f$ if the particle center lies in the support of the kernel function of node A, else 0.
     *
     * @param[in,out] mMatrix_ACD_RowMajor Moment matrices @f$ M_{ACD} @f$, row-major.
     */
    virtual void vci_compute_MMatrix( double* mMatrix_ACD_RowMajor ) override
    {
      for ( int A = 0; A < _nNodes; A++ ) {
        const double R_A = _assignedKernelFunctions[A]->isInSupport( _centerReferenceIntermediate.data() ) ? 1.0 : 0.0;

        for ( int C = 0; C < _nVCIConstraints; C++ )
          for ( int D = 0; D < _nVCIConstraints; D++ )
            mMatrix_ACD_RowMajor[A * ( _nVCIConstraints * _nVCIConstraints ) + C * _nVCIConstraints +
                                 D] += R_A * _P( C ) * _P( D ) * _volReferenceIntermediate;
      }
    };

    /**
     * @brief Correct the test function gradients, @f$ \partial T_A/\partial Y_i \mathrel{+}= \chi_A\,\eta_{AiC}\,P_C
     * @f$.
     * @param[in] eta_AiC_RowMajor Correction coefficients @f$ \eta_{AiC} @f$, row-major.
     */
    virtual void vci_assignTestFunctionCorrectionTerms( const double* eta_AiC_RowMajor ) override
    {
      for ( int A = 0; A < _nNodes; A++ ) {
        const double R_A = _assignedKernelFunctions[A]->isInSupport( _centerReferenceIntermediate.data() ) ? 1.0 : 0.0;
        for ( int i = 0; i < nDim; i++ ) {
          for ( int C = 0; C < _nVCIConstraints; C++ ) {
            _dT_dY( i, A ) += eta_AiC_RowMajor[A * ( nDim * _nVCIConstraints ) + i * _nVCIConstraints + C] * R_A *
                              _P( C );
          }
        }
      }
    };

    /**
     * @brief Undeformed volume.
     * @return @f$ V_0 @f$ of the material point.
     */
    virtual double getVolumeUndeformed() const { return _mp.getVolumeUndeformed(); };

    /**
     * @brief Apply an initial condition; `geostaticstress` is forwarded to the material point.
     * @param[in] conditionName Name of the condition.
     * @param[in] value         Its value(s).
     * @throws std::invalid_argument for any other condition.
     */
    virtual void setInitialCondition( const std::string& conditionName, const double* value ) override
    {
      if ( conditionName == "geostaticstress" ) {
        _mp.setInitialCondition( conditionName, value );
      }
      else {
        throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid initial condition" );
      }
    };

  private:
    /// @brief Set the center in the intermediate reference configuration to the accepted material point position.
    virtual void updateParticlePositionToReferenceIntermediate()
    {
      _mp.getVertexCoordinates( _centerReferenceIntermediate.data() );
    };

    /// @brief Set @f$ V_Y = V_0\det\boldsymbol{F}_n @f$.
    virtual void updateVolumeToReferenceIntermediate()
    {
      _volReferenceIntermediate = _mp.getVolumeUndeformed() * determinant( _mp.dY_dX() );
    };

    /**
     * @brief Set the VCI order, size the monomial basis and evaluate it at the current center.
     * @param[in] order Polynomial order of the VCI basis.
     */
    void setVCIOrder( int order )
    {
      _vciOrder = order;
      // number of monomials of degree <= order in nDim variables
      _nVCIConstraints = nDim == 2 ? ( order + 1 ) * ( order + 2 ) / 2
                                   : ( order + 1 ) * ( order + 2 ) * ( order + 3 ) / 6;
      _P.resize( _nVCIConstraints );
      _P_Gradient.resize( _nVCIConstraints, nDim );
      // evaluate the basis right away: VCI may run before the first accepted increment updates it
      Math::computeMonomialBasis( order, _centerReferenceIntermediate, _P );
      Math::computeMonomialBasisGradient( order, _centerReferenceIntermediate, _P_Gradient );
    };

    /**
     * @brief Coordinates of the evaluation points: the particle center.
     * @param[out] coordinates Coordinates (nDim values).
     */
    virtual void getEvaluationCoordinates( double* coordinates ) const { getVertexCoordinates( coordinates ); }

    /**
     * @brief Number of evaluation points.
     * @return 1.
     */
    virtual int getNumberOfEvaluationPoints() const
    {
      return 1; // only one evaluation point at the center of the particle
    };
  };

  template < int nDim >
  StateView GradientEnhancedFiniteStrainParticle< nDim >::getStateView( const std::string& stateName, int qp ) const
  {
    if ( stateName == "vertex displacements" )
      // the point particle's single vertex is the material point itself
      return _mp.getStateView( "displacement" );
    return _mp.getStateView( stateName );
  }

  template < int nDim >
  GradientEnhancedFiniteStrainParticle< nDim >::GradientEnhancedFiniteStrainParticle(
    int                                                  elementID,
    const double*                                        centerCoordinates0,
    int                                                  sizeCenterCoordinates0,
    double                                               volume,
    const std::string&                                   materialName,
    const double*                                        materialProperties,
    int                                                  nMaterialProperties,
    const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
    : _centerCoordinatesUndeformed( Eigen::Map< const Eigen::Matrix< double, nDim, 1 > >( centerCoordinates0 ) ),
      _centerReferenceIntermediate( _centerCoordinatesUndeformed ),
      _volReferenceIntermediate( volume ), // the reference is the undeformed configuration until the first accept
      __mp( []( int elementID_, const double* coordinates_, int nCoordinates_, double volume_ )
              -> MaterialPoints::GradientEnhancedFiniteStrainMaterialPoint< nDim >* {
        if constexpr ( nDim == 2 )
          return new MaterialPoints::GradientEnhancedFiniteStrainMaterialPoint2D( elementID_,
                                                                                  coordinates_,
                                                                                  nCoordinates_,
                                                                                  volume_ );
        else
          return new MaterialPoints::GradientEnhancedFiniteStrainMaterialPoint3D( elementID_,
                                                                                  coordinates_,
                                                                                  nCoordinates_,
                                                                                  volume_ );
      }( elementID, _centerCoordinatesUndeformed.data(), _centerCoordinatesUndeformed.size(), volume ) ),
      _mp( *__mp ),
      _meshfreeApproximation( approximation ),
      _newmark_beta( 0. ),
      _newmark_gamma( 0. ),
      _vciOrder( 0 )
  {
    if ( sizeCenterCoordinates0 != nDim ) {
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": size of center coordinates must be "
                                                << nDim << ", but got " << sizeCenterCoordinates0 );
    }

    MarmotMaterialSection section( materialName, materialProperties, nMaterialProperties );

    _mp.assignMaterial( section );

    this->setVCIOrder( _vciOrder );
  }

  template < int nDim >
  void GradientEnhancedFiniteStrainParticle< nDim >::computePhysicsKernels( const double* dQ,
                                                                            double*       fInt,
                                                                            double*       dFInt_ddQ,
                                                                            double        timeNew,
                                                                            double        dT )
  {
    using namespace Marmot::FastorIndices;
    using namespace Fastor;
    using to_jk = Fastor::OIndex< j_, k_ >;

    constexpr int nodeBlockSize = nDim + 1;

    Tensor< double, nDim >       du( 0.0 );
    Tensor< double, nDim, nDim > du_dY( 0.0 );

    double                 dn = 0.0;
    Tensor< double, nDim > dn_dY( 0.0 );

    for ( int B = 0; B < _nNodes; B++ ) {

      const int idxB_u = nodeBlockSize * B;
      const int idxB_n = nodeBlockSize * B + nDim;

      const double N_B     = _N( B );
      const auto   dN_B_dY = Tensor< double, nDim >( _dN_dY.col( B ).data() ); // works because ColumnMajor of Eigen

      const auto dQU = Tensor< double, nDim >( dQ + idxB_u );
      const auto dQN = dQ[idxB_n];

      du += N_B * dQU;
      dn += N_B * dQN;

      du_dY += einsum< i, j >( dQU, dN_B_dY );
      dn_dY += ( dQN * dN_B_dY );
    }

    _mp.prepareYourself( timeNew, dT );
    _mp.incrementDeformation( du, du_dY, dn );
    _mp.computeYourself( timeNew, dT );

    const double density0 = _mp.getDensityUndeformed();

    auto v = _mp.getVelocity();
    auto a = _mp.getAcceleration();

    Tensor< double, nDim, nDim > da_ddu( 0.0 );
    Marmot::TimeIntegration::newmarkBetaIntegration< nDim >( du.data(),
                                                             v.data(),
                                                             a.data(),
                                                             dT,
                                                             this->_newmark_beta,
                                                             this->_newmark_gamma,
                                                             da_ddu.data() );
    _mp.setVelocity( v );
    _mp.setAcceleration( a );

    Tensor< double, nDim > r_U( 0.0 );
    double                 r_N( 0.0 );

    Tensor< double, nDim, nDim > k_UU( 0.0 );
    Tensor< double, nDim >       k_UN( 0.0 );
    Tensor< double, nDim >       k_NU( 0.0 );
    double                       k_NN( 0.0 );

    const auto&  S           = _mp.response.S;
    const auto&  dLocalField = _mp.response.dL;
    const double c           = _mp.response.nonLocalRadius * _mp.response.nonLocalRadius;

    const double V0 = getVolumeUndeformed();

    const auto& t = _mp.tangents;

    Eigen::Map< Eigen::VectorXd > P( fInt, _nNodes * nodeBlockSize );
    Eigen::Map< Eigen::MatrixXd > K( dFInt_ddQ, _nNodes * nodeBlockSize, _nNodes * nodeBlockSize );

    // clang-format off
    for ( int A = 0; A < _nNodes; A++ ) {

      const double T_A = _T( A );
      const auto                   dT_A_dY = TensorMap< const double, nDim >( _dT_dY.col( A ).data() );
      const Tensor< double, nDim > dT_A_dx = einsum< ji, j >( inv( _mp.dx_dY() ), dT_A_dY );
      const Tensor< double, nDim > dT_A_dX = einsum< ji, j >( _mp.dY_dX(), dT_A_dY );

      const Tensor< double, nDim > dn_dX = einsum< ji, j >( _mp.dY_dX(), dn_dY );

        const int idxA_u = nodeBlockSize * A;
        const int idxA_n = nodeBlockSize * A + nDim;

        r_U = ( +einsum< i, ij >( dT_A_dx, S ) ) * V0;
        r_N = evaluate( ( T_A * dn + c * einsum< i, i >( dT_A_dX, dn_dX ) - T_A * dLocalField ) * V0 ).toscalar();

        // add inertia
        r_U += density0 * a * T_A * V0;

        {
          using namespace Eigen;
          P.template segment< nDim >( idxA_u ) += Map< Matrix< double, nDim, 1 > >( r_U.data() );
          P( idxA_n ) += r_N;
        }

        for ( int B = 0; B < _nNodes; B++ ) {

          const int idxB_u = nodeBlockSize * B;
          const int idxB_n = nodeBlockSize * B + nDim;

          const double N_B     = _N( B );
          const auto dN_B_dY = TensorMap< const double, nDim >( _dN_dY.col(B).data() );
          const auto dN_B_dx = evaluate( einsum< ji, j >( inv( _mp.dx_dY() ), dN_B_dY ) );
          const auto dN_B_dX = evaluate( einsum< ji, j >( _mp.dY_dX(), dN_B_dY ) ); // no dependence on current deformations!

          // aux stiffness tensors
          const auto dS_dqU_B = evaluate ( + einsum < ijkl, l > ( t.dS_dDeltaF, dN_B_dY ) );
          const auto dS_dqN_B = evaluate (                      ( t.dS_dN *      N_B    ) );
          const auto dL_dqU_B = evaluate ( + einsum <   kl, l > ( t.dL_dDeltaF, dN_B_dY ) );

          k_UU  = ( + einsum< i, ijk > ( dT_A_dx, dS_dqU_B )                              ) * V0;
          k_UN  = ( + einsum< i,  ij > ( dT_A_dx, dS_dqN_B )                              ) * V0;

          k_NU  = (                     - ( T_A * dL_dqU_B )                              ) * V0;
          k_NN  = ( + T_A * N_B + inner( dT_A_dX, dN_B_dX ) * c - T_A * N_B * t.dL_dN     ) * V0;

          k_UU += ( - einsum< k, ij, i, to_jk >( dT_A_dx, S, dN_B_dx ) ) * V0;

          k_UU += density0 * da_ddu * T_A * N_B * V0;

          {
              using namespace Eigen;
              K.template block< nDim, nDim >( idxA_u, idxB_u ) += Map< Matrix< double, nDim, nDim > >( torowmajor( k_UU ).data() );
              K.template block< nDim,    1 >( idxA_u, idxB_n ) += Map< Matrix< double, nDim,    1 > >( torowmajor( k_UN ).data() );
              K.template block<    1, nDim >( idxA_n, idxB_u ) += Map< Matrix< double,    1, nDim > >( torowmajor( k_NU ).data() );
              K                             ( idxA_n, idxB_n ) +=                                                  k_NN           ;
          }
      }
    }
    // clang-format on
  }

  template < int nDim >
  void GradientEnhancedFiniteStrainParticle< nDim >::computeDistributedLoad( int           type,
                                                                             int           surfaceID,
                                                                             const double* load,
                                                                             double*       fExt,
                                                                             double*       dExt_dQ,
                                                                             double        timeNew,
                                                                             double        dT ) const
  {
    throw std::invalid_argument( MakeString()
                                 << __PRETTY_FUNCTION__ << ": a point particle has no faces for a distributed load" );
  }

  template < int nDim >
  void GradientEnhancedFiniteStrainParticle< nDim >::computeBodyLoad( int           type,
                                                                      const double* load,
                                                                      double*       fExt,
                                                                      double*       dExt_dQ,
                                                                      double        timeNew,
                                                                      double        dT ) const
  {
    switch ( type ) {
    case BodyForce: {
      constexpr int nodeBlockSize = nDofPerNodeU + nDofPerNodeN;
      const double  V0            = getVolumeUndeformed();
      for ( int A = 0; A < this->_nNodes; A++ )
        for ( int i = 0; i < nDofPerNodeU; i++ )
          fExt[nodeBlockSize * A + i] -= this->_T( A ) * load[i] * V0;
      break;
    }
    default: throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid body load type" );
    }
  }

} // namespace Marmot::Meshfree
