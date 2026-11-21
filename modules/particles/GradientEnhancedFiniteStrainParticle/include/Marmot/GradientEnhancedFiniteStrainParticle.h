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
#include "Marmot/MarmotParticleTangent.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotUtils.h"
#include "Marmot/NewmarkBetaIntegrator.h"
#include <array>
#include <vector>

namespace Marmot::Discretizations::Meshfree {

  /**
   * @class Marmot::Discretizations::Meshfree::GradientEnhancedFiniteStrainParticle
   * @brief Nodally integrated point particle for gradient-enhanced (implicit-gradient) finite-strain materials
   *        without a micropolar continuum.
   *
   * The particle couples the displacement field @f$ \boldsymbol{u} @f$ to @p nNonlocal nonlocal fields
   * @f$ \bar{N}_K @f$ (nodal fields `displacement` and, for one nonlocal variable, `nonlocal damage`, for six,
   * `strain symmetric`; see Marmot::Discretizations::MaterialPoints::GradientEnhancedFiniteStrainNames), node block
   * @f$ n_\mathrm{dim}+n_\mathrm{nonlocal} @f$, dofs ordered node by node as
   * @f$ \{u_1,\dots,u_{n_\mathrm{dim}},\bar{N}_1,\dots,\bar{N}_{n_\mathrm{nonlocal}}\} @f$. The weak forms below
   * are written for one nonlocal variable; with several, each @f$ \bar{N}_K @f$ has its own balance with its own
   * @f$ c_K @f$, and they couple only through the material (@f$ \partial L_K/\partial\bar{N}_M @f$). It owns
   * a single Marmot::Discretizations::MaterialPoints::GradientEnhancedFiniteStrainMaterialPoint, which drives a
   * Marmot::Materials::MaterialGradientEnhancedFiniteStrain, and is integrated with that single point (direct nodal
   * integration).
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
   * (GradientEnhancedFiniteStrainMaterialPoint::response) as its source: an implicit increment starts from the
   * balanced state of the last one, so the increment form leaves out only a residual that is zero. The explicit scheme
   * writes the balance in total form instead (computePhysicsKernelsExplicit()), because an explicit increment does not
   * start from a balanced state; the material point accumulates @f$ \nabla_X\bar{N} @f$ for it. The acceleration
   * @f$ \boldsymbol{a} @f$ follows from the Newmark-beta update of @f$ \Delta\boldsymbol{u} @f$ (@f$ \beta = 0 @f$
   * switches the inertia off).
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
   * **VCI.** The particle implements the variationally consistent integration hooks of Particle: the test
   * function gradients are corrected by @f$ \partial T_A/\partial Y_i \mathrel{+}= \eta_{AiC}\,P_C(\boldsymbol{Y}) @f$
   * with a monomial basis @f$ P_C @f$ of order `VCI order` (Chen, Hillman, Rüter, 2013). The correction applies to
   * both weak forms.
   *
   * The point particle has no faces and supports no distributed loads; use GradientEnhancedFiniteStrainParticleSQCNI
   * for boundary loads. computeBodyLoad() applies a body force per unit undeformed volume, also for the derived
   * particles. The loads exist in the implicit form only; the explicit variants are those of Particle, which throw.
   *
   * @tparam nDim      Spatial dimension (2: plane strain, 3: 3D).
   * @tparam nNonlocal Number of nonlocal variables (1: a scalar nonlocal damage driver, 6: a nonlocal strain).
   */
  template < int nDim, int nNonlocal >
  class GradientEnhancedFiniteStrainParticle : public Marmot::Discretizations::Meshfree::Particle {

  protected:
    Eigen::Matrix< double, nDim, 1 > _centerCoordinatesUndeformed; ///< center @f$ \boldsymbol{X}_c @f$, undeformed
    Eigen::Matrix< double, nDim, 1 > _centerReferenceIntermediate; ///< center in the intermediate reference
                                                                   ///< configuration @f$ \boldsymbol{Y} @f$
    double _volReferenceIntermediate; ///< volume @f$ V_Y = V_0 \det\boldsymbol{F}_n @f$ (updated at each accept)

    /// the material point type
    using MaterialPoint = MaterialPoints::GradientEnhancedFiniteStrainMaterialPoint< nDim, nNonlocal >;

    MaterialPoint*                       __mp;                     ///< the owned material point
    MaterialPoint&                       _mp;                      ///< reference to the material point
    const Approximation&                 _meshfreeApproximation;   ///< the meshfree approximation
    std::vector< const KernelFunction* > _assignedKernelFunctions; ///< kernel functions of the nodes

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

    static constexpr int nDofPerNodeU = nDim;      ///< dofs per node of the displacement field U
    static constexpr int nDofPerNodeN = nNonlocal; ///< dofs per node of the nonlocal fields N

    using Material = Marmot::Materials::MaterialGradientEnhancedFiniteStrain< nNonlocal >; ///< the consumed material
    using NonlocalVector = typename Material::NonlocalVector; ///< a value per nonlocal variable

    using ForceSized = Eigen::Matrix< double, nDim, 1 >;      ///< force vector of size nDim

    /**
     * @brief Set all properties, in the order of getPropertyNames().
     * @param[in] properties  `newmark-beta beta`, `newmark-beta gamma`, `VCI order`.
     * @param[in] nProperties Number of properties; must be 3.
     * @throws std::runtime_error for a wrong number of properties.
     */
    void setProperties( const double* properties, int nProperties ) override
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
    void setProperty( const std::string& propertyName, const double* property ) override
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
    std::vector< std::string > getPropertyNames() const override { return _validProperties; };

    /**
     * @brief Number of state variables: those of the material point (including the material).
     * @return Required size of the state variable vector.
     */
    int getNumberOfRequiredStateVars() const override { return _mp.getNumberOfRequiredStateVars(); };

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
    void assignMeshfreeKernelFunctions( const std::vector< const KernelFunction* >& kernelFunctions ) override
    {
      _assignedKernelFunctions = kernelFunctions;

      _nNodes = static_cast< int >( _assignedKernelFunctions.size() );

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
     * @return @f$ n_\mathrm{dim} + n_\mathrm{nonlocal} @f$.
     */
    int getNBaseDof() const override { return nDofPerNodeU + nDofPerNodeN; }

    /**
     * @brief Fields per node.
     * @return `displacement`, then the nonlocal fields: `nonlocal damage` for one nonlocal variable, `strain
     * symmetric` for six (Marmot::Discretizations::MaterialPoints::GradientEnhancedFiniteStrainNames).
     */
    const std::vector< std::string >& getFields() const override
    {
      static const std::vector< std::string > nodeFields = [] {
        std::vector< std::string > fields = { "displacement" };
        for ( const auto& field : MaterialPoints::GradientEnhancedFiniteStrainNames< nNonlocal >::nonlocalNodeFields() )
          fields.push_back( field );
        return fields;
      }();
      return nodeFields;
    };

    /**
     * @brief Coordinates of the single vertex: the material point position of the last accepted state.
     * @param[out] coordinates Coordinates (nDim values).
     */
    void getVertexCoordinates( double* coordinates ) const override { _mp.getVertexCoordinates( coordinates ); }

    /**
     * @brief Not available: a point particle has no faces.
     * @param[in]  faceID      Face id.
     * @param[out] coordinates Unused.
     * @throws std::runtime_error always.
     */
    void getFaceCoordinates( int faceID, double* coordinates ) const override
    {
      throw std::runtime_error( "Error: GradientEnhancedFiniteStrainParticle::getFaceCoordinates not implemented." );
    }

    /**
     * @brief Center coordinates, identical to getVertexCoordinates().
     * @param[out] coordinates Coordinates (nDim values).
     */
    void getCenterCoordinates( double* coordinates ) const override { getVertexCoordinates( coordinates ); }

    /**
     * @brief Vertex coordinates for visualization, identical to getVertexCoordinates().
     * @param[out] coordinates Coordinates (nDim values).
     */
    void getVisualizationVertexCoordinates( double* coordinates ) const override
    {
      getVertexCoordinates( coordinates );
    };

    /**
     * @brief Number of vertices.
     * @return 1.
     */
    int getNumberOfVertices() const override { return 1; };

    /**
     * @brief Shape of the particle.
     * @return "point".
     */
    std::string getParticleShape() const override { return "point"; }

    /**
     * @brief Spatial dimension.
     * @return nDim.
     */
    int getDimension() const override { return nDim; };

    /**
     * @brief Construct a point particle and its material point, and assign the material.
     * @param[in] elementID           Label of the particle (also of its material point).
     * @param[in] nodeCoordinates     Undeformed center coordinates (nDim values).
     * @param[in] nNodeCoordiantes    Number of coordinates; must equal nDim.
     * @param[in] volume              Undeformed volume @f$ V_0 @f$.
     * @param[in] materialName        Name of a Marmot::Materials::MaterialGradientEnhancedFiniteStrain material with
     *                                @p nNonlocal nonlocal variables.
     * @param[in] materialProperties  Material properties.
     * @param[in] nMaterialProperties Number of material properties.
     * @param[in] approximation       Meshfree approximation used for the shape functions.
     * @throws std::invalid_argument for a wrong number of coordinates or an unsuitable material.
     */
    GradientEnhancedFiniteStrainParticle( int                  elementID,
                                          const double*        nodeCoordinates,
                                          int                  nNodeCoordiantes,
                                          double               volume,
                                          const std::string&   materialName,
                                          const double*        materialProperties,
                                          int                  nMaterialProperties,
                                          const Approximation& approximation );

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
    void acceptStateAndPosition() override
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
     * to @p fInt and @p dFInt_ddQ (dofs node by node, @f$ \{u_1,\dots,u_{n_\mathrm{dim}},\bar{N}_1,\dots\} @f$).
     *
     * @param[in]     dQ        Increment of the nodal dofs since the last accepted state.
     * @param[in,out] fInt      Internal force vector.
     * @param[in,out] dFInt_ddQ Tangent, column-major.
     * @param[in]     timeNew   Time at the end of the increment.
     * @param[in]     dT        Time increment.
     */
    void computePhysicsKernels( const double* dQ, double* fInt, double* dFInt_ddQ, double timeNew, double dT ) override;

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
    void computeBodyLoad( int type, const double* load, double* fExt, double* dExt_dQ, double timeNew, double dT )
      const override;

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
    void computeDistributedLoad( int           type,
                                 int           surfaceID,
                                 const double* load,
                                 double*       fExt,
                                 double*       dExt_dQ,
                                 double        timeNew,
                                 double        dT ) const override;

    /**
     * @brief Explicit scheme: advance the material point by the increment @p dQ.
     *
     * Interpolates @f$ \Delta\boldsymbol{u} @f$, @f$ \Delta\boldsymbol{F} @f$, @f$ \Delta\bar{N}_K @f$ and their
     * gradients as computePhysicsKernels() does, evaluates the material point and stores the rates of the increment,
     * @f$ \boldsymbol{v} = \Delta\boldsymbol{u}/\Delta t @f$, @f$ \boldsymbol{a} = (\boldsymbol{v} -
     * \boldsymbol{v}_n)/\Delta t @f$ and @f$ \dot{\bar{N}}_K = \Delta\bar{N}_K/\Delta t @f$ (not for a zero-length
     * increment). The forces follow from computePhysicsKernelsExplicit(), after the host has accepted the state and
     * assigned the kernel functions of the new position.
     *
     * @param[in] dQ      Increment of the nodal dofs, node by node.
     * @param[in] timeNew Time at the end of the increment.
     * @param[in] dT      Time increment.
     */
    void updatePhysicsExplicit( const double* dQ, double timeNew, double dT ) override;

    /**
     * @brief Explicit scheme: the internal force of the state of the last updatePhysicsExplicit().
     *
     * The displacement rows are @f$ r^U_{Aj} = \partial T_A/\partial x_i\,\tau_{ij}\,V_0 @f$, without inertia (the
     * solver divides by the lumped inertia). The nonlocal rows are the balance in TOTAL form,
     * @f[
     *   r^N_{AK} = \Bigl( T_A\,\bar{N}_K + c_K\,\frac{\partial T_A}{\partial X_i}\frac{\partial\bar{N}_K}{\partial X_i}
     *     - T_A\,L_K \Bigr) V_0 ,
     * @f]
     * with the totals the material point keeps: an explicit increment does not start from a balanced state, so the
     * increment form of computePhysicsKernels() would drop the residual of the last one. Its time derivatives,
     * @f$ m_K\ddot{\bar{N}}_K + \eta_K\dot{\bar{N}}_K @f$, are the solver's, from computeLumpedInertia() and
     * computeLumpedDamping().
     *
     * @param[in,out] fInt Internal force vector, node by node; the contribution is added.
     */
    void computePhysicsKernelsExplicit( double* fInt ) override;

    /**
     * @brief Lumped inertia, the coefficients of the second time derivatives: @f$ \rho_0\,T_A\,V_0 @f$ on the
     * displacement rows and @f$ m_K\,T_A\,V_0 @f$ on the nonlocal rows, with the micro-inertia @f$ m_K @f$ of the
     * material (zero for a nonlocal field that is first order in time).
     * @param[in,out] mLumped Lumped inertia vector, node by node; the contribution is added.
     */
    void computeLumpedInertia( double* mLumped ) const override;

    /**
     * @brief Lumped damping, the coefficients of the first time derivatives: @f$ \eta_K\,T_A\,V_0 @f$ on the nonlocal
     * rows, with the nonlocal viscosity @f$ \eta_K @f$ of the material; the displacement rows are left untouched.
     * @param[in,out] cLumped Lumped damping vector, node by node; the contribution is added.
     */
    void computeLumpedDamping( double* cLumped ) const override;

    /**
     * @brief Lumped momentum: @f$ \rho_0\,T_A\,V_0\,v_i @f$ on the displacement rows and
     * @f$ m_K\,T_A\,V_0\,\dot{\bar{N}}_K @f$ on the nonlocal rows, i.e. the lumped inertia times the rates of the last
     * explicit increment, so that a velocity projected from it keeps the rate of a second-order nonlocal field.
     * @param[in,out] mLumped Lumped momentum vector, node by node; the contribution is added.
     */
    void computeLumpedMomentum( double* mLumped ) const override;

    /**
     * @brief Add an upper bound of @f$ |\boldsymbol{K}|\boldsymbol{W} @f$ for the stiffness at the current state,
     * see Particle::computeBoundOfAbsoluteStiffnessTimes(): the material point is evaluated at a zero increment, and
     * the bound follows from addStiffnessRowsOfNode() node by node (AbsoluteTangentBound), without the stiffness.
     * Used by all particles derived from this one.
     * @param[in] columnWeights The nonnegative weights, column-major, nDof x nColumns.
     * @param[in] nColumns The number of columns of @p columnWeights.
     * @param[in,out] bound Column-major, nDof x nColumns: the bound is added.
     * @param[in] timeNew The time of the evaluation.
     * @param[in] dT The time increment of the evaluation.
     */
    void computeBoundOfAbsoluteStiffnessTimes( const double* columnWeights,
                                               int           nColumns,
                                               double*       bound,
                                               double        timeNew,
                                               double        dT ) override;
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
    Marmot::Utilities::StateView getStateView( const std::string& stateName, int qp ) const override;

    /**
     * @brief Meshfree shape functions of the assigned nodes at a point.
     * @param[out] vec         Shape functions (one per assigned node).
     * @param[in]  coordinates Coordinates of the point.
     */
    void getInterpolationVector( double* vec, const double* coordinates ) const override
    {
      _meshfreeApproximation.computeShapeFunctions( coordinates, _assignedKernelFunctions, vec );
    };

    // VCI:

    /**
     * @brief Number of VCI constraints (monomials) per node and direction.
     * @return _nVCIConstraints.
     */
    int vci_getNumberOfConstraints() override { return _nVCIConstraints; }

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
    void vci_compute_Test_P_BoundaryIntegral( double*       R_AiC_RowMajor,
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
    void vci_compute_TestGradient_P_Integral( double* R_AiC_RowMajor ) override
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
    void vci_compute_Test_PGradient_Integral( double* R_AiC_RowMajor ) override
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
    void vci_compute_MMatrix( double* mMatrix_ACD_RowMajor ) override
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
    void vci_assignTestFunctionCorrectionTerms( const double* eta_AiC_RowMajor ) override
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
    double getVolumeUndeformed() const override { return _mp.getVolumeUndeformed(); };

    /**
     * @brief Apply an initial condition; `geostaticstress` is forwarded to the material point.
     * @param[in] conditionName Name of the condition.
     * @param[in] value         Its value(s).
     * @throws std::invalid_argument for any other condition.
     */
    void setInitialCondition( const std::string& conditionName, const double* value ) override
    {
      if ( conditionName == "geostaticstress" ) {
        _mp.setInitialCondition( conditionName, value );
      }
      else {
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__ << ": invalid initial condition" );
      }
    };

  protected:
    /// The rows of one node of the tangent, factored by the chain rule (see TangentRowsOfNode).
    using TangentRows = TangentRowsOfNode< nDim, nDim + nNonlocal >;

    /**
     * @brief The rows of node @p A of the stiffness at the current response of the material point: the tangent of
     * computePhysicsKernels() without the inertia, as derivatives of the force of the node with respect to the
     * interpolated increment (TangentRowsOfNode), in the dof order of a node, @f$ [\Delta\boldsymbol{u},
     * \Delta\bar{\boldsymbol{N}}] @f$.
     * @details The forces of node @f$ A @f$ are @f$ f^U_{Aj} = \partial T_A/\partial x_i\, \tau_{ij} V_0 @f$ and
     * @f$ f^N_{AK} = \big( T_A (\bar{N}_K - L_K) + c_K\, \partial T_A/\partial X_i\, \partial \bar{N}_K/\partial X_i
     * \big) V_0 @f$, with @f$ c_K @f$ the square of the nonlocal radius. Their derivatives:
     * @f[
     *   \frac{\partial f^U_{Aj}}{\partial (\partial \Delta u_k / \partial Y_l)} = V_0 \Big(
     *   \frac{\partial T_A}{\partial x_i} \frac{\partial \tau_{ij}}{\partial \Delta F_{kl}}
     *   - \frac{\partial T_A}{\partial x_k}\, \tau_{ij} \frac{\partial Y_l}{\partial x_i} \Big), \qquad
     *   \frac{\partial f^U_{Aj}}{\partial \Delta\bar{N}_M} = V_0 \frac{\partial T_A}{\partial x_i}
     *   \frac{\partial \tau_{ij}}{\partial \bar{N}_M},
     * @f]
     * @f[
     *   \frac{\partial f^N_{AK}}{\partial (\partial \Delta u_k / \partial Y_l)} = -V_0\, T_A
     *   \frac{\partial L_K}{\partial \Delta F_{kl}}, \qquad
     *   \frac{\partial f^N_{AK}}{\partial \Delta\bar{N}_M} = V_0\, T_A \Big( \delta_{KM} -
     *   \frac{\partial L_K}{\partial \bar{N}_M} \Big), \qquad
     *   \frac{\partial f^N_{AK}}{\partial (\partial \Delta\bar{N}_K / \partial Y_l)} = V_0\, c_K
     *   \frac{\partial T_A}{\partial X_i} \frac{\partial Y_l}{\partial X_i}.
     * @f]
     * Derived particles add their own terms (e.g. a stabilization).
     * @param[in] A The node.
     * @param[in] dY_dx The inverse of the gradient of the motion of the particle, @f$ \partial Y / \partial x @f$,
     * the same for all its nodes: computed once per particle by the caller.
     * @param[in,out] rows The derivatives; the stiffness of the node is added.
     */
    virtual void addStiffnessRowsOfNode( int                                         A,
                                         const Fastor::Tensor< double, nDim, nDim >& dY_dx,
                                         TangentRows&                                rows ) const;

    /**
     * @brief The second gradients of the trial shape functions, in the layout of addTangentRowsOfNode().
     * @return An empty matrix: the point particle has none; a stabilized particle overrides it.
     */
    virtual Eigen::MatrixXd secondGradientsOfTrialFunctions() const { return Eigen::MatrixXd(); }

    /**
     * @brief Interpolate the increment of the current step from the nodal increments.
     * @param[in]  dQ    Increment of the nodal dofs, node by node.
     * @param[out] du    @f$ \Delta\boldsymbol{u} = N_B\,\Delta\boldsymbol{q}^U_B @f$.
     * @param[out] du_dY @f$ \partial\Delta u_i/\partial Y_j = \Delta q^U_{Bi}\,\partial N_B/\partial Y_j @f$.
     * @param[out] dn    @f$ \Delta\bar{N}_K = N_B\,\Delta q^N_{BK} @f$.
     * @param[out] dn_dY @f$ \partial\Delta\bar{N}_K/\partial Y_j = \Delta q^N_{BK}\,\partial N_B/\partial Y_j @f$.
     */
    void interpolateIncrement( const double*                                            dQ,
                               Fastor::Tensor< double, nDim >&                          du,
                               Fastor::Tensor< double, nDim, nDim >&                    du_dY,
                               NonlocalVector&                                          dn,
                               std::array< Fastor::Tensor< double, nDim >, nNonlocal >& dn_dY ) const;

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
    void getEvaluationCoordinates( double* coordinates ) const override { getVertexCoordinates( coordinates ); }

    /**
     * @brief Number of evaluation points.
     * @return 1.
     */
    int getNumberOfEvaluationPoints() const override
    {
      return 1; // only one evaluation point at the center of the particle
    };
  };

  template < int nDim, int nNonlocal >
  Marmot::Utilities::StateView GradientEnhancedFiniteStrainParticle< nDim, nNonlocal >::getStateView(
    const std::string& stateName,
    int                qp ) const
  {
    if ( stateName == "vertex displacements" )
      // the point particle's single vertex is the material point itself
      return _mp.getStateView( "displacement" );
    return _mp.getStateView( stateName );
  }

  template < int nDim, int nNonlocal >
  GradientEnhancedFiniteStrainParticle< nDim, nNonlocal >::GradientEnhancedFiniteStrainParticle(
    int                                                     elementID,
    const double*                                           centerCoordinates0,
    int                                                     sizeCenterCoordinates0,
    double                                                  volume,
    const std::string&                                      materialName,
    const double*                                           materialProperties,
    int                                                     nMaterialProperties,
    const Marmot::Discretizations::Meshfree::Approximation& approximation )
    : _centerCoordinatesUndeformed( Eigen::Map< const Eigen::Matrix< double, nDim, 1 > >( centerCoordinates0 ) ),
      _centerReferenceIntermediate( _centerCoordinatesUndeformed ),
      _volReferenceIntermediate( volume ), // the reference is the undeformed configuration until the first accept
      __mp( new MaterialPoint( elementID,
                               _centerCoordinatesUndeformed.data(),
                               static_cast< int >( _centerCoordinatesUndeformed.size() ),
                               volume ) ),
      _mp( *__mp ),
      _meshfreeApproximation( approximation ),
      _newmark_beta( 0. ),
      _newmark_gamma( 0. ),
      _vciOrder( 0 )
  {
    if ( sizeCenterCoordinates0 != nDim ) {
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__ << ": size of center coordinates must be " << nDim
                                   << ", but got " << sizeCenterCoordinates0 );
    }

    Marmot::Discretizations::MaterialSection section( materialName, materialProperties, nMaterialProperties );

    _mp.assignMaterial( section );

    this->setVCIOrder( _vciOrder );
  }

  template < int nDim, int nNonlocal >
  void GradientEnhancedFiniteStrainParticle< nDim, nNonlocal >::computePhysicsKernels( const double* dQ,
                                                                                       double*       fInt,
                                                                                       double*       dFInt_ddQ,
                                                                                       double        timeNew,
                                                                                       double        dT )
  {
    using namespace Marmot::Math::Tensors::FastorIndices;
    using namespace Fastor;

    constexpr int nodeBlockSize = nDim + nNonlocal;

    Tensor< double, nDim >                          du;
    Tensor< double, nDim, nDim >                    du_dY;
    NonlocalVector                                  dn;
    std::array< Tensor< double, nDim >, nNonlocal > dn_dY;
    interpolateIncrement( dQ, du, du_dY, dn, dn_dY );

    _mp.prepareYourself( timeNew, dT );
    _mp.incrementDeformation( du, du_dY, dn, dn_dY );
    _mp.computeYourself( timeNew, dT );

    const double density0 = _mp.getDensityUndeformed();

    auto v = _mp.getVelocity();
    auto a = _mp.getAcceleration();

    Tensor< double, nDim, nDim > da_ddu( 0.0 );
    Marmot::ContinuumMechanics::TimeIntegration::newmarkBetaIntegration< nDim >( du.data(),
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

    const auto& S           = _mp.response.S;
    const auto& dLocalField = _mp.response.dL;
    const auto& R           = _mp.response.nonLocalRadius;

    NonlocalVector c;
    for ( int kNl = 0; kNl < nNonlocal; kNl++ )
      c( kNl ) = R( kNl ) * R( kNl );

    const double V0 = getVolumeUndeformed();

    Eigen::Map< Eigen::VectorXd > P( fInt, _nNodes * nodeBlockSize );
    Eigen::Map< Eigen::MatrixXd > K( dFInt_ddQ, _nNodes * nodeBlockSize, _nNodes * nodeBlockSize );

    const Tensor< double, nDim, nDim > dY_dx    = inv( _mp.dx_dY() );
    const Eigen::MatrixXd              d2N_dYdY = secondGradientsOfTrialFunctions();

    TangentRows rows;
    for ( int A = 0; A < _nNodes; A++ ) {

      const double                 T_A     = _T( A );
      const auto                   dT_A_dY = TensorMap< const double, nDim >( _dT_dY.col( A ).data() );
      const Tensor< double, nDim > dT_A_dx = einsum< ji, j >( dY_dx, dT_A_dY );
      const Tensor< double, nDim > dT_A_dX = einsum< ji, j >( _mp.dY_dX(), dT_A_dY );

      const int idxA_u = nodeBlockSize * A;
      const int idxA_n = nodeBlockSize * A + nDim;

      r_U = ( +einsum< i, ij >( dT_A_dx, S ) ) * V0;

      // add inertia
      r_U += density0 * a * T_A * V0;

      P.template segment< nDim >( idxA_u ) += Eigen::Map< Eigen::Matrix< double, nDim, 1 > >( r_U.data() );

      for ( int kNl = 0; kNl < nNonlocal; kNl++ ) {
        const Tensor< double, nDim > dn_dX = einsum< ji, j >( _mp.dY_dX(), dn_dY[kNl] );

        r_N = evaluate( ( T_A * dn( kNl ) + c( kNl ) * einsum< i, i >( dT_A_dX, dn_dX ) - T_A * dLocalField( kNl ) ) *
                        V0 )
                .toscalar();

        P( idxA_n + kNl ) += r_N;
      }

      rows.setZero();
      addStiffnessRowsOfNode( A, dY_dx, rows );
      rows.df_dIncrement
        .template topLeftCorner< nDim,
                                 nDim >() += density0 * T_A * V0 *
                                             Eigen::Map< const Eigen::Matrix< double, nDim, nDim, Eigen::RowMajor > >(
                                               da_ddu.data() );
      addTangentRowsOfNode( rows, A, _N, _dN_dY, d2N_dYdY, K );
    }
  }

  template < int nDim, int nNonlocal >
  void GradientEnhancedFiniteStrainParticle< nDim, nNonlocal >::addStiffnessRowsOfNode(
    int                                         A,
    const Fastor::Tensor< double, nDim, nDim >& dY_dx,
    TangentRows&                                rows ) const
  {
    using namespace Marmot::Math::Tensors::FastorIndices;
    using namespace Fastor;

    const auto&                        S       = _mp.response.S;
    const auto&                        R       = _mp.response.nonLocalRadius;
    const auto&                        t       = _mp.tangents;
    const double                       V0      = getVolumeUndeformed();
    const double                       T_A     = _T( A );
    const Tensor< double, nDim, nDim > dY_dX   = _mp.dY_dX();
    const auto                         dT_A_dY = TensorMap< const double, nDim >( _dT_dY.col( A ).data() );
    const Tensor< double, nDim >       dT_A_dx = einsum< ji, j >( dY_dx, dT_A_dY );
    const Tensor< double, nDim >       dT_A_dX = einsum< ji, j >( dY_dX, dT_A_dY );

    // displacement rows, displacement columns: the material and the geometric part
    for ( int l = 0; l < nDim; l++ )
      for ( int j = 0; j < nDim; j++ ) {
        double geometric = 0.0;
        for ( int i = 0; i < nDim; i++ )
          geometric += S( i, j ) * dY_dx( l, i );
        for ( int k = 0; k < nDim; k++ ) {
          double material = 0.0;
          for ( int i = 0; i < nDim; i++ )
            material += dT_A_dx( i ) * t.dS_dDeltaF( i, j, k, l );
          rows.df_dIncrementGradient[l]( j, k ) += V0 * ( material - dT_A_dx( k ) * geometric );
        }
      }

    for ( int K = 0; K < nNonlocal; K++ ) {
      const int nK = nDim + K;
      // displacement rows, nonlocal column K
      for ( int j = 0; j < nDim; j++ ) {
        double dfU_dN = 0.0;
        for ( int i = 0; i < nDim; i++ )
          dfU_dN += dT_A_dx( i ) * t.dS_dN[K]( i, j );
        rows.df_dIncrement( j, nK ) += V0 * dfU_dN;
      }
      // nonlocal row K: displacement columns, and the gradient term of its own nonlocal column
      for ( int l = 0; l < nDim; l++ ) {
        for ( int k = 0; k < nDim; k++ )
          rows.df_dIncrementGradient[l]( nK, k ) -= V0 * T_A * t.dL_dDeltaF[K]( k, l );
        double dT_A_dX_dot_dYl_dX = 0.0;
        for ( int i = 0; i < nDim; i++ )
          dT_A_dX_dot_dYl_dX += dT_A_dX( i ) * dY_dX( l, i );
        rows.df_dIncrementGradient[l]( nK, nK ) += V0 * R( K ) * R( K ) * dT_A_dX_dot_dYl_dX;
      }
      // nonlocal row K, nonlocal columns M
      for ( int M = 0; M < nNonlocal; M++ )
        rows.df_dIncrement( nK, nDim + M ) += V0 * T_A * ( ( K == M ? 1.0 : 0.0 ) - t.dL_dN( K, M ) );
    }
  }

  template < int nDim, int nNonlocal >
  void GradientEnhancedFiniteStrainParticle< nDim, nNonlocal >::computeBoundOfAbsoluteStiffnessTimes(
    const double* columnWeights,
    int           nColumns,
    double*       bound,
    double        timeNew,
    double        dT )
  {
    // the material point at the current state: a zero increment
    std::array< Fastor::Tensor< double, nDim >, nNonlocal > dn_dY;
    dn_dY.fill( Fastor::Tensor< double, nDim >( 0.0 ) );
    _mp.prepareYourself( timeNew, dT );
    _mp.incrementDeformation( Fastor::Tensor< double, nDim >( 0.0 ),
                              Fastor::Tensor< double, nDim, nDim >( 0.0 ),
                              NonlocalVector::Zero(),
                              dn_dY );
    _mp.computeYourself( timeNew, dT );

    const Fastor::Tensor< double, nDim, nDim > dY_dx = inv( _mp.dx_dY() );
    addBoundOfAbsoluteTangentTimesWeights< nDim, nDim + nNonlocal >(
      columnWeights,
      nColumns,
      _N,
      _dN_dY,
      secondGradientsOfTrialFunctions(),
      [&]( int A, TangentRows& rows ) { addStiffnessRowsOfNode( A, dY_dx, rows ); },
      bound );
  }

  template < int nDim, int nNonlocal >
  void GradientEnhancedFiniteStrainParticle< nDim, nNonlocal >::interpolateIncrement(
    const double*                                            dQ,
    Fastor::Tensor< double, nDim >&                          du,
    Fastor::Tensor< double, nDim, nDim >&                    du_dY,
    NonlocalVector&                                          dn,
    std::array< Fastor::Tensor< double, nDim >, nNonlocal >& dn_dY ) const
  {
    using namespace Marmot::Math::Tensors::FastorIndices;
    using namespace Fastor;

    constexpr int nodeBlockSize = nDim + nNonlocal;

    du.zeros();
    du_dY.zeros();
    dn.setZero();
    dn_dY.fill( Tensor< double, nDim >( 0.0 ) );

    for ( int B = 0; B < _nNodes; B++ ) {

      const int idxB_u = nodeBlockSize * B;
      const int idxB_n = nodeBlockSize * B + nDim;

      const double N_B     = _N( B );
      const auto   dN_B_dY = Tensor< double, nDim >( _dN_dY.col( B ).data() ); // works because ColumnMajor of Eigen

      const auto dQU = Tensor< double, nDim >( dQ + idxB_u );

      du += N_B * dQU;
      du_dY += einsum< i, j >( dQU, dN_B_dY );

      for ( int kNl = 0; kNl < nNonlocal; kNl++ ) {
        const double dQN = dQ[idxB_n + kNl];
        dn( kNl ) += N_B * dQN;
        dn_dY[kNl] += ( dQN * dN_B_dY );
      }
    }
  }

  template < int nDim, int nNonlocal >
  void GradientEnhancedFiniteStrainParticle< nDim, nNonlocal >::updatePhysicsExplicit( const double* dQ,
                                                                                       double        timeNew,
                                                                                       double        dT )
  {
    Fastor::Tensor< double, nDim >                          du;
    Fastor::Tensor< double, nDim, nDim >                    du_dY;
    NonlocalVector                                          dn;
    std::array< Fastor::Tensor< double, nDim >, nNonlocal > dn_dY;
    interpolateIncrement( dQ, du, du_dY, dn, dn_dY );

    _mp.prepareYourself( timeNew, dT );
    _mp.incrementDeformation( du, du_dY, dn, dn_dY );
    _mp.computeYourself( timeNew, dT );

    // a zero-length increment (the start of a step) leaves the rates as they are
    if ( dT <= 0.0 )
      return;

    const Fastor::Tensor< double, nDim > v_n  = _mp.getVelocity();
    const Fastor::Tensor< double, nDim > v_np = du / dT;
    _mp.setVelocity( v_np );
    _mp.setAcceleration( Fastor::Tensor< double, nDim >( ( v_np - v_n ) / dT ) );
    _mp.setNonlocalFieldRate( dn / dT );
  }

  template < int nDim, int nNonlocal >
  void GradientEnhancedFiniteStrainParticle< nDim, nNonlocal >::computePhysicsKernelsExplicit( double* fInt )
  {
    using namespace Marmot::Math::Tensors::FastorIndices;
    using namespace Fastor;

    constexpr int nodeBlockSize = nDim + nNonlocal;

    const auto&                                           S     = _mp.response.S;
    const auto&                                           R     = _mp.response.nonLocalRadius;
    const NonlocalVector                                  N     = _mp.getNonlocalField();
    const NonlocalVector                                  L     = _mp.getLocalField();
    const std::array< Tensor< double, nDim >, nNonlocal > dN_dX = _mp.getNonlocalFieldGradient();
    const Tensor< double, nDim, nDim >                    dY_dx = inv( _mp.dx_dY() );
    const Tensor< double, nDim, nDim >                    dY_dX = _mp.dY_dX();
    const double                                          V0    = getVolumeUndeformed();

    Eigen::Map< Eigen::VectorXd > P( fInt, _nNodes * nodeBlockSize );

    for ( int A = 0; A < _nNodes; A++ ) {

      const double                 T_A     = _T( A );
      const auto                   dT_A_dY = TensorMap< const double, nDim >( _dT_dY.col( A ).data() );
      const Tensor< double, nDim > dT_A_dx = einsum< ji, j >( dY_dx, dT_A_dY );
      const Tensor< double, nDim > dT_A_dX = einsum< ji, j >( dY_dX, dT_A_dY );

      const int idxA_u = nodeBlockSize * A;
      const int idxA_n = nodeBlockSize * A + nDim;

      const Tensor< double, nDim > r_U = einsum< i, ij >( dT_A_dx, S ) * V0;
      P.template segment< nDim >( idxA_u ) += Eigen::Map< const Eigen::Matrix< double, nDim, 1 > >( r_U.data() );

      for ( int kNl = 0; kNl < nNonlocal; kNl++ ) {
        const double c = R( kNl ) * R( kNl );
        P( idxA_n + kNl ) += ( T_A * N( kNl ) + c * inner( dT_A_dX, dN_dX[kNl] ) - T_A * L( kNl ) ) * V0;
      }
    }
  }

  template < int nDim, int nNonlocal >
  void GradientEnhancedFiniteStrainParticle< nDim, nNonlocal >::computeLumpedInertia( double* mLumped ) const
  {
    constexpr int nodeBlockSize = nDim + nNonlocal;

    const double         density0 = _mp.getDensityUndeformed();
    const NonlocalVector m        = _mp.getNonlocalMicroInertia();
    const double         V0       = getVolumeUndeformed();

    for ( int A = 0; A < _nNodes; A++ ) {
      const double T_A = _T( A );
      for ( int i = 0; i < nDim; i++ )
        mLumped[nodeBlockSize * A + i] += density0 * T_A * V0;
      for ( int kNl = 0; kNl < nNonlocal; kNl++ )
        mLumped[nodeBlockSize * A + nDim + kNl] += m( kNl ) * T_A * V0;
    }
  }

  template < int nDim, int nNonlocal >
  void GradientEnhancedFiniteStrainParticle< nDim, nNonlocal >::computeLumpedDamping( double* cLumped ) const
  {
    constexpr int nodeBlockSize = nDim + nNonlocal;

    const NonlocalVector eta = _mp.getNonlocalViscosity();
    const double         V0  = getVolumeUndeformed();

    for ( int A = 0; A < _nNodes; A++ ) {
      const double T_A = _T( A );
      for ( int kNl = 0; kNl < nNonlocal; kNl++ )
        cLumped[nodeBlockSize * A + nDim + kNl] += eta( kNl ) * T_A * V0;
    }
  }

  template < int nDim, int nNonlocal >
  void GradientEnhancedFiniteStrainParticle< nDim, nNonlocal >::computeLumpedMomentum( double* mLumped ) const
  {
    constexpr int nodeBlockSize = nDim + nNonlocal;

    const double                         density0 = _mp.getDensityUndeformed();
    const NonlocalVector                 m        = _mp.getNonlocalMicroInertia();
    const Fastor::Tensor< double, nDim > v        = _mp.getVelocity();
    const NonlocalVector                 dN_dt    = _mp.getNonlocalFieldRate();
    const double                         V0       = getVolumeUndeformed();

    for ( int A = 0; A < _nNodes; A++ ) {
      const double T_A = _T( A );
      for ( int i = 0; i < nDim; i++ )
        mLumped[nodeBlockSize * A + i] += density0 * T_A * V0 * v( i );
      for ( int kNl = 0; kNl < nNonlocal; kNl++ )
        mLumped[nodeBlockSize * A + nDim + kNl] += m( kNl ) * T_A * V0 * dN_dt( kNl );
    }
  }

  template < int nDim, int nNonlocal >
  void GradientEnhancedFiniteStrainParticle< nDim, nNonlocal >::computeDistributedLoad( int           type,
                                                                                        int           surfaceID,
                                                                                        const double* load,
                                                                                        double*       fExt,
                                                                                        double*       dExt_dQ,
                                                                                        double        timeNew,
                                                                                        double        dT ) const
  {
    throw std::invalid_argument( Marmot::Utilities::MakeString()
                                 << __PRETTY_FUNCTION__ << ": a point particle has no faces for a distributed load" );
  }

  template < int nDim, int nNonlocal >
  void GradientEnhancedFiniteStrainParticle< nDim, nNonlocal >::computeBodyLoad( int           type,
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
    default:
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__ << ": invalid body load type" );
    }
  }

} // namespace Marmot::Discretizations::Meshfree
