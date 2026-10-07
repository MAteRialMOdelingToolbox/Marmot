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

#include "Marmot/DisplacementMaterialPoint.h"
#include "Marmot/GenericParticle.h" // New base class
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotMeshfreeApproximation.h"
#include "Marmot/MarmotMonomialBasisFunctions.h"
#include "Marmot/MarmotParticle.h"
#include "Marmot/MarmotUtils.h"
#include "Marmot/NewmarkBetaIntegrator.h"
#include <vector>

/// the displacement material point of dimension nDim: DisplacementMaterialPoint2D (plane strain) or
/// DisplacementMaterialPoint3D
template < int nDim >
using MaterialPointType = std::conditional_t< nDim == 2,
                                              Marmot::MaterialPoints::DisplacementMaterialPoint2D,
                                              Marmot::MaterialPoints::DisplacementMaterialPoint3D >;

namespace Marmot::Meshfree {

  /**
   * @class Marmot::Meshfree::DisplacementParticle
   * @brief RKPM particle of the finite-strain displacement formulation, integrated at its center (direct nodal
   * integration); registered as "Displacement/PlaneStrain/Point".
   *
   * @details The particle owns one DisplacementMaterialPoint (see MaterialPointType), which carries the kinematics and
   * the state of a MarmotMaterialFiniteStrain; configurations and deformation gradients are those of the material
   * point:
   * @f$ \boldsymbol{X} @f$ undeformed, @f$ \boldsymbol{Y} @f$ intermediate (last accepted state), @f$ \boldsymbol{x}
   * @f$ current, @f$ \boldsymbol{F} = \Delta\boldsymbol{F}\,\boldsymbol{F}_n @f$. The dofs @f$ \Delta q_{Bk} @f$ are
   * the increments of the current step at the nodes (kernel functions) of the particle. With the trial functions
   * @f$ N_B @f$ and the test functions @f$ T_A @f$ (identical, unless VCI corrects the test gradients) and their
   * gradients with respect to @f$ \boldsymbol{Y} @f$, the residual and the tangent are
   * @f[
   *   \Delta F_{iJ} = \delta_{iJ} + \Delta q_{Bi}\,\frac{\partial N_B}{\partial Y_J}, \qquad
   *   r_{Aj} = \frac{\partial T_A}{\partial x_i}\,\tau_{ij}\,V_0 + \rho_0\,a_j\,T_A\,V_0,
   * @f]
   * @f[
   *   \frac{\partial r_{Aj}}{\partial \Delta q_{Bk}} =
   *     \left( \frac{\partial T_A}{\partial x_i}\,\frac{\partial \tau_{ij}}{\partial \Delta F_{kL}}\,
   *     \frac{\partial N_B}{\partial Y_L}
   *     - \frac{\partial T_A}{\partial x_k}\,\tau_{ij}\,\frac{\partial N_B}{\partial x_i} \right) V_0
   *     + \rho_0\,\frac{\partial a_j}{\partial \Delta u_k}\,T_A\,N_B\,V_0,
   * @f]
   * with @f$ \partial(\bullet)/\partial x_i = \Delta F^{-1}_{Ji}\,\partial(\bullet)/\partial Y_J @f$, the undeformed
   * volume @f$ V_0 @f$, the undeformed density @f$ \rho_0 @f$ and the acceleration @f$ \boldsymbol{a} @f$ of a
   * Newmark-beta integration of the center displacement increment (no inertia for the default
   * @f$ \beta = 0 @f$). The second term of the tangent is the geometric stiffness.
   *
   * The derived classes DisplacementParticleSQCNI (SQCNI, SNNI) and DisplacementParticleSQCNIxNSNI (NSNI) change
   * where and how the shape functions and their gradients are evaluated; the subdomain-integrated
   * DisplacementParticleSQCNIxSDI is not derived from this class, but uses the same material point and kernel.
   *
   * Properties (in this order): "VCI order" (from GenericParticle), "newmark-beta beta", "newmark-beta gamma".
   *
   * @tparam nDim Spatial dimension (2: plane strain, 3: 3D).
   */
  template < int nDim >
  class DisplacementParticle : public Marmot::Meshfree::GenericParticle< nDim > {

    using TensorD  = Fastor::Tensor< double, nDim >;       ///< vector of size nDim
    using TensorDD = Fastor::Tensor< double, nDim, nDim >; ///< second-order tensor of size nDim

  protected:
    std::unique_ptr< MaterialPointType< nDim > > _mp; ///< the material point at the particle center

    double _newmark_beta;  ///< Newmark parameter @f$ \beta @f$ (property "newmark-beta beta", default 0)
    double _newmark_gamma; ///< Newmark parameter @f$ \gamma @f$ (property "newmark-beta gamma", default 0)

    /// static vector of valid properties
    inline static const std::vector< std::string > _validProperties = {
      "newmark-beta beta",
      "newmark-beta gamma",
    };

    /**
     * @brief Deformation gradient of the last accepted state.
     * @return @f$ \boldsymbol{F}_n = \partial\boldsymbol{Y}/\partial\boldsymbol{X} @f$ of the material point.
     */
    virtual TensorDD dY_dX() const { return ( this->_mp->dY_dX() ); }

    /**
     * @brief Incremental deformation gradient.
     * @return @f$ \Delta\boldsymbol{F} = \partial\boldsymbol{x}/\partial\boldsymbol{Y} @f$ of the material point.
     */
    virtual TensorDD dx_dY() const { return ( this->_mp->dx_dY() ); }

    /**
     * @brief Total displacement of the particle center at the last accepted state.
     * @return The displacement @f$ \boldsymbol{u} @f$.
     */
    virtual TensorD getDisplacementAtCenter() const
    {
      TensorD u( 0.0 );
      this->_mp->getCenterDisplacement( u.data() );
      return u;
    }

  public:
    /// Body load types.
    enum BodyLoadTypes {
      BodyForce, ///< body force ("BODYFORCE")
    };

    /// Distributed load types.
    enum DistributedLoadTypes {
      Pressure,     ///< follower pressure on a particle face ("PRESSURE")
      CWFCorrection ///< consistent weak form boundary correction on a particle face ("CWFCORRECTION")
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
     * @brief Supported distributed loads: none, a point particle has no faces. The particles with a domain
     * (DisplacementParticleSQCNI and its derivatives) support "PRESSURE" and "CWFCORRECTION".
     * @return An empty map.
     */
    const std::unordered_map< std::string, int >& getSupportedDistributedLoadTypes() const override
    {
      static const std::unordered_map< std::string, int > _supportedDistributedLoadTypes = {};
      return _supportedDistributedLoadTypes;
    };

    static constexpr int nDofPerNodeU = nDim;            ///< dofs per node of the displacement field

    using Material = MarmotMaterialFiniteStrain;         ///< material interface consumed by the material point

    using ForceSized = Eigen::Matrix< double, nDim, 1 >; ///< a force vector

    /**
     * @brief Sets all properties, in the order of getPropertyNames().
     * @param[in] properties The property values.
     * @param[in] nProperties Number of values; must equal the number of property names.
     * @throws std::runtime_error if nProperties does not match.
     */
    virtual void setProperties( const double* properties, int nProperties ) override
    {
      // Combine property names from base and derived
      std::vector< std::string > allPropertyNames = Marmot::Meshfree::GenericParticle< nDim >::getPropertyNames();
      allPropertyNames.insert( allPropertyNames.end(), _validProperties.begin(), _validProperties.end() );

      if ( nProperties != static_cast< int >( allPropertyNames.size() ) ) {
        std::ostringstream oss;
        oss << "Error in " << __PRETTY_FUNCTION__ << ": ";
        oss << "Expected " << allPropertyNames.size() << " properties, but got " << nProperties << ". ";
        oss << "Valid properties are: ";
        for ( const auto& prop : allPropertyNames ) {
          oss << prop << ", ";
        }
        throw std::runtime_error( oss.str() );
      }

      for ( int i = 0; i < nProperties; i++ ) {
        setProperty( allPropertyNames[i], &properties[i] );
      }
    };

    /**
     * @brief Sets a single property; unknown names are passed to GenericParticle::setProperty().
     * @param[in] propertyName Name of the property.
     * @param[in] property Its value.
     * @throws std::runtime_error (from GenericParticle) for an unknown name.
     */
    virtual void setProperty( const std::string& propertyName, const double* property ) override
    {
      if ( propertyName == "newmark-beta beta" ) {
        _newmark_beta = property[0];
      }
      else if ( propertyName == "newmark-beta gamma" ) {
        _newmark_gamma = property[0];
      }
      else {
        // If not a DisplacementParticle specific property, try the base class
        Marmot::Meshfree::GenericParticle< nDim >::setProperty( propertyName, property );
      }
    };

    /**
     * @brief Names of the properties.
     * @return "VCI order", "newmark-beta beta", "newmark-beta gamma".
     */
    virtual std::vector< std::string > getPropertyNames() const override
    {
      std::vector< std::string > names = Marmot::Meshfree::GenericParticle< nDim >::getPropertyNames();
      names.insert( names.end(), _validProperties.begin(), _validProperties.end() );
      return names;
    };

    /**
     * @brief Number of state variables: those of the material point, including the material state.
     * @return The number of state variables.
     */
    virtual int getNumberOfRequiredStateVars() const override { return _mp->getNumberOfRequiredStateVars(); };

    /**
     * @brief Assigns the state vector to the material point.
     * @param[in] stateVars State vector, owned by the host.
     * @param[in] nStateVars Its length.
     */
    void assignStateVars( double* stateVars, int nStateVars ) override
    {
      _mp->assignStateVars( stateVars, nStateVars );
    }

    /**
     * @brief Number of dofs per node.
     * @return nDim.
     */
    virtual int getNBaseDof() const { return nDofPerNodeU; }

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
     * @brief Constructs the particle and its material point.
     * @param[in] elementID Label of the particle (also the label of the material point).
     * @param[in] nodeCoordinates Coordinates of the particle center in the undeformed configuration (nDim values).
     * @param[in] nNodeCoordiantes Number of coordinates; must be nDim.
     * @param[in] volume Volume of the particle in the undeformed configuration.
     * @param[in] materialName Name of the finite-strain material.
     * @param[in] materialProperties Material properties.
     * @param[in] nMaterialProperties Number of material properties.
     * @param[in] approximation The meshfree approximation for the shape functions.
     * @throws std::invalid_argument for a wrong number of coordinates or an unknown material.
     */
    DisplacementParticle( int                                elementID,
                          const double*                      nodeCoordinates,
                          int                                nNodeCoordiantes,
                          double                             volume,
                          const std::string&                 materialName,
                          const double*                      materialProperties,
                          int                                nMaterialProperties,
                          const MarmotMeshfreeApproximation& approximation );

    /**
     * @brief Initializes the material point.
     */
    void initializeYourself() override { _mp->initializeYourself(); };

    /**
     * @brief Accepts the increment: updates the material point, the intermediate volume
     * @f$ V_Y = V_0 \det\boldsymbol{F}_n @f$ and position, and the VCI monomial basis at the new position.
     */
    virtual void acceptStateAndPosition() override
    {
      _mp->acceptStateAndPosition();
      _mp->prepareYourself( 0, 0 );

      this->updateVolumeToReferenceIntermediate();
      this->updateParticlePositionToReferenceIntermediate();

      // Use base class members for VCI
      Math::computeMonomialBasis( this->_vciOrder, this->_centerReferenceIntermediate, this->_P );
      Math::computeMonomialBasisGradient( this->_vciOrder, this->_centerReferenceIntermediate, this->_P_Gradient );
    };

    /**
     * @brief Updates the material point with the increment dQ and assembles the residual and the tangent of the class
     * description. Also updates velocity and acceleration of the material point by the Newmark-beta scheme.
     * @param[in] dQ Nodal displacement increments of the current step (nDim values per node).
     * @param[in,out] fInt Residual (internal and inertia forces), the contribution is added.
     * @param[in,out] dFInt_ddQ Tangent @f$ \partial r/\partial\Delta q @f$, the contribution is added.
     * @param[in] timeNew Time at the end of the increment.
     * @param[in] dT Time increment.
     */
    virtual void computePhysicsKernels( const double* dQ,
                                        double*       fInt,
                                        double*       dFInt_ddQ,
                                        double        timeNew,
                                        double        dT ) override;

    /**
     * @brief Body load: a body force @f$ \boldsymbol{b} @f$ per unit undeformed volume (dead load),
     * @f$ P_{Ai} \mathrel{-}= T_A\,b_i\,V_0 @f$ (the host's sign convention for external loads, as in the cells);
     * the tangent is zero. Used by all particles derived from this one.
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
                                  double        dT ) const override;

    /**
     * @brief Distributed load: a point particle has no faces, so there is none (see
     * getSupportedDistributedLoadTypes()); DisplacementParticleSQCNI implements them.
     * @throws std::invalid_argument always.
     * @param[in] type The distributed load type.
     * @param[in] surfaceID The face of the particle.
     * @param[in] load The load values.
     * @param[in,out] fExt Load vector (not modified).
     * @param[in,out] dExt_dQ Tangent (not modified).
     * @param[in] timeNew Time at the end of the increment.
     * @param[in] dT Time increment.
     */
    virtual void computeDistributedLoad( int           type,
                                         int           surfaceID,
                                         const double* load,
                                         double*       fExt,
                                         double*       dExt_dQ,
                                         double        timeNew,
                                         double        dT ) const override;

    /**
     * @brief Lumped mass, @f$ m_{Ai} \mathrel{+}= \rho_0\,T_A\,V_0 @f$.
     * @param[in,out] mLumped Lumped mass vector (nDim values per node), the contribution is added.
     */
    virtual void computeLumpedInertia( double* mLumped ) const override;

    /**
     * @brief Lumped momentum, @f$ p_{Ai} \mathrel{+}= \rho_0\,T_A\,V_0\,v_i @f$.
     * @param[in,out] mLumped Lumped momentum vector (nDim values per node), the contribution is added.
     */
    virtual void computeLumpedMomentum( double* mLumped ) const override;

    /**
     * @brief State of the material point (see DisplacementMaterialPoint::getStateView()).
     * @param[in] stateName Name of the state.
     * @param[in] qp Evaluation point (not used; the particle has a single material point).
     * @return The view on the state.
     */
    virtual StateView getStateView( const std::string& stateName, int qp ) const override;

    // VCI methods are now in GenericParticle, but vci_compute_Test_P_BoundaryIntegral needs override
    // because it depends on dY_dX() which is physics-specific.
    /**
     * @brief VCI boundary term @f$ R_{AiC} \mathrel{+}= T_A\,P_C\,(N\,dA_Y)_i @f$, with the boundary vector
     * transformed from the undeformed to the intermediate configuration by Nanson's formula,
     * @f$ \boldsymbol{N}\,dA_Y = J_n\,\boldsymbol{F}_n^{-\mathsf T}\,\boldsymbol{N}\,dA_0 @f$.
     * @param[in,out] R_AiC_RowMajor VCI matrix (nNodes x nDim x nVCIConstraints, row major), the contribution is
     * added.
     * @param[in] boundarySurfaceVector Boundary surface vector @f$ \boldsymbol{N}\,dA_0 @f$ (nDim values).
     * @param[in] boundaryFaceID Face ID (not used).
     */
    virtual void vci_compute_Test_P_BoundaryIntegral( double*       R_AiC_RowMajor,
                                                      const double* boundarySurfaceVector,
                                                      int           boundaryFaceID ) override
    {
      using namespace Fastor;

      Tensor< double, nDim > n_dA0( boundarySurfaceVector ); // undeformed load vector p * N_I * dA_0

      // apply Nanson's formula
      const Tensor< double, nDim, nDim > FInv = inverse( dY_dX() );
      const double                       J    = determinant( dY_dX() );

      const Tensor< double, nDim > n_dAY = J * transpose( FInv ) % n_dA0;

      for ( int A = 0; A < this->_nNodes; A++ )                            // Use base class _nNodes
        for ( int i = 0; i < nDim; i++ )
          for ( int C = 0; C < this->_nVCIConstraints; C++ )               // Use base class _nVCIConstraints
            R_AiC_RowMajor[A * ( nDim * this->_nVCIConstraints ) + i * this->_nVCIConstraints +
                           C] += this->_T( A ) * this->_P( C ) * n_dAY[i]; // Use base class _T, _P
    };

    /**
     * @brief Volume in the undeformed configuration.
     * @return @f$ V_0 @f$ of the material point.
     */
    virtual double getVolumeUndeformed() const { return _mp->getVolumeUndeformed(); };

    /**
     * @brief Initial conditions are forwarded to the material point, which supports none.
     * @param[in] conditionName Name of the initial condition.
     * @param[in] value Its values.
     * @throws std::invalid_argument always (from DisplacementMaterialPoint::setInitialCondition()).
     */
    virtual void setInitialCondition( const std::string& conditionName, const double* value ) override
    {
      _mp->setInitialCondition( conditionName, value );
    };

  private:
    /**
     * @brief Sets the particle center of GenericParticle to the intermediate position
     * @f$ \boldsymbol{Y} = \boldsymbol{X} + \boldsymbol{u} @f$ of the material point.
     */
    virtual void updateParticlePositionToReferenceIntermediate()
    {
      _mp->getVertexCoordinates(
        this->_centerReferenceIntermediate.data() ); // Use base class _centerReferenceIntermediate
    };

    /**
     * @brief Sets the intermediate volume of GenericParticle, @f$ V_Y = V_0 \det\boldsymbol{F}_n @f$.
     */
    virtual void updateVolumeToReferenceIntermediate()
    {
      this->_volReferenceIntermediate = getVolumeUndeformed() *
                                        determinant( dY_dX() ); // Use base class _volReferenceIntermediate
    };
  };

  template < int nDim >
  StateView DisplacementParticle< nDim >::getStateView( const std::string& stateName, int qp ) const
  {
    return _mp->getStateView( stateName );
  }

  template < int nDim >
  DisplacementParticle< nDim >::DisplacementParticle(
    int                                                  elementID,
    const double*                                        centerCoordinates0,
    int                                                  sizeCenterCoordinates0,
    double                                               volume,
    const std::string&                                   materialName,
    const double*                                        materialProperties,
    int                                                  nMaterialProperties,
    const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
    : Marmot::Meshfree::GenericParticle< nDim >( elementID,
                                                 centerCoordinates0,
                                                 sizeCenterCoordinates0,
                                                 volume,
                                                 approximation ), // Call virtual base constructor
      _mp( std::make_unique< MaterialPointType< nDim > >( elementID,
                                                          Eigen::Map< const Eigen::Matrix< double, nDim, 1 > >(
                                                            centerCoordinates0 )
                                                            .data(),
                                                          sizeCenterCoordinates0,
                                                          volume ) ),
      _newmark_beta( 0. ),
      _newmark_gamma( 0. )
  {
    MarmotMaterialSection section( materialName, materialProperties, nMaterialProperties );

    _mp->assignMaterial( section );
  }

  template < int nDim >
  void DisplacementParticle< nDim >::computePhysicsKernels( const double* dQ,
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

    Tensor< double, nDim > du( 0.0 );

    Tensor< double, nDim, nDim > du_dY( 0.0 );

    for ( int B = 0; B < this->_nNodes; B++ ) { // Use base class _nNodes

      const int idxB_u = nodeBlockSize * B;

      const double N_B     = this->_N( B );                                          // Use base class _N
      const auto   dN_B_dY = Tensor< double, nDim >( this->_dN_dY.col( B ).data() ); // Use base class _dN_dY

      const auto dQU = Tensor< double, nDim >( dQ + idxB_u );

      du += N_B * dQU;

      du_dY += einsum< i, j >( dQU, dN_B_dY );
    }

    _mp->prepareYourself( timeNew, dT );
    _mp->incrementDeformation( du, du_dY );
    _mp->computeYourself( timeNew, dT );

    const double density0 = _mp->getDensityUndeformed();

    auto v = _mp->getVelocity();
    auto a = _mp->getAcceleration();

    Tensor< double, nDim, nDim > da_ddu( 0.0 );
    Marmot::TimeIntegration::newmarkBetaIntegration< nDim >( du.data(),
                                                             v.data(),
                                                             a.data(),
                                                             dT,
                                                             this->_newmark_beta,
                                                             this->_newmark_gamma,
                                                             da_ddu.data() );
    _mp->setVelocity( v );
    _mp->setAcceleration( a );

    Tensor< double, nDim > r_U( 0.0 );

    Tensor< double, nDim, nDim > k_UU( 0.0 );

    const auto& S = _mp->response.S;

    const double V0 = getVolumeUndeformed();

    const auto& t = _mp->tangents;

    Eigen::Map< Eigen::VectorXd > P( fInt, this->_nNodes * nodeBlockSize ); // Use base class _nNodes
    Eigen::Map< Eigen::MatrixXd > K( dFInt_ddQ,
                                     this->_nNodes * nodeBlockSize,
                                     this->_nNodes * nodeBlockSize ); // Use base class _nNodes

    // clang-format off
    for ( int A = 0; A < this->_nNodes; A++ ) { // Use base class _nNodes

      const double T_A = this->_T( A ); // Use base class _T
      const auto                   dT_A_dY = TensorMap< const double, nDim >( this->_dT_dY.col( A ).data() ); // Use base class _dT_dY
      const Tensor< double, nDim > dT_A_dx = einsum< ji, j >( inv( _mp->dx_dY() ), dT_A_dY );

        const int idxA_u = nodeBlockSize * A;

        r_U = ( +einsum< i, ij >( dT_A_dx, S ) ) * V0;

        // add inertia
        r_U += density0 * a * T_A * V0;

        {
          using namespace Eigen;
          P.template segment< nDim >( idxA_u ) += Map< Matrix< double, nDim, 1 > >( r_U.data() );
        }

        for ( int B = 0; B < this->_nNodes; B++ ) { // Use base class _nNodes

          const int idxB_u = nodeBlockSize * B;

          const double                 N_B     = this->_N( B ); // Use base class _N
          const auto dN_B_dY = TensorMap< const double, nDim >( this->_dN_dY.col(B).data() ); // Use base class _dN_dY
          const auto dN_B_dx = evaluate( einsum< ji, j >( inv( _mp->dx_dY() ), dN_B_dY ) );

          // aux stiffness tensors
          const auto dS_dqU_B = evaluate ( + einsum < ijkl, l > ( t.dS_dDeltaF, dN_B_dY )                                            );

          k_UU  = ( + einsum< i, ijk        >  ( dT_A_dx, dS_dqU_B )   ) * V0;
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

  template < int nDim >
  void DisplacementParticle< nDim >::computeDistributedLoad( int           type,
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
  void DisplacementParticle< nDim >::computeBodyLoad( int           type,
                                                      const double* load,
                                                      double*       fExt,
                                                      double*       dExt_dQ,
                                                      double        timeNew,
                                                      double        dT ) const
  {
    switch ( type ) {
    case BodyForce: {
      const double V0 = getVolumeUndeformed();
      for ( int A = 0; A < this->_nNodes; A++ )
        for ( int i = 0; i < nDofPerNodeU; i++ )
          fExt[nDofPerNodeU * A + i] -= this->_T( A ) * load[i] * V0;
      break;
    }
    default: throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid body load type" );
    }
  }

  template < int nDim >
  void DisplacementParticle< nDim >::computeLumpedInertia( double* mLumped ) const
  {
    const double density0 = _mp->getDensityUndeformed();
    const double V0       = getVolumeUndeformed();

    for ( int A = 0; A < this->_nNodes; A++ ) { // Use base class _nNodes
      const double T_A    = this->_T( A );      // Use base class _T
      const int    idxA_u = nDofPerNodeU * A;

      for ( int i = 0; i < nDofPerNodeU; i++ ) {
        mLumped[idxA_u + i] += density0 * T_A * V0;
      }
    }
  }

  template < int nDim >
  void DisplacementParticle< nDim >::computeLumpedMomentum( double* mLumped ) const
  {
    const double density0 = _mp->getDensityUndeformed();
    const double V0       = getVolumeUndeformed();
    const auto   v        = _mp->getVelocity();

    for ( int A = 0; A < this->_nNodes; A++ ) { // Use base class _nNodes
      const double T_A    = this->_T( A );      // Use base class _T
      const int    idxA_u = nDofPerNodeU * A;

      for ( int i = 0; i < nDofPerNodeU; i++ ) {
        mLumped[idxA_u + i] += density0 * T_A * V0 * v[i];
      }
    }
  }

} // namespace Marmot::Meshfree
