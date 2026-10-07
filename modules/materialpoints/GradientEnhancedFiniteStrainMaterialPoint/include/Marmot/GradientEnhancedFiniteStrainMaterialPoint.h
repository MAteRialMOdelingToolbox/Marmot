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
#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotGeostaticStress.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotLowerDimensionalStress.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrain.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrainFactory.h"
#include "Marmot/MarmotMaterialPoint.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

#include <memory>
#include <vector>

namespace Marmot::MaterialPoints {

  /**
   * @class Marmot::MaterialPoints::GradientEnhancedFiniteStrainMaterialPoint
   * @brief Material point for gradient-enhanced (implicit-gradient) finite-strain materials without a micropolar
   *        continuum.
   *
   * The material point carries the displacement field @f$ \boldsymbol{u} @f$ and a single scalar nonlocal field
   * @f$ \bar{N} @f$ and drives a MarmotMaterialGradientEnhancedFiniteStrain, which returns the Kirchhoff stress
   * @f$ \boldsymbol{\tau} @f$, the local driving force @f$ L @f$ (the source of @f$ \bar{N} - c\,\nabla^2\bar{N} = L
   * @f$, @f$ c = R^2 @f$ with the nonlocal radius @f$ R @f$) and the four algorithmic tangents. It is the material
   * point both of the MPM cell GradientEnhancedFiniteStrainCell and of the RKPM particle
   * GradientEnhancedFiniteStrainParticle (and its SQCNI / NSNI variants). There is no micro-rotation field and hence
   * no couple stress.
   *
   * **Kinematics (semi-Lagrangian).** The state variable `deformation gradient` (@c dY_dX) holds the deformation
   * gradient @f$ \boldsymbol{F}_n = \partial\boldsymbol{Y}/\partial\boldsymbol{X} @f$ of the last accepted state,
   * i.e. of the intermediate reference configuration @f$ \boldsymbol{Y} @f$ with respect to the undeformed
   * configuration @f$ \boldsymbol{X} @f$, and `delta deformation gradient` (@c dx_dY) the increment
   * @f$ \Delta\boldsymbol{F} = \partial\boldsymbol{x}/\partial\boldsymbol{Y} @f$ accumulated since. The material is
   * evaluated at
   * @f[
   *   \boldsymbol{F} = \Delta\boldsymbol{F}\,\boldsymbol{F}_n, \qquad \bar{N} = \bar{N}_n + \Delta\bar{N},
   * @f]
   * where @f$ \Delta F_{ij} = \delta_{ij} + \partial\Delta u_i/\partial Y_j @f$ and @f$ \Delta\bar{N} @f$ are
   * accumulated by incrementDeformation() from the increments interpolated by the cell or particle. Because the
   * unknowns of an increment are
   * @f$ \Delta\boldsymbol{F} @f$ and @f$ \Delta\bar{N} @f$, computeYourself() converts the stress and driving-force
   * tangents with respect to @f$ \boldsymbol{F} @f$ to tangents with respect to @f$ \Delta\boldsymbol{F} @f$ by the
   * chain rule
   * @f[
   *   \frac{\partial F_{iI}}{\partial \Delta F_{jJ}} = \delta_{ij}\,F_{n,JI},
   * @f]
   * and stores them in #tangents together with the response in #response.
   *
   * In plane strain (@c nDim = 2) all tensors are stored in 3D, the material is evaluated through
   * MarmotMaterialGradientEnhancedFiniteStrain::computePlaneStrain, and the in-plane parts are handed out.
   *
   * @tparam nDim Spatial dimension (2: plane strain, 3: 3D).
   */
  template < int nDim >
  class GradientEnhancedFiniteStrainMaterialPoint : public MarmotMaterialPoint {

  protected:
    constexpr static int _nVertices = 1;                     ///< a material point is represented by a single vertex

    using TensorD    = Fastor::Tensor< double, nDim >;       ///< vector of size nDim
    using TensorDD   = Fastor::Tensor< double, nDim, nDim >; ///< second-order tensor of size nDim
    using TensorDDDD = Fastor::Tensor< double, nDim, nDim, nDim, nDim >; ///< fourth-order tensor of size nDim

    int _mpNumber;   ///< label of the material point, also passed to the material as its number

    TensorD _x0;     ///< coordinates @f$ \boldsymbol{X} @f$ in the undeformed configuration

    double _vol0;    ///< volume @f$ V_0 @f$ in the undeformed configuration
    double _density; ///< mass density in the reference configuration, as reported by the material

    using Material = MarmotMaterialGradientEnhancedFiniteStrain;            ///< the consumed material interface

    std::unique_ptr< MarmotMaterialGradientEnhancedFiniteStrain > material; ///< the assigned material

    /**
     * @class Marmot::MaterialPoints::GradientEnhancedFiniteStrainMaterialPoint::MPStateVarManager
     * @brief State variable layout of the material point, followed by the state of the material.
     *
     * All kinematic quantities are stored with 3 (vectors) or 9 (tensors) entries, also in plane strain.
     */
    class MPStateVarManager : public MarmotStateVarVectorManager {

      /// names and lengths of the state variables; the material state starts at `begin of material state`
      inline const static auto layout = makeLayout( {
        { .name = "displacement", .length = 3 },
        { .name = "velocity", .length = 3 },
        { .name = "acceleration", .length = 3 },
        { .name = "delta displacement", .length = 3 },
        { .name = "delta deformation gradient", .length = 9 },
        { .name = "deformation gradient", .length = 9 },
        { .name = "nonlocal damage", .length = 1 },
        { .name = "local damage", .length = 1 },
        { .name = "stress", .length = 9 },
        { .name = "F0 XX", .length = 1 },
        { .name = "F0 YY", .length = 1 },
        { .name = "F0 ZZ", .length = 1 },
        { .name = "begin of material state", .length = 0 },
      } );

    public:
      FastorStandardTensors::TensorMap3d  u;     ///< displacement of the last accepted state
      FastorStandardTensors::TensorMap3d  v;     ///< velocity (Newmark-beta, set by the particle)
      FastorStandardTensors::TensorMap3d  a;     ///< acceleration (Newmark-beta, set by the particle)
      FastorStandardTensors::TensorMap3d  du;    ///< displacement increment of the current step
      FastorStandardTensors::TensorMap33d dx_dY; ///< deformation gradient increment @f$ \Delta\boldsymbol{F} @f$
      FastorStandardTensors::TensorMap33d dY_dX; ///< deformation gradient @f$ \boldsymbol{F}_n @f$, last accepted
      double& nonLocalDamage; ///< nonlocal field @f$ \bar{N} @f$ (total; incrementDeformation() adds the increment of
                              ///< the step to the committed value, see the state contract of MarmotMaterialPoint)
      double& localDamage;    ///< local driving force @f$ L @f$ of the last material evaluation
      FastorStandardTensors::TensorMap33d stress;        ///< Kirchhoff stress @f$ \boldsymbol{\tau} @f$ (3x3)
      double&                             F0_XX;         ///< eigen deformation (geostatic stress), XX component
      double&                             F0_YY;         ///< eigen deformation (geostatic stress), YY component
      double&                             F0_ZZ;         ///< eigen deformation (geostatic stress), ZZ component
      Eigen::Map< Eigen::VectorXd >       materialState; ///< state variables of the material

      /**
       * @brief Number of state variables of the material point itself (without the material).
       * @return Size of the fixed part of the layout.
       */
      static int getNumberOfRequiredStateVars() { return layout.nRequiredStateVars; };

      /**
       * @brief Map the layout onto a state variable vector.
       * @param[in,out] theStateVarVector State variable vector of the material point.
       * @param[in]     nStateVars        Total size of the vector (material point + material).
       */
      MPStateVarManager( double* theStateVarVector, int nStateVars )
        : MarmotStateVarVectorManager( theStateVarVector, layout ),
          u( &find( "displacement" ) ),
          v( &find( "velocity" ) ),
          a( &find( "acceleration" ) ),
          du( &find( "delta displacement" ) ),
          dx_dY( &find( "delta deformation gradient" ) ),
          dY_dX( &find( "deformation gradient" ) ),
          nonLocalDamage( find( "nonlocal damage" ) ),
          localDamage( find( "local damage" ) ),
          stress( &find( "stress" ) ),
          F0_XX( find( "F0 XX" ) ),
          F0_YY( find( "F0 YY" ) ),
          F0_ZZ( find( "F0 ZZ" ) ),
          materialState( &find( "begin of material state" ), nStateVars - getNumberOfRequiredStateVars() ){};
    };

    std::unique_ptr< MPStateVarManager > state; ///< view on the assigned state variable vector

  public:
    bool hasEigenDeformation; ///< true once a geostatic initial condition has set an eigen deformation

    /**
     * @brief Construct a material point.
     * @param[in] mpNumber           Label of the material point.
     * @param[in] vertexCoordinates  Undeformed coordinates @f$ \boldsymbol{X} @f$ (nDim values).
     * @param[in] nVertexCoordinates Number of coordinates (not checked).
     * @param[in] volume             Undeformed volume @f$ V_0 @f$.
     */
    GradientEnhancedFiniteStrainMaterialPoint( int           mpNumber,
                                               const double* vertexCoordinates,
                                               int           nVertexCoordinates,
                                               double        volume )
      : _mpNumber( mpNumber ), hasEigenDeformation( false )
    {
      assignVertexCoordinates( vertexCoordinates );
      assignVolume( volume );
    };

    /**
     * @brief Assign the state variable vector (material point layout followed by the material state).
     * @param[in,out] stateVars  State variable vector, owned by the host.
     * @param[in]     nStateVars Its size, getNumberOfRequiredStateVars().
     */
    void assignStateVars( double* stateVars, int nStateVars );

    /**
     * @brief Access a state variable of the material point or, if it has none of that name, of the material.
     * @param[in] stateName Name of the state variable (see MPStateVarManager).
     * @return View on the state variable.
     */
    StateView getStateView( const std::string& stateName ) const;

    /**
     * @brief Shape of the material point.
     * @return Always "point".
     */
    std::string getMaterialPointShape() const { return "point"; };

    /**
     * @brief Create the material from the MarmotMaterialGradientEnhancedFiniteStrainFactory.
     * @param[in] property Material section (name and properties).
     * @throws std::invalid_argument if the material does not implement MarmotMaterialGradientEnhancedFiniteStrain.
     */
    void assignMaterial( const MarmotMaterialSection& property );

    /**
     * @brief Initialize the state: @f$ \boldsymbol{F}_n = \boldsymbol{I} @f$, unit eigen deformation, a zero
     * increment, the material state, and the density.
     *
     * The density is queried from the material here already, so that the inertia can be assembled before the first
     * computation. Requires assignMaterial() and assignStateVars() to have been called.
     */
    void initializeYourself();

    /**
     * @brief Label of the material point.
     * @return The label passed to the constructor.
     */
    int getMaterialPointNumber() const { return _mpNumber; }

    /**
     * @brief Spatial dimension.
     * @return nDim.
     */
    int getDimension() const { return nDim; }

    /**
     * @brief Number of vertices.
     * @return 1.
     */
    int getNumberOfVertices() const { return _nVertices; };

    /**
     * @brief Number of state variables: material point layout plus material.
     * @return Required size of the state variable vector; requires an assigned material.
     */
    int getNumberOfRequiredStateVars() const
    {
      return MPStateVarManager::getNumberOfRequiredStateVars() + material->getNumberOfRequiredStateVars();
    };

    /**
     * @brief Set the undeformed volume.
     * @param[in] volume Undeformed volume @f$ V_0 @f$.
     */
    void assignVolume( double volume ) { _vol0 = volume; };

    /**
     * @brief Undeformed volume.
     * @return @f$ V_0 @f$.
     */
    double getVolumeUndeformed() const { return _vol0; }

    /**
     * @brief Set the undeformed coordinates.
     * @param[in] coordinates Undeformed coordinates @f$ \boldsymbol{X} @f$ (nDim values).
     */
    void assignVertexCoordinates( const double* coordinates ) { _x0 = TensorD( coordinates ); };

    /**
     * @brief Coordinates of the single vertex, identical to getCoordinatesAtCenter().
     * @param[out] coordinates Coordinates @f$ \boldsymbol{X} + \boldsymbol{u} @f$ (nDim values).
     */
    void getVertexCoordinates( double* coordinates ) const { return getCoordinatesAtCenter( coordinates ); };

    /**
     * @brief Coordinates of the last accepted state, @f$ \boldsymbol{X} + \boldsymbol{u} @f$.
     *
     * The displacement increment of the current step is not included, so this is the position in the intermediate
     * reference configuration @f$ \boldsymbol{Y} @f$.
     *
     * @param[out] coordinates Coordinates (nDim values).
     */
    void getCoordinatesAtCenter( double* coordinates ) const
    {
      Eigen::Map< const Eigen::Matrix< double, nDim, 1 > > x0( _x0.data() );
      Eigen::Map< Eigen::Matrix< double, nDim, 1 > >       newCoords( coordinates );
      Eigen::Map< Eigen::Matrix< double, nDim, 1 > >       u( state->u.data() );

      newCoords = x0 + u;
    };

    /**
     * @brief Displacement of the last accepted state.
     * @param[out] displacement Displacement @f$ \boldsymbol{u} @f$ (nDim values).
     */
    void getCenterDisplacement( double* displacement ) const
    {
      for ( int i = 0; i < nDim; i++ )
        displacement[i] = state->u( i );
    };

    /**
     * @brief Mass density in the reference configuration.
     * @return Density reported by the material at initializeYourself() and after each computeYourself().
     */
    double getDensityUndeformed() const { return _density; };

    /**
     * @brief Undeformed coordinates.
     * @return @f$ \boldsymbol{X} @f$.
     */
    const TensorD& coordinates() const { return _x0; };

    /**
     * @brief Reset the increment of the current step: @f$ \Delta\boldsymbol{u} = \boldsymbol{0} @f$,
     * @f$ \Delta\boldsymbol{F} = \boldsymbol{I} @f$.
     *
     * The nonlocal field is not reset; incrementDeformation() adds to the stored total @f$ \bar{N} @f$.
     *
     * @param[in] timeNew Time at the end of the increment (unused).
     * @param[in] dT      Time increment (unused).
     */
    virtual void prepareYourself( double timeNew, double dT );

    /**
     * @brief Evaluate the material at @f$ \boldsymbol{F} = \Delta\boldsymbol{F}\,\boldsymbol{F}_n @f$ and
     * @f$ \bar{N} @f$ and fill #response and #tangents.
     * @param[in] timeNew Time at the end of the increment.
     * @param[in] dT      Time increment.
     */
    virtual void computeYourself( double timeNew, double dT ) = 0;

    /**
     * @brief Accept the current increment: @f$ \boldsymbol{u} \leftarrow \boldsymbol{u} + \Delta\boldsymbol{u} @f$,
     * @f$ \boldsymbol{F}_n \leftarrow \Delta\boldsymbol{F}\,\boldsymbol{F}_n @f$.
     *
     * The increment itself is not reset here; the next prepareYourself() does that.
     */
    virtual void acceptStateAndPosition()
    {
      const auto&                           u_n  = state->u;
      const FastorStandardTensors::Tensor3d u_np = u_n + state->du;

      const FastorStandardTensors::Tensor33d dx_dX_n  = state->dY_dX;
      const FastorStandardTensors::Tensor33d dx_dX_np = state->dx_dY % dx_dX_n;

      mapEigenToFastor( state->u )     = mapEigenToFastor( u_np );
      mapEigenToFastor( state->dY_dX ) = mapEigenToFastor( dx_dX_np );
    };

    /**
     * @brief Add an interpolated increment to the state.
     *
     * @f$ \Delta\boldsymbol{u} \mathrel{+}= \delta\boldsymbol{u} @f$,
     * @f$ \Delta\boldsymbol{F} \mathrel{+}= \partial\delta\boldsymbol{u}/\partial\boldsymbol{Y} @f$ and
     * @f$ \bar{N} \mathrel{+}= \delta\bar{N} @f$; in plane strain the inputs are expanded to 3D.
     *
     * @param[in] displacementIncrement         Displacement increment @f$ \delta\boldsymbol{u} @f$.
     * @param[in] displacementGradientIncrement Its gradient with respect to the intermediate reference
     *                                          configuration, @f$ \partial\delta u_i/\partial Y_j @f$.
     * @param[in] nonLocalDamage                Increment @f$ \delta\bar{N} @f$ of the nonlocal field.
     */
    virtual void incrementDeformation( const TensorD&  displacementIncrement,
                                       const TensorDD& displacementGradientIncrement,
                                       double          nonLocalDamage ) = 0;

    /// @brief Response of the last computeYourself(), in nDim components.
    struct {
      Fastor::Tensor< double, nDim, nDim > S;  ///< Kirchhoff stress @f$ \boldsymbol{\tau} @f$
      double                               dL; ///< @f$ \Delta L @f$: the local driving force of this evaluation
                                               ///< minus the value stored in the state variable `local damage`
                                               ///< (which is then overwritten by the new value)
      double nonLocalRadius;                   ///< nonlocal radius @f$ R @f$ of the material, @f$ c = R^2 @f$
    } response;

    /// @brief Algorithmic tangents of the last computeYourself(), with respect to @f$ \Delta\boldsymbol{F} @f$ and
    /// @f$ \bar{N} @f$, in nDim components.
    struct {
      Fastor::Tensor< double, nDim, nDim, nDim, nDim > dS_dDeltaF; ///< @f$ \partial\tau_{ij}/\partial\Delta F_{kL} @f$
      Fastor::Tensor< double, nDim, nDim >             dS_dN;      ///< @f$ \partial\tau_{ij}/\partial\bar{N} @f$
      Fastor::Tensor< double, nDim, nDim >             dL_dDeltaF; ///< @f$ \partial L/\partial\Delta F_{kL} @f$
      double                                           dL_dN;      ///< @f$ \partial L/\partial\bar{N} @f$
    } tangents;

    /**
     * @brief Deformation gradient increment.
     * @return @f$ \Delta\boldsymbol{F} = \partial\boldsymbol{x}/\partial\boldsymbol{Y} @f$ (nDim x nDim).
     */
    TensorDD dx_dY() const { return state->dx_dY( Fastor::seq( 0, nDim ), Fastor::seq( 0, nDim ) ); };

    /**
     * @brief Deformation gradient of the last accepted state.
     * @return @f$ \boldsymbol{F}_n = \partial\boldsymbol{Y}/\partial\boldsymbol{X} @f$ (nDim x nDim).
     */
    TensorDD dY_dX() const { return state->dY_dX( Fastor::seq( 0, nDim ), Fastor::seq( 0, nDim ) ); };

    /**
     * @brief Velocity.
     * @return Stored velocity (nDim components).
     */
    TensorD getVelocity() const { return state->v( Fastor::seq( 0, nDim ) ); };

    /**
     * @brief Acceleration.
     * @return Stored acceleration (nDim components).
     */
    TensorD getAcceleration() const { return state->a( Fastor::seq( 0, nDim ) ); };

    /**
     * @brief Store the velocity.
     * @param[in] velocity Velocity (nDim components).
     */
    void setVelocity( const TensorD& velocity )
    {
      for ( int i = 0; i < nDim; i++ )
        state->v( i ) = velocity( i );
    };

    /**
     * @brief Store the acceleration.
     * @param[in] acceleration Acceleration (nDim components).
     */
    void setAcceleration( const TensorD& acceleration )
    {
      for ( int i = 0; i < nDim; i++ )
        state->a( i ) = acceleration( i );
    };

    /**
     * @brief Apply an initial condition.
     *
     * Only `geostaticstress` is supported: the eigen deformation that produces the hydrostatic stress
     * @f$ \tau_{XX} = \tau_{YY} = \tau_{ZZ} = \mathrm{value}[0] @f$ is found by
     * MarmotMaterialGradientEnhancedFiniteStrain::findEigenDeformationForEigenStress and applied in every subsequent
     * material evaluation.
     *
     * @param[in] conditionName Name of the condition.
     * @param[in] value         Its value(s); `value[0]` is the geostatic normal stress.
     * @throws std::invalid_argument for any other condition.
     */
    virtual void setInitialCondition( const std::string& conditionName, const double* value ) override
    {
      if ( conditionName == "geostaticstress" ) {
        std::tuple< double, double, double > geostaticNormalStressComponents = { value[0], value[0], value[0] };
        const auto [F0_XX,
                    F0_YY,
                    F0_ZZ] = material->findEigenDeformationForEigenStress( { state->F0_XX, state->F0_YY, state->F0_ZZ },
                                                                           geostaticNormalStressComponents,
                                                                           state->materialState.data() );

        state->F0_XX = F0_XX;
        state->F0_YY = F0_YY;
        state->F0_ZZ = F0_ZZ;

        hasEigenDeformation = true;
      }
      else {
        throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid initial condition" );
      }
    }
  };

  template < int nDim >
  void GradientEnhancedFiniteStrainMaterialPoint< nDim >::assignStateVars( double* stateVars, int nStateVars )
  {
    state = std::make_unique< MPStateVarManager >( stateVars, nStateVars );
  }

  template < int nDim >
  StateView GradientEnhancedFiniteStrainMaterialPoint< nDim >::getStateView( const std::string& stateName ) const
  {
    if ( state->contains( stateName ) )
      return state->getStateView( stateName );
    else
      return material->getStateView( stateName, state->materialState.data() );
  }

  template < int nDim >
  void GradientEnhancedFiniteStrainMaterialPoint< nDim >::assignMaterial( const MarmotMaterialSection& section )
  {
    material = std::unique_ptr< MarmotMaterialGradientEnhancedFiniteStrain >(
      dynamic_cast< MarmotMaterialGradientEnhancedFiniteStrain* >(
        MarmotLibrary::MarmotMaterialGradientEnhancedFiniteStrainFactory::createMaterial( section.materialName,
                                                                                          section.materialProperties,
                                                                                          section.nMaterialProperties,
                                                                                          _mpNumber ) ) );

    if ( !material )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__
                                                << ": invalid material assigned; cannot cast to "
                                                   "MarmotMaterialGradientEnhancedFiniteStrain!" );
  }

  template < int nDim >
  void GradientEnhancedFiniteStrainMaterialPoint< nDim >::initializeYourself()
  {
    state->dY_dX.eye();
    state->F0_XX = 1.0;
    state->F0_YY = 1.0;
    state->F0_ZZ = 1.0;
    this->prepareYourself( 0, 0 );
    material->initializeYourself( state->materialState.data(), state->materialState.size() );
    // known from the start, so that inertia can be assembled before the first computation
    _density = material->getDensity( state->materialState.data() );
  }

  template < int nDim >
  void GradientEnhancedFiniteStrainMaterialPoint< nDim >::prepareYourself( double timeNew, double dT )
  {
    state->du.zeros();
    state->dx_dY.eye();
  }

  /**
   * @class Marmot::MaterialPoints::GradientEnhancedFiniteStrainMaterialPoint2D
   * @brief Plane-strain gradient-enhanced finite-strain material point (registered as
   *        `GradientEnhancedFiniteStrain/PlaneStrain`).
   *
   * The material is evaluated through MarmotMaterialGradientEnhancedFiniteStrain::computePlaneStrain with the 3D
   * deformation gradient whose out-of-plane entries are those of the identity; the response and tangents are reduced
   * to their in-plane components.
   */
  class GradientEnhancedFiniteStrainMaterialPoint2D : public GradientEnhancedFiniteStrainMaterialPoint< 2 > {

  public:
    using GradientEnhancedFiniteStrainMaterialPoint::GradientEnhancedFiniteStrainMaterialPoint;

    /**
     * @brief Evaluate the material in plane strain and fill #response and #tangents.
     * @param[in] timeNew Time at the end of the increment.
     * @param[in] dT      Time increment.
     */
    void computeYourself( double timeNew, double dT );

    /**
     * @brief Add an interpolated in-plane increment (expanded to 3D) to the state.
     * @param[in] displacementIncrement         Displacement increment.
     * @param[in] displacementGradientIncrement Its gradient with respect to @f$ \boldsymbol{Y} @f$.
     * @param[in] nonLocalDamage                Increment of the nonlocal field.
     */
    void incrementDeformation( const TensorD&  displacementIncrement,
                               const TensorDD& displacementGradientIncrement,
                               double          nonLocalDamage );
  };

  /**
   * @class Marmot::MaterialPoints::GradientEnhancedFiniteStrainMaterialPoint3D
   * @brief Three-dimensional gradient-enhanced finite-strain material point (registered as
   *        `GradientEnhancedFiniteStrain/3D`).
   */
  class GradientEnhancedFiniteStrainMaterialPoint3D : public GradientEnhancedFiniteStrainMaterialPoint< 3 > {

  public:
    using GradientEnhancedFiniteStrainMaterialPoint::GradientEnhancedFiniteStrainMaterialPoint;

    /**
     * @brief Evaluate the material in 3D and fill #response and #tangents.
     * @param[in] timeNew Time at the end of the increment.
     * @param[in] dT      Time increment.
     */
    void computeYourself( double timeNew, double dT );

    /**
     * @brief Add an interpolated increment to the state.
     * @param[in] displacementIncrement         Displacement increment.
     * @param[in] displacementGradientIncrement Its gradient with respect to @f$ \boldsymbol{Y} @f$.
     * @param[in] nonLocalDamage                Increment of the nonlocal field.
     */
    void incrementDeformation( const TensorD&  displacementIncrement,
                               const TensorDD& displacementGradientIncrement,
                               double          nonLocalDamage );
  };

} // namespace Marmot::MaterialPoints
