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
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
#include "Marmot/MarmotMaterialPoint.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotTensor.h"

#include <memory>
#include <stdexcept>
#include <vector>

namespace Marmot::MaterialPoints {

  /**
   * @class Marmot::MaterialPoints::DisplacementMaterialPoint
   * @brief MPM material point of the finite-strain displacement formulation.
   *
   * @details The material point carries the state of a MarmotMaterialFiniteStrain and the kinematics of a total
   * deformation gradient that is updated multiplicatively, increment by increment. Three configurations are
   * distinguished: the undeformed configuration @f$ \boldsymbol{X} @f$, the intermediate configuration
   * @f$ \boldsymbol{Y} @f$ (the last accepted state, i.e., the configuration at the beginning of the current
   * increment) and the current configuration @f$ \boldsymbol{x} @f$. The deformation gradient is
   * @f[
   *   \boldsymbol{F} = \frac{\partial \boldsymbol{x}}{\partial \boldsymbol{X}}
   *     = \Delta\boldsymbol{F}\,\boldsymbol{F}_n, \qquad
   *   \Delta\boldsymbol{F} = \frac{\partial \boldsymbol{x}}{\partial \boldsymbol{Y}}
   *     = \boldsymbol{I} + \frac{\partial \Delta\boldsymbol{u}}{\partial \boldsymbol{Y}}, \qquad
   *   \boldsymbol{F}_n = \frac{\partial \boldsymbol{Y}}{\partial \boldsymbol{X}},
   * @f]
   * where the host (a cell or a particle) supplies the displacement increment @f$ \Delta\boldsymbol{u} @f$ and its
   * gradient with respect to @f$ \boldsymbol{Y} @f$ through incrementDeformation(). computeYourself() evaluates the
   * material with @f$ \boldsymbol{F} @f$ and provides the Kirchhoff stress @f$ \boldsymbol{\tau} @f$ (response.S)
   * and its derivative with respect to the increment @f$ \Delta\boldsymbol{F} @f$ (tangents.dS_dDeltaF),
   * @f[
   *   \frac{\partial \tau_{ij}}{\partial \Delta F_{kL}}
   *     = \frac{\partial \tau_{ij}}{\partial F_{kN}}\,F_{n,LN}.
   * @f]
   * acceptStateAndPosition() then sets @f$ \boldsymbol{u} \leftarrow \boldsymbol{u} + \Delta\boldsymbol{u} @f$
   * and @f$ \boldsymbol{F}_n \leftarrow \Delta\boldsymbol{F}\,\boldsymbol{F}_n @f$.
   *
   * Internally, all kinematic quantities are stored in 3D; in 2D (plane strain) the out-of-plane components of
   * @f$ \Delta\boldsymbol{F} @f$ and @f$ \boldsymbol{F}_n @f$ remain those of the identity. The volume and the
   * density are those of the undeformed configuration.
   *
   * @tparam nDim Spatial dimension (2: plane strain, 3: 3D).
   */
  template < int nDim >
  class DisplacementMaterialPoint : public MarmotMaterialPoint {

  protected:
    constexpr static int _nVertices = 1; ///< a material point has a single vertex, its center

    /// number of rotational dofs in nDim dimensions (not used by this class)
    static constexpr int nRot = Marmot::ContinuumMechanics::CommonTensors::getNumberOfDofForRotation( nDim );

    using TensorD    = Fastor::Tensor< double, nDim >;                   ///< vector of size nDim
    using TensorDD   = Fastor::Tensor< double, nDim, nDim >;             ///< second-order tensor of size nDim
    using TensorDDDD = Fastor::Tensor< double, nDim, nDim, nDim, nDim >; ///< fourth-order tensor of size nDim

    int _mpNumber;   ///< label of the material point, also passed to the material

    TensorD _x0;     ///< coordinates in the undeformed configuration @f$ \boldsymbol{X} @f$

    double _vol0;    ///< volume in the undeformed configuration @f$ V_0 @f$
    double _density; ///< mass density in the undeformed configuration, as provided by the material

    using Material = MarmotMaterialFiniteStrain; ///< material interface consumed by the material point

    std::unique_ptr< Material > material;        ///< the finite-strain material

    /**
     * @class Marmot::MaterialPoints::DisplacementMaterialPoint::MPStateVarManager
     * @brief State of the material point in the state vector provided by the host.
     *
     * @details Layout (all quantities 3D, regardless of nDim): total displacement @f$ \boldsymbol{u} @f$
     * ("displacement"), velocity, acceleration, the displacement increment of the current increment
     * @f$ \Delta\boldsymbol{u} @f$ ("delta displacement"), the incremental deformation gradient
     * @f$ \Delta\boldsymbol{F} @f$ ("delta deformation gradient"), the deformation gradient of the last accepted
     * state @f$ \boldsymbol{F}_n @f$ ("deformation gradient"), the Kirchhoff stress ("stress"), followed by the
     * state of the material ("begin of material state").
     */
    class MPStateVarManager : public MarmotStateVarVectorManager {

      /// the layout of the state vector (names and lengths)
      inline const static auto layout = makeLayout( {
        { .name = "displacement", .length = 3 },
        { .name = "velocity", .length = 3 },
        { .name = "acceleration", .length = 3 },
        { .name = "delta displacement", .length = 3 },
        { .name = "delta deformation gradient", .length = 9 },
        { .name = "deformation gradient", .length = 9 },
        { .name = "stress", .length = 9 },
        { .name = "begin of material state", .length = 0 },
      } );

    public:
      FastorStandardTensors::TensorMap3d  u;     ///< total displacement of the last accepted state
      FastorStandardTensors::TensorMap3d  v;     ///< velocity
      FastorStandardTensors::TensorMap3d  a;     ///< acceleration
      FastorStandardTensors::TensorMap3d  du;    ///< displacement increment @f$ \Delta\boldsymbol{u} @f$
      FastorStandardTensors::TensorMap33d dx_dY; ///< incremental deformation gradient @f$ \Delta\boldsymbol{F} @f$
      FastorStandardTensors::TensorMap33d
        dY_dX; ///< deformation gradient of the last accepted state @f$ \boldsymbol{F}_n @f$
      FastorStandardTensors::TensorMap33d stress;        ///< Kirchhoff stress (see computeYourself())
      Eigen::Map< Eigen::VectorXd >       materialState; ///< state variables of the material

      /**
       * @brief Number of state variables required by the material point itself, without the material state.
       * @return The length of the layout.
       */
      static int getNumberOfRequiredStateVars() { return layout.nRequiredStateVars; };

      /**
       * @brief Maps the state layout onto a state vector.
       * @param[in] theStateVarVector State vector of the material point (owned by the host).
       * @param[in] nStateVars Total length of the state vector, including the material state.
       */
      MPStateVarManager( double* theStateVarVector, int nStateVars )
        : MarmotStateVarVectorManager( theStateVarVector, layout ),
          u( &find( "displacement" ) ),
          v( &find( "velocity" ) ),
          a( &find( "acceleration" ) ),
          du( &find( "delta displacement" ) ),
          dx_dY( &find( "delta deformation gradient" ) ),
          dY_dX( &find( "deformation gradient" ) ),
          stress( &find( "stress" ) ),
          materialState( &find( "begin of material state" ), nStateVars - getNumberOfRequiredStateVars() ){};
    };

    std::unique_ptr< MPStateVarManager > state; ///< view on the state vector, set by assignStateVars()

  public:
    /**
     * @brief Constructs a material point.
     * @param[in] mpNumber Label of the material point.
     * @param[in] vertexCoordinates Coordinates of the material point in the undeformed configuration (nDim values).
     * @param[in] nVertexCoordinates Number of coordinates (not used).
     * @param[in] volume Volume in the undeformed configuration.
     */
    DisplacementMaterialPoint( int mpNumber, const double* vertexCoordinates, int nVertexCoordinates, double volume )
      : _mpNumber( mpNumber )
    {

      assignVertexCoordinates( vertexCoordinates );
      assignVolume( volume );
    };

    /**
     * @brief Assigns the state vector (see MPStateVarManager for its layout).
     * @param[in] stateVars State vector, owned by the host.
     * @param[in] nStateVars Length of the state vector, see getNumberOfRequiredStateVars().
     */
    void assignStateVars( double* stateVars, int nStateVars );

    /**
     * @brief Returns a view on a state of the material point or, if not found, of the material.
     * @param[in] stateName Name of the state, e.g. "deformation gradient", "stress" or a material state.
     * @return The view on the state.
     */
    StateView getStateView( const std::string& stateName ) const;

    /**
     * @brief Shape of the material point.
     * @return Always "point".
     */
    std::string getMaterialPointShape() const { return "point"; };

    /**
     * @brief Creates the material through the MarmotMaterialFiniteStrainFactory.
     * @param[in] property Section with the material name and properties.
     * @throws std::invalid_argument if no finite-strain material of that name exists.
     */
    void assignMaterial( const MarmotMaterialSection& property );

    /**
     * @brief Initializes the state: @f$ \boldsymbol{F}_n = \boldsymbol{I} @f$, a zero increment, the material state
     * and the density (so that the inertia can be assembled before the first computation).
     */
    void initializeYourself();

    /**
     * @brief Label of the material point.
     * @return The label.
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
     * @brief Number of state variables required by the material point, including those of the material.
     * @return The number of state variables.
     */
    int getNumberOfRequiredStateVars() const
    {
      return MPStateVarManager::getNumberOfRequiredStateVars() + material->getNumberOfRequiredStateVars();
    };

    /**
     * @brief Sets the volume in the undeformed configuration.
     * @param[in] volume The volume @f$ V_0 @f$.
     */
    void assignVolume( double volume ) { _vol0 = volume; };

    /**
     * @brief Volume in the undeformed configuration.
     * @return @f$ V_0 @f$.
     */
    double getVolumeUndeformed() const { return _vol0; }

    /**
     * @brief Sets the coordinates in the undeformed configuration.
     * @param[in] coordinates The coordinates @f$ \boldsymbol{X} @f$ (nDim values).
     */
    void assignVertexCoordinates( const double* coordinates ) { _x0 = TensorD( coordinates ); };

    /**
     * @brief Coordinates of the single vertex, identical to getCoordinatesAtCenter().
     * @param[out] coordinates The coordinates (nDim values).
     */
    void getVertexCoordinates( double* coordinates ) const { return getCoordinatesAtCenter( coordinates ); };

    /**
     * @brief Coordinates in the intermediate configuration, @f$ \boldsymbol{Y} = \boldsymbol{X} + \boldsymbol{u} @f$,
     * with the displacement of the last accepted state (the increment of the current step is not included).
     * @param[out] coordinates The coordinates (nDim values).
     */
    void getCoordinatesAtCenter( double* coordinates ) const
    {
      Eigen::Map< const Eigen::Matrix< double, nDim, 1 > > x0( _x0.data() );
      Eigen::Map< Eigen::Matrix< double, nDim, 1 > >       newCoords( coordinates );
      Eigen::Map< Eigen::Matrix< double, nDim, 1 > >       u( state->u.data() );

      newCoords = x0 + u;
    };

    /**
     * @brief Total displacement of the last accepted state.
     * @param[out] displacement The displacement @f$ \boldsymbol{u} @f$ (nDim values).
     */
    void getCenterDisplacement( double* displacement ) const
    {
      for ( int i = 0; i < nDim; i++ )
        displacement[i] = state->u( i );
    };

    /**
     * @brief Mass density in the undeformed configuration, as provided by the material.
     * @return The density.
     */
    double getDensityUndeformed() const { return _density; };

    /**
     * @brief Coordinates in the undeformed configuration.
     * @return @f$ \boldsymbol{X} @f$.
     */
    const TensorD& getCoordinatesUndeformed() const { return _x0; };

    /**
     * @brief Starts a new evaluation of the increment: resets @f$ \Delta\boldsymbol{u} = \boldsymbol{0} @f$ and
     * @f$ \Delta\boldsymbol{F} = \boldsymbol{I} @f$, so that incrementDeformation() can accumulate the full
     * increment again.
     * @param[in] timeNew Time at the end of the increment (not used).
     * @param[in] dT Time increment (not used).
     */
    virtual void prepareYourself( double timeNew, double dT );

    /**
     * @brief Evaluates the material with @f$ \boldsymbol{F} = \Delta\boldsymbol{F}\,\boldsymbol{F}_n @f$ and fills
     * response and tangents.
     * @param[in] timeNew Time at the end of the increment.
     * @param[in] dT Time increment.
     */
    virtual void computeYourself( double timeNew, double dT ) = 0;

    /**
     * @brief Accepts the increment: @f$ \boldsymbol{u} \leftarrow \boldsymbol{u} + \Delta\boldsymbol{u} @f$ and
     * @f$ \boldsymbol{F}_n \leftarrow \Delta\boldsymbol{F}\,\boldsymbol{F}_n @f$. The increment itself is reset
     * by the next prepareYourself().
     */
    virtual void acceptStateAndPosition()
    {

      const auto&                           u_n  = state->u;
      const FastorStandardTensors::Tensor3d u_np = u_n + state->du;

      // TODO: make auto&
      const FastorStandardTensors::Tensor33d dx_dX_n  = state->dY_dX;
      const FastorStandardTensors::Tensor33d dx_dX_np = state->dx_dY % dx_dX_n;

      mapEigenToFastor( state->u )     = mapEigenToFastor( u_np );
      mapEigenToFastor( state->dY_dX ) = mapEigenToFastor( dx_dX_np );
    };

    /**
     * @brief Adds a contribution to the displacement increment and to its gradient,
     * @f$ \Delta\boldsymbol{u} \mathrel{+}= \delta\boldsymbol{u} @f$ and
     * @f$ \Delta\boldsymbol{F} \mathrel{+}= \partial\,\delta\boldsymbol{u}/\partial\boldsymbol{Y} @f$.
     * @param[in] displacementIncrement Contribution @f$ \delta\boldsymbol{u} @f$ to the displacement increment.
     * @param[in] displacementGradientIncrement Its gradient with respect to the intermediate configuration
     * @f$ \boldsymbol{Y} @f$.
     */
    virtual void incrementDeformation( const TensorD&  displacementIncrement,
                                       const TensorDD& displacementGradientIncrement ) = 0;

    /// response of the last computeYourself()
    struct {
      Fastor::Tensor< double, nDim, nDim > S; ///< Kirchhoff stress @f$ \boldsymbol{\tau} @f$ (in-plane part in 2D)
    } response;

    /// algorithmic tangents of the last computeYourself()
    struct {
      /// @f$ \partial \tau_{ij} / \partial \Delta F_{kL} @f$ (in-plane part in 2D)
      Fastor::Tensor< double, nDim, nDim, nDim, nDim > dS_dDeltaF;
    } tangents;

    /**
     * @brief Incremental deformation gradient.
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
     * @return The velocity (nDim values).
     */
    TensorD getVelocity() const { return state->v( Fastor::seq( 0, nDim ) ); };

    /**
     * @brief Acceleration.
     * @return The acceleration (nDim values).
     */
    TensorD getAcceleration() const { return state->a( Fastor::seq( 0, nDim ) ); };

    /**
     * @brief Sets the velocity (typically after a time integration in the host).
     * @param[in] velocity The velocity (nDim values).
     */
    void setVelocity( const TensorD& velocity )
    {
      for ( int i = 0; i < nDim; i++ )
        state->v( i ) = velocity( i );
    };

    /**
     * @brief Sets the acceleration (typically after a time integration in the host).
     * @param[in] acceleration The acceleration (nDim values).
     */
    void setAcceleration( const TensorD& acceleration )
    {
      for ( int i = 0; i < nDim; i++ )
        state->a( i ) = acceleration( i );
    };

    /**
     * @brief Initial conditions are not supported (in particular not `geostaticstress`, which would need an eigen
     *        deformation as in GradientEnhancedFiniteStrainMaterialPoint).
     * @param[in] conditionName Name of the initial condition.
     * @param[in] value Values of the initial condition (unused).
     * @throws std::invalid_argument always, so that an input with an initial condition does not silently start from
     *         the unloaded state.
     */
    virtual void setInitialCondition( const std::string& conditionName, const double* value ) override
    {
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": initial condition '" << conditionName
                                                << "' is not supported" );
    };
  };

  template < int nDim >
  void DisplacementMaterialPoint< nDim >::assignStateVars( double* stateVars, int nStateVars )
  {
    state = std::make_unique< MPStateVarManager >( stateVars, nStateVars );
  }

  template < int nDim >
  StateView DisplacementMaterialPoint< nDim >::getStateView( const std::string& stateName ) const
  {

    if ( state->contains( stateName ) )
      return state->getStateView( stateName );
    else
      return material->getStateView( stateName, state->materialState.data() );
  }

  template < int nDim >
  void DisplacementMaterialPoint< nDim >::assignMaterial( const MarmotMaterialSection& section )
  {
    material = std::unique_ptr< Material >(
      MarmotLibrary::MarmotMaterialFiniteStrainFactory::createMaterial( section.materialName,
                                                                        section.materialProperties,
                                                                        section.nMaterialProperties,
                                                                        _mpNumber ) );

    if ( !material )
      throw std::invalid_argument( MakeString()
                                   << __PRETTY_FUNCTION__ << ": invalid finite strain material assigned!" );
  }

  template < int nDim >
  void DisplacementMaterialPoint< nDim >::initializeYourself()
  {
    state->dY_dX.eye();
    /* state->dx_dY.eye(); */
    this->prepareYourself( 0, 0 );
    material->initializeYourself( state->materialState.data(), state->materialState.size() );
    // known from the start, so that inertia can be assembled before the first computation
    _density = material->getDensity( state->materialState.data() );
  }

  template < int nDim >
  void DisplacementMaterialPoint< nDim >::prepareYourself( double timeNew, double dT )
  {
    state->du.zeros();
    state->dx_dY.eye();
  }

  /**
   * @class Marmot::MaterialPoints::DisplacementMaterialPoint2D
   * @brief Plane-strain displacement material point (registered as "Displacement/PlaneStrain").
   *
   * @details The in-plane increment is expanded to 3D (out-of-plane components of the identity), the material is
   * evaluated by MarmotMaterialFiniteStrain::computePlaneStrain() with the 3D deformation gradient, and the in-plane
   * parts of @f$ \boldsymbol{\tau} @f$ and @f$ \partial\boldsymbol{\tau}/\partial\Delta\boldsymbol{F} @f$ are
   * kept. The full 3D Kirchhoff stress is also written to the state "stress".
   */
  class DisplacementMaterialPoint2D : public DisplacementMaterialPoint< 2 > {

  public:
    using DisplacementMaterialPoint::DisplacementMaterialPoint;

    /**
     * @brief Evaluates the material in plane strain, see DisplacementMaterialPoint::computeYourself().
     * @param[in] timeNew Time at the end of the increment.
     * @param[in] dT Time increment.
     */
    void computeYourself( double timeNew, double dT );

    /**
     * @brief Adds an in-plane contribution to the increment, see DisplacementMaterialPoint::incrementDeformation().
     * @param[in] displacementIncrement In-plane contribution to the displacement increment.
     * @param[in] displacementGradientIncrement Its in-plane gradient with respect to @f$ \boldsymbol{Y} @f$.
     */
    void incrementDeformation( const TensorD& displacementIncrement, const TensorDD& displacementGradientIncrement );
  };

  /**
   * @class Marmot::MaterialPoints::DisplacementMaterialPoint3D
   * @brief 3D displacement material point (registered as "Displacement/3D").
   *
   * @details The material is evaluated with the 3D deformation gradient by
   * MarmotMaterialFiniteStrain::computeStress(); the Kirchhoff stress is written to the state "stress".
   */
  class DisplacementMaterialPoint3D : public DisplacementMaterialPoint< 3 > {

  public:
    using DisplacementMaterialPoint::DisplacementMaterialPoint;

    /**
     * @brief Evaluates the material in 3D, see DisplacementMaterialPoint::computeYourself().
     * @param[in] timeNew Time at the end of the increment.
     * @param[in] dT Time increment.
     */
    void computeYourself( double timeNew, double dT );

    /**
     * @brief Adds a contribution to the increment, see DisplacementMaterialPoint::incrementDeformation().
     * @param[in] displacementIncrement Contribution to the displacement increment.
     * @param[in] displacementGradientIncrement Its gradient with respect to @f$ \boldsymbol{Y} @f$.
     */
    void incrementDeformation( const TensorD& displacementIncrement, const TensorDD& displacementGradientIncrement );
  };

} // namespace Marmot::MaterialPoints
