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
#include "Marmot/Marmot.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
#include "Marmot/MarmotMaterialPoint.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotTensor.h"

#include <memory>
#include <vector>

namespace Marmot::MaterialPoints {

  template < int nDim >
  class DisplacementMaterialPoint : public MarmotMaterialPoint {

  protected:
    constexpr static int _nVertices = 1;

    static constexpr int nRot = Marmot::ContinuumMechanics::CommonTensors::getNumberOfDofForRotation( nDim );

    using TensorD    = Fastor::Tensor< double, nDim >;
    using TensorDD   = Fastor::Tensor< double, nDim, nDim >;
    using TensorDDDD = Fastor::Tensor< double, nDim, nDim, nDim, nDim >;

    int _mpNumber;

    TensorD _x0;

    double _vol0;
    double _density;

    using Material = MarmotMaterialFiniteStrain;

    std::unique_ptr< Material > material;

    class MPStateVarManager : public MarmotStateVarVectorManager {

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
      FastorStandardTensors::TensorMap3d  u;
      FastorStandardTensors::TensorMap3d  v;
      FastorStandardTensors::TensorMap3d  a;
      FastorStandardTensors::TensorMap3d  du;
      FastorStandardTensors::TensorMap33d dx_dY;
      FastorStandardTensors::TensorMap33d dY_dX;
      FastorStandardTensors::TensorMap33d stress;
      Eigen::Map< Eigen::VectorXd >       materialState;

      static int getNumberOfRequiredStateVars() { return layout.nRequiredStateVars; };

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

    std::unique_ptr< MPStateVarManager > state;

  public:
    DisplacementMaterialPoint( int mpNumber, const double* vertexCoordinates, int nVertexCoordinates, double volume )
      : _mpNumber( mpNumber )
    {

      assignVertexCoordinates( vertexCoordinates );
      assignVolume( volume );
    };

    void assignStateVars( double* stateVars, int nStateVars );

    StateView getStateView( const std::string& stateName ) const;

    std::string getMaterialPointShape() const { return "point"; };

    void assignMaterial( const MarmotMaterialSection& property );

    void initializeYourself();

    int getMaterialPointNumber() const { return _mpNumber; }

    int getDimension() const { return nDim; }

    int getNumberOfVertices() const { return _nVertices; };

    int getNumberOfRequiredStateVars() const
    {
      return MPStateVarManager::getNumberOfRequiredStateVars() + material->getNumberOfRequiredStateVars();
    };

    void assignVolume( double volume ) { _vol0 = volume; };

    double getVolumeUndeformed() const { return _vol0; }

    void assignVertexCoordinates( const double* coordinates ) { _x0 = TensorD( coordinates ); };

    void getVertexCoordinates( double* coordinates ) const { return getCoordinatesAtCenter( coordinates ); };

    void getCoordinatesAtCenter( double* coordinates ) const
    {
      Eigen::Map< const Eigen::Matrix< double, nDim, 1 > > x0( _x0.data() );
      Eigen::Map< Eigen::Matrix< double, nDim, 1 > >       newCoords( coordinates );
      Eigen::Map< Eigen::Matrix< double, nDim, 1 > >       u( state->u.data() );

      newCoords = x0 + u;
    };

    void getCenterDisplacement( double* displacement ) const
    {
      for ( int i = 0; i < nDim; i++ )
        displacement[i] = state->u( i );
    };

    double getDensityUndeformed() const { return _density; };

    const TensorD& getCoordinatesUndeformed() const { return _x0; };

    virtual void prepareYourself( double timeNew, double dT );

    virtual void computeYourself( double timeNew, double dT ) = 0;

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

    virtual void incrementDeformation( const TensorD&  displacementIncrement,
                                       const TensorDD& displacementGradientIncrement ) = 0;

    struct {
      Fastor::Tensor< double, nDim, nDim > S; // Kirchhoff stress
    } response;

    struct {
      Fastor::Tensor< double, nDim, nDim, nDim, nDim > dS_dDeltaF;
    } tangents;

    TensorDD dx_dY() const { return state->dx_dY( Fastor::seq( 0, nDim ), Fastor::seq( 0, nDim ) ); };

    TensorDD dY_dX() const { return state->dY_dX( Fastor::seq( 0, nDim ), Fastor::seq( 0, nDim ) ); };

    TensorD getVelocity() const { return state->v( Fastor::seq( 0, nDim ) ); };

    TensorD getAcceleration() const { return state->a( Fastor::seq( 0, nDim ) ); };

    void setVelocity( const TensorD& velocity )
    {
      for ( int i = 0; i < nDim; i++ )
        state->v( i ) = velocity( i );
    };

    void setAcceleration( const TensorD& acceleration )
    {
      for ( int i = 0; i < nDim; i++ )
        state->a( i ) = acceleration( i );
    };

    virtual void setInitialCondition( const std::string& conditionName, const double* value ) override{};
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
    material = std::unique_ptr< Material >( MarmotLibrary::MarmotMaterialFiniteStrainFactory::createMaterial(
      section.materialName,
      section.materialProperties,
      section.nMaterialProperties,
      _mpNumber ) );

    if ( !material )
      throw std::invalid_argument( MakeString()
                                   << __PRETTY_FUNCTION__
                                   << ": invalid finite strain material assigned!" );
  }

  template < int nDim >
  void DisplacementMaterialPoint< nDim >::initializeYourself()
  {
    state->dY_dX.eye();
    /* state->dx_dY.eye(); */
    this->prepareYourself( 0, 0 );
    material->initializeYourself( state->materialState.data(), state->materialState.size() );
  }

  template < int nDim >
  void DisplacementMaterialPoint< nDim >::prepareYourself( double timeNew, double dT )
  {
    state->du.zeros();
    state->dx_dY.eye();
  }

  class DisplacementMaterialPoint2D : public DisplacementMaterialPoint< 2 > {

  public:
    using DisplacementMaterialPoint::DisplacementMaterialPoint;

    void computeYourself( double timeNew, double dT );

    void incrementDeformation( const TensorD& displacementIncrement, const TensorDD& displacementGradientIncrement );
  };

  class DisplacementMaterialPoint3D : public DisplacementMaterialPoint< 3 > {

  public:
    using DisplacementMaterialPoint::DisplacementMaterialPoint;

    void computeYourself( double timeNew, double dT );

    void incrementDeformation( const TensorD& displacementIncrement, const TensorDD& displacementGradientIncrement );
  };

} // namespace Marmot::MaterialPoints
