/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck.
 *
 * Thomas Mader thomas.mader@boku.ac.at
 *
 * This file is part of the MAteRialMOdellingToolbox (marmot).
 * LGPL v2.1+, see LICENSE.md at the top level directory of marmot.
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
   * @brief Material point for gradient-enhanced (implicit-gradient) finite-strain materials
   *        WITHOUT a micropolar continuum.
   *
   * The non-micropolar sibling of @ref GradientEnhancedMicropolarMaterialPoint: it carries the
   * displacement field and a single scalar nonlocal field @f$\bar N@f$, and drives a
   * @ref MarmotMaterialGradientEnhancedFiniteStrain, which returns the Kirchhoff stress
   * @f$\boldsymbol\tau@f$, the local driving force @f$L@f$ (the source term of
   * @f$\bar N - c\nabla^2\bar N = L@f$) and the four algorithmic tangents.  There is no
   * microrotation field and hence no couple stress.
   *
   * Kinematics follow the semi-Lagrangian convention of the meshfree framework: `dY_dX` holds
   * the total deformation gradient of the last accepted state (reference X -> intermediate
   * reference Y) and `dx_dY` the increment accumulated since, so the deformation gradient handed
   * to the material is @f$\boldsymbol F = \mathrm{dx\_dY}\cdot\mathrm{dY\_dX}@f$.
   */
  template < int nDim >
  class GradientEnhancedFiniteStrainMaterialPoint : public MarmotMaterialPoint {

  protected:
    constexpr static int _nVertices = 1;

    using TensorD    = Fastor::Tensor< double, nDim >;
    using TensorDD   = Fastor::Tensor< double, nDim, nDim >;
    using TensorDDDD = Fastor::Tensor< double, nDim, nDim, nDim, nDim >;

    int _mpNumber;

    TensorD _x0;

    double _vol0;
    double _density;

    using Material = MarmotMaterialGradientEnhancedFiniteStrain;

    std::unique_ptr< MarmotMaterialGradientEnhancedFiniteStrain > material;

    class MPStateVarManager : public MarmotStateVarVectorManager {

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
      FastorStandardTensors::TensorMap3d  u;
      FastorStandardTensors::TensorMap3d  v;
      FastorStandardTensors::TensorMap3d  a;
      FastorStandardTensors::TensorMap3d  du;
      FastorStandardTensors::TensorMap33d dx_dY;
      FastorStandardTensors::TensorMap33d dY_dX;
      double&                             nonLocalDamage;
      double&                             localDamage;
      FastorStandardTensors::TensorMap33d stress;
      double&                             F0_XX;
      double&                             F0_YY;
      double&                             F0_ZZ;
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
          nonLocalDamage( find( "nonlocal damage" ) ),
          localDamage( find( "local damage" ) ),
          stress( &find( "stress" ) ),
          F0_XX( find( "F0 XX" ) ),
          F0_YY( find( "F0 YY" ) ),
          F0_ZZ( find( "F0 ZZ" ) ),
          materialState( &find( "begin of material state" ), nStateVars - getNumberOfRequiredStateVars() ){};
    };

    std::unique_ptr< MPStateVarManager > state;

  public:
    bool hasEigenDeformation;

    GradientEnhancedFiniteStrainMaterialPoint( int           mpNumber,
                                               const double* vertexCoordinates,
                                               int           nVertexCoordinates,
                                               double        volume )
      : _mpNumber( mpNumber ), hasEigenDeformation( false )
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

    const TensorD& coordinates() const { return _x0; };

    virtual void prepareYourself( double timeNew, double dT );

    virtual void computeYourself( double timeNew, double dT ) = 0;

    virtual void acceptStateAndPosition()
    {
      const auto&                           u_n  = state->u;
      const FastorStandardTensors::Tensor3d u_np = u_n + state->du;

      const FastorStandardTensors::Tensor33d dx_dX_n  = state->dY_dX;
      const FastorStandardTensors::Tensor33d dx_dX_np = state->dx_dY % dx_dX_n;

      mapEigenToFastor( state->u )     = mapEigenToFastor( u_np );
      mapEigenToFastor( state->dY_dX ) = mapEigenToFastor( dx_dX_np );
    };

    virtual void incrementDeformation( const TensorD&  displacementIncrement,
                                       const TensorDD& displacementGradientIncrement,
                                       double          nonLocalDamage ) = 0;

    struct {
      Fastor::Tensor< double, nDim, nDim > S;  ///< Kirchhoff stress tau
      double                               dL; ///< increment of the local damage driving force
      double                               nonLocalRadius;
    } response;

    struct {
      Fastor::Tensor< double, nDim, nDim, nDim, nDim > dS_dDeltaF;
      Fastor::Tensor< double, nDim, nDim >             dS_dN;
      Fastor::Tensor< double, nDim, nDim >             dL_dDeltaF;
      double                                           dL_dN;
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

  class GradientEnhancedFiniteStrainMaterialPoint2D : public GradientEnhancedFiniteStrainMaterialPoint< 2 > {

  public:
    using GradientEnhancedFiniteStrainMaterialPoint::GradientEnhancedFiniteStrainMaterialPoint;

    void computeYourself( double timeNew, double dT );

    void incrementDeformation( const TensorD&  displacementIncrement,
                               const TensorDD& displacementGradientIncrement,
                               double          nonLocalDamage );
  };

  class GradientEnhancedFiniteStrainMaterialPoint3D : public GradientEnhancedFiniteStrainMaterialPoint< 3 > {

  public:
    using GradientEnhancedFiniteStrainMaterialPoint::GradientEnhancedFiniteStrainMaterialPoint;

    void computeYourself( double timeNew, double dT );

    void incrementDeformation( const TensorD&  displacementIncrement,
                               const TensorDD& displacementGradientIncrement,
                               double          nonLocalDamage );
  };

} // namespace Marmot::MaterialPoints
