#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotLowerDimensionalStress.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotVoigt.h"
#include <array>

using namespace Eigen;

void Marmot::Materials::MaterialHypoElastic::setCharacteristicElementLength( double length )
{
  characteristicElementLength = length;
}

namespace {

  /**
   * A reduced stress state of a MaterialHypoElastic: ReducedStress::computeReducedStress() around
   * MaterialHypoElastic::computeStress(), each evaluation from the start of the increment. On failure, the state
   * variables are reset to the start of the increment.
   */
  template < std::size_t nP, std::size_t nF >
  auto computeReducedStress( const Marmot::Materials::MaterialHypoElastic&           material,
                             Marmot::Materials::MaterialHypoElastic::state3D&        state,
                             Marmot::Math::Vector6d&                                 dStrain,
                             const std::array< int, nP >&                            p,
                             const std::array< int, nF >&                            f,
                             const Marmot::Materials::MaterialHypoElastic::timeInfo& timeInfo )
  {
    Map< VectorXd >              stateVars( state.stateVars, material.getNumberOfRequiredStateVars() );
    const VectorXd               stateVarsOld            = stateVars;
    const Marmot::Math::Vector6d stressOld               = state.stress;
    const double                 elasticEnergyDensityOld = state.elasticEnergyDensity;
    const double                 dissipationOld          = state.dissipation;

    auto computeStress3D = [&]( const Marmot::Math::Vector6d& dStrain_,
                                Marmot::Math::Vector6d&       stress,
                                Marmot::Math::Matrix6d&       dStress_dStrain ) {
      stateVars                  = stateVarsOld;
      state.stress               = stressOld;
      state.elasticEnergyDensity = elasticEnergyDensityOld;
      state.dissipation          = dissipationOld;
      material.computeStress( state, dStress_dStrain, dStrain_, timeInfo );
      stress = state.stress;
    };

    try {
      return Marmot::ContinuumMechanics::ReducedStress::computeReducedStress( computeStress3D, dStrain, p, f );
    }
    catch ( const Marmot::Utilities::StressUpdateFailed& ) {
      stateVars = stateVarsOld;
      throw;
    }
  }

} // namespace

void Marmot::Materials::MaterialHypoElastic::computePlaneStress( state2D&                      state2D,
                                                                 Marmot::Math::Matrix3d&       dStress_dStrain2D,
                                                                 const Marmot::Math::Vector3d& dStrain2D,
                                                                 const timeInfo&               timeInfo ) const
{
  using namespace Marmot::ContinuumMechanics;
  using namespace Marmot::ContinuumMechanics::Voigt;

  Marmot::Math::Vector6d dStrain = make3DVoigt< VoigtSize::TwoD >( dStrain2D );
  // initial guess: isochoric, gamma13 = gamma23 = 0
  dStrain( 2 ) = -dStrain2D( 0 ) - dStrain2D( 1 );

  state3D state( make3DVoigt< VoigtSize::TwoD >( state2D.stress ),
                 state2D.elasticEnergyDensity,
                 state2D.dissipation,
                 state2D.stateVars );

  const auto condensation = computeReducedStress( *this,
                                                  state,
                                                  dStrain,
                                                  ReducedStress::planeStressPrescribed,
                                                  ReducedStress::planeStressFree,
                                                  timeInfo );

  state2D.stress               = reduce3DVoigt< VoigtSize::TwoD >( state.stress );
  state2D.elasticEnergyDensity = state.elasticEnergyDensity;
  state2D.dissipation          = state.dissipation;
  dStress_dStrain2D            = condensation.dStressP_dStrainP;
}

void Marmot::Materials::MaterialHypoElastic::computeUniaxialStress( state1D&        state1D,
                                                                    double&         dStress_dStrain1D,
                                                                    const double    dStrain1D,
                                                                    const timeInfo& timeInfo ) const
{
  using namespace Marmot::ContinuumMechanics;

  // initial guess: no lateral and shear strains
  Marmot::Math::Vector6d dStrain = Marmot::Math::Vector6d::Zero();
  dStrain( 0 )                   = dStrain1D;
  Marmot::Math::Vector6d stress  = Marmot::Math::Vector6d::Zero();
  stress( 0 )                    = state1D.stress;

  state3D state( stress, state1D.elasticEnergyDensity, state1D.dissipation, state1D.stateVars );

  const auto condensation = computeReducedStress( *this,
                                                  state,
                                                  dStrain,
                                                  ReducedStress::uniaxialStressPrescribed,
                                                  ReducedStress::uniaxialStressFree,
                                                  timeInfo );

  state1D.stress               = state.stress( 0 );
  state1D.elasticEnergyDensity = state.elasticEnergyDensity;
  state1D.dissipation          = state.dissipation;
  dStress_dStrain1D            = condensation.dStressP_dStrainP( 0, 0 );
}

double Marmot::Materials::MaterialHypoElastic::getMaximumWaveSpeed( const state3D& state ) const
{
  const int nStateVars = getNumberOfRequiredStateVars();

  std::vector< double > stateVarsCopy( nStateVars, 0.0 );
  if ( state.stateVars != nullptr && nStateVars > 0 ) {
    std::copy_n( state.stateVars, nStateVars, stateVarsCopy.begin() );
  }

  state3D stateCopy( state.stress, state.elasticEnergyDensity, state.dissipation, stateVarsCopy.data() );

  Marmot::Math::Matrix6d dStress_dStrain = Marmot::Math::Matrix6d::Zero();
  const auto             dStrain         = Marmot::Math::Vector6d::Zero();
  const timeInfo         timeInfo{ 0.0, 1.0 };

  computeStress( stateCopy, dStress_dStrain, dStrain, timeInfo );

  const double maxStiffnessDiagonal = std::max( { dStress_dStrain( 0, 0 ),
                                                  dStress_dStrain( 1, 1 ),
                                                  dStress_dStrain( 2, 2 ),
                                                  dStress_dStrain( 3, 3 ),
                                                  dStress_dStrain( 4, 4 ),
                                                  dStress_dStrain( 5, 5 ) } );
  const double density              = getDensity( stateCopy.stateVars );

  return density > 0.0 ? std::sqrt( std::max( 0.0, maxStiffnessDiagonal ) / density ) : 0.0;
}
