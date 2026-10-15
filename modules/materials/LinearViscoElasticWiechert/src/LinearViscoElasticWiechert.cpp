#include "Marmot/LinearViscoElasticWiechert.h"

#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotUtils.h"
#include "Marmot/MarmotViscoelasticity.h"

#include <Eigen/Core>

#include <stdexcept>

namespace Marmot::Materials {

  LinearViscoElasticWiechert::LinearViscoElasticWiechert( const double* materialProperties,
                                                          int           nMaterialProperties,
                                                          int           materialNumber )
    : Marmot::Materials::MaterialHypoElastic( materialProperties, nMaterialProperties, materialNumber ),
      E( Marmot::Utilities::checkedMaterialProperty( materialProperties, nMaterialProperties, 0 ) ),
      nu( Marmot::Utilities::checkedMaterialProperty( materialProperties, nMaterialProperties, 1 ) ),
      m( Marmot::Utilities::checkedMaterialProperty( materialProperties, nMaterialProperties, 2 ) ),
      n( Marmot::Utilities::checkedMaterialProperty( materialProperties, nMaterialProperties, 3 ) ),
      nMaxwell( static_cast< size_t >(
        Marmot::Utilities::checkedMaterialProperty( materialProperties, nMaterialProperties, 4 ) ) ),
      minTau( Marmot::Utilities::checkedMaterialProperty( materialProperties, nMaterialProperties, 5 ) ),
      timeToDays( Marmot::Utilities::checkedMaterialProperty( materialProperties, nMaterialProperties, 6 ) ),
      zerothWiechertStiffness( 0.0 )
  {
    if ( materialProperties[4] < 1.0 ) {
      throw std::invalid_argument( "LinearViscoElasticWiechert requires at least one Maxwell element." );
    }
    if ( minTau <= 0.0 ) {
      throw std::invalid_argument( "LinearViscoElasticWiechert requires a positive minimum relaxation time." );
    }
    if ( m < 0.0 || n <= 0.0 ) {
      throw std::invalid_argument( "LinearViscoElasticWiechert requires m >= 0 and n > 0." );
    }

    stateLayout.add( "maxwellStateVars", 6 * nMaxwell );
    stateLayout.finalize();

    unitStiffness = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( 1.0, nu );

    const double spacing = std::sqrt( 10. );
    relaxationTimes      = ContinuumMechanics::Viscoelasticity::Wiechert::generateRelaxationTimes( static_cast< int >(
                                                                                                nMaxwell ),
                                                                                              minTau,
                                                                                              spacing );

    using namespace Marmot::ContinuumMechanics::Viscoelasticity;
    auto relaxationFunction = [&]( autodiff::Real< powerLawApproximationOrder, double > time ) {
      return RelaxationFunctions::powerLaw( time, m, n );
    };

    elasticModuli = ContinuumMechanics::Viscoelasticity::Wiechert::computeElasticModuli<
      powerLawApproximationOrder >( relaxationFunction, relaxationTimes, spacing );

    const double spectrumCoefficient     = m * n * ( n + 1. ) * std::pow( 2., -n );
    const double slowestResolvedBoundary = relaxationTimes( relaxationTimes.size() - 1 ) * std::sqrt( spacing );
    zerothWiechertStiffness              = spectrumCoefficient / n * std::pow( slowestResolvedBoundary, -n );
  }

  void LinearViscoElasticWiechert::computeStress( state3D&                      state,
                                                  Marmot::Math::Matrix6d&       dStressDDStrain,
                                                  const Marmot::Math::Vector6d& dStrain,
                                                  const timeInfo&               timeInfo ) const
  {
    Marmot::Math::mVector6d                    nominalStress( state.stress.data() );
    Marmot::Math::mMatrix6d                    tangent( dStressDDStrain.data() );
    Eigen::Map< const Marmot::Math::Vector6d > strainIncrement( dStrain.data() );

    if ( ( strainIncrement.array() == 0.0 ).all() && timeInfo.dT == 0.0 ) {
      tangent = E * unitStiffness;
      return;
    }

    Eigen::Map< ContinuumMechanics::Viscoelasticity::Wiechert::StateVarMatrix > maxwellStateVars( state.stateVars,
                                                                                                  6,
                                                                                                  nMaxwell );
    const double                                                                dTimeDays = timeInfo.dT * timeToDays;

    Marmot::Math::Vector6d maxwellStressIncrement = Marmot::Math::Vector6d::Zero();
    double                 maxwellStiffness       = 0.0;

    ContinuumMechanics::Viscoelasticity::Wiechert::evaluateWiechert( dTimeDays,
                                                                     elasticModuli,
                                                                     relaxationTimes,
                                                                     maxwellStateVars,
                                                                     maxwellStiffness,
                                                                     maxwellStressIncrement,
                                                                     1.0 );

    const double effectiveStiffness = E + zerothWiechertStiffness + maxwellStiffness;

    // stiffnessTensor is exactly linear in the modulus: it builds C from nu alone and scales it by
    // E / ((1+nu)(1-2nu)). Scaling the cached unit tensor is therefore the same quantity, and removes
    // the last 6x6 assembly from the per-iteration path.
    tangent                                      = effectiveStiffness * unitStiffness;
    const Marmot::Math::Vector6d stressIncrement = tangent * strainIncrement - maxwellStressIncrement;
    nominalStress += stressIncrement;

    ContinuumMechanics::Viscoelasticity::Wiechert::updateStateVarMatrix( dTimeDays,
                                                                         elasticModuli,
                                                                         relaxationTimes,
                                                                         maxwellStateVars,
                                                                         strainIncrement,
                                                                         unitStiffness );
  }

  double LinearViscoElasticWiechert::getDensity( const double* ) const
  {
    if ( nMaterialProperties < 8 ) {
      throw std::runtime_error( Marmot::Utilities::MakeString()
                                << __PRETTY_FUNCTION__ << ": Density not provided in material properties array!" );
    }

    return materialProperties[7];
  }

} // namespace Marmot::Materials
