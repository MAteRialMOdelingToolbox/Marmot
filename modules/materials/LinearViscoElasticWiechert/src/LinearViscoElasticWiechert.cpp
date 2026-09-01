#include "Marmot/LinearViscoElasticWiechert.h"

#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotViscoelasticity.h"

#include <Eigen/Core>

#include <stdexcept>

namespace Marmot::Materials {

  namespace {

    const double& checkedMaterialProperty( const double* materialProperties, int nMaterialProperties, int index )
    {
      if ( nMaterialProperties < 7 ) {
        throw std::invalid_argument( "LinearViscoElasticWiechert requires at least 7 material properties." );
      }
      if ( materialProperties == nullptr ) {
        throw std::invalid_argument( "LinearViscoElasticWiechert requires a valid material property array." );
      }

      return materialProperties[index];
    }

  } // namespace

  LinearViscoElasticWiechert::LinearViscoElasticWiechert( const double* materialProperties,
                                                          int           nMaterialProperties,
                                                          int           materialNumber )
    : MarmotMaterialHypoElastic( materialProperties, nMaterialProperties, materialNumber ),
      E( checkedMaterialProperty( materialProperties, nMaterialProperties, 0 ) ),
      nu( checkedMaterialProperty( materialProperties, nMaterialProperties, 1 ) ),
      m( checkedMaterialProperty( materialProperties, nMaterialProperties, 2 ) ),
      n( checkedMaterialProperty( materialProperties, nMaterialProperties, 3 ) ),
      nMaxwell( static_cast< size_t >( checkedMaterialProperty( materialProperties, nMaterialProperties, 4 ) ) ),
      minTau( checkedMaterialProperty( materialProperties, nMaterialProperties, 5 ) ),
      timeToDays( checkedMaterialProperty( materialProperties, nMaterialProperties, 6 ) ),
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
    relaxationTimes      = Wiechert::generateRelaxationTimes( static_cast< int >( nMaxwell ), minTau, spacing );

    using namespace Marmot::ContinuumMechanics::Viscoelasticity;
    auto relaxationFunction = [&]( autodiff::Real< powerLawApproximationOrder, double > time ) {
      return RelaxationFunctions::powerLaw( time, m, n );
    };

    elasticModuli = Wiechert::computeElasticModuli< powerLawApproximationOrder >( relaxationFunction,
                                                                                  relaxationTimes,
                                                                                  spacing );

    const double spectrumCoefficient     = m * n * ( n + 1. ) * std::pow( 2., -n );
    const double slowestResolvedBoundary = relaxationTimes( relaxationTimes.size() - 1 ) * std::sqrt( spacing );
    zerothWiechertStiffness              = spectrumCoefficient / n * std::pow( slowestResolvedBoundary, -n );
  }

  void LinearViscoElasticWiechert::computeStress( state3D&                state,
                                                  Marmot::Matrix6d&       dStressDDStrain,
                                                  const Marmot::Vector6d& dStrain,
                                                  const timeInfo&         timeInfo ) const
  {
    mVector6d                    nominalStress( state.stress.data() );
    mMatrix6d                    tangent( dStressDDStrain.data() );
    Eigen::Map< const Vector6d > strainIncrement( dStrain.data() );

    if ( ( strainIncrement.array() == 0.0 ).all() && timeInfo.dT == 0.0 ) {
      tangent = E * unitStiffness;
      return;
    }

    Eigen::Map< Wiechert::StateVarMatrix > maxwellStateVars( state.stateVars, 6, nMaxwell );
    const double                           dTimeDays = timeInfo.dT * timeToDays;

    Vector6d maxwellStressIncrement = Vector6d::Zero();
    double   maxwellStiffness       = 0.0;

    Wiechert::evaluateWiechert( dTimeDays,
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
    tangent                        = effectiveStiffness * unitStiffness;
    const Vector6d stressIncrement = tangent * strainIncrement - maxwellStressIncrement;
    nominalStress += stressIncrement;

    Wiechert::updateStateVarMatrix( dTimeDays,
                                    elasticModuli,
                                    relaxationTimes,
                                    maxwellStateVars,
                                    strainIncrement,
                                    unitStiffness );
  }

  double LinearViscoElasticWiechert::getDensity( const double* ) const
  {
    if ( nMaterialProperties < 8 ) {
      throw std::runtime_error( MakeString()
                                << __PRETTY_FUNCTION__ << ": Density not provided in material properties array!" );
    }

    return materialProperties[7];
  }

} // namespace Marmot::Materials
