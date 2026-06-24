#include "Marmot/LinearViscoElasticWiechert.h"

#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotViscoelasticity.h"

#include <Eigen/Core>

#include <stdexcept>

namespace Marmot::Materials {

  LinearViscoElasticWiechert::LinearViscoElasticWiechert( const double* materialProperties,
                                                          int           nMaterialProperties,
                                                          int           materialNumber )
    : MarmotMaterialHypoElastic( materialProperties, nMaterialProperties, materialNumber ),
      E( materialProperties[0] ),
      nu( materialProperties[1] ),
      m( materialProperties[2] ),
      n( materialProperties[3] ),
      nMaxwell( static_cast< size_t >( materialProperties[4] ) ),
      minTau( materialProperties[5] ),
      timeToDays( materialProperties[6] ),
      zerothWiechertStiffness( 0.0 )
  {
    if ( nMaterialProperties < 7 ) {
      throw std::invalid_argument( "LinearViscoElasticWiechert requires at least 7 material properties." );
    }
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
      tangent = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );
      return;
    }

    Eigen::Map< Wiechert::StateVarMatrix > maxwellStateVars( state.stateVars, 6, nMaxwell );
    const double                           dTimeDays = timeInfo.dT * timeToDays;

    const Matrix6d unitStiffness = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( 1.0, nu );

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

    tangent = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( effectiveStiffness, nu );
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
