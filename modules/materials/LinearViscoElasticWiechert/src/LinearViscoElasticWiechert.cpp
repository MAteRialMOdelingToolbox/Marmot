#include "Marmot/LinearViscoElasticWiechert.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotViscoelasticity.h"
#include "Marmot/MarmotWiechert.h"

#include <Eigen/Core>

#include <stdexcept>

namespace Marmot::Materials {

  using namespace Eigen;
  using namespace Marmot;

  LinearViscoElasticWiechert::LinearViscoElasticWiechert( const double* materialProperties,
                                                          int           nMaterialProperties,
                                                          int           materialNumber )
    : MarmotMaterialHypoElastic( materialProperties, nMaterialProperties, materialNumber ),
      // clang-format off
      // elasticity parameters
      E( materialProperties[0] ),
      nu( materialProperties[1] ),
      // power-law relaxation parameters
      m( materialProperties[2] ),
      n( materialProperties[3] ),
      nMaxwell( static_cast< size_t >( materialProperties[4] ) ),
      minTau( materialProperties[5] ),
      timeToDays( materialProperties[6] ),
      zerothWiechertStiffness( 0.0 )
  // clang-format on
  {
    if ( materialProperties[4] < 1.0 )
      throw std::invalid_argument( "LinearViscoElasticWiechert requires at least one Maxwell element." );
    if ( minTau <= 0.0 )
      throw std::invalid_argument( "LinearViscoElasticWiechert requires a positive minimum relaxation time." );
    if ( m < 0.0 || n <= 0.0 )
      throw std::invalid_argument( "LinearViscoElasticWiechert requires m >= 0 and n > 0." );

    initializeStateLayout();

    const double spacing = std::sqrt( 10. );
    relaxationTimes      = Wiechert::generateRelaxationTimes( nMaxwell, minTau, spacing );

    using namespace Marmot::ContinuumMechanics::Viscoelasticity;
    auto relaxationFunction = [&]( autodiff::Real< powerLawApproximationOrder, double > time ) {
      return RelaxationFunctions::powerLaw( time, m, n );
    };
    elasticModuli = Wiechert::computeElasticModuli< powerLawApproximationOrder >( relaxationFunction,
                                                                                  relaxationTimes,
                                                                                  spacing );

    // Maxwell branches slower than the resolved range behave elastically over the represented time horizon.
    const double spectrumCoefficient     = m * n * ( n + 1. ) * std::pow( 2., -n );
    const double slowestResolvedBoundary = relaxationTimes( relaxationTimes.size() - 1 ) * std::sqrt( spacing );
    zerothWiechertStiffness              = spectrumCoefficient / n * std::pow( slowestResolvedBoundary, -n );
  }

  void LinearViscoElasticWiechert::computeStress( state3D&        state,
                                                  Matrix6d&       dStressDDStrain,
                                                  const Vector6d& dStrain,
                                                  const timeInfo& timeInfo ) const
  {
    mVector6d     nomStress( state.stress.data() );
    mMatrix6d     D( dStressDDStrain.data() );
    const auto    dE = Map< const Vector6d >( dStrain.data() );
    const double& dT = timeInfo.dT;

    if ( ( dE.array() == 0 ).all() && dT == 0 ) {
      D = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );
      return;
    }

    Eigen::Map< Wiechert::StateVarMatrix > creepStateVars( state.stateVars, 6, nMaxwell );
    const double                           dTimeDays = dT * timeToDays;

    Matrix6d CelUnit = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( 1.0, nu );

    Vector6d creepStressIncrement = Vector6d::Zero();
    double   creepStiffness       = 0;

    Wiechert::evaluateWiechert( dTimeDays,
                                elasticModuli,
                                relaxationTimes,
                                creepStateVars,
                                creepStiffness,
                                creepStressIncrement,
                                1.0 );

    double effectiveStiffness = E + zerothWiechertStiffness + creepStiffness;

    D                    = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( effectiveStiffness, nu );
    Vector6d deltaStress = D * dE - creepStressIncrement;
    nomStress            = nomStress + deltaStress;

    Wiechert::updateStateVarMatrix( dTimeDays, elasticModuli, relaxationTimes, creepStateVars, dE, CelUnit );

    return;
  }

  double LinearViscoElasticWiechert::getDensity( const double* stateVars ) const
  {
    if ( nMaterialProperties < 8 ) {
      throw std::runtime_error(
        std::string( MakeString() << __PRETTY_FUNCTION__ << ": Density not specified for this material." ) );
    }
    return this->materialProperties[7];
  }
} // namespace Marmot::Materials
