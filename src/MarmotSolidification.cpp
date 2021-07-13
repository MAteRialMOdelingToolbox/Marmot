#include "Marmot/MarmotSolidification.h"

namespace Marmot::Materials {

  namespace SolidificationTheory {

    double solidifiedVolume( double timeInDays, Parameters params )
    {
      return 1. / ( pow( params.lambda0 / timeInDays, params.m ) + params.q3 / params.q2 );
    }

    double computeZerothElasticModul( double tauMin, double n, int order )
    {

      const double tau0 = tauMin / sqrt( 10. );
      double       E0;
      switch ( order ) {
      case 1: E0 = 1. / ( log( 1. + pow( tau0, n ) ) ); break;
      case 2:
        E0 = 1. / ( log( 1. + pow( 2. * tau0, n ) ) - n * pow( 2. * tau0, n ) / ( 1. + pow( 2. * tau0, n ) ) );
        break;
      case 3: E0 = 1. / ( ( 1. - n ) * log( 1. + pow( 3. * tau0, n ) ) ); break;
      default: throw std::invalid_argument( "Invalid approximation order requested" ); break;
      }
      return E0;
    }

    Result computeCreepStrainIncrementAndComplianceComponents(
      double                                                 tStartDays,
      double                                                 dTimeDays,
      double                                                 amplificationFactor,
      const Marmot::Matrix6d&                                flowUnitCompliance,
      const Marmot::Vector6d&                                stressOld,
      const Parameters&                                      parameters,
      const KelvinChainProperties&                           kelvinChainProperties,
      const Eigen::Ref< const KelvinChain::StateVarMatrix >& kelvinStateVars )
    {

      Vector6d                     dECreep = Vector6d::Zero();
      UniaxialComplianceComponents complianceComponents;

      const double v = solidifiedVolume( tStartDays + dTimeDays / 2., parameters );

      complianceComponents.elastic = parameters.q1 * 1e-6;

      if ( dTimeDays <= 1e-14 )
        return { dECreep, complianceComponents };

      complianceComponents.flow = amplificationFactor * parameters.q4 / ( tStartDays + dTimeDays / 2. ) / 2. *
                                  dTimeDays * 1e-6;
      dECreep += complianceComponents.flow * 2. * flowUnitCompliance * stressOld;

      complianceComponents.viscoelastic = parameters.q2 / ( v * kelvinChainProperties.E0 ) * amplificationFactor * 1e-6;

      KelvinChain::evaluateKelvinChain( dTimeDays,
                                        kelvinChainProperties.elasticModuli,
                                        kelvinChainProperties.retardationTimes,
                                        kelvinStateVars,
                                        complianceComponents.viscoelastic,
                                        dECreep,
                                        1e-6 * amplificationFactor * parameters.q2 / v );

      return { dECreep, complianceComponents };
    }

  } // namespace SolidificationTheory
} // namespace Marmot::Materials
