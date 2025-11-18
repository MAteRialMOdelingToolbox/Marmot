#include "Marmot/LinearViscoelasticPowerLaw.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotViscoelasticity.h"

using namespace Marmot;
using namespace Eigen;

namespace Marmot::Materials {

  LinearViscoelasticPowerLaw::LinearViscoelasticPowerLaw( const double* materialProperties,
                                                          int           nMaterialProperties,
                                                          int           materialLabel )
    : MarmotMaterialHypoElastic( materialProperties, nMaterialProperties, materialLabel ),
      // clang-format off
      // elastic parameters
      E                                 ( materialProperties[0] ),
      nu                                ( materialProperties[1] ),
      // viscoelastic parameters
      m                                 ( materialProperties[2] ),
      n                                 ( materialProperties[3] ),
      nKelvin                           ( static_cast< size_t > ( materialProperties[4] ) ),
      minTau                            ( materialProperties[5] ),
      timeToDays                        ( materialProperties[6] )
  // clang-format on
  {
    initializeStateLayout();

    // assume sqrt( 10 ) spacing between retardation times
    retardationTimes = KelvinChain::generateRetardationTimes( nKelvin, minTau, sqrt( 10. ) );

    using namespace Marmot::ContinuumMechanics::Viscoelasticity;
    auto phi_ = [&]( autodiff::Real< powerLawApproximationOrder, double > tau ) {
      return ComplianceFunctions::powerLaw( tau, m, n );
    };

    elasticModuli = KelvinChain::computeElasticModuli< powerLawApproximationOrder >( phi_, retardationTimes );

    // for 2nd order approximations
    zerothKelvinChainCompliance = m * ( 1. - n ) * pow( 2., n ) * pow( minTau / sqrt( sqrt( 10. ) ), n );
  }

  void LinearViscoelasticPowerLaw::computeStress( state3D&        state,
                                                  double*         dStressDDStrain,
                                                  const double*   dStrain,
                                                  const timeInfo& timeInfo )

  {
    mVector6d nomStress( state.stress.data() );
    Vector6d  dE( dStrain );
    mMatrix6d C( dStressDDStrain );

    const double& dT = timeInfo.dT;

    if ( ( dE.array() == 0 ).all() && dT == 0 ) {
      C = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );
      return;
    }

    Eigen::Map< KelvinChain::StateVarMatrix > creepStateVars( state.stateVars, 6, nKelvin );

    const double dTimeDays = dT * timeToDays;

    Matrix6d CelUnitInv = ContinuumMechanics::Elasticity::Isotropic::complianceTensor( 1.0, nu );

    Vector6d creepStrainIncrement = Vector6d::Zero();
    double   creepCompliance      = 0;

    KelvinChain::evaluateKelvinChain( dTimeDays,
                                      elasticModuli,
                                      retardationTimes,
                                      creepStateVars,
                                      creepCompliance,
                                      creepStrainIncrement,
                                      1.0 );

    using namespace Marmot::ContinuumMechanics::Viscoelasticity;
    double effectiveCompliance = 1. / E + zerothKelvinChainCompliance + creepCompliance;

    C                    = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( 1. / effectiveCompliance, nu );
    Vector6d deltaStress = C * ( dE - creepStrainIncrement );
    nomStress            = nomStress + deltaStress;

    KelvinChain::updateStateVarMatrix( dTimeDays,
                                       elasticModuli,
                                       retardationTimes,
                                       creepStateVars,
                                       deltaStress,
                                       CelUnitInv );

    return;
  }
} // namespace Marmot::Materials
