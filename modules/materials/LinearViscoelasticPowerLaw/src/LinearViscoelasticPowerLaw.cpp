#include "Marmot/LinearViscoelasticPowerLaw.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotViscoelasticity.h"

using namespace Marmot;
using namespace Eigen;

namespace Marmot::Materials {

  double LinearViscoelasticPowerLaw::getDensity( [[maybe_unused]] const double* stateVars ) const
  {
    if ( nMaterialProperties < 8 ) {
      throw std::runtime_error( Marmot::Utilities::MakeString()
                                << __PRETTY_FUNCTION__ << ": Density not provided in material properties array!" );
    }
    else {
      return materialProperties[7];
    }
  }
  LinearViscoelasticPowerLaw::LinearViscoelasticPowerLaw( const double* materialProperties,
                                                          int           nMaterialProperties,
                                                          int           materialLabel )
    : Marmot::Materials::MaterialHypoElastic( materialProperties, nMaterialProperties, materialLabel ),
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
    stateLayout.add( "kelvinStateVars", 6 * nKelvin );
    stateLayout.finalize();

    // assume sqrt( 10 ) spacing between retardation times
    retardationTimes = ContinuumMechanics::Viscoelasticity::KelvinChain::generateRetardationTimes( static_cast< int >(
                                                                                                     nKelvin ),
                                                                                                   minTau,
                                                                                                   sqrt( 10. ) );

    using namespace Marmot::ContinuumMechanics::Viscoelasticity;
    auto phi_ = [&]( autodiff::Real< powerLawApproximationOrder, double > tau ) {
      return ComplianceFunctions::powerLaw( tau, m, n );
    };

    elasticModuli = ContinuumMechanics::Viscoelasticity::KelvinChain::computeElasticModuli<
      powerLawApproximationOrder >( phi_, retardationTimes );

    // for 2nd order approximations
    zerothKelvinChainCompliance = m * ( 1. - n ) * pow( 2., n ) * pow( minTau / sqrt( sqrt( 10. ) ), n );
  }

  void LinearViscoelasticPowerLaw::computeStress( state3D&                      state,
                                                  Marmot::Math::Matrix6d&       dStressDDStrain,
                                                  const Marmot::Math::Vector6d& dStrain,
                                                  const timeInfo&               timeInfo ) const

  {
    Marmot::Math::mVector6d             nomStress( state.stress.data() );
    Map< const Marmot::Math::Vector6d > dE( dStrain.data() );
    Marmot::Math::mMatrix6d             C( dStressDDStrain.data() );

    const double& dT = timeInfo.dT;

    if ( ( dE.array() == 0 ).all() && dT == 0 ) {
      C = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );
      return;
    }

    Eigen::Map< ContinuumMechanics::Viscoelasticity::KelvinChain::StateVarMatrix > creepStateVars( state.stateVars,
                                                                                                   6,
                                                                                                   nKelvin );

    const double dTimeDays = dT * timeToDays;

    Marmot::Math::Matrix6d CelUnitInv = ContinuumMechanics::Elasticity::Isotropic::complianceTensor( 1.0, nu );

    Marmot::Math::Vector6d creepStrainIncrement = Marmot::Math::Vector6d::Zero();
    double                 creepCompliance      = 0;

    ContinuumMechanics::Viscoelasticity::KelvinChain::evaluateKelvinChain( dTimeDays,
                                                                           elasticModuli,
                                                                           retardationTimes,
                                                                           creepStateVars,
                                                                           creepCompliance,
                                                                           creepStrainIncrement,
                                                                           1.0 );

    using namespace Marmot::ContinuumMechanics::Viscoelasticity;
    double effectiveCompliance = 1. / E + zerothKelvinChainCompliance + creepCompliance;

    C = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( 1. / effectiveCompliance, nu );
    Marmot::Math::Vector6d deltaStress = C * ( dE - creepStrainIncrement );
    nomStress                          = nomStress + deltaStress;

    ContinuumMechanics::Viscoelasticity::KelvinChain::updateStateVarMatrix( dTimeDays,
                                                                            elasticModuli,
                                                                            retardationTimes,
                                                                            creepStateVars,
                                                                            deltaStress,
                                                                            CelUnitInv );

    return;
  }
} // namespace Marmot::Materials
