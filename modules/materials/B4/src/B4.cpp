#include "Marmot/B4.h"
#include "Marmot/B4Shrinkage.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotTypedefs.h"
#include <string>

using namespace Marmot;
using namespace Eigen;

namespace Marmot::Materials {

  B4::B4( const double* materialProperties, int nMaterialProperties, int materialLabel )
    : Marmot::Materials::MaterialHypoElastic( materialProperties, nMaterialProperties, materialLabel ),
      // clang-format off
      // elastic parameters
      nu                                ( materialProperties[0] ),
      q1                                ( materialProperties[1] ),
      // basic creep parameters
      q2                                ( materialProperties[2] ),
      q3                                ( materialProperties[3] ),
      q4                                ( materialProperties[4] ),
      n                                 ( materialProperties[5] ),
      m                                 ( materialProperties[6] ),
      nKelvinBasic                      ( static_cast< size_t > ( materialProperties[7] ) ),
      minTauBasic                       ( materialProperties[8] ),
      // autogenous shrinkage parameters
      ultimateAutogenousShrinkageStrain ( materialProperties[9] ),
      autogenousShrinkageHalfTime       ( materialProperties[10] ),
      alpha                             ( materialProperties[11] ),
      rt                                ( materialProperties[12] ),
      // drying shrinkage parameters
      ultimateDryingShrinkageStrain     ( materialProperties[13] ),
      dryingShrinkageHalfTime           ( materialProperties[14] ),
      dryingStart                       ( materialProperties[15] ),
      hEnv                              ( materialProperties[16] ),
      // drying creep paramters
      q5                                ( materialProperties[17] ),
      nKelvinDrying                     ( static_cast< size_t > ( materialProperties[18] ) ),
      minTauDrying                      ( materialProperties[19] ),
      // additional parameters
      castTime                          ( materialProperties[20] ),
      timeToDays                        ( materialProperties[21] ),
      solidificationParameters          ( { q1, q2, q3, q4, n, m } )
  // clang-format on
  {
    stateLayout.add( "basicCreepStateVars", nKelvinBasic * 6 );
    stateLayout.add( "dryingCreepStateVars", nKelvinDrying * 6 );
    stateLayout.finalize();
    solidificationKelvinProperties
      .retardationTimes = ContinuumMechanics::Viscoelasticity::KelvinChain::generateRetardationTimes( nKelvinBasic,
                                                                                                      minTauBasic,
                                                                                                      10. );

    auto phiBasic = [&]( autodiff::Real< basicCreepComplianceApproximationOrder, double > tau ) {
      return ContinuumMechanics::Viscoelasticity::SolidificationTheory::phi( tau, solidificationParameters );
    };

    solidificationKelvinProperties
      .elasticModuli = ContinuumMechanics::Viscoelasticity::KelvinChain::computeElasticModuli<
      basicCreepComplianceApproximationOrder >( phiBasic, solidificationKelvinProperties.retardationTimes );

    solidificationKelvinProperties.E0 = ContinuumMechanics::Viscoelasticity::SolidificationTheory::
      computeZerothElasticModul( minTauBasic, n, basicCreepComplianceApproximationOrder );
  }

  void B4::computeStress( state3D&                      state,
                          Marmot::Math::Matrix6d&       dStressDDStrain,
                          const Marmot::Math::Vector6d& dStrain,
                          const timeInfo&               timeInfo ) const

  {
    Marmot::Math::mVector6d             nomStress( state.stress.data() );
    Map< const Marmot::Math::Vector6d > dE( dStrain.data() );
    Marmot::Math::mMatrix6d             C( dStressDDStrain.data() );

    const double& dT   = timeInfo.dT;
    const double& time = timeInfo.time;

    if ( ( dE.array() == 0 ).all() && dT == 0 ) {
      C = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( 1e6 / q1, nu );
      return;
    }

    // map state variables of basic creep Kelvin chain
    auto basicCreepStateVars = stateLayout.getAs< Eigen::Map< Eigen::MatrixXd > >( state.stateVars,
                                                                                   "basicCreepStateVars",
                                                                                   6,
                                                                                   nKelvinBasic );

    // map state variables of drying creep Kelvin chain
    auto dryingCreepStateVars = stateLayout.getAs< Eigen::Map< Eigen::MatrixXd > >( state.stateVars,
                                                                                    "dryingCreepStateVars",
                                                                                    6,
                                                                                    nKelvinDrying );

    const double dTimeDays  = dT * timeToDays;
    const double tStartDays = ( time - dT - castTime ) * timeToDays;

    Marmot::Math::Matrix6d CelUnitInv = ContinuumMechanics::Elasticity::Isotropic::complianceTensor( 1.0, nu );

    // compute basic creep strains and compliance
    auto [basicCreepStrainIncrement, basicCreepUniaxialComplianceComponents] = ContinuumMechanics::Viscoelasticity::
      SolidificationTheory::computeCreepStrainIncrementAndComplianceComponents( tStartDays,
                                                                                dTimeDays,
                                                                                1.0,
                                                                                CelUnitInv,
                                                                                nomStress,
                                                                                solidificationParameters,
                                                                                solidificationKelvinProperties,
                                                                                basicCreepStateVars );

    // compute drying creep strains and compliancei
    ContinuumMechanics::Viscoelasticity::KelvinChain::Properties dryingCreepRetardationTimes = ContinuumMechanics::
      Viscoelasticity::KelvinChain::generateRetardationTimes( nKelvinDrying, minTauDrying, sqrt( 10. ) );
    double b      = 8. * ( 1. - hEnv );
    double xiZero = std::min( -1e-16, ( dryingStart - tStartDays - dTimeDays / 2. ) / dryingShrinkageHalfTime );

    auto phiDrying = [&]( autodiff::Real< dryingCreepComplianceApproximationOrder, double > tau ) {
      return phi( tau, b, xiZero );
    };

    ContinuumMechanics::Viscoelasticity::KelvinChain::Properties
      dryingCreepElasticModuli = ContinuumMechanics::Viscoelasticity::KelvinChain::computeElasticModuli<
        dryingCreepComplianceApproximationOrder >( phiDrying, dryingCreepRetardationTimes );

    Marmot::Math::Vector6d dryingCreepStrainIncrement = Marmot::Math::Vector6d::Zero();
    double                 dryingCreepCompliance      = 0;

    ContinuumMechanics::Viscoelasticity::KelvinChain::evaluateKelvinChain( dTimeDays / dryingShrinkageHalfTime,
                                                                           dryingCreepElasticModuli,
                                                                           dryingCreepRetardationTimes,
                                                                           dryingCreepStateVars,
                                                                           dryingCreepCompliance,
                                                                           dryingCreepStrainIncrement,
                                                                           1e-6 * q5 / exp( 4 ) );

    double effectiveCompliance = basicCreepUniaxialComplianceComponents.elastic +
                                 basicCreepUniaxialComplianceComponents.viscoelastic +
                                 basicCreepUniaxialComplianceComponents.flow + dryingCreepCompliance;

    // compute shrinkage strain increment
    Marmot::Math::Vector6d
      shrinkageStrainIncrement = B4Shrinkage::computeShrinkageStrainIncrement( tStartDays,
                                                                               dTimeDays,
                                                                               ultimateAutogenousShrinkageStrain,
                                                                               autogenousShrinkageHalfTime,
                                                                               alpha,
                                                                               rt,
                                                                               ultimateDryingShrinkageStrain,
                                                                               dryingShrinkageHalfTime,
                                                                               1. - hEnv * hEnv * hEnv,
                                                                               dryingStart );

    C = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( 1. / effectiveCompliance, nu );
    Marmot::Math::Vector6d deltaStress = C * ( dE - basicCreepStrainIncrement - dryingCreepStrainIncrement -
                                               shrinkageStrainIncrement );
    nomStress                          = nomStress + deltaStress;

    ContinuumMechanics::Viscoelasticity::KelvinChain::updateStateVarMatrix( dTimeDays / dryingShrinkageHalfTime,
                                                                            dryingCreepElasticModuli,
                                                                            dryingCreepRetardationTimes,
                                                                            dryingCreepStateVars,
                                                                            deltaStress,
                                                                            CelUnitInv );

    ContinuumMechanics::Viscoelasticity::KelvinChain::updateStateVarMatrix( dTimeDays,
                                                                            solidificationKelvinProperties
                                                                              .elasticModuli,
                                                                            solidificationKelvinProperties
                                                                              .retardationTimes,
                                                                            basicCreepStateVars,
                                                                            deltaStress,
                                                                            CelUnitInv );
    return;
  }

  double B4::getDensity( const double* stateVars ) const
  {
    if ( nMaterialProperties >= 23 )
      return materialProperties[22];
    throw std::runtime_error(
      std::string( Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << ": Density not specified for B4." ) );
  }

} // namespace Marmot::Materials
