#include "Marmot/B3.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotUtility.h"
#include "Marmot/MarmotVoigt.h"
#include "autodiff/forward/real.hpp"
#include <iostream>
#include <map>
#include <string>

using namespace Marmot;
using namespace Eigen;

namespace Marmot::Materials {

  B3::B3( const double* materialProperties, int nMaterialProperties, int materialLabel )
    : MarmotMaterialHypoElastic( materialProperties, nMaterialProperties, materialLabel ),
      // clang-format off
      nu                      ( materialProperties[0] ),
      q1                      ( materialProperties[1] ),
      q2                      ( materialProperties[2] ),
      q3                      ( materialProperties[3] ),
      q4                      ( materialProperties[4] ),
      q5                      ( materialProperties[5] ),
      hEnv                    ( materialProperties[6] ),
      shrinkageHalfTime       ( materialProperties[7] ),
      ultimateShrinkageStrain ( materialProperties[8] ),
      n                       ( materialProperties[9] ),
      m                       ( materialProperties[10] ),
      nKelvinBasic            ( static_cast< size_t > ( materialProperties[11] ) ),
      minTauBasic             ( materialProperties[12] ),
      nKelvinDrying           ( static_cast< size_t > ( materialProperties[13] ) ),
      minTauDrying            ( materialProperties[14] ),
      dryingStart             ( materialProperties[15] ),
      dTStatic                ( materialProperties[16] ),
      timeToDays              ( materialProperties[17] ),
      castTime                ( materialProperties[18] ),
      solidificationParameters ( { q1, q2, q3, q4, n, m } )
  // clang-format on
  {
    solidificationKelvinProperties.retardationTimes = KelvinChain::generateRetardationTimes( nKelvinBasic,
                                                                                             minTauBasic,
                                                                                             10. );

    auto phiBasic = [&]( autodiff::Real< basicCreepComplianceApproximationOrder, double > tau ) {
      return SolidificationTheory::phi( tau, solidificationParameters );
    };

    solidificationKelvinProperties.elasticModuli = KelvinChain::computeElasticModuli<
      basicCreepComplianceApproximationOrder >( phiBasic, solidificationKelvinProperties.retardationTimes );

    solidificationKelvinProperties
      .E0 = SolidificationTheory::computeZerothElasticModul( minTauBasic, n, basicCreepComplianceApproximationOrder );
  }

  void B3::computeStress( double*       stress,
                          double*       dStressDDStrain,
                          const double* dStrain,
                          const double* timeOld,
                          const double  dT,
                          double&       pNewDT )

  {
    mVector6d nomStress( stress );
    Vector6d  dE( dStrain );
    mMatrix6d C( dStressDDStrain );

    if ( ( dE.array() == 0 ).all() && dT == 0 ) {
      C = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( 1e6 / q1, nu );
      return;
    }

    Eigen::Ref< KelvinChain::mapStateVarMatrix > basicCreepStateVars(
      ( stateVarManager->kelvinStateVars ).leftCols( nKelvinBasic ) );
    Eigen::Ref< KelvinChain::mapStateVarMatrix > dryingCreepStateVars(
      ( stateVarManager->kelvinStateVars ).rightCols( nKelvinDrying ) );

    const double dTimeDays  = dT * timeToDays;
    const double tStartDays = ( timeOld[1] - castTime ) * timeToDays;

    Matrix6d CelUnitInv = ContinuumMechanics::Elasticity::Isotropic::complianceTensor( 1.0, nu );

    // compute basic creep strains and compliance
    auto [basicCreepStrainIncrement, basicCreepUniaxialComplianceComponents] = SolidificationTheory::
      computeCreepStrainIncrementAndComplianceComponents( tStartDays,
                                                          dTimeDays,
                                                          1.0,
                                                          CelUnitInv,
                                                          nomStress,
                                                          solidificationParameters,
                                                          solidificationKelvinProperties,
                                                          basicCreepStateVars );

    // compute drying creep strains and compliance
    KelvinChain::Properties dryingCreepRetardationTimes = KelvinChain::generateRetardationTimes( nKelvinDrying,
                                                                                                 minTauDrying,
                                                                                                 sqrt( 10. ) );
    double                  b                           = 8. * ( 1. - hEnv );
    double                  xiZero = std::min( -1e-16, ( dryingStart - timeOld[1] - dT / 2. ) / shrinkageHalfTime );

    auto phiDrying = [&]( autodiff::Real< dryingCreepComplianceApproximationOrder, double > tau ) {
      return phi( tau, b, xiZero );
    };

    KelvinChain::Properties dryingCreepElasticModuli = KelvinChain::computeElasticModuli<
      dryingCreepComplianceApproximationOrder >( phiDrying, dryingCreepRetardationTimes );

    Vector6d dryingCreepStrainIncrement = Vector6d::Zero();
    double   dryingCreepCompliance      = 0;

    KelvinChain::evaluateKelvinChain( dTimeDays,
                                      dryingCreepElasticModuli,
                                      dryingCreepRetardationTimes,
                                      dryingCreepStateVars,
                                      dryingCreepCompliance,
                                      dryingCreepStrainIncrement,
                                      1e-6 * q5 / exp( 4 ) );

    double effectiveCompliance = basicCreepUniaxialComplianceComponents.elastic +
                                 basicCreepUniaxialComplianceComponents.viscoelastic +
                                 basicCreepUniaxialComplianceComponents.flow + dryingCreepCompliance;

    C                    = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( 1. / effectiveCompliance, nu );
    Vector6d deltaStress = C * ( dE - basicCreepStrainIncrement - dryingCreepStrainIncrement );
    nomStress            = nomStress + deltaStress;

    KelvinChain::updateStateVarMatrix( dTimeDays,
                                       dryingCreepElasticModuli,
                                       dryingCreepRetardationTimes,
                                       dryingCreepStateVars,
                                       deltaStress,
                                       CelUnitInv );

    KelvinChain::updateStateVarMatrix( dTimeDays,
                                       solidificationKelvinProperties.elasticModuli,
                                       solidificationKelvinProperties.retardationTimes,
                                       basicCreepStateVars,
                                       deltaStress,
                                       CelUnitInv );
    return;
  }

  void B3::assignStateVars( double* stateVars_, int nStateVars )
  {
    if ( nStateVars < getNumberOfRequiredStateVars() )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": Not sufficient stateVars!" );

    this->stateVarManager = std::make_unique< B3StateVarManager >( stateVars_, nKelvinBasic + nKelvinDrying );

    MarmotMaterial::assignStateVars( stateVars_, nStateVars );
  }

  StateView B3::getStateView( const std::string& stateName ) { return stateVarManager->getStateView( stateName ); }

  int B3::getNumberOfRequiredStateVars()
  {
    return B3StateVarManager::layout.nRequiredStateVars + ( nKelvinBasic + nKelvinDrying ) * 6;
  }
} // namespace Marmot::Materials
