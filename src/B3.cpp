#include "Marmot/B3.h"
#include "Marmot/DryingCreepComplianceSpectrumApproximation.h"
#include "Marmot/LogPowerLawSpectrumApproximation.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotUtility.h"
#include "Marmot/MarmotVoigt.h"
#include "Marmot/Solidification.h"
#include <iostream>
#include <locale>
#include <map>
#include <stdexcept>
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
      numberOfKelvinUnits     ( static_cast< size_t > ( materialProperties[11] ) ),
      minimalRetardationTime  ( materialProperties[12] ),
      dTStatic                ( materialProperties[13] ),
      timeToDays              ( materialProperties[14] ),
      castTime                ( materialProperties[15] ),
      basicCreepMaterialParameters { q1, q2, q3, q4, n, m }
  // clang-format on
  {
    double                           kelvinE0;
    Solidification::KelvinProperties kelvinElasticModuli( numberOfKelvinUnits );
    Solidification::KelvinProperties
      kelvinRetardationTimes = SpectrumApproximation::generateLogSpacedRetardationTimes( minimalRetardationTime,
                                                                                         numberOfKelvinUnits );
    SpectrumApproximation::computeApproximation( kelvinElasticModuli,
                                                 kelvinRetardationTimes,
                                                 kelvinE0,
                                                 n,
                                                 SpectrumApproximation::ApproximationOrder::secondOrder );

    kelvinChainProperties = { kelvinE0, kelvinElasticModuli, kelvinRetardationTimes };
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

    double&                               EStatic = stateVarManager->EStatic;
    Eigen::Ref< Solidification::mKelvinStateVarMatrix > kelvinStateVars( ( stateVarManager->kelvinStateVars ).leftCols( numberOfKelvinUnits ) );
    Eigen::Ref< Solidification::mKelvinStateVarMatrix > dryingCreepStateVars( ( stateVarManager->kelvinStateVars ).rightCols( numberOfKelvinUnits )   );

    const double dTimeDays  = dT * timeToDays;
    const double tStartDays = ( timeOld[1] - castTime ) * timeToDays;

    Vector6d dECreep;
    Matrix6d CelUnitInv = ContinuumMechanics::Elasticity::Isotropic::complianceTensor( 1.0, nu );

    Solidification solidification( basicCreepMaterialParameters, kelvinChainProperties, dTStatic );

    EStatic = solidification.getStaticModulus( tStartDays, kelvinStateVars );
    /// 2. Calculate effective stiffness and creep strains according to %Solidification Theory.
    double E = solidification.getStiffnessAndCreepStrains( tStartDays,
                                                           dTimeDays,
                                                           CelUnitInv,
                                                           nomStress,
                                                           kelvinStateVars,
                                                           dECreep,
                                                           1.0 );

    // compute drying creep strains and compliance
    Solidification::KelvinProperties dryingCreepElasticModuli( numberOfKelvinUnits );
    Solidification::KelvinProperties dryingCreepRetardationTimes = SpectrumApproximation::DryingCreepB3::
      generateRetardationTimes( minimalRetardationTime, numberOfKelvinUnits );
    SpectrumApproximation::DryingCreepB3::
      computeApproximation( dryingCreepElasticModuli,
                            dryingCreepRetardationTimes,
                            hEnv,
                            dT / 2. / shrinkageHalfTime,
                            ( castTime - ( timeOld[1] + dT / 2 ) ) / shrinkageHalfTime,
                            SpectrumApproximation::DryingCreepB3::ApproximationOrder::firstOrder );

    Vector6d dEDryingCreep;
    dEDryingCreep.setZero();
    double invE = 0;

    for ( int i = 0; i < dryingCreepRetardationTimes.size(); i++ ) {

      double&       tau = dryingCreepRetardationTimes( i );
      const double& D   = dryingCreepElasticModuli( i );

      double lambda, beta;
      solidification.computeLambdaAndBeta( dT, tau, lambda, beta );
      invE +=  q5 / exp( 4 ) * ( 1 - lambda ) / D;

      dEDryingCreep += q5 / exp( 4 ) * ( 1 - beta ) * dryingCreepStateVars.col( i );
    }
    invE *= 1e-6;
    dEDryingCreep *= 1e-6;

    //                               std::cout << "dECreep = " << dECreep << std::endl;
    /// 3. Calculate the increment of the stress tensor as \f$ \Delta\boldsymbol{\sigma} = \mathbb{C} :
    /// \Delta\boldsymbol{\varepsilon}^{\scriptsize\mbox{el}} \f$.
    C                    = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( 1. / ( 1. / E + invE ), nu );
    Vector6d deltaStress = C * ( dE - dECreep - dEDryingCreep );
    nomStress            = nomStress + deltaStress;

    //    std::cout << "nominal stress = " << nomStress << std::endl;

    /// 4. Update the state variables of the solidifying Kelvin chain.
    solidification.updateCreepStateVars( tStartDays, dTimeDays, CelUnitInv, deltaStress, kelvinStateVars );
    solidification.updateCreepStateVars( tStartDays, dTimeDays, CelUnitInv, deltaStress, dryingCreepStateVars );
    return;
  }

  void B3::assignStateVars( double* stateVars_, int nStateVars )
  {
    if ( nStateVars < getNumberOfRequiredStateVars() )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": Not sufficient stateVars!" );

    this->stateVarManager = std::make_unique< B3StateVarManager >( stateVars_, 2 * numberOfKelvinUnits );

    MarmotMaterial::assignStateVars( stateVars_, nStateVars );
  }

  StateView B3::getStateView( const std::string& stateName ) { return stateVarManager->getStateView( stateName ); }

  int B3::getNumberOfRequiredStateVars()
  {
    return B3StateVarManager::layout.nRequiredStateVars + numberOfKelvinUnits * 6 * 2;
  }
} // namespace Marmot::Materials
