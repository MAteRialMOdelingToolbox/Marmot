#include "Marmot/LinearViscoelasticOrthotropicPowerLaw.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotViscoelasticity.h"
#include "Marmot/MarmotVoigt.h"
#include <iostream>
#include <map>
#include <string>

using namespace Marmot;
using namespace Eigen;

namespace Marmot::Materials {

  LinearViscoelasticOrthotropicPowerLaw::LinearViscoelasticOrthotropicPowerLaw( const double* materialProperties,
                                                                                int           nMaterialProperties,
                                                                                int           materialLabel )
    : MarmotMaterialHypoElastic( materialProperties, nMaterialProperties, materialLabel ),
      // clang-format off
      // elastic parameters
      stiffnessScaleFactor              ( materialProperties[0] ),
      E1                                ( materialProperties[1] ),
      E2                                ( materialProperties[2] ),
      E3                                ( materialProperties[3] ),
      nu12                              ( materialProperties[4] ),
      nu23                              ( materialProperties[5] ),
      nu13                              ( materialProperties[6] ),
      G12                               ( materialProperties[7] ),
      G23                               ( materialProperties[8] ),
      G13                               ( materialProperties[9] ),
      // viscoelastic parameters
      m                                 ( materialProperties[10] ),
      n                                 ( materialProperties[11] ),
      powerLawApproximationOrder        ( static_cast< size_t > ( materialProperties[12] ) ),
      nKelvin                           ( static_cast< size_t > ( materialProperties[13] ) ),
      minTau                            ( materialProperties[14] ),
      spacing                           ( materialProperties[15] ),
      timeToDays                        ( materialProperties[16] ),
      // material coordinate system
      direction1                        ( { materialProperties[17], materialProperties[18], materialProperties[19] } ),
      direction2                        ( { materialProperties[20], materialProperties[21], materialProperties[22] } )
  // clang-format on
  {
    retardationTimes = KelvinChain::generateRetardationTimes( nKelvin, minTau, spacing );

    auto computeZerothKelvinChainCompliance = [&]( const int order, double tau ) {
      const int& k   = order;
      double     fac = 1;
      for ( int a = 1; a < k; a++ )
        fac *= ( n - a );
      return -pow( -k, k ) / Math::factorial( k - 1 ) * m * fac * pow( k, n - k ) * pow( tau, n );
    };

    // compute moduli for a given approximation order in [2, 3, 4, 7]
    switch ( powerLawApproximationOrder ) {
    case 2: {
      using namespace Marmot::ContinuumMechanics::Viscoelasticity;
      auto phi_2nd  = [&]( autodiff::Real< 2, double > tau ) { return ComplianceFunctions::powerLaw( tau, m, n ); };
      elasticModuli = KelvinChain::computeElasticModuli< 2 >( phi_2nd, retardationTimes );
      zerothKelvinChainCompliance = computeZerothKelvinChainCompliance( 2, minTau / sqrt( spacing ) );
      break;
    }
    case 3: {
      using namespace Marmot::ContinuumMechanics::Viscoelasticity;
      auto phi_3rd  = [&]( autodiff::Real< 3, double > tau ) { return ComplianceFunctions::powerLaw( tau, m, n ); };
      elasticModuli = KelvinChain::computeElasticModuli< 3 >( phi_3rd, retardationTimes );
      zerothKelvinChainCompliance = computeZerothKelvinChainCompliance( 3, minTau / sqrt( spacing ) );
      break;
    }
    case 4: {
      using namespace Marmot::ContinuumMechanics::Viscoelasticity;
      auto phi_4th  = [&]( autodiff::Real< 4, double > tau ) { return ComplianceFunctions::powerLaw( tau, m, n ); };
      elasticModuli = KelvinChain::computeElasticModuli< 4 >( phi_4th, retardationTimes );
      zerothKelvinChainCompliance = computeZerothKelvinChainCompliance( 4, minTau / sqrt( spacing ) );
      break;
    }
    case 7: {
      using namespace Marmot::ContinuumMechanics::Viscoelasticity;
      auto phi_7th  = [&]( autodiff::Real< 7, double > tau ) { return ComplianceFunctions::powerLaw( tau, m, n ); };
      elasticModuli = KelvinChain::computeElasticModuli< 7 >( phi_7th, retardationTimes );
      zerothKelvinChainCompliance = computeZerothKelvinChainCompliance( 7, minTau / sqrt( spacing ) );
      break;
    }

    default: throw std::invalid_argument( "powerLawApproximationOrder must be in [2,3,4,7]" );
    }

    // local normalized stiffness and compliance tensors
    CInv = 1. / stiffnessScaleFactor *
           ContinuumMechanics::Elasticity::Orthotropic::complianceTensor( E1, E2, E3, nu12, nu23, nu13, G12, G23, G13 );
    Cel = stiffnessScaleFactor *
          ContinuumMechanics::Elasticity::Orthotropic::stiffnessTensor( E1, E2, E3, nu12, nu23, nu13, G12, G23, G13 );

    CelUnitInv = E1 * stiffnessScaleFactor * CInv;
    CelUnit    = 1. / E1 / stiffnessScaleFactor * Cel;

    // material coordinate system
    localCoordinateSystem = Marmot::Math::orthonormalCoordinateSystem( direction1, direction2 );

    // unit stiffness matrix in global coordinate system
    CelUnitGlobal = ContinuumMechanics::VoigtNotation::Transformations::
      transformStiffnessToGlobalSystem( CelUnit, localCoordinateSystem );
  }

  void LinearViscoelasticOrthotropicPowerLaw::computeStress( double*       stress,
                                                             double*       dStressDDStrain,
                                                             const double* dStrain,
                                                             const double* timeOld,
                                                             const double  dT,
                                                             double&       pNewDT )

  {
    mVector6d nomStress( stress );
    Vector6d  dE( dStrain );
    mMatrix6d C( dStressDDStrain );

    using namespace Marmot::ContinuumMechanics::VoigtNotation;
    Vector6d dELocal = Transformations::transformStrainToLocalSystem( dE, localCoordinateSystem );

    if ( ( dE.array() == 0 ).all() && dT == 0 ) {
      C = E1 * CelUnitGlobal;
      return;
    }

    Eigen::Ref< KelvinChain::mapStateVarMatrix > creepStateVars( stateVarManager->kelvinStateVars );

    const double dTimeDays = dT * timeToDays;

    Vector6d creepStrainIncrement = Vector6d::Zero();
    double   creepCompliance      = 0;

    // evaluate Kelvin-Chain
    KelvinChain::evaluateKelvinChain( dTimeDays,
                                      elasticModuli,
                                      retardationTimes,
                                      creepStateVars,
                                      creepCompliance,
                                      creepStrainIncrement,
                                      1.0 );

    // compute total 1D effective compliance
    const double effectiveCompliance = 1. / E1 / stiffnessScaleFactor + zerothKelvinChainCompliance + creepCompliance;

    // compute local effective stiffness tensor
    Matrix6d localEffectiveStiffness = 1. / effectiveCompliance * CelUnit;

    // compute stress increment in local coordinate system
    Vector6d deltaStressLocal = localEffectiveStiffness * ( dELocal - creepStrainIncrement );

    // update stress in global coordinate system
    nomStress += Transformations::transformStressToGlobalSystem( deltaStressLocal, localCoordinateSystem );

    // transform local stiffness to global coordinate system
    C = 1. / effectiveCompliance * CelUnitGlobal;

    // update internal state variables
    KelvinChain::updateStateVarMatrix( dTimeDays,
                                       elasticModuli,
                                       retardationTimes,
                                       creepStateVars,
                                       deltaStressLocal,
                                       CelUnitInv );
  }

  void LinearViscoelasticOrthotropicPowerLaw::assignStateVars( double* stateVars_, int nStateVars )
  {
    if ( nStateVars < getNumberOfRequiredStateVars() )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": Not sufficient stateVars!" );

    this->stateVarManager = std::make_unique< LinearViscoelasticOrthotropicPowerLawStateVarManager >( stateVars_,
                                                                                                      nKelvin );

    MarmotMaterial::assignStateVars( stateVars_, nStateVars );
  }

  StateView LinearViscoelasticOrthotropicPowerLaw::getStateView( const std::string& stateName )
  {
    return stateVarManager->getStateView( stateName );
  }

  int LinearViscoelasticOrthotropicPowerLaw::getNumberOfRequiredStateVars()
  {
    return LinearViscoelasticOrthotropicPowerLawStateVarManager::layout.nRequiredStateVars + nKelvin * 6;
  }
} // namespace Marmot::Materials
