#include "Marmot/RankineDamagePlasticity.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTypedefs.h"

namespace Marmot::Materials {

  using namespace Eigen;
  using namespace Marmot;

  double RankineDamagePlasticity::getDensity( const double* stateVars ) const
  {
    if ( this->nMaterialProperties < 5 ) {
      throw std::runtime_error(
        std::string( MakeString() << __PRETTY_FUNCTION__ << ": No density given! nMaterialProperties < 5." ) );
    }
    return this->materialProperties[4];
  }

  RankineDamagePlasticity::RankineDamagePlasticity( const double* materialProperties,
                                                    const int     nMaterialProperties,
                                                    const int     materialLabel )
    : MarmotMaterialHypoElastic( materialProperties, nMaterialProperties, materialLabel ),
      effectiveProperties( [&] {
        if ( nMaterialProperties < 4 )
          throw std::runtime_error( MakeString() << __PRETTY_FUNCTION__ << ": Insufficient material properties." );
        if ( !( materialProperties[3] > 0. ) )
          throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": epsilonF must be positive." );
        // E, nu, ft, no hardening
        return std::array< double, 6 >{ materialProperties[0],
                                        materialProperties[1],
                                        materialProperties[2],
                                        0.,
                                        0.,
                                        0. };
      }() ),
      effectiveMaterial( effectiveProperties.data(), effectiveProperties.size(), materialLabel )
  {
    stateLayout.add( "kappaD", 1 );
    stateLayout.add( "acYield", 1 );
    stateLayout.add( "omega", 1 );
    stateLayout.finalize();
  }

  void RankineDamagePlasticity::computeStress( state3D&        state,
                                               Matrix6d&       dStress_dStrain,
                                               const Vector6d& dStrain,
                                               const timeInfo& timeInfo ) const
  {
    const double& E        = this->materialProperties[0];
    const double& nu       = this->materialProperties[1];
    const double& epsilonF = this->materialProperties[3];

    mVector6d  S( state.stress.data() );
    mMatrix6d  dS_dE( dStress_dStrain.data() );
    const auto dE = Map< const Vector6d >( dStrain.data() );

    double& kappaD  = stateLayout.getAs< double& >( state.stateVars, "kappaD" );
    double& acYield = stateLayout.getAs< double& >( state.stateVars, "acYield" );
    double& omega   = stateLayout.getAs< double& >( state.stateVars, "omega" );

    const Matrix6d Cel    = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );
    const Matrix6d CelInv = ContinuumMechanics::Elasticity::Isotropic::complianceTensor( E, nu );

    // integrity 1 - omega, evaluated directly to retain its relative accuracy for omega -> 1
    const double integrityOld = std::exp( -kappaD / epsilonF );
    if ( !( integrityOld > 0. ) )
      throw StressUpdateFailed( MakeString()
                                << __PRETTY_FUNCTION__ << ": complete damage, effective stress undefined" );

    const Vector6d effStressOld = S / integrityOld;
    const double   psiEffOld    = 0.5 * effStressOld.dot( CelInv * effStressOld );

    // handle zero strain increment
    if ( dE.isZero( 1e-14 ) ) {
      dS_dE                               = integrityOld * Cel;
      state.dElasticEnergyDensity_dStrain = integrityOld * effStressOld;
      state.dDissipation_dStrain.setZero();
      return;
    }

    // perfectly plastic Rankine return mapping of the effective stress; its hardening variable is kappaD
    Vector2d effectiveStateVars( kappaD, acYield );
    state3D  effectiveState( effStressOld, 0., 0., effectiveStateVars.data() );
    Matrix6d dEffStress_dE;
    effectiveMaterial.computeStress( effectiveState, dEffStress_dE, dE, timeInfo );
    const Vector6d& effStress = effectiveState.stress;

    const bool plastic = effectiveStateVars( 1 ) == 1.;
    kappaD             = effectiveStateVars( 0 );
    acYield            = effectiveStateVars( 1 );

    const double integrity      = std::exp( -kappaD / epsilonF );
    const double dOmega_dKappaD = integrity / epsilonF;
    omega                       = 1. - integrity;

    // dKappaD = tr( dE - Cel^-1 dEffStress ), which vanishes for elastic steps
    const Vector6d I          = ( Vector6d() << 1., 1., 1., 0., 0., 0. ).finished();
    const Vector6d dKappaD_dE = plastic ? Vector6d( I - dEffStress_dE.transpose() * ( CelInv * I ) ) : Vector6d::Zero();
    const Vector6d dOmega_dE  = dOmega_dKappaD * dKappaD_dE;
    const Vector6d effElastStrain = CelInv * effStress;
    const double   psiEff         = 0.5 * effStress.dot( effElastStrain );

    S     = integrity * effStress;
    dS_dE = integrity * dEffStress_dE - effStress * dOmega_dE.transpose();

    // free energy (1 - omega) psiEff, dissipation (1 - omega) effStress : dEp + psiEff dOmega
    state.elasticEnergyDensity += integrity * psiEff - integrityOld * psiEffOld;
    state.dElasticEnergyDensity_dStrain = integrity * dEffStress_dE.transpose() * effElastStrain - psiEff * dOmega_dE;

    if ( !plastic ) {
      state.dDissipation_dStrain.setZero();
      return;
    }

    const double   dOmega         = integrityOld - integrity;
    const Vector6d dStrainPlastic = dE - CelInv * ( effStress - effStressOld );
    const Matrix6d dEp_dE         = Matrix6d::Identity() - CelInv * dEffStress_dE;
    const double   plasticWorkEff = effStress.dot( dStrainPlastic );
    state.dissipation += integrity * plasticWorkEff + psiEff * dOmega;
    state.dDissipation_dStrain = -plasticWorkEff * dOmega_dE +
                                 integrity *
                                   ( dEffStress_dE.transpose() * dStrainPlastic + dEp_dE.transpose() * effStress ) +
                                 dOmega * dEffStress_dE.transpose() * effElastStrain + psiEff * dOmega_dE;
  }

} // namespace Marmot::Materials
