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
    if ( this->nMaterialProperties < 6 ) {
      throw std::runtime_error( std::string(
        Utilities::MakeString() << __PRETTY_FUNCTION__ << ": No density given! nMaterialProperties < 6." ) );
    }
    return this->materialProperties[5];
  }

  RankineDamagePlasticity::RankineDamagePlasticity( const double* materialProperties,
                                                    const int     nMaterialProperties,
                                                    const int     materialLabel )
    : MaterialHypoElastic( materialProperties, nMaterialProperties, materialLabel ),
      effectiveProperties( [&] {
        if ( nMaterialProperties < 5 )
          throw std::runtime_error( Utilities::MakeString()
                                    << __PRETTY_FUNCTION__ << ": Insufficient material properties." );
        if ( !( materialProperties[3] > 0. ) )
          throw std::invalid_argument( Utilities::MakeString()
                                       << __PRETTY_FUNCTION__ << ": epsilonF must be positive." );
        if ( !( materialProperties[4] > 0. && materialProperties[4] <= 1. ) )
          throw std::invalid_argument( Utilities::MakeString()
                                       << __PRETTY_FUNCTION__ << ": omegaMax must be in (0, 1]." );
        // E, nu, ft, no hardening
        return std::array< double, 6 >{ materialProperties[0],
                                        materialProperties[1],
                                        materialProperties[2],
                                        0.,
                                        0.,
                                        0. };
      }() ),
      effectiveMaterial( effectiveProperties.data(), static_cast< int >( effectiveProperties.size() ), materialLabel )
  {
    stateLayout.add( "kappaD", 1 );
    stateLayout.add( "acYield", 1 );
    stateLayout.add( "omega", 1 );
    stateLayout.finalize();
  }

  void RankineDamagePlasticity::computeStress( state3D&              state,
                                               Math::Matrix6d&       dStress_dStrain,
                                               const Math::Vector6d& dStrain,
                                               const timeInfo&       timeInfo ) const
  {
    computeStressImplicit( state, dStress_dStrain, dStrain, timeInfo, nullptr, nullptr );
  }

  void RankineDamagePlasticity::computeStressWithEnergyTangents( state3D&              state,
                                                                 Math::Matrix6d&       dStress_dStrain,
                                                                 const Math::Vector6d& dStrain,
                                                                 const timeInfo&       timeInfo,
                                                                 Math::Vector6d&       dElasticEnergyDensity_dStrain,
                                                                 Math::Vector6d&       dDissipation_dStrain ) const
  {
    computeStressImplicit( state,
                           dStress_dStrain,
                           dStrain,
                           timeInfo,
                           &dElasticEnergyDensity_dStrain,
                           &dDissipation_dStrain );
  }

  void RankineDamagePlasticity::computeStressImplicit( state3D&              state,
                                                       Math::Matrix6d&       dStress_dStrain,
                                                       const Math::Vector6d& dStrain,
                                                       const timeInfo&       timeInfo,
                                                       Math::Vector6d*       dPsi_dE,
                                                       Math::Vector6d*       dDissipation_dE ) const
  {
    updateStress( state, dStress_dStrain, dStrain, timeInfo, dPsi_dE, dDissipation_dE );

    // a diverging iterate (e.g., of an arc-length iteration) may overflow the results, e.g., NaN energy tangents from
    // 0 * inf: request a cutback instead of passing non-finite values to the solver
    const bool energyTangentsFinite = !( dPsi_dE && dDissipation_dE ) ||
                                      ( dPsi_dE->allFinite() && dDissipation_dE->allFinite() );
    if ( !( state.stress.allFinite() && dStress_dStrain.allFinite() && std::isfinite( state.elasticEnergyDensity ) &&
            std::isfinite( state.dissipation ) && energyTangentsFinite ) )
      throw Utilities::StressUpdateFailed(
        Utilities::MakeString() << __PRETTY_FUNCTION__ << ": non-finite results, e.g., from a diverging iterate" );
  }

  void RankineDamagePlasticity::updateStress( state3D&              state,
                                              Math::Matrix6d&       dStress_dStrain,
                                              const Math::Vector6d& dStrain,
                                              const timeInfo&       timeInfo,
                                              Math::Vector6d*       dPsi_dE,
                                              Math::Vector6d*       dDissipation_dE ) const
  {
    const bool energyTangents = dPsi_dE && dDissipation_dE;

    const double& E        = this->materialProperties[0];
    const double& nu       = this->materialProperties[1];
    const double& epsilonF = this->materialProperties[3];
    const double& omegaMax = this->materialProperties[4];

    Math::mVector6d S( state.stress.data() );
    Math::mMatrix6d dS_dE( dStress_dStrain.data() );
    const auto      dE = Map< const Math::Vector6d >( dStrain.data() );

    double& kappaD  = stateLayout.getAs< double& >( state.stateVars, "kappaD" );
    double& acYield = stateLayout.getAs< double& >( state.stateVars, "acYield" );
    double& omega   = stateLayout.getAs< double& >( state.stateVars, "omega" );

    const Math::Matrix6d Cel    = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );
    const Math::Matrix6d CelInv = ContinuumMechanics::Elasticity::Isotropic::complianceTensor( E, nu );

    // integrity 1 - omega = 1 - omegaMax + omegaMax exp(-kappaD / (omegaMax epsilonF)), evaluated directly to retain
    // its relative accuracy for omega -> 1
    const double integrityOld = ( 1. - omegaMax ) + omegaMax * std::exp( -kappaD / ( omegaMax * epsilonF ) );
    if ( !( integrityOld > 0. ) )
      throw Utilities::StressUpdateFailed( Utilities::MakeString()
                                           << __PRETTY_FUNCTION__ << ": complete damage, effective stress undefined" );

    const Math::Vector6d effStressOld = S / integrityOld;
    const double         psiEffOld    = 0.5 * effStressOld.dot( CelInv * effStressOld );

    // handle zero strain increment
    if ( dE.isZero( 1e-14 ) ) {
      dS_dE   = integrityOld * Cel;
      acYield = 0;
      if ( energyTangents ) {
        *dPsi_dE = integrityOld * effStressOld;
        dDissipation_dE->setZero();
      }
      return;
    }

    // perfectly plastic Rankine return mapping of the effective stress; its hardening variable is kappaD
    Vector2d       effectiveStateVars( kappaD, acYield );
    state3D        effectiveState( effStressOld, 0., 0., effectiveStateVars.data() );
    Math::Matrix6d dEffStress_dE;
    effectiveMaterial.computeStress( effectiveState, dEffStress_dE, dE, timeInfo );
    const Math::Vector6d& effStress = effectiveState.stress;

    const bool plastic = effectiveStateVars( 1 ) == 1.;
    kappaD             = effectiveStateVars( 0 );
    acYield            = effectiveStateVars( 1 );

    const double decay          = std::exp( -kappaD / ( omegaMax * epsilonF ) );
    const double integrity      = ( 1. - omegaMax ) + omegaMax * decay;
    const double dOmega_dKappaD = decay / epsilonF;
    omega                       = 1. - integrity;
    // complete damage only from a diverging iterate: request a cutback
    if ( !( integrity > 0. ) )
      throw Utilities::StressUpdateFailed( Utilities::MakeString()
                                           << __PRETTY_FUNCTION__ << ": complete damage, kappaD = " << kappaD );

    // dKappaD = tr( dE - Cel^-1 dEffStress ), which vanishes for elastic steps
    const Math::Vector6d I              = ( Math::Vector6d() << 1., 1., 1., 0., 0., 0. ).finished();
    const Math::Vector6d dKappaD_dE     = plastic ? Math::Vector6d( I - dEffStress_dE.transpose() * ( CelInv * I ) )
                                                  : Math::Vector6d::Zero();
    const Math::Vector6d dOmega_dE      = dOmega_dKappaD * dKappaD_dE;
    const Math::Vector6d effElastStrain = CelInv * effStress;
    const double         psiEff         = 0.5 * effStress.dot( effElastStrain );

    S     = integrity * effStress;
    dS_dE = integrity * dEffStress_dE - effStress * dOmega_dE.transpose();

    // free energy (1 - omega) psiEff, dissipation (1 - omega) effStress : dEp + psiEff dOmega
    state.elasticEnergyDensity += integrity * psiEff - integrityOld * psiEffOld;
    if ( energyTangents )
      *dPsi_dE = integrity * dEffStress_dE.transpose() * effElastStrain - psiEff * dOmega_dE;

    if ( !plastic ) {
      if ( energyTangents )
        dDissipation_dE->setZero();
      return;
    }

    const double         dOmega         = integrityOld - integrity;
    const Math::Vector6d dStrainPlastic = dE - CelInv * ( effStress - effStressOld );
    const double         plasticWorkEff = effStress.dot( dStrainPlastic );
    state.dissipation += integrity * plasticWorkEff + psiEff * dOmega;

    if ( energyTangents ) {
      const Math::Matrix6d dEp_dE = Math::Matrix6d::Identity() - CelInv * dEffStress_dE;
      *dDissipation_dE            = -plasticWorkEff * dOmega_dE +
                         integrity * ( dEffStress_dE.transpose() * dStrainPlastic + dEp_dE.transpose() * effStress ) +
                         dOmega * dEffStress_dE.transpose() * effElastStrain + psiEff * dOmega_dE;
    }
  }

} // namespace Marmot::Materials
