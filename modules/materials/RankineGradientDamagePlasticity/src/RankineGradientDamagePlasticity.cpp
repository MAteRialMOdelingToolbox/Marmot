#include "Marmot/RankineGradientDamagePlasticity.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTypedefs.h"
#include <algorithm>
#include <cmath>

namespace Marmot::Materials {

  using namespace Eigen;
  using namespace Marmot;

  namespace {
    // E, nu, ft, epsilonF, damage threshold omegaMax, nonlocal length, over-nonlocal parameter m
    constexpr int nRequiredProperties  = 7;
    constexpr int idxDensity           = 7;
    constexpr int idxNonlocalViscosity = 8;
  } // namespace

  RankineGradientDamagePlasticity::RankineGradientDamagePlasticity( const double* materialProperties,
                                                                    int           nMaterialProperties,
                                                                    int           materialLabel )
    : Base( materialProperties, nMaterialProperties, materialLabel ),
      effectiveProperties( [&] {
        if ( nMaterialProperties < nRequiredProperties )
          throw std::runtime_error( Utilities::MakeString()
                                    << __PRETTY_FUNCTION__ << ": Insufficient material properties." );
        if ( !( materialProperties[3] > 0. ) )
          throw std::invalid_argument( Utilities::MakeString()
                                       << __PRETTY_FUNCTION__ << ": epsilonF must be positive." );
        if ( !( materialProperties[4] > 0. && materialProperties[4] <= 1. ) )
          throw std::invalid_argument( Utilities::MakeString()
                                       << __PRETTY_FUNCTION__ << ": omegaMax must be in (0, 1]." );
        if ( !( materialProperties[5] >= 0. ) )
          throw std::invalid_argument( Utilities::MakeString()
                                       << __PRETTY_FUNCTION__ << ": nonlocal length must be >= 0." );
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
    stateLayout.add( "kappaDNonlocal", 1 );
    stateLayout.add( "kappaDWeightedMax", 1 );
    stateLayout.finalize();
  }

  double RankineGradientDamagePlasticity::getDensity( const double* ) const
  {
    if ( nMaterialProperties <= idxDensity )
      throw std::runtime_error( Utilities::MakeString()
                                << __PRETTY_FUNCTION__ << ": No density given! nMaterialProperties < 8." );
    return materialProperties[idxDensity];
  }

  std::vector< double > RankineGradientDamagePlasticity::getNonlocalViscosity( const double* ) const
  {
    return { nMaterialProperties > idxNonlocalViscosity ? materialProperties[idxNonlocalViscosity] : 0.0 };
  }

  void RankineGradientDamagePlasticity::computeStress( response& res, tangents& tan, const increment& inc ) const
  {
    // a non-finite nonlocal field would not enter the checked results, since the history maximum ignores it
    if ( !inc.K.allFinite() )
      throw Utilities::StressUpdateFailed( Utilities::MakeString()
                                           << __PRETTY_FUNCTION__ << ": non-finite nonlocal field" );

    updateStress( res, tan, inc );

    // a diverging iterate may overflow the results (e.g., the energies to inf or 0 * inf = NaN): request a cutback
    // instead of passing non-finite values to the solver
    if ( !( res.stress.allFinite() && res.KLocal.allFinite() && std::isfinite( res.elasticEnergyDensity ) &&
            std::isfinite( res.dissipation ) && tan.dStressddStrain.allFinite() && tan.dStressddK.allFinite() &&
            tan.dKLocalddStrain.allFinite() ) )
      throw Utilities::StressUpdateFailed(
        Utilities::MakeString() << __PRETTY_FUNCTION__ << ": non-finite results, e.g., from a diverging iterate" );
  }

  void RankineGradientDamagePlasticity::updateStress( response& res, tangents& tan, const increment& inc ) const
  {
    const double& E              = materialProperties[0];
    const double& nu             = materialProperties[1];
    const double& epsilonF       = materialProperties[3];
    const double& omegaMax       = materialProperties[4];
    const double& nonlocalLength = materialProperties[5];
    const double& m              = materialProperties[6];

    Math::mVector6d S( res.stress.data() );
    const auto      dE = Map< const Math::Vector6d >( inc.dStrain.data() );
    const double    K  = inc.K( 0 );

    double& kappaD            = stateLayout.getAs< double& >( res.stateVars, "kappaD" );
    double& acYield           = stateLayout.getAs< double& >( res.stateVars, "acYield" );
    double& omega             = stateLayout.getAs< double& >( res.stateVars, "omega" );
    double& kappaDNonlocal    = stateLayout.getAs< double& >( res.stateVars, "kappaDNonlocal" );
    double& kappaDWeightedMax = stateLayout.getAs< double& >( res.stateVars, "kappaDWeightedMax" );

    res.c( 0 ) = nonlocalLength * nonlocalLength;
    tan.dcddK.setZero();
    tan.d2cddK2.setZero();
    tan.dKLocalddK.setZero();

    const Math::Matrix6d Cel    = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );
    const Math::Matrix6d CelInv = ContinuumMechanics::Elasticity::Isotropic::complianceTensor( E, nu );

    // integrity 1 - omega = 1 - omegaMax + omegaMax exp(-kappa / (omegaMax epsilonF)), evaluated from the history
    // variable to retain its relative accuracy for omega -> 1
    const double integrityOld = ( 1. - omegaMax ) + omegaMax * std::exp( -kappaDWeightedMax / ( omegaMax * epsilonF ) );
    if ( !( integrityOld > 0. ) )
      throw Utilities::StressUpdateFailed( Utilities::MakeString()
                                           << __PRETTY_FUNCTION__ << ": complete damage, effective stress undefined" );

    const Math::Vector6d effStressOld = S / integrityOld;
    const double         psiEffOld    = 0.5 * effStressOld.dot( CelInv * effStressOld );

    // perfectly plastic Rankine return mapping of the effective stress (independent of the damage); its hardening
    // variable is kappaD
    Math::Vector6d effStress     = effStressOld;
    Math::Matrix6d dEffStress_dE = Cel;
    bool           plastic       = false;
    acYield                      = 0;
    if ( !dE.isZero( 1e-14 ) ) {
      Vector2d                     effectiveStateVars( kappaD, acYield );
      MaterialHypoElastic::state3D effectiveState( effStressOld, 0., 0., effectiveStateVars.data() );
      effectiveMaterial.computeStress( effectiveState, dEffStress_dE, dE, { .time = inc.time, .dT = inc.dT } );
      effStress = effectiveState.stress;
      plastic   = effectiveStateVars( 1 ) == 1.;
      kappaD    = effectiveStateVars( 0 );
      acYield   = effectiveStateVars( 1 );
    }

    // dKappaD = tr( dE - Cel^-1 dEffStress ), which vanishes for elastic steps
    const Math::Vector6d I          = ( Math::Vector6d() << 1., 1., 1., 0., 0., 0. ).finished();
    const Math::Vector6d dKappaD_dE = plastic ? Math::Vector6d( I - dEffStress_dE.transpose() * ( CelInv * I ) )
                                              : Math::Vector6d::Zero();

    // over-nonlocal damage driving variable and its history maximum
    const double kappaDWeighted     = m * K + ( 1. - m ) * kappaD;
    const bool   loading            = kappaDWeighted > kappaDWeightedMax;
    kappaDWeightedMax               = std::max( kappaDWeightedMax, kappaDWeighted );
    kappaDNonlocal                  = K;
    const double decay              = std::exp( -kappaDWeightedMax / ( omegaMax * epsilonF ) );
    const double integrity          = ( 1. - omegaMax ) + omegaMax * decay;
    const double dOmega_dKappaDWeig = loading ? decay / epsilonF : 0.;
    omega                           = 1. - integrity;
    // complete damage only from a diverging iterate (e.g., at a snap-through): request a cutback
    if ( !( integrity > 0. ) )
      throw Utilities::StressUpdateFailed(
        Utilities::MakeString() << __PRETTY_FUNCTION__ << ": complete damage, kappaDWeighted = " << kappaDWeighted );

    const Math::Vector6d dOmega_dE = dOmega_dKappaDWeig * ( 1. - m ) * dKappaD_dE;
    const double         dOmega_dK = dOmega_dKappaDWeig * m;

    S                   = integrity * effStress;
    tan.dStressddStrain = integrity * dEffStress_dE - effStress * dOmega_dE.transpose();
    tan.dStressddK      = -effStress * dOmega_dK;

    res.KLocal( 0 )              = kappaD;
    tan.dKLocalddStrain.row( 0 ) = dKappaD_dE.transpose();

    // free energy (1 - omega) psiEff, dissipation (1 - omega) effStress : dEp + psiEff dOmega
    const double psiEff = 0.5 * effStress.dot( CelInv * effStress );
    res.elasticEnergyDensity += integrity * psiEff - integrityOld * psiEffOld;
    res.dissipation += psiEff * ( integrityOld - integrity );
    if ( plastic ) {
      const Math::Vector6d dStrainPlastic = dE - CelInv * ( effStress - effStressOld );
      res.dissipation += integrity * effStress.dot( dStrainPlastic );
    }
  }

} // namespace Marmot::Materials
