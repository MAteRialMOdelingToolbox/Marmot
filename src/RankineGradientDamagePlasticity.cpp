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
    // E, nu, ft, epsilonF, nonlocal length, over-nonlocal parameter m
    constexpr int nRequiredProperties  = 6;
    constexpr int idxDensity           = 6;
    constexpr int idxNonlocalViscosity = 7;
  } // namespace

  RankineGradientDamagePlasticity::RankineGradientDamagePlasticity( const double* materialProperties,
                                                                    int           nMaterialProperties,
                                                                    int           materialLabel )
    : Base( materialProperties, nMaterialProperties, materialLabel ),
      effectiveProperties( [&] {
        if ( nMaterialProperties < nRequiredProperties )
          throw std::runtime_error( MakeString() << __PRETTY_FUNCTION__ << ": Insufficient material properties." );
        if ( !( materialProperties[3] > 0. ) )
          throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": epsilonF must be positive." );
        if ( !( materialProperties[4] >= 0. ) )
          throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": nonlocal length must be >= 0." );
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
    stateLayout.add( "kappaDNonlocal", 1 );
    stateLayout.add( "kappaDWeightedMax", 1 );
    stateLayout.finalize();
  }

  double RankineGradientDamagePlasticity::getDensity( const double* ) const
  {
    if ( nMaterialProperties <= idxDensity )
      throw std::runtime_error( MakeString() << __PRETTY_FUNCTION__ << ": No density given! nMaterialProperties < 7." );
    return materialProperties[idxDensity];
  }

  std::vector< double > RankineGradientDamagePlasticity::getNonlocalViscosity( const double* ) const
  {
    return { nMaterialProperties > idxNonlocalViscosity ? materialProperties[idxNonlocalViscosity] : 0.0 };
  }

  void RankineGradientDamagePlasticity::computeStress( response& res, tangents& tan, const increment& inc ) const
  {
    const double& E              = materialProperties[0];
    const double& nu             = materialProperties[1];
    const double& epsilonF       = materialProperties[3];
    const double& nonlocalLength = materialProperties[4];
    const double& m              = materialProperties[5];

    mVector6d    S( res.stress.data() );
    const auto   dE = Map< const Vector6d >( inc.dStrain.data() );
    const double K  = inc.K( 0 );

    double& kappaD            = stateLayout.getAs< double& >( res.stateVars, "kappaD" );
    double& acYield           = stateLayout.getAs< double& >( res.stateVars, "acYield" );
    double& omega             = stateLayout.getAs< double& >( res.stateVars, "omega" );
    double& kappaDNonlocal    = stateLayout.getAs< double& >( res.stateVars, "kappaDNonlocal" );
    double& kappaDWeightedMax = stateLayout.getAs< double& >( res.stateVars, "kappaDWeightedMax" );

    res.c( 0 ) = nonlocalLength * nonlocalLength;
    tan.dcddK.setZero();
    tan.d2cddK2.setZero();
    tan.dKLocalddK.setZero();

    const Matrix6d Cel    = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );
    const Matrix6d CelInv = ContinuumMechanics::Elasticity::Isotropic::complianceTensor( E, nu );

    // integrity 1 - omega, evaluated from the history variable to retain its relative accuracy for omega -> 1
    const double integrityOld = std::exp( -kappaDWeightedMax / epsilonF );
    if ( !( integrityOld > 0. ) )
      throw StressUpdateFailed( MakeString()
                                << __PRETTY_FUNCTION__ << ": complete damage, effective stress undefined" );

    const Vector6d effStressOld = S / integrityOld;
    const double   psiEffOld    = 0.5 * effStressOld.dot( CelInv * effStressOld );

    // perfectly plastic Rankine return mapping of the effective stress (independent of the damage); its hardening
    // variable is kappaD
    Vector6d effStress     = effStressOld;
    Matrix6d dEffStress_dE = Cel;
    bool     plastic       = false;
    if ( !dE.isZero( 1e-14 ) ) {
      Vector2d                           effectiveStateVars( kappaD, acYield );
      MarmotMaterialHypoElastic::state3D effectiveState( effStressOld, 0., 0., effectiveStateVars.data() );
      effectiveMaterial.computeStress( effectiveState, dEffStress_dE, dE, { .time = inc.time, .dT = inc.dT } );
      effStress = effectiveState.stress;
      plastic   = effectiveStateVars( 1 ) == 1.;
      kappaD    = effectiveStateVars( 0 );
      acYield   = effectiveStateVars( 1 );
    }

    // dKappaD = tr( dE - Cel^-1 dEffStress ), which vanishes for elastic steps
    const Vector6d I          = ( Vector6d() << 1., 1., 1., 0., 0., 0. ).finished();
    const Vector6d dKappaD_dE = plastic ? Vector6d( I - dEffStress_dE.transpose() * ( CelInv * I ) ) : Vector6d::Zero();

    // over-nonlocal damage driving variable and its history maximum
    const double kappaDWeighted     = m * K + ( 1. - m ) * kappaD;
    const bool   loading            = kappaDWeighted > kappaDWeightedMax;
    kappaDWeightedMax               = std::max( kappaDWeightedMax, kappaDWeighted );
    kappaDNonlocal                  = K;
    const double integrity          = std::exp( -kappaDWeightedMax / epsilonF );
    const double dOmega_dKappaDWeig = loading ? integrity / epsilonF : 0.;
    omega                           = 1. - integrity;

    const Vector6d dOmega_dE = dOmega_dKappaDWeig * ( 1. - m ) * dKappaD_dE;
    const double   dOmega_dK = dOmega_dKappaDWeig * m;

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
      const Vector6d dStrainPlastic = dE - CelInv * ( effStress - effStressOld );
      res.dissipation += integrity * effStress.dot( dStrainPlastic );
    }
  }

} // namespace Marmot::Materials
