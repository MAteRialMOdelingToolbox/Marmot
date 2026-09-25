#include "Marmot/GradientEnhancedCompressibleNeoHookeDamage.h"
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotEnergyDensityFunctions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotStressMeasures.h"
#include "Marmot/MarmotUtils.h"
#include <cmath>
#include <stdexcept>

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Fastor;
  using namespace FastorIndices;
  using namespace FastorStandardTensors;

  GradientEnhancedCompressibleNeoHookeDamage::GradientEnhancedCompressibleNeoHookeDamage(
    const double* materialProperties,
    int           nMaterialProperties,
    int           materialNumber )
    : MarmotMaterialGradientEnhancedFiniteStrain( materialProperties, nMaterialProperties, materialNumber ),
      K( materialProperties[0] ),
      G( materialProperties[1] ),
      kappa0( materialProperties[2] ),
      kappaF( materialProperties[3] ),
      nonLocalRadius( materialProperties[4] )
  {
    if ( nMaterialProperties < 5 )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__
                                                << ": expected at least 5 material properties (K, G, kappa0, "
                                                   "kappaF, l), got "
                                                << nMaterialProperties );
    if ( !( kappaF > kappa0 ) )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": kappaF must be greater than kappa0" );

    stateLayout.add( "kappa", 1 );
    stateLayout.finalize();
  }

  double GradientEnhancedCompressibleNeoHookeDamage::getDensity( const double* stateVars ) const
  {
    if ( nMaterialProperties <= 5 )
      throw std::runtime_error( MakeString() << __PRETTY_FUNCTION__ << ": density not provided (property 5)" );
    return materialProperties[5];
  }

  std::pair< double, double > GradientEnhancedCompressibleNeoHookeDamage::damage( double kappa ) const
  {
    if ( kappa <= kappa0 )
      return { 0.0, 0.0 };

    const double e = std::exp( -( kappa - kappa0 ) / ( kappaF - kappa0 ) );

    const double D      = 1.0 - kappa0 / kappa * e;
    const double dD_dKa = kappa0 / ( kappa * kappa ) * e + kappa0 / kappa * e / ( kappaF - kappa0 );

    return { D, dD_dKa };
  }

  void GradientEnhancedCompressibleNeoHookeDamage::computeStress( ConstitutiveResponse< 3 >& response,
                                                                  AlgorithmicModuli< 3 >&    tangents,
                                                                  const Deformation< 3 >&    deformation,
                                                                  const TimeIncrement&       timeIncrement ) const
  {
    using namespace ContinuumMechanics;

    const auto& F = deformation.F;

    // undamaged neo-Hookean response, as in CompressibleNeoHooke
    const auto [C, dC_dF]                  = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( F );
    const auto [psi0, dPsi_dC, d2Psi_dCdC] = EnergyDensityFunctions::SecondOrderDerived::PenceGouPotentialB( C, K, G );

    const Tensor33d PK2                     = 2. * dPsi_dC;
    const auto [tau0, dTau_dPK2, dTau0_dFe] = StressMeasures::FirstOrderDerived::KirchhoffStressFromPK2( PK2, F );
    const Tensor3333d dTau0_dF = 2.0 * einsum< ijKL, KLMN >( einsum< ijKL, IJKL >( dTau_dPK2, d2Psi_dCdC ), dC_dF ) +
                                 dTau0_dFe;

    // local driving force: energy-equivalent strain
    const double E = 9. * K * G / ( 3. * K + G );
    const double L = std::sqrt( 2. * std::max( psi0, 0.0 ) / E );

    // damage from the history of the nonlocal field
    double&      kappa    = stateLayout.getAs< double& >( response.stateVars, "kappa" );
    const double kappaOld = std::max( kappa, kappa0 );
    const bool   loading  = deformation.N > kappaOld;
    const double kappaNew = loading ? deformation.N : kappaOld;

    const auto [D, dD_dKappa] = damage( kappaNew );

    response.tau                  = ( 1. - D ) * tau0;
    response.L                    = L;
    response.nonLocalRadius       = nonLocalRadius;
    response.elasticEnergyDensity = ( 1. - D ) * psi0;
    response.dissipation          = 0.0;

    tangents.dTau_dF = ( 1. - D ) * dTau0_dF;
    tangents.dTau_dN = loading ? Tensor33d( -dD_dKappa * tau0 ) : Tensor33d( 0.0 );

    // dL/dF = (dpsi0/dF) / (E L), with dpsi0/dF = P0 = F S0; zero in the undeformed state, where L has a kink
    if ( L > 1e-14 ) {
      const Tensor33d P0 = einsum< iK, KJ >( F, PK2 );
      tangents.dL_dF     = ( 1. / ( E * L ) ) * P0;
    }
    else
      tangents.dL_dF = 0.0;

    tangents.dL_dN = 0.0;

    kappa = kappaNew;
  }

} // namespace Marmot::Materials
