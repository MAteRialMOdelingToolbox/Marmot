#include "Marmot/AT2PhaseField.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotPhaseFieldEnergyDegradation.h"
#include "Marmot/MarmotTypedefs.h"

namespace Marmot::Materials {

  using namespace Eigen;
  using namespace Marmot;

  AT2PhaseField::AT2PhaseField( const double* materialProperties, int nMaterialProperties, int materialNumber )
    : MarmotMaterialGeneralGradientEnhancedHypoElastic< 1 >( materialProperties, nMaterialProperties, materialNumber ),
      C( ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( materialProperties[0], materialProperties[1] ) )
  {
    initializeStateLayout();
  }

  double AT2PhaseField::getDensity( const double* stateVars ) const
  {
    if ( nMaterialProperties <= 4 )
      throw std::runtime_error( "Density not provided in material properties for AT2PhaseField." );
    return materialProperties[4];
  }

  std::vector< double > AT2PhaseField::getNonlocalViscosity( const double* stateVars ) const
  {
    if ( nMaterialProperties <= 5 )
      throw std::runtime_error( "Viscosity not provided in material properties for AT2PhaseField." );
    return { materialProperties[5] };
  }

  void AT2PhaseField::computeStress( response& res, tangents& tan, const increment& inc ) const
  {
    // material properties
    const double& Gc = materialProperties[2];
    const double& l  = materialProperties[3];

    // access state variables
    double& H      = stateLayout.getAs< double& >( res.stateVars, "maxCrackDrivingForce" );
    auto    strain = stateLayout.getAs< Map< Vector6d > >( res.stateVars, "strain" );

    // accumulate total strain
    const Vector6d eps = strain + inc.dStrain;

    // phase-field value at the current Gauss point (nonlocal variable K)
    const double& phi = inc.K( 0 );

    // quadratic degradation function: g(phi) = (1-phi)^2
    const auto [g, dg_dphi, d2g_dphi2] = PhaseField::EnergyDegradationFunctions::SecondOrderDerived::quadratic( phi );

    // positive elastic strain energy density (no tension-compression split)
    const double psiPlus = 0.5 * eps.dot( C * eps );

    // irreversibility: update crack driving force history variable
    const bool   loading = psiPlus > H;
    const double H_new   = loading ? psiPlus : H;

    // compute degraded stress
    res.stress = g * ( C * eps );

    // local crack driving force (right-hand side of the phase-field equation)
    //   phi - l^2 * Delta(phi) = KLocal = 2*l/Gc * (1-phi) * H
    res.KLocal( 0 ) = ( 2. * l / Gc ) * ( 1. - phi ) * H_new;

    // gradient enhancement coefficient c = l^2
    res.c( 0 ) = l * l;

    // --- tangent moduli ---

    // d(stress) / d(dStrain) = g(phi) * C
    tan.dStressddStrain = g * C;

    // d(stress) / d(phi) = g'(phi) * C * eps  (stored as 6x1 column)
    tan.dStressddK.col( 0 ) = dg_dphi * ( C * eps );

    // d(KLocal) / d(dStrain):  2l/Gc * (1-phi) * dH/d(eps)
    //   Active loading: dH/d(eps) = d(psi+)/d(eps) = C * eps
    //   Unloading / elastic: dH/d(eps) = 0
    if ( loading )
      tan.dKLocalddStrain.row( 0 ) = ( 2. * l / Gc ) * ( 1. - phi ) * ( C * eps ).transpose();
    else
      tan.dKLocalddStrain.setZero();

    // d(KLocal) / d(phi) = -2*l/Gc * H_new
    tan.dKLocalddK( 0, 0 ) = -( 2. * l / Gc ) * H_new;

    // dcddK and d2cddK2 remain zero (c = l^2 is constant)
    tan.dcddK.setZero();
    tan.d2cddK2.setZero();

    // update state variables
    strain = eps;
    H      = H_new;
  }

} // namespace Marmot::Materials
