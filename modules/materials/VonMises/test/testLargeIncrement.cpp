#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include <cmath>
#include <functional>
#include <memory>
#include <string>
#include <vector>

using namespace Marmot::Testing;

namespace {

  using MHE = Marmot::Materials::MaterialHypoElastic;
  using Marmot::Utilities::MakeString;

  /// one single, very large strain increment of order 1 from the virgin state: the trial stress is
  /// rho_trial = ||s_trial|| = 1.7e5 MPa = 887 sqrt(2/3) f_y, whose round-off eps * rho_trial = 3.8e-11 MPa lies above
  /// the absolute tolerance 1e-12 the return mapping used before; with the tolerance relative to ||s_trial|| the
  /// return mapping converges, and the stress lies on the yield surface of the updated hardening variable,
  /// sigma_eq = f_y(kappa)
  Marmot::Math::Vector6d largeStrainIncrement()
  {
    Marmot::Math::Vector6d dStrain;
    dStrain << 1.0, -0.3, 0.2, 0.5, -0.4, 0.3;
    return dStrain;
  }

  struct Hardening {
    std::string name;
    double      HLin;
    double      deltaYieldStress;
    double      delta;
  };

  /// perfect plasticity and linear hardening: g(dKappa) is linear, Newton converges in one step; nonlinear hardening:
  /// several iterations under the relative criterion
  const std::vector< Hardening > hardeningLaws = { { "perfectly plastic", 0., 0., 0. },
                                                   { "linear hardening", 2100., 0., 0. },
                                                   { "nonlinear hardening", 2100., 100., 20. } };

  void checkOnYieldSurface( const std::string&            material,
                            const Hardening&              h,
                            bool                          explicitReturnMapping,
                            const Marmot::Math::Vector6d& dStrain )
  {
    const double fy0 = 235.;
    //                                         E,        nu,  fy,  HLin,   deltaYieldStress,   delta,   density
    const std::vector< double > properties = { 210000., 0.3, fy0, h.HLin, h.deltaYieldStress, h.delta, 7.85e-9 };

    auto mat = std::unique_ptr< MHE >(
      Marmot::Registry::MaterialHypoElasticFactory::createMaterial( material,
                                                                    properties.data(),
                                                                    static_cast< int >( properties.size() ),
                                                                    1 ) );
    std::vector< double > stateVars( mat->getNumberOfRequiredStateVars(), 0.0 );
    mat->initializeYourself( stateVars.data(), static_cast< int >( stateVars.size() ) );

    MHE::state3D state;
    state.stateVars = stateVars.data();
    if ( explicitReturnMapping )
      mat->computeStressExplicit( state, dStrain, { 0.0, 1.0 } );
    else {
      Marmot::Math::Matrix6d dStress_dStrain;
      mat->computeStress( state, dStress_dStrain, dStrain, { 0.0, 1.0 } );
    }

    const double kappa = mat->getStateView( "kappa", stateVars.data() ).stateLocation[0];
    const double fy    = fy0 + h.HLin * kappa + h.deltaYieldStress * ( 1. - std::exp( -h.delta * kappa ) );
    // the deviator first: J2 = I1^2 / 3 - I2 of the full stress cancels the hydrostatic part (here ~ 1e3 f_y) and
    // loses ~ 1e-9 relative accuracy; J2 of the deviator (I1 = 0) does not
    const Marmot::Math::Vector6d s       = Marmot::ContinuumMechanics::Voigt::IDev * state.stress;
    const double                 sigmaEq = Marmot::ContinuumMechanics::Invariants::vonMisesEquivalentStress( s );

    throwExceptionOnFailure( kappa > 0 && std::abs( sigmaEq / fy - 1 ) < 1e-10,
                             MakeString()
                               << material << ", " << h.name << ( explicitReturnMapping ? ", explicit" : "" )
                               << ": sigma_eq / f_y(kappa) - 1 = " << sigmaEq / fy - 1 << ", kappa = " << kappa );
  }

  void testLargeIncrement()
  {
    for ( const auto& h : hardeningLaws ) {
      checkOnYieldSurface( "VONMISES", h, false, largeStrainIncrement() );
      checkOnYieldSurface( "VONMISES", h, true, largeStrainIncrement() );
      checkOnYieldSurface( "ADVONMISES", h, false, largeStrainIncrement() );
    }
  }

} // namespace

int main()
{
  std::vector< std::function< void( void ) > > tests = { testLargeIncrement };

  executeTestsAndCollectExceptions( tests );
  return 0;
}
