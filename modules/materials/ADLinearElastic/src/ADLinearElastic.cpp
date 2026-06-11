#include "Marmot/ADLinearElastic.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotTypedefs.h"
#include "autodiff/forward/dual/eigen.hpp"

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace autodiff;
  using namespace Eigen;
  using namespace ContinuumMechanics::Elasticity;

  ADLinearElastic::ADLinearElastic( const double* materialProperties, int nMaterialProperties, int materialNumber )
    : MarmotMaterialHypoElasticAD::MarmotMaterialHypoElasticAD( materialProperties,
                                                                nMaterialProperties,
                                                                materialNumber ),
      C( Isotropic::stiffnessTensor( materialProperties[0], materialProperties[1] ) )
  {
    stateLayout.finalize();
  }

  double ADLinearElastic::getDensity( const double* stateVars ) const
  {
    if ( nMaterialProperties >= 3 )
      return materialProperties[2];
    throw std::runtime_error(
      std::string( MakeString() << __PRETTY_FUNCTION__ << ": Density not specified for ADLinearElastic." ) );
  }

  void ADLinearElastic::computeStressAD( state3DAD&                 state,
                                         const Marmot::Vector6dual& dStrain,
                                         const timeInfo&            timeInfo ) const
  {

    const Vector6dual dStress = C * dStrain;
    state.stress += dStress;
    state.elasticEnergyDensity += 0.5 * Math::makeReal( dStrain.dot( dStress ) );
    state.dissipation = 0.0;
  }
} // namespace Marmot::Materials
