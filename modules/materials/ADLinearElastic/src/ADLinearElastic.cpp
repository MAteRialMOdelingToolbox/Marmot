#include "Marmot/ADLinearElastic.h"
#include "Marmot/MarmotAutomaticDifferentiationTypes.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotTypedefs.h"

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Marmot::Math::Differentiation::Automatic;
  using namespace Eigen;
  using namespace ContinuumMechanics::Elasticity;

  ADLinearElastic::ADLinearElastic( const double* materialProperties, int nMaterialProperties, int materialNumber )
    : Marmot::Materials::MaterialHypoElasticAD::MaterialHypoElasticAD( materialProperties,
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
    throw std::runtime_error( std::string( Marmot::Utilities::MakeString()
                                           << __PRETTY_FUNCTION__ << ": Density not specified for ADLinearElastic." ) );
  }

  void ADLinearElastic::computeStressAD( state3DAD& state, const Vector6AD& dStrain, const timeInfo& timeInfo ) const
  {

    const Vector6AD dStress = C * dStrain;

    state.elasticEnergyDensity += Math::Tensors::makeReal( dStrain.dot( state.stress + 0.5 * dStress ) );
    state.stress += dStress;
    state.dissipation = 0.0;
  }
} // namespace Marmot::Materials
