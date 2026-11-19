#include "Marmot/MarmotMaterialHypoElasticAD.h"
#include "Marmot/MarmotAutomaticDifferentiation.h"
#include "Marmot/MarmotAutomaticDifferentiationTypes.h"
#include "Marmot/MarmotTypedefs.h"

using namespace Eigen;
using namespace Marmot::Math::Differentiation::Automatic;

void Marmot::Materials::MaterialHypoElasticAD::computeStress( state3D&                      state,
                                                              Marmot::Math::Matrix6d&       dStressDDStrain,
                                                              const Marmot::Math::Vector6d& dStrain,
                                                              const timeInfo&               timeInfo

) const
{

  using namespace Marmot;
  Marmot::Math::mVector6d      S( state.stress.data() );
  const Marmot::Math::Vector6d dEps = Map< const Marmot::Math::Vector6d >( dStrain.data() );
  Map< VectorXd >              stateVars( state.stateVars, this->getNumberOfRequiredStateVars() );

  // remember old state
  const VectorXd stateVarsOld = stateVars;

  double elasticEnergyDensity;
  double dissipation;

  Marmot::Math::mMatrix6d C( dStressDDStrain.data() );
  // ----------------------------------------
  // autodiff part
  // ----------------------------------------
  // compute stress and tangent with autodiff
  std::tie( S, C ) = dF_dX< 6 >(
    [&]( const Vector6AD& dE_ ) {
      // reset stateVars to old state
      stateVars = stateVarsOld;

      // construct AD state
      state3DAD stateAD = state3DAD( state );
      // compute stress
      computeStressAD( stateAD, dE_, timeInfo );

      elasticEnergyDensity = stateAD.elasticEnergyDensity;
      dissipation          = stateAD.dissipation;

      Vector6AD res = stateAD.stress;
      return res;
    },
    dEps );
  // ----------------------------------------

  state.elasticEnergyDensity = elasticEnergyDensity;
  state.dissipation          = dissipation;
}
