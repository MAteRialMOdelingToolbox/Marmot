#include "Marmot/MarmotMaterialHypoElasticAD.h"
#include "Marmot/MarmotAutomaticDifferentiation.h"
#include "Marmot/MarmotTypedefs.h"
#include <autodiff/forward/dual/eigen.hpp>

using namespace Eigen;
using namespace autodiff;

void MarmotMaterialHypoElasticAD::computeStress( state3D&                state,
                                                 Marmot::Math::Matrix6d&       dStressDDStrain,
                                                 const Marmot::Math::Vector6d& dStrain,
                                                 const timeInfo&         timeInfo

) const
{

  using namespace Marmot;
  Marmot::Math::mVector6d       S( state.stress.data() );
  const Marmot::Math::Vector6d  dEps = Map< const Marmot::Math::Vector6d >( dStrain.data() );
  Map< VectorXd > stateVars( state.stateVars, this->getNumberOfRequiredStateVars() );

  // remember old state
  const VectorXd stateVarsOld = stateVars;
  const Marmot::Math::Vector6d SOld         = S;

  Marmot::Math::mMatrix6d C( dStressDDStrain.data() );
  // ----------------------------------------
  // autodiff part
  // ----------------------------------------
  // compute stress and tangent with autodiff
  std::tie( S, C ) = Marmot::Math::AutomaticDifferentiation::dF_dX(
    [&]( const Marmot::Math::Vector6dual dE_ ) {
      // reset stateVars to old state
      stateVars = stateVarsOld;

      Marmot::Math::Vector6dual s( SOld );

      // construct AD state
      state3DAD stateAD = state3DAD( state );
      // compute stress
      computeStressAD( stateAD, dE_, timeInfo );

      Marmot::Math::Vector6dual res = stateAD.stress;
      return res;
    },
    dEps );
  // ----------------------------------------
}
