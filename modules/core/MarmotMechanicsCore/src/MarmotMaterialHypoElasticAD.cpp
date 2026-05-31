#include "Marmot/MarmotMaterialHypoElasticAD.h"
#include "Marmot/MarmotAutomaticDifferentiation.h"
#include "Marmot/MarmotTypedefs.h"
#include <autodiff/forward/dual/eigen.hpp>

using namespace Eigen;
using namespace autodiff;

void MarmotMaterialHypoElasticAD::computeStress( state3D&                state,
                                                 Marmot::Matrix6d&       dStressDDStrain,
                                                 const Marmot::Vector6d& dStrain,
                                                 const timeInfo&         timeInfo

) const
{

  using namespace Marmot;
  mVector6d       S( state.stress.data() );
  const Vector6d  dEps = Map< const Vector6d >( dStrain.data() );
  Map< VectorXd > stateVars( state.stateVars, this->getNumberOfRequiredStateVars() );

  // remember old state
  const VectorXd stateVarsOld = stateVars;
  const Vector6d SOld         = S;

  mMatrix6d C( dStressDDStrain.data() );
  // ----------------------------------------
  // autodiff part
  // ----------------------------------------
  // compute stress and tangent with autodiff
  std::tie( S, C ) = Marmot::AutomaticDifferentiation::dF_dX(
    [&]( const Marmot::Vector6dual dE_ ) {
      // reset stateVars to old state
      stateVars = stateVarsOld;

      Marmot::Vector6dual s( SOld );

      // construct AD state
      state3DAD stateAD = state3DAD( state );
      // compute stress
      computeStressAD( stateAD, dE_, timeInfo );

      Marmot::Vector6dual res = stateAD.stress;
      return res;
    },
    dEps );
  // ----------------------------------------
}
