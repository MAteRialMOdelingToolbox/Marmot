#include "Marmot/MarmotMaterialHypoElasticAD.h"
#include "Marmot/MarmotAutomaticDifferentiation.h"
#include "Marmot/MarmotTypedefs.h"

using namespace Eigen;
using namespace autodiff;

void MarmotMaterialHypoElasticAD::computeStress( state3D&        state,
                                                 double*         dStressDDStrain,
                                                 const double*   dStrain,
                                                 const timeInfo& timeInfo

) const
{

  using namespace Marmot;
  mVector6d       S( state.stress.data() );
  const Vector6d  dEps = Map< const Vector6d >( dStrain );
  Map< VectorXd > stateVars( state.stateVars, this->getNumberOfRequiredStateVars() );

  // remember old state
  const VectorXd stateVarsOld = stateVars;
  const Vector6d SOld         = S;

  mMatrix6d C( dStressDDStrain );
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
      state3DAD stateAD;
      stateAD.stress              = s.data();
      stateAD.strainEnergyDensity = state.strainEnergyDensity;
      stateAD.stateVars           = stateVars.data();

      // compute stress
      computeStressAD( stateAD, dE_.data(), timeInfo );

      Marmot::Vector6dual stressAD( stateAD.stress );

      return stressAD;
    },
    dEps );
  // ----------------------------------------
}
