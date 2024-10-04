#include "Marmot/MarmotMaterialHypoElasticAD.h"
#include "Marmot/MarmotAutomaticDifferentiation.h"
#include "Marmot/MarmotTypedefs.h"

using namespace Eigen;
using namespace autodiff;

void MarmotMaterialHypoElasticAD::computeStress( double*       stress,
                                                 double*       dStressDDStrain,
                                                 const double* dStrain,
                                                 const double* time,
                                                 const double  dT,
                                                 double&       pNewDT )
{

  using namespace Marmot;
  mVector6d       S( stress );
  const Vector6d  dEps = Map< const Vector6d >( dStrain );
  Map< VectorXd > stateVars( this->stateVars, this->nStateVars );

  // remember old state
  const VectorXd stateVarsOld = stateVars;
  const Vector6d SOld         = S;

  mMatrix6d C( dStressDDStrain );
  // ----------------------------------------
  // autodiff part
  // ----------------------------------------
  // compute stress and tangent with autodiff
  std::tie( S, C ) = Marmot::AutomaticDifferentiation::jacobian(
    [&]( const autodiff::VectorXdual& dE_ ) {
      // reset stateVars to old state
      stateVars = stateVarsOld;

      // set up dual vector for stress
      autodiff::VectorXdual s( SOld );

      // compute stress
      computeStressAD( s.data(), dE_.data(), time, dT, pNewDT );

      return s;
    },
    dEps );
  // ----------------------------------------
}
