#include "LinearElastic.h"
#include "MarmotAbaqusUtility.h"
#include "MarmotFunctions.h"
#include "MarmotTypedefs.h"
#include "MarmotVoigt.h"

using namespace marmot;
using namespace Eigen;

void LinearElastic::computeStress( double* stress,
                                   double* dStressDDStrain,
                                   const double* dStrain,
                                   const double* timeOld,
                                   const double  dT,
                                   double&       pNewDT )
{
    // material properties
    const double& E  = materialProperties[0];
    const double& nu = materialProperties[1];

    mVector6           S( stress );
    Map<const Vector6> dE( dStrain );
    mMatrix6           C( dStressDDStrain );

    C = mechanics::Cel( E, nu );

    // Zero strain  increment check
    if ( ( dE.array() == 0 ).all() )
        return;

    S = S + C * dE;
}
