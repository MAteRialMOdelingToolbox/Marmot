#include "LinearElastic.h"
#include "bftAbaqusUtility.h"
#include "bftFunctions.h"
#include "bftTypedefs.h"
#include "bftVoigt.h"

using namespace bft;
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

    const int nTensor = 6;

    mVector6           S( stress );
    Map<const Vector6> dE( dStrain );
    mMatrix6           C( dStressDDStrain );

    C = mechanics::Cel( E, nu );

    // Zero strain  increment check
    if ( ( dE.array() == 0 ).all() )
        return;

    S = S + C * dE;
}
