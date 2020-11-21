#include "Marmot/LinearElastic.h"
#include "Marmot/MarmotAbaqusUtility.h"
#include "Marmot/MarmotFunctions.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

using namespace Marmot;
using namespace Eigen;

LinearElastic::LinearElastic( const double* materialProperties, int nMaterialProperties, int materialNumber )
    : MarmotMaterialHypoElastic::MarmotMaterialHypoElastic( materialProperties, nMaterialProperties, materialNumber ),
      E( materialProperties[0] ),
      nu( materialProperties[1] )
{
    assert( nMaterialProperties >= 2 );
}

void LinearElastic::computeStress( double*       stress,
                                   double*       dStressDDStrain,
                                   const double* dStrain,
                                   const double* timeOld,
                                   const double  dT,
                                   double&       pNewDT )
{
    mVector6             S( stress );
    Map< const Vector6 > dE( dStrain );
    mMatrix6             C( dStressDDStrain );

    C = mechanics::Cel( E, nu );

    // Zero strain  increment check
    if ( ( dE.array() == 0 ).all() )
        return;

    S = S + C * dE;
}
