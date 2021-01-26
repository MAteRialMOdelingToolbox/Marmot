#include "Marmot/LinearElastic.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotUtility.h"
#include "Marmot/MarmotVoigt.h"

namespace Marmot::Materials {

    using namespace Marmot;
    using namespace Eigen;

    LinearElastic::LinearElastic( const double* materialProperties, int nMaterialProperties, int materialNumber )
        : MarmotMaterialHypoElastic::MarmotMaterialHypoElastic( materialProperties,
                                                                nMaterialProperties,
                                                                materialNumber ),
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
        mVector6d             S( stress );
        Map< const Vector6d > dE( dStrain );
        mMatrix6d             C( dStressDDStrain );

        C = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );

        // Zero strain  increment check
        if ( ( dE.array() == 0 ).all() )
            return;

        S = S + C * dE;
    }
} // namespace Marmot::Materials
