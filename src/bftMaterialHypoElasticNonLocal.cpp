#include "bftMaterialHypoElasticNonLocal.h"
#include "bftFunctions.h"
#include "bftVoigt.h"

using namespace Eigen;

void BftMaterialHypoElasticNonLocal::computePlaneStress( double*       stress_,
                                                         double&       K_local,
                                                         double&       nonLocalRadius,
                                                         double*       dStressDDStrain_,
                                                         double*       dK_localDDStrain,
                                                         double*       dStressDK,
                                                         const double* strainOld,
                                                         double*       dStrain_,
                                                         double        KOld,
                                                         double        dK,
                                                         const double* timeOld,
                                                         const double  dT,
                                                         double&       pNewDT )
{
    using namespace bft;

    Map<Vector6>  stress( stress_ );
    Map<Matrix6>  dStressDDStrain( dStressDDStrain_ );
    Map<Vector6>  dStrain( dStrain_ );
    Map<VectorXd> stateVars( this->stateVars, this->nStateVars );

    Vector6  stressTemp;
    VectorXd stateVarsOld = stateVars;
    Vector6  dStrainTemp  = dStrain;

    // assumption of isochoric deformation for initial guess
    dStrainTemp( 2 ) = ( -dStrain( 0 ) - dStrain( 1 ) );

    int planeStressCount = 1;
    while ( true ) {
        stressTemp = stress;
        stateVars  = stateVarsOld;

        computeStress( stressTemp.data(),
                       K_local,
                       nonLocalRadius,
                       dStressDDStrain.data(),
                       dK_localDDStrain,
                       dStressDK,
                       strainOld,
                       dStrainTemp.data(),
                       KOld,
                       dK,
                       timeOld,
                       dT,
                       pNewDT );

        if ( pNewDT < 1.0 ) {
            return;
        }

        double residual = stressTemp.array().abs()[2];

        if ( residual < 1.e-10 || ( planeStressCount > 7 && residual < 1e-5 ) ) {
            break;
        }

        double tangentCompliance = 1. / dStressDDStrain( 2, 2 );
        if ( isNaN( tangentCompliance ) || std::abs( tangentCompliance ) > 1e10 )
            tangentCompliance = 1e10;

        dStrainTemp[2] -= tangentCompliance * stressTemp[2];

        planeStressCount += 1;
        if ( planeStressCount > 10 ) {
            pNewDT = 0.25;
            warningToMSG( "PlaneStressWrapper requires cutback" );
            return;
        }
    }

    dStrain = dStrainTemp;
    stress  = stressTemp;
}
