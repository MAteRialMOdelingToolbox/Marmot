#include "Marmot/HughesWinget.h"
#include "Marmot/MarmotFunctions.h"
#include "Marmot/MarmotKinematics.h"
#include "Marmot/MarmotMaterialHypoElasticNonLocal.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotVoigt.h"

using namespace Eigen;

void MarmotMaterialHypoElasticNonLocal::computeStress( double*       stress_,
                                                    double&       K_local,
                                                    double&       nonLocalRadius,
                                                    double*       dStressDDDeformationGradient_,
                                                    double*       dK_localDDeformationGradient_,
                                                    double*       dStressDK,
                                                    const double* FOld_,
                                                    const double* FNew_,
                                                    const double  KOld,
                                                    const double  dK,
                                                    const double* timeOld,
                                                    const double  dT,
                                                    double&       pNewDT )
{
    // Standard implemenation of the Abaqus like Hughes-Winget algorithm
    // Approximation of the algorithmic tangent in order to
    // facilitate the dCauchy_dStrain tangent provided by
    // small strain material models

    using namespace Marmot;
    const Map<const Matrix3d> FNew( FNew_ );
    const Map<const Matrix3d> FOld( FOld_ );
    Marmot::mVector6d             stress( stress_ );

    Matrix6 CJaumann                = Matrix6::Zero();
    Vector6d dK_LocalDStretchingRate = Vector6d::Zero();

    HughesWinget hughesWingetIntegrator( FOld, FNew, HughesWinget::Formulation::AbaqusLike );

    auto dEps = hughesWingetIntegrator.getStrainIncrement();
    stress    = hughesWingetIntegrator.rotateTensor( stress );

    computeStress( stress.data(),
                   K_local,
                   nonLocalRadius,
                   CJaumann.data(),
                   dK_LocalDStretchingRate.data(),
                   dStressDK,
                   dEps.data(),
                   KOld,
                   dK,
                   timeOld,
                   dT,
                   pNewDT );

    TensorMap<Eigen::Tensor<double, 3>> dS_dF( dStressDDDeformationGradient_, 6, 3, 3 );
    Map<Matrix3d>                       dKLocal_dF( dK_localDDeformationGradient_ );

    Matrix3d FInv = FNew.inverse();
    dS_dF         = hughesWingetIntegrator.compute_dS_dF( stress, FInv, CJaumann );
    dKLocal_dF    = hughesWingetIntegrator.compute_dScalar_dF( FInv, dK_LocalDStretchingRate );
}

void MarmotMaterialHypoElasticNonLocal::computePlaneStress( double*       stress_,
                                                         double&       K_local,
                                                         double&       nonLocalRadius,
                                                         double*       dStressDDStrain_,
                                                         double*       dK_localDDStrain,
                                                         double*       dStressDK,
                                                         double*       dStrain_,
                                                         double        KOld,
                                                         double        dK,
                                                         const double* timeOld,
                                                         const double  dT,
                                                         double&       pNewDT )
{
    using namespace Marmot;

    Map<Vector6d>  stress( stress_ );
    Map<Matrix6>  dStressDDStrain( dStressDDStrain_ );
    Map<Vector6d>  dStrain( dStrain_ );
    Map<VectorXd> stateVars( this->stateVars, this->nStateVars );

    Vector6d  stressTemp;
    VectorXd stateVarsOld = stateVars;
    Vector6d  dStrainTemp  = dStrain;

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
            MarmotJournal::warningToMSG( "PlaneStressWrapper requires cutback" );
            return;
        }
    }

    dStrain = dStrainTemp;
    stress  = stressTemp;
}
