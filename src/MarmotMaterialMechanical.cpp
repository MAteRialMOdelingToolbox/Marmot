#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialMechanical.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotVoigt.h"
#include <iostream>

using namespace Eigen;

//void MarmotMaterialHypoElastic::computeStress( double*       stress_,
                                            //double*       dStressDDStrain_,
                                            //const double* FOld_,
                                            //const double* FNew_,
                                            //const double* timeOld,
                                            //const double  dT,
                                            //double&       pNewDT )
//{
    //const Map<const Matrix3d> FNew( FNew_ );
    //const Map<const Matrix3d> FOld( FOld_ );

     ////Marmot::Vector6d dEps = 1./2 * ( Marmot::ContinuumMechanics::VoigtNotation::StrainToVoigt( H + H.tranpose() ) );

     ////computeStress (stress_, dStressDDStrain_, dEps.data(), timeOld, dT, pNewDT);
//}

void MarmotMaterialMechanical::computePlaneStress( double*       stress_,
                                             double*       dStressDDDeformationGradient_,
                                             const double* FOld_,
                                             double*       FNew_,
                                             const double* timeOld,
                                             const double  dT,
                                             double&       pNewDT )
{

    using namespace Marmot;

    Map<Vector6d>  stress( stress_ );
    Map<Matrix<double, 6,9>>  dStressDDDeformationGradient( dStressDDDeformationGradient_);
    Map<Matrix3d>  FNew( FNew_);
    Map<VectorXd> stateVars( this->stateVars, this->nStateVars );

    Vector6d  stressTemp;
    VectorXd stateVarsOld = stateVars;
    Matrix3d FNewTemp = FNew;

    // assumption of isochoric deformation for initial guess
    FNewTemp( 2,2 ) = 1./ ( FNew( 0,0 ) * FNew( 1,1 ) );

    int planeStressCount = 1;
    while ( true ) {
        stressTemp = stress;
        stateVars  = stateVarsOld;

        computeStress( stressTemp.data(), dStressDDDeformationGradient.data(), FOld_, FNewTemp.data(), timeOld, dT, pNewDT );

        if ( pNewDT < 1.0 ) {
            return;
        }

        double residual = stressTemp.array().abs()[2];

        if ( residual < 1.e-10 || ( planeStressCount > 7 && residual < 1e-8 ) ) {
            break;
        }

        const double dS33_dF33 = dStressDDDeformationGradient( 2, 8 );

        double tangentCompliance = 1. / dS33_dF33;
        if ( Math::isNaN( tangentCompliance ) || std::abs( tangentCompliance ) > 1e10 )
            tangentCompliance = 1e10;

        FNewTemp(2,2) -= tangentCompliance * stressTemp(2);

        planeStressCount += 1;
        if ( planeStressCount > 13 ) {
            pNewDT = 0.25;
            MarmotJournal::warningToMSG( "PlaneStressWrapper requires cutback" );
            return;
        }
    }

    FNew = FNewTemp;
    stress  = stressTemp;
}

void MarmotMaterialMechanical::computeUniaxialStress( double*       stress_,
                                                double*       dStressDDStrain_,
                                                const double* FOld,
                                                double*       FNew,
                                                const double* timeOld,
                                                const double  dT,
                                                double&       pNewDT )
{
    // using namespace Marmot;

    // Map<Vector6d>  stress( stress_ );
    // Map<Matrix6d>  dStressDDStrain( dStressDDStrain_ );
    // Map<Vector6d>  dStrain( dStrain_ );
    // Map<VectorXd> stateVars( this->stateVars, this->nStateVars );

    // Vector6d  stressTemp;
    // VectorXd stateVarsOld = stateVars;
    // Vector6d  dStrainTemp  = dStrain;

    // dStrainTemp( 2 ) = 0.0;
    // dStrainTemp( 1 ) = 0.0;

    // int count = 1;
    // while ( true ) {
    // stressTemp = stress;
    // stateVars  = stateVarsOld;

    // computeStress( stressTemp.data(), dStressDDStrain.data(), dStrainTemp.data(), timeOld, dT, pNewDT );

    // if ( pNewDT < 1.0 ) {
    // return;
    //}

    // double residual = stressTemp.array().abs()[1] + stressTemp.array().abs()[2];

    // if ( residual < 1.e-10 || ( count > 7 && residual < 1e-8 ) ) {
    // break;
    //}

    // dStrainTemp.segment<2>( 1 ) -= dStressDDStrain.block<2, 2>( 1, 1 ).colPivHouseholderQr().solve(
    // stressTemp.segment<2>( 1 ) );

    // count += 1;
    // if ( count > 13 ) {
    // pNewDT = 0.25;
    // warningToMSG( "UniaxialStressWrapper requires cutback" );
    // return;
    //}
    //}

    // dStrain = dStrainTemp;
    // stress  = stressTemp;
}
