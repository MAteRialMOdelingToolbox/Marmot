#include "bftMaterialHyperElastic.h"
#include "bftFunctions.h"
#include "bftVoigt.h"
#include <iostream>

using namespace Eigen;

void BftMaterialHyperElastic::computeStress( double*       S,
                                                  double*       dSdE,
                                                  const double* F_n,
                                                  const double* F_np,
                                                  const double* timeOld,
                                                  const double  dT,
                                                  double&       pNewDT )
{



}

void BftMaterialHyperElastic::computePlaneStressPK2( double*       S,
                                                  double*       dSdE,
                                                  double* E,
                                                  //const double* F,
                                                  const double* timeOld,
                                                  const double  dT,
                                                  double&       pNewDT )

{
    // using namespace bft;

    // Map<Vector6>  stress( stress_ );
    // Map<Matrix6>  dStressDDStrain( dStressDDStrain_ );
    // Map<Vector6>  dStrain( dStrain_ );
    // Map<VectorXd> stateVars( this->stateVars, this->nStateVars );

    // Vector6  stressTemp;
    // VectorXd stateVarsOld = stateVars;
    // Vector6  dStrainTemp  = dStrain;

    //// assumption of isochoric deformation for initial guess
    // dStrainTemp( 2 ) = ( -dStrain( 0 ) - dStrain( 1 ) );

    // int planeStressCount = 1;
    // while ( true ) {
    // stressTemp = stress;
    // stateVars  = stateVarsOld;


    // if ( pNewDT < 1.0 ) {
    // return;
    //}

    // double residual = stressTemp.array().abs()[2];

    // if ( residual < 1.e-10 || ( planeStressCount > 7 && residual < 1e-8 ) ) {
    // break;
    //}

    // double tangentCompliance = 1. / dStressDDStrain( 2, 2 );
    // if ( isNaN( tangentCompliance ) || std::abs( tangentCompliance ) > 1e10 )
    // tangentCompliance = 1e10;

    // dStrainTemp[2] -= tangentCompliance * stressTemp[2];

    // planeStressCount += 1;
    // if ( planeStressCount > 13 ) {
    // pNewDT = 0.25;
    // BftJournal::warningToMSG( "PlaneStressWrapper requires cutback" );
    // return;
    //}
    //}

    // dStrain = dStrainTemp;
    // stress  = stressTemp;
}

void BftMaterialHyperElastic::computeUniaxialStressPK2( double*       S,
                                                     double*       dSdE,
                                                     double* E,
                                                     //const double* F,
                                                     const double* timeOld,
                                                     const double  dT,
                                                     double&       pNewDT )
{
    // using namespace bft;

    // Map<Vector6>  stress( stress_ );
    // Map<Matrix6>  dStressDDStrain( dStressDDStrain_ );
    // Map<Vector6>  dStrain( dStrain_ );
    // Map<VectorXd> stateVars( this->stateVars, this->nStateVars );

    // Vector6  stressTemp;
    // VectorXd stateVarsOld = stateVars;
    // Vector6  dStrainTemp  = dStrain;

    // dStrainTemp( 2 ) = 0.0;
    // dStrainTemp( 1 ) = 0.0;

    // int count = 1;
    // while ( true ) {
    // stressTemp = stress;
    // stateVars  = stateVarsOld;


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
    // BftJournal::warningToMSG( "UniaxialStressWrapper requires cutback" );
    // return;
    //}
    //}

    // dStrain = dStrainTemp;
    // stress  = stressTemp;
}
