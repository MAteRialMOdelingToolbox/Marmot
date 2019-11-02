#include "bftMaterialNonLocalCosseratHypoElastic.h"
#include "bftTypedefs.h"
#include "bftCosserat2DConversionTools.h"
#include <iostream>

using namespace Eigen;
using namespace bft;


void BftMaterialNonLocalCosseratHypoElastic::computePlaneStrain( 
                                ConstitutiveResponse<2>& response,
                                AlgorithmicModuli<2>& algorithmicModuli,
                                const DeformationIncrement<2>& deformationIncrement,
                                const TimeIncrement& timeIncrement,
                                double& pNewDT )
{
    ConstitutiveResponse<3> response3D
    {
        bft::Cosserat::to3D( response.stress ),
        bft::Cosserat::momentTo3D( response.coupleStress),
        response.localField,
        response.nonLocalRadius
    };

    AlgorithmicModuli<3> algorithmicModuli3D;

    DeformationIncrement<3> deformationIncrement3D
    {
        bft::Cosserat::to3D( deformationIncrement.dStrain),
        bft::Cosserat::momentTo3D( deformationIncrement.dCurvature),
        deformationIncrement.nonLocalField
    };

    this->computeStress(response3D, algorithmicModuli3D, deformationIncrement3D, timeIncrement, pNewDT);

    response = {
        bft::Cosserat::S_to2D( response3D.stress ),
        bft::Cosserat::M_to2D( response3D.coupleStress),
        response3D.localField,
        response3D.nonLocalRadius
    };

    algorithmicModuli = {

        bft::Cosserat::dSdE_to2D ( algorithmicModuli3D.dStressDDStrain),
        bft::Cosserat::dSdK_to2D ( algorithmicModuli3D.dStressDDCurvature),
        bft::Cosserat::S_to2D (   algorithmicModuli3D.dStressDNonLocalField),

        bft::Cosserat::dMdE_to2D ( algorithmicModuli3D.dCoupleStressDDStrain),
        bft::Cosserat::dMdK_to2D ( algorithmicModuli3D.dCoupleStressDDCurvature),
        bft::Cosserat::M_to2D (   algorithmicModuli3D.dCoupleStressDNonLocalField),

        bft::Cosserat::S_to2D (   algorithmicModuli3D.dLocalFieldDDStrain),
        bft::Cosserat::M_to2D (   algorithmicModuli3D.dLocalFieldDDCurvature),

        algorithmicModuli3D.dLocalFieldDNonLocalField
    };

}

//void BftMaterialNonLocalCosseratHypoElastic::computePlaneStress( ConstitutiveResponse* response,
                                                                 //AlgorithmicModuli*    algorithmicModuli,
                                                                 //DeformationIncrement* deformationIncrement,
                                                                 //const TimeIncrement*  timeIncrement,
                                                                 //double&               pNewDT )
//{
    //using mVector9 = Map< Matrix<double, 9, 1>> ;
    //using mMatrix9 = Map< Matrix<double, 9, 9>> ;

    //Map<VectorXd> stateVars( this->stateVars, this->nStateVars );
    //VectorXd stateVarsOld (  stateVars ) ;


    //constexpr int nS = 5;
    //constexpr int nM = 7;

                                               ////13, 23, 31,32,33
    //const static Array<double, nS, 1> sIndices {  2,   5, 6, 7, 8};
                                               ////   11, 12, xx, 21, 22, xx, 31, 32, 33 
    //const static Array<double, nM, 1> mIndices {  0,   1,      3,  4,      6,  7,  8};

    //constexpr int sizeEq = nS + nM;

    //Vector9d stressOld ( mVector9 ( response-> stress ));
    //Vector9d coupleStressOld ( mVector9 ( response-> coupleStress ));

    //using XSized = Matrix<double, sizeEq, 1>;
    //using XXSized = Matrix<double, sizeEq, sizeEq>;

    //XSized dX = XSized::Zero();
    //XSized residual;

    //XXSized  ReducedStiffness;

    //int planeStressCount = 0;
    //while ( true ) {
        //mVector9 ( response->stress ) = stressOld;
        //mVector9 ( response->coupleStress ) = coupleStressOld;
        //stateVars  = stateVarsOld;

        //computeStress ( response, 
                        //algorithmicModuli,
                        //deformationIncrement,
                        //timeIncrement,
                        //pNewDT);

        //if ( pNewDT < 1.0 ) {
            //return;
        //}

        //residual.head( nS ) = mVector9 ( response->stress ) ( sIndices );
        //residual.tail( nM ) = mVector9 ( response->coupleStress ) ( mIndices );

        //const double resNorm = residual.array().abs().maxCoeff();

        //if ( resNorm < 1.e-14 ){
            //std::cout << resNorm << std::endl; 
            //return;
        //} 

        //// clang-format off
        //ReducedStiffness.topLeftCorner      (nS, nS) = mMatrix9 ( algorithmicModuli->dStressDDStrain )          ( sIndices, sIndices );
        //ReducedStiffness.topRightCorner     (nS, nM) = mMatrix9 ( algorithmicModuli->dStressDDCurvature)        ( sIndices, mIndices );
        //ReducedStiffness.bottomLeftCorner   (nM, nS) = mMatrix9 ( algorithmicModuli->dCoupleStressDDStrain )    ( mIndices, sIndices );
        //ReducedStiffness.bottomRightCorner  (nM, nM) = mMatrix9 ( algorithmicModuli->dCoupleStressDDCurvature)  ( mIndices, mIndices );
        //// clang-format on

        //dX -= ReducedStiffness.fullPivHouseholderQr().solve( residual );

        //mVector9 ( deformationIncrement-> dStrain ) (sIndices) = dX.head(nS);
        //mVector9 ( deformationIncrement-> dCurvature ) (mIndices ) = dX.tail ( nM );

        //planeStressCount += 1;
        //if ( planeStressCount > 5 ) {
            //pNewDT = 0.25;
            //return;
        //}
    //}
//}
