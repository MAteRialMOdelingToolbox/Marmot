#include "bftCosserat2D.h"
#include "bftCosserat2DConversionTools.h"
#include "bftMaterialNonLocalCosseratHypoElastic.h"
#include "bftTypedefs.h"
#include <iostream>

using namespace Eigen;
using namespace bft;

void BftMaterialNonLocalCosseratHypoElastic::computePlaneStrain( ConstitutiveResponse<2>&       response,
                                                                 Eigen::Map<Eigen::Matrix3d>&   stress3D,
                                                                 Eigen::Map<Eigen::Matrix3d>&   coupleStress3D,
                                                                 AlgorithmicModuli<2>&          algorithmicModuli,
                                                                 const DeformationIncrement<2>& deformationIncrement,
                                                                 const TimeIncrement&           timeIncrement,
                                                                 double&                        pNewDT )
{
    ConstitutiveResponse<3> response3D{stress3D,
                                       coupleStress3D,
                                       response.localField,
                                       response.nonLocalRadius};

    AlgorithmicModuli<3> algorithmicModuli3D;

    DeformationIncrement<3> deformationIncrement3D{bft::Cosserat::to3D( deformationIncrement.dStrain ),
                                                   bft::Cosserat::momentTo3D( deformationIncrement.dCurvature ),
                                                   deformationIncrement.nonLocalField};

    this->computeStress( response3D, algorithmicModuli3D, deformationIncrement3D, timeIncrement, pNewDT );

    stress3D       = response3D.stress.reshaped( 3, 3 );
    coupleStress3D = response3D.coupleStress.reshaped( 3, 3 );

    response = {bft::Cosserat::S_to2D( response3D.stress ),
                bft::Cosserat::M_to2D( response3D.coupleStress ),
                response3D.localField,
                response3D.nonLocalRadius};

    algorithmicModuli = {

        bft::Cosserat::dSdE_to2D_( algorithmicModuli3D.dStressDDStrain ),
        bft::Cosserat::dSdK_to2D_( algorithmicModuli3D.dStressDDCurvature ),
        bft::Cosserat::S_to2D( algorithmicModuli3D.dStressDNonLocalField ),

        bft::Cosserat::dMdE_to2D_( algorithmicModuli3D.dCoupleStressDDStrain ),
        bft::Cosserat::dMdK_to2D_( algorithmicModuli3D.dCoupleStressDDCurvature ),
        bft::Cosserat::M_to2D( algorithmicModuli3D.dCoupleStressDNonLocalField ),

        bft::Cosserat::S_to2D( algorithmicModuli3D.dLocalFieldDDStrain ),
        bft::Cosserat::M_to2D( algorithmicModuli3D.dLocalFieldDDCurvature ),
        algorithmicModuli3D.dLocalFieldDNonLocalField};
}

//void BftMaterialNonLocalCosseratHypoElastic::computePlaneStress( ConstitutiveResponse<2>&       response,
                                                                 //AlgorithmicModuli<2>&          algorithmicModuli,
                                                                 //const DeformationIncrement<2>& deformationIncrement,
                                                                 //const TimeIncrement&           timeIncrement,
                                                                 //double&                        pNewDT )
//{

    //Map<VectorXd> stateVars( this->stateVars, this->nStateVars );
    //VectorXd      stateVarsOld( stateVars );

    //constexpr int nCompS = 5;
    //constexpr int nCompM = 7;

    //// 13, 23, 31,32,33
    //const static Array<double, nCompS, 1> compS{2, 5, 6, 7, 8};
    ////   11, 12, xx, 21, 22, xx, 31, 32, 33
    //const static Array<double, nCompM, 1> compM{0, 1, 2, 3, 4, 5, 8};

    //constexpr int sizeEq = nCompS + nCompM;

    //using XSized  = Matrix<double, sizeEq, 1>;
    //using XXSized = Matrix<double, sizeEq, sizeEq>;

    //XSized dX = XSized::Zero();
    //XSized residual;

    //XXSized ReducedStiffness;

    //ConstitutiveResponse<3> response3D;
    //ConstitutiveResponse<3> response3DOld{bft::Cosserat::to3D_( response.stress ),
                                          //bft::Cosserat::to3D_( response.coupleStress ),
                                          //response.localField,
                                          //response.nonLocalRadius};
    //AlgorithmicModuli<3>    algorithmicModuli3D;

    //DeformationIncrement<3> deformationIncrement3D = {bft::Cosserat::to3D_( deformationIncrement.dStrain ),
                                                      //bft::Cosserat::to3D_( deformationIncrement.dCurvature ),
                                                      //deformationIncrement.nonLocalField};

    //int planeStressCount = 0;
    //while ( true ) {
        //response3D = response3DOld;
        //stateVars  = stateVarsOld;

        //computeStress( response3D, algorithmicModuli3D, deformationIncrement3D, timeIncrement, pNewDT );

        //if ( pNewDT < 1.0 ) {
            //return;
        //}

        //residual.head( nCompS ) = response3D.stress( compS );
        //residual.tail( nCompM ) = response3D.coupleStress( compM );

        //const double resNorm = residual.array().abs().maxCoeff();

        //if ( resNorm < 1.e-12 ) {
            //break;
        //}

        //// clang-format off
         //ReducedStiffness.topLeftCorner      (nCompS, nCompS) = algorithmicModuli3D.dStressDDStrain           ( compS, compS ); 
         //ReducedStiffness.topRightCorner     (nCompS, nCompM) = algorithmicModuli3D.dStressDDCurvature        ( compS, compM ); 
         //ReducedStiffness.bottomLeftCorner   (nCompM, nCompS) = algorithmicModuli3D.dCoupleStressDDStrain     ( compM, compS ); 
         //ReducedStiffness.bottomRightCorner  (nCompM, nCompM) = algorithmicModuli3D.dCoupleStressDDCurvature  ( compM, compM );
        //// clang-format on

        //dX -= ReducedStiffness.colPivHouseholderQr().solve( residual );

        //deformationIncrement3D.dStrain( compS )    = dX.head( nCompS );
        //deformationIncrement3D.dCurvature( compM ) = dX.tail( nCompM );

        //planeStressCount += 1;
        //if ( planeStressCount > 10 ) {
            //pNewDT = 0.25;
            //return;
        //}
    //}

    //response = {bft::Cosserat::S_to2D_( response3D.stress ),
                //bft::Cosserat::M_to2D_( response3D.coupleStress ),
                //response3D.localField,
                //response3D.nonLocalRadius};

    //const auto psAlgorithmicModuli = bft::Cosserat::PlaneStress::
        //computePlaneStressTangents( algorithmicModuli3D.dStressDDStrain,
                                    //algorithmicModuli3D.dStressDDCurvature,
                                    //algorithmicModuli3D.dCoupleStressDDStrain,
                                    //algorithmicModuli3D.dCoupleStressDDCurvature );

    //algorithmicModuli = {psAlgorithmicModuli.dStressdStrain,
                         //psAlgorithmicModuli.dStressdCurvature,
                         //bft::Cosserat::S_to2D_( algorithmicModuli3D.dStressDNonLocalField ),

                         //psAlgorithmicModuli.dCoupleStressdStrain,
                         //psAlgorithmicModuli.dCoupleStressdCurvature,
                         //bft::Cosserat::M_to2D_( algorithmicModuli3D.dCoupleStressDNonLocalField ),

                         //bft::Cosserat::S_to2D_( algorithmicModuli3D.dLocalFieldDDStrain ),
                         //bft::Cosserat::M_to2D_( algorithmicModuli3D.dLocalFieldDDCurvature ),
                         //algorithmicModuli3D.dLocalFieldDNonLocalField};
//}
