#include "bftAbaqusUmatWrapper.h"
#include "bftVoigt.h"
#include "bftTypedefs.h"
#include <iostream>
#include "bftFunctions.h"

namespace bft{

    void simpleUmat(
		        Ref<MatrixXd> abqC,
		        Ref<VectorXd> abqStress,
                Ref<VectorXd> abqStateVars,
		        const Ref<const VectorXd>& abqStrain,
		        const Ref<const VectorXd>& abqDStrain,
                const Ref<const VectorXd>& matProps,
		        const int nProps,
                double& pNewdT,
                const double charElemlen,
                const double time[2],
                const double dT,
                const int nDirect,
                const int nShear,
		        const int nTensor,
		        const int noEl,
                const int materialID,
                pUmatType umatPointer)
        {

        // dummy values for UMAT which are not needed for stress update
        double sSE = 0.0;
        double sPD = 0.0;
        double sCD = 0.0;
        double rpl = 0.0;
        double ddSigmaddTemp = 0.0;
        double dRpldE = 0.0;
        double dRpldT = 0.0;
        const double temp = 0.0;
        const double dTemp = 0.0;
        const double preDefFieldVarArr = 0.0;
        const double dpreDefFieldVarArr = 0.0;
        const double coords[3] = {0};
        const double dRot[9] = {0};
        const double dfGrd0[9] = {0};
        const double dfGrd1[9] = {0};
        const int noPt = 0;
        const int layer = 0;
        const int kSecPt = 0;
        const int jStep[4] = {0};
        const int kInc = 0;
        const int matNameLength = 80;
        const char matName[80] = "UMAT";

        (*umatPointer)( &abqStress(0), &abqStateVars(0), &abqC(0,0), sSE, sPD, sCD, rpl, 
					            &ddSigmaddTemp, &dRpldE, dRpldT, abqStrain.data(), abqDStrain.data(), 
					            &time[0], dT, temp, dTemp, &preDefFieldVarArr, 	&dpreDefFieldVarArr, 
					            matName, nDirect, nShear, nTensor, abqStateVars.size(), matProps.data(), nProps, 
					            &coords[0], &dRot[0], pNewdT, charElemlen, &dfGrd0[0], &dfGrd1[0],
                				noEl, noPt, layer, kSecPt, &jStep[0], kInc, matNameLength);

        }

	void simpleUmatNonLocal(Ref<MatrixXd>                       dStressdStrain,
                            Ref<VectorXd>                       stress,
                            double&                             intParamLocal,
                            Ref<Vector6>                        dStressDIntParamNonlocal,
                            Ref<Vector6>                        dIntParamLocalDStrain,
                            double&                             nonLocalRadius,
                            Ref<VectorXd>                       stateVars,
                            const  Ref<const VectorXd>&         strainOld,
                            const  Ref<const VectorXd>&         dStrain,
                            double                              intParamNonLocalOld,
                            double                              dIntParamNonLocal,
                            const  Ref<const VectorXd>&         matProps,
                            double&                             pNewdT,
                            double                              time[2],
                            double                              dT,
                            int                                 nDirect,
                            int                                 nShear,
                            int                                 nTensor,
                            int                                 noEl,
                            int                                 noGaussPt,
                            pUmatType                           umatPointer)
    {
        // dummy values for UMAT which are not needed for stress update
        double sSE = 0.0;
        double sPD = 0.0;
        double sCD = 0.0;
        const double preDefFieldVarArr = 0.0;
        const double dpreDefFieldVarArr = 0.0;
        const double coords[3] = {0};
        const double dRot[9] = {0};
        const double dfGrd0[9] = {0};
        const double dfGrd1[9] = {0};
        const int noPt = 0;
        const int layer = 0;
        const int kSecPt = 0;
        const int jStep[4] = {0};
        const int kInc = 0;
        const int matNameLength = 80;
        const char matName[80] = "UMAT";

        (*umatPointer)( &stress(0), &stateVars(0), &dStressdStrain(0,0), sSE, sPD, sCD, intParamLocal, 
					            &dStressDIntParamNonlocal(0), &dIntParamLocalDStrain(0), nonLocalRadius, strainOld.data(), dStrain.data(), 
					            &time[0], dT, intParamNonLocalOld, dIntParamNonLocal, &preDefFieldVarArr, 	&dpreDefFieldVarArr, 
					            matName, nDirect, nShear, nTensor, stateVars.size(), matProps.data(), matProps.size(), 
					            &coords[0], &dRot[0], pNewdT, 0.0 , &dfGrd0[0], &dfGrd1[0],
                				noEl, noPt, layer, kSecPt, &jStep[0], kInc, matNameLength);
    }


    void umatPlaneStress(	Ref<Matrix6> Cep,
	                        Ref<Vector6> stress,
                            Ref<VectorXd> stateVars,
	                        const Ref<const Vector6>& strain,
	                        Ref<Vector6> dStrain,
                            const Ref<const VectorXd>& matProps,
	                        const int nProps,
                            double& pNewdT,
                            const double charElemlen,
                            const double time[2],
                            const double dT,
	                        const int noEl,
                            const int materialID,
                            pUmatType umatPointer)
    {
        int planeStressCount =          1;

        Vector6 stressTemp ;
        VectorXd stateVarsTemp;
        Vector6 dStrainTemp = dStrain;
        
        //assumption of isochoric deformation for initial guess
        dStrainTemp(2) = (- dStrain(0) - dStrain(1));

        while (true)
            {
                stressTemp =        stress;
                stateVarsTemp =     stateVars;

                simpleUmat(Cep, stressTemp, stateVarsTemp, strain, dStrainTemp, 
                            matProps, nProps, pNewdT, charElemlen, time, 
                            dT, 3, 3, 6, noEl, materialID, umatPointer);

                if(pNewdT < 1.0){
                 //   pNewdT = 1e36;
                    //return umatPlaneStressBisectionMethod( Cep, stress, stateVars, strain, 
                    //        dStrain, matProps, nProps, pNewdT, charElemlen, 
                    //        time, dT, noEl, materialID, umatPointer);
                    {pNewdT = 0.25; return;}
                }

                double residual = stressTemp.array().abs()[2];
                

                if (residual <1.e-10 || (planeStressCount > 7 && residual < 1e-5) ) {
                    break;}

                double tangentCompliance = 1./Cep(2,2);
                if(isNaN(tangentCompliance) || std::abs(tangentCompliance) > 1e10)
                    tangentCompliance=1e10;

                dStrainTemp[2] -= tangentCompliance *  stressTemp[2];
                 
                planeStressCount += 1;
                if (planeStressCount > 10)
                  //  return umatPlaneStressBisectionMethod( Cep, stress, stateVars, strain, 
                   //         dStrain, matProps, nProps, pNewdT, charElemlen, 
                    //        time, dT, noEl, materialID, umatPointer);
                     {pNewdT = 0.25;
                         warningToMSG("PlaneStressWrapper requires cutback"); return;}
            }

        dStrain =   dStrainTemp;
        stress =    stressTemp;
        stateVars = stateVarsTemp;
    }

	void simpleUmatPlaneStressNonLocal(
                            Ref<MatrixXd>                       dStressdStrain,
                            Ref<VectorXd>                       stress,
                            double&                             intParamLocal,
                            Ref<Vector6>                        dStressDIntParamNonlocal,
                            Ref<Vector6>                        dIntParamLocalDStrain,
                            double&                             nonLocalRadius,
                            Ref<VectorXd>                       stateVars,
                            const  Ref<const VectorXd>&         strainOld,
                            Ref<VectorXd>                       dStrain,
                            double                              intParamNonLocalOld,
                            double                              dIntParamNonLocal,
                            const  Ref<const VectorXd>&         matProps,
                            double&                             pNewdT,
                            double                              time[2],
                            double                              dT,
                            int                                 nDirect,
                            int                                 nShear,
                            int                                 nTensor,
                            int                                 noEl,
                            int                                 noGaussPt, 
                            pUmatType                           umatPointer) {

        int planeStressCount =          1;

        Vector6 stressTemp ;
        VectorXd stateVarsTemp;
        Vector6 dStrainTemp = dStrain;
        //assumption of isochoric deformation for initial guess
        dStrainTemp(2) = (- dStrain(0) - dStrain(1));

        while (true)
            {
                stressTemp =        stress;
                stateVarsTemp =     stateVars;
               

                simpleUmatNonLocal( dStressdStrain, stressTemp, intParamLocal, dStressDIntParamNonlocal, dIntParamLocalDStrain, 
                        nonLocalRadius, stateVarsTemp, strainOld, dStrainTemp, intParamNonLocalOld, dIntParamNonLocal, 
                        matProps, pNewdT, &time[2], dT, 3, 3, 6, noEl, noGaussPt, umatPointer);

                if(pNewdT < 1.0){
                   // pNewdT = 1e36;
                   // return umatPlaneStressBisectionMethodNonLocal( dStressdStrain, stressTemp, intParamLocal, dStressDIntParamNonlocal, dIntParamLocalDStrain, 
                   //     nonLocalRadius, stateVarsTemp, strainOld, dStrainTemp, intParamNonLocalOld, dIntParamNonLocal, 
                   //     matProps, pNewdT, &time[2], dT, 3, 3, 6, noEl, noGaussPt, umatPointer);
                    {pNewdT = 0.25; return;}
                }

                double residual = stressTemp.array().abs()[2];
                
                if (residual <1.e-10 || (planeStressCount > 7 && residual < 1e-5) ) {
                    break;}

                double tangentCompliance = 1./dStressdStrain(2,2);
                if(isNaN(tangentCompliance) || std::abs(tangentCompliance) > 1e10)
                    tangentCompliance=1e10;

                dStrainTemp[2] -= tangentCompliance *  stressTemp[2];
                 
                planeStressCount += 1;
                if (planeStressCount > 10)
                 //   return umatPlaneStressBisectionMethodNonLocal( dStressdStrain, stressTemp, intParamLocal, dStressDIntParamNonlocal, dIntParamLocalDStrain, 
                  //      nonLocalRadius, stateVarsTemp, strainOld, dStrainTemp, intParamNonLocalOld, dIntParamNonLocal, 
                   //     matProps, pNewdT, &time[2], dT, 3, 3, 6, noEl, noGaussPt, umatPointer);
                     {  
                         warningToMSG("PlaneStressWrapper requires cutback"); 
                         pNewdT = 0.25; return;}
            }

        dStrain =   dStrainTemp;
        stress =    stressTemp;
        stateVars = stateVarsTemp;

    }

    void umatPlaneStressBisectionMethod(	Ref<Matrix6> Cep,
	                        Ref<Vector6> stress,
                            Ref<VectorXd> stateVars,
	                        const Ref<const Vector6>& strain,
	                        Ref<Vector6> dStrain,
                            const Ref<const VectorXd>& matProps,
	                        const int nProps,
                            double& pNewdT,
                            const double charElemlen,
                            const double time[2],
                            const double dT,
	                        const int noEl,
                            const int materialID,
                            pUmatType umatPointer)
    {
        int planeStressCount = 1;

        Vector6 stressTemp ;
        VectorXd stateVarsTemp;
        Vector6 dStrainTemp = dStrain;


        // factor 2 is arbitrary, maybe it can be improved
        double dStrainOutside = 2 * std::abs(dStrain(0)) > std::abs(dStrain(1)) ? -dStrain(0) : -dStrain(1);
        double dStrainInside = -dStrainOutside;
        double sgn = dStrainOutside > 0 ? 1 : -1;

        double dStrainCenter;
        double s33;

        std::ostringstream str;

        while (true){
            planeStressCount +=1;
            stressTemp =        stress;
            stateVarsTemp =     stateVars;

            dStrainCenter  = (dStrainInside + dStrainOutside) / 2;
            dStrainTemp(2) = dStrainCenter;
            simpleUmat(Cep, stressTemp, stateVarsTemp, strain, dStrainTemp, 
                        matProps, nProps, pNewdT, charElemlen, time, 
                        dT, 3, 3, 6, noEl, materialID, umatPointer);

            if (pNewdT < 1.0)
                return ;

            if(planeStressCount > 20) {
                pNewdT = 0.25;
                return; }

            s33 = stressTemp(2);

            if (std::abs(s33)<1.e-10 || (planeStressCount > 15 && std::abs(s33)<1.e-2 ) ) 
                break;

            if( (s33 < 0  && sgn > 0) || (s33 >= 0 && sgn < 0) )
                dStrainInside = dStrainCenter;
            else
                dStrainOutside = dStrainCenter;
        }

        dStrain =   dStrainTemp;
        stress =    stressTemp;
        stateVars = stateVarsTemp;
    }

void umatPlaneStressBisectionMethodNonLocal(
                            Ref<MatrixXd>                       dStressdStrain,
                            Ref<VectorXd>                       stress,
                            double&                             intParamLocal,
                            Ref<Vector6>                        dStressDIntParamNonlocal,
                            Ref<Vector6>                        dIntParamLocalDStrain,
                            double&                             nonLocalRadius,
                            Ref<VectorXd>                       stateVars,
                            const  Ref<const VectorXd>&         strainOld,
                            Ref<VectorXd>                       dStrain,
                            double                              intParamNonLocalOld,
                            double                              dIntParamNonLocal,
                            const  Ref<const VectorXd>&         matProps,
                            double&                             pNewdT,
                            double                              time[2],
                            double                              dT,
                            int                                 nDirect,
                            int                                 nShear,
                            int                                 nTensor,
                            int                                 noEl,
                            int                                 noGaussPt, 
                            pUmatType                           umatPointer) {
        int planeStressCount = 1;

        Vector6 stressTemp ;
        VectorXd stateVarsTemp;
        Vector6 dStrainTemp = dStrain;


        // factor 2 is arbitrary, maybe it can be improved
        double dStrainOutside = 2 * std::abs(dStrain(0)) > std::abs(dStrain(1)) ? -dStrain(0) : -dStrain(1);
        double dStrainInside = -dStrainOutside;
        double sgn = dStrainOutside > 0 ? 1 : -1;

        double dStrainCenter;
        double s33;

        std::ostringstream str;

        while (true){
            planeStressCount +=1;
            stressTemp =        stress;
            stateVarsTemp =     stateVars;

            dStrainCenter  = (dStrainInside + dStrainOutside) / 2;
            dStrainTemp(2) = dStrainCenter;
            simpleUmatNonLocal( dStressdStrain, stressTemp, intParamLocal, dStressDIntParamNonlocal, dIntParamLocalDStrain, 
                        nonLocalRadius, stateVarsTemp, strainOld, dStrainTemp, intParamNonLocalOld, dIntParamNonLocal, 
                        matProps, pNewdT, &time[2], dT, 3, 3, 6, noEl, noGaussPt, umatPointer);

            if (pNewdT < 1.0)
                return ;

            if(planeStressCount > 20) {
                pNewdT = 0.25;
                return; }

            s33 = stressTemp(2);

            if (std::abs(s33)<1.e-10 || (planeStressCount > 15 && std::abs(s33)<1.e-2 ) ) 
                break;

            if( (s33 < 0  && sgn > 0) || (s33 >= 0 && sgn < 0) )
                dStrainInside = dStrainCenter;
            else
                dStrainOutside = dStrainCenter;
        }

        dStrain =   dStrainTemp;
        stress =    stressTemp;
        stateVars = stateVarsTemp;
    }
}
