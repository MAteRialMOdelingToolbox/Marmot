#include "bftVoigt.h"
#include "bftTypedefs.h"
#include <iostream>

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


    void umatPlaneStress(	Ref<Matrix6> Cep,
	                        Ref<Vector6> stress,
                            Ref<VectorXd> stateVarLocal,
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
        if (isNaN(stress))
            std::cout << "PLANE STRESS WRAPPER RECEIVED NAN STRESS" << std::endl;

        bool planeStressConvergence = false;
        int planeStressCount = 1;
        bool flagDebug = true;

        const int nDirect = 3;
        const int nShear = 3;
        const int nTensor = 6;

        Vector6d stressTemp ;
        VectorXd stateVarTemp ;
        Vector6d dStrainTemp = dStrain;

        while (planeStressConvergence == false)
            {
                stressTemp      = stress;
                stateVarTemp    = stateVarLocal;

	            simpleUmat(Cep, stressTemp, stateVarTemp, strain, dStrainTemp, 
                            matProps, nProps, pNewdT, charElemlen, time, 
                            dT, nDirect, nShear, nTensor, noEl, materialID, umatPointer);

                if (stressTemp.array().abs()[2]<1.e-8)
                {
                    planeStressConvergence = true;
                    break;
                }

                if (Cep(2,2) < 1.e-10)
                {   
                    std::cout << "division by 0 Cep: " << Cep(2,2) << std::endl;
                    Cep(2,2) = 1.e-10;
                } 
                if (isNaN(dStrainTemp))
                    std::cout <<  planeStressCount << "    PlaneStressLoop dStrain\n" << dStrain << std::endl;
                               
                dStrainTemp[2] -= 1./Cep(2,2) * ( stressTemp[2] );
                 
                planeStressCount += 1;

                if (planeStressCount > 30)
                {
                    std::cout << "PlaneStressWrapper: too many iterations; R[i] was "<< stressTemp.array().abs()[2] <<"ELEMENT " << noEl <<  "time" << time[2] << std::endl;
                    pNewdT = 0.25;
                    std::cout << "dT was " << dT << " pNewDt is " << pNewdT << std::endl;
                    return; 
                }
            }

        dStrain = dStrainTemp;
        stress = stressTemp;
        stateVarLocal = stateVarTemp;
    }

}//end of namespace bft
