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
        bool planeStressConvergence =   false;
        int planeStressCount =          1;

        const int nDirect =             3;
        const int nShear =              3;
        const int nTensor =             6;

        Vector6 stressTemp ;
        VectorXd stateVarsTemp;
        Vector6 dStrainTemp = dStrain;

        while (true)
            {
                stressTemp =        stress;
                stateVarsTemp =     stateVars;

	            simpleUmat(Cep, stressTemp, stateVarsTemp, strain, dStrainTemp, 
                            matProps, nProps, pNewdT, charElemlen, time, 
                            dT, nDirect, nShear, nTensor, noEl, materialID, umatPointer);

                if (stressTemp.array().abs()[2]<1.e-8) {
                    break;}

                if (Cep(2,2) < 1.e-12) 
                    Cep(2,2) = 1.e-12; 

                dStrainTemp[2] -= 1./Cep(2,2) *  stressTemp[2];
                 
                planeStressCount += 1;
                if (planeStressCount > 10) {
                    pNewdT = 0.25;
                    return; }
            }

        dStrain =   dStrainTemp;
        stress =    stressTemp;
        stateVars = stateVarsTemp;
    }
}
