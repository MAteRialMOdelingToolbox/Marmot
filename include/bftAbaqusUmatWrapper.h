#pragma once 
#include "bftTypedefs.h"

namespace bft{
	void simpleUmat( Ref<MatrixXd> abqC,
						 Ref<VectorXd> abqStress,
            			 Ref<VectorXd> abqStateVars,
						const  Ref<const VectorXd>& abqStrain,
						const  Ref<const VectorXd>& abqDStrain,
            			const  Ref<const Eigen::VectorXd>& matProps,
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
                        pUmatType umatPointer);

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
                            pUmatType umatPointer);
}//end of namespace bft
