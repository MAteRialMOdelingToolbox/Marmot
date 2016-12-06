#include "bftAbaqusUtility.h"
#include "bftFunctions.h"

/*
 * convinience functions for umats in Abaqus
 * */

namespace bft{

    void backToAbaqus(const Matrix6& jacobian, Map<MatrixXd>& ABQJacobian, const Vector6& stress, Map<VectorXd>& ABQStress, int nTensor)
    {
        ABQStress =     stress.head(nTensor);
        ABQJacobian =   jacobian.topLeftCorner(nTensor, nTensor);
        return;
    }

    void backToAbaqusPlaneStress(const Matrix6& jacobian, Map<MatrixXd>& ABQJacobian, const Vector6& stress, 
                                Map<VectorXd>& ABQStress)
    {
        const int nDirect = 2;
        const int nShear  = 1;
        const int nTensor = 3;

        for (int i=0; i<nDirect; i++){
            ABQStress(i) = stress(i);
        }
        for (int j=0; j<nShear; j++)
            ABQStress(j+nDirect) = stress(j+3);
        ABQStress(0) = stress(0);
        ABQStress(1) = stress(1);
        ABQStress(2) = stress(3);
        ABQJacobian.topLeftCorner(nDirect, nDirect) = jacobian.topLeftCorner(nDirect, nDirect);
        ABQJacobian.block<nShear,nDirect>(nDirect,0) =  jacobian.block<nShear,nDirect>(3,0);
		ABQJacobian.block<nDirect,nShear>(0,nDirect) = jacobian.block<nDirect,nShear>(0,3);
        ABQJacobian(nDirect,nDirect) = jacobian(3,3);
        return;
    }

    void discardIncrementAndBackToAbaqus(double& pNewDT, double value, const std::string& message)
    {
        pNewDT = value;
        warningToMSG(message);
        return;
    }

}
