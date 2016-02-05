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

    void discardIncrementAndBackToAbaqus(double& pNewDT, double value, const std::string& message)
    {
        pNewDT = value;
        warningToMSG(message);
        return;
    }

}
