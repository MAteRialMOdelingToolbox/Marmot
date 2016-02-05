#pragma once 
#include "bftTypedefs.h"

/*
 * convinience functions for umats in Abaqus
 * */

namespace bft{
    void backToAbaqus(const Matrix6& jacobian, Map<MatrixXd>& ABQJacobian, const Vector6& stress, Map<VectorXd>& ABQStress, int nTensor);
    void discardIncrementAndBackToAbaqus(double& pNewDT, double value, const std::string& message);
}
