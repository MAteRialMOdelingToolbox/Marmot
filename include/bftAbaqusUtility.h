#pragma once 
#include "bftTypedefs.h"

/*
 * convinience functions for umats in Abaqus
 * */

namespace bft{
    void backToAbaqus(const Matrix6& jacobian, Map<MatrixXd>& ABQJacobian, const Vector6& stress, Map<VectorXd>& ABQStress, int nTensor);
    void backToAbaqusPlaneStress(const Matrix6& jacobian, Map<MatrixXd>& ABQJacobian, const Vector6& stress, Map<VectorXd>& ABQStress);
    void backToAbaqusNonLocal(  const Matrix6& dStressdStrain,              Ref<MatrixXd>  ABQdStressDStrain, 
                                const Vector6& stress,                      Ref<VectorXd>  ABQStress,
                                double intParameterLocal,                   double&         ABQParameterLocal,
                                const Vector6& dStressDIntParamNonLocal,    Ref<VectorXd>   ABQDStressDIntParamNonLocal,
                                const Vector6& dIntParamLocalDStrain,       Ref<VectorXd>  ABQDIntParameterLocalDStrain,
                                double nonLocalRadius,                      double&         ABQNonLocalRadius,
                                int nTensor);
    void discardIncrementAndBackToAbaqus(double& pNewDT, double value, const std::string& message);
}
