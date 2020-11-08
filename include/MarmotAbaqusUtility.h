#pragma once
#include "MarmotTypedefs.h"

/*
 * convinience functions for umats in Abaqus
 * */

namespace marmot {
    [[deprecated("MarmotMaterial input is always Voigt6 - sized -> these aux. functions are not necessary anymore")]]
    void backToAbaqus( const Matrix6&               jacobian,
                       Eigen::Map<Eigen::MatrixXd>& ABQJacobian,
                       const marmot::Vector6&          stress,
                       Eigen::Map<Eigen::VectorXd>& ABQStress,
                       int                          nTensor );

    [[deprecated("MarmotMaterial input is always Voigt6 - sized -> these aux. functions are not necessary anymore")]]
    void backToAbaqusPlaneStress( const Matrix6&               jacobian,
                                  Eigen::Map<Eigen::MatrixXd>& ABQJacobian,
                                  const marmot::Vector6&          stress,
                                  Eigen::Map<Eigen::VectorXd>& ABQStress );

    [[deprecated("MarmotMaterial input is always Voigt6 - sized -> these aux. functions are not necessary anymore")]]
    void backToAbaqusNonLocal( const Matrix6&              dStressdStrain,
                               Eigen::Ref<Eigen::MatrixXd> ABQdStressDStrain,
                               const marmot::Vector6&         stress,
                               Eigen::Ref<Eigen::VectorXd> ABQStress,
                               double                      intParameterLocal,
                               double&                     ABQParameterLocal,
                               const marmot::Vector6&         dStressDIntParamNonLocal,
                               Eigen::Ref<Eigen::VectorXd> ABQDStressDIntParamNonLocal,
                               const marmot::Vector6&         dIntParamLocalDStrain,
                               Eigen::Ref<Eigen::VectorXd> ABQDIntParameterLocalDStrain,
                               double                      nonLocalRadius,
                               double&                     ABQNonLocalRadius,
                               int                         nTensor );

    void discardTheIncrement( double& pNewDT, double value, const std::string& message );
} // namespace marmot
