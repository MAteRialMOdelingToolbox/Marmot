#pragma once
#include "Marmot/MarmotTypedefs.h"

/*
 * convinience functions for umats in Abaqus
 * */

namespace Marmot {
    [[deprecated("MarmotMaterial input is always Voigt6 - sized -> these aux. functions are not necessary anymore")]]
    void backToAbaqus( const Matrix6d&               jacobian,
                       Eigen::Map<Eigen::MatrixXd>& ABQJacobian,
                       const Marmot::Vector6d&          stress,
                       Eigen::Map<Eigen::VectorXd>& ABQStress,
                       int                          nTensor );

    [[deprecated("MarmotMaterial input is always Voigt6 - sized -> these aux. functions are not necessary anymore")]]
    void backToAbaqusPlaneStress( const Matrix6d&               jacobian,
                                  Eigen::Map<Eigen::MatrixXd>& ABQJacobian,
                                  const Marmot::Vector6d&          stress,
                                  Eigen::Map<Eigen::VectorXd>& ABQStress );

    [[deprecated("MarmotMaterial input is always Voigt6 - sized -> these aux. functions are not necessary anymore")]]
    void backToAbaqusNonLocal( const Matrix6d&              dStressdStrain,
                               Eigen::Ref<Eigen::MatrixXd> ABQdStressDStrain,
                               const Marmot::Vector6d&         stress,
                               Eigen::Ref<Eigen::VectorXd> ABQStress,
                               double                      intParameterLocal,
                               double&                     ABQParameterLocal,
                               const Marmot::Vector6d&         dStressDIntParamNonLocal,
                               Eigen::Ref<Eigen::VectorXd> ABQDStressDIntParamNonLocal,
                               const Marmot::Vector6d&         dIntParamLocalDStrain,
                               Eigen::Ref<Eigen::VectorXd> ABQDIntParameterLocalDStrain,
                               double                      nonLocalRadius,
                               double&                     ABQNonLocalRadius,
                               int                         nTensor );

    void discardTheIncrement( double& pNewDT, double value, const std::string& message );
} // namespace Marmot
