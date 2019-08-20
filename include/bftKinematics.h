#pragma once
#include "bftFunctions.h"
#include "bftTypedefs.h"
#include "bftVoigt.h"

namespace bft {
    namespace kinematics {
        namespace strain {
            bft::Vector6 GreenLagrange( const Eigen::Matrix3d& F );
            bft::Tensor633d dGreenLagrangedDeformationGradient ( const Eigen::Matrix3d& F);
        }

        namespace velocityGradient {
            extern const Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3, 3, 3>> dOmega_dVelocityGradient;

            extern const Eigen::TensorFixedSize<double, Eigen::Sizes<6, 3, 3>> dStretchingRate_dVelocityGradient;
        } // namespace velocityGradient

        namespace deformationGradient {
            template <int nDim>
            Eigen::Matrix3d make3D( const Eigen::Ref<const Eigen::Matrix<double, nDim, nDim>>& tensor );
        }
    } // namespace kinematics
} // namespace bft
