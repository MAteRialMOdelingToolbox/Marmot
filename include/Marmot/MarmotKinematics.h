#pragma once
#include "MarmotFunctions.h"
#include "MarmotTypedefs.h"
#include "MarmotVoigt.h"

namespace Marmot {
    namespace kinematics {
        namespace strain {
            Marmot::Vector6 GreenLagrange( const Eigen::Matrix3d& F );
            Marmot::Tensor633d dGreenLagrangedDeformationGradient ( const Eigen::Matrix3d& F);
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
} // namespace Marmot
