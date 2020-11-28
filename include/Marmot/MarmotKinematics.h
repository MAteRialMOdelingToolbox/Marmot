#pragma once
#include "Marmot/MarmotFunctions.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

namespace Marmot {
    namespace ContinuumMechanics::Kinematics {
        namespace strain {
            Marmot::Vector6d GreenLagrange( const Eigen::Matrix3d& F );
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
    } // namespace ContinuumMechanics::Kinematics
} // namespace Marmot
