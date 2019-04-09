#pragma once
#include "bftFunctions.h"
#include "bftTypedefs.h"
#include "bftVoigt.h"

namespace bft {
    namespace kinematics {

        template <int nDim>
        Eigen::Matrix<double, VOIGTFROMDIM( nDim ), 1> GreenLagrangeStrain(
            const Eigen::Matrix<double, nDim, nDim>& DeformationGradient )
        {
            Eigen::Matrix<double, nDim, nDim> H = DeformationGradient - Eigen::Matrix<double, nDim, nDim>::Identity();
            return bft::Vgt::VoigtFromStrainMatrix<nDim>( 0.5 * ( H + H.transpose() + H.transpose() * H ) );
        }

        bft::Vector6 infinitesimalStrainIncrement( const Eigen::Ref<const Eigen::Matrix3d>& F_n,
                                                   const Eigen::Ref<const Eigen::Matrix3d>& F_np,
                                                   double                                   alpha = 0.5 );

        template <int nDim>
        Eigen::Matrix3d makeDeformationGradient3D( const Eigen::Ref<const Eigen::Matrix<double, nDim, nDim>>& F );

        Vector6 transformVoigtLeft(const Eigen::Matrix3d& Q, const Vector6& voigt);
        Vector6 transformVoigtRight(const Vector6& voigt, const Eigen::Matrix3d& Q);


    } // namespace kinematics
} // namespace bft
