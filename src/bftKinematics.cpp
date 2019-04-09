#include "bftKinematics.h"

using namespace Eigen;

namespace bft {
    namespace kinematics {

        bft::Vector6 infinitesimalStrainIncrement( const Eigen::Ref<const Eigen::Matrix3d>& F_n,
                                                   const Eigen::Ref<const Eigen::Matrix3d>& F_np,
                                                   double                                   alpha )
        {
            Matrix3d e;

            const Matrix3d dF = F_np - F_n;
            const Matrix3d H  = dF * ( F_n + alpha * dF ).inverse();

            alpha == 0.5 ? e = 0.5 * ( H + H.transpose() )
                         : e = 0.5 * ( H + H.transpose() + ( 1 - 2 * alpha ) * H.transpose() * H );

            return bft::Vgt::VoigtFromStrainMatrix( e );
        }

        template <int nDim>
        Eigen::Matrix3d makeDeformationGradient3D( const Eigen::Ref<const Eigen::Matrix<double, 3, 3>>& F )
        {
            return F;
        }

        template <int nDim>
        Eigen::Matrix3d makeDeformationGradient3D( const Eigen::Ref<const Eigen::Matrix<double, 2, 2>>& F )
        {
            Matrix3d F3D              = Matrix3d::Zero();
            F3D.topLeftCorner( 2, 2 ) = F;
            return F3D;
        }

        template <int nDim>
        Eigen::Matrix3d makeDeformationGradient3D( const Eigen::Ref<const Eigen::Matrix<double, 1, 1>>& F )
        {
            Matrix3d F3D = Matrix3d::Zero();
            F3D( 0, 0 )  = F( 0, 0 );
            return F3D;
        }

    } // namespace kinematics
}     // namespace bft
