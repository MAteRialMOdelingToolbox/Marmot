#include "Marmot/MarmotKinematics.h"
#include "Marmot/MarmotTensor.h"

using namespace Eigen;

namespace Marmot::ContinuumMechanics::Kinematics {

  namespace Strain {

    Marmot::Math::Vector6d GreenLagrange( const Eigen::Matrix3d& F )
    {
      Eigen::Matrix3d H = F - Eigen::Matrix3d::Identity();
      return Marmot::ContinuumMechanics::Voigt::voigtFromStrainMatrix< 3 >( 0.5 *
                                                                            ( H + H.transpose() + H.transpose() * H ) );
    }

    Marmot::Math::Tensors::Tensor633d dGreenLagrangedDeformationGradient( const Eigen::Matrix3d& F )
    {
      Marmot::Math::Tensors::Tensor633d dEdF;
      auto                              kron = Marmot::Math::Matrix3d::Identity();

      for ( int IJ = 0; IJ < 6; IJ++ ) {
        auto [I, J] = Marmot::Math::Tensors::IndexNotation::fromVoigt< 3 >( IJ );
        for ( int k = 0; k < 3; k++ )
          for ( int L = 0; L < 3; L++ )

            dEdF( IJ, k, L ) = 0.5 * ( kron( I, L ) * F( k, J ) + kron( J, L ) * F( k, I ) ) * ( I == J ? 1 : 2 );
      }

      return dEdF;
    }

  } // namespace Strain
  namespace DeformationGradient {
    template <>
    Eigen::Matrix3d make3D( const Eigen::Ref< const Eigen::Matrix< double, 1, 1 > >& tensor )
    {
      Marmot::Math::Matrix3d tensor3D = Marmot::Math::Matrix3d::Identity();
      tensor3D( 0, 0 )                = tensor( 0, 0 );
      return tensor3D;
    }

    template <>
    Eigen::Matrix3d make3D( const Eigen::Ref< const Eigen::Matrix< double, 2, 2 > >& tensor )
    {
      Marmot::Math::Matrix3d tensor3D = Marmot::Math::Matrix3d::Identity();
      tensor3D.topLeftCorner( 2, 2 )  = tensor;
      return tensor3D;
    }

    template <>
    Eigen::Matrix3d make3D( const Eigen::Ref< const Eigen::Matrix< double, 3, 3 > >& tensor )
    {
      return tensor;
    }
  } // namespace DeformationGradient

} // namespace Marmot::ContinuumMechanics::Kinematics
