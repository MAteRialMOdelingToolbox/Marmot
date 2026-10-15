#include "Marmot/MarmotLowerDimensionalStress.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotVoigt.h"

using namespace Eigen;

namespace Marmot::ContinuumMechanics {

  namespace UniaxialStress {
    double getUniaxialStressTangent( const Ref< const Marmot::Math::Matrix6d >& C )
    {
      Marmot::Math::Vector6d b;
      b << 1, 0, 0, 0, 0, 0;
      Marmot::Math::Matrix6d A = C;
      A.row( 0 ) << 1, 0, 0, 0, 0, 0;
      Marmot::Math::Vector6d dEdEUniaxial = A.colPivHouseholderQr().solve( b );
      return C.row( 0 ) * dEdEUniaxial;
    }

  } // namespace UniaxialStress

  namespace PlaneStrain {

    Marmot::Math::Tensors::Tensor322d reduce3D_dStress_dDeformationGradient(
      const Marmot::Math::Tensors::Tensor633d& dStressdDeformationGradient3D )
    {
      static constexpr int     planeVoigtIndices[] = { 0, 1, 3 };
      Marmot::Math::Tensors::Tensor322d tangent2D;
      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 2; j++ )
          for ( int k = 0; k < 2; k++ )
            tangent2D( i, j, k ) = dStressdDeformationGradient3D( planeVoigtIndices[i], j, k );

      return tangent2D;
    }

    Marmot::Math::Matrix3d getPlaneStrainTangent( const Marmot::Math::Matrix6d& C )
    {
      Marmot::Math::Matrix3d CPlaneStrain              = Marmot::Math::Matrix3d::Zero();
      CPlaneStrain.topLeftCorner( 2, 2 ) = C.topLeftCorner( 2, 2 );
      CPlaneStrain( 2, 2 )               = C( 3, 3 );
      CPlaneStrain.block< 1, 2 >( 2, 0 ) = C.block< 1, 2 >( 3, 0 );
      CPlaneStrain.block< 2, 1 >( 0, 2 ) = C.block< 2, 1 >( 0, 3 );

      return CPlaneStrain;
    }

    Matrix< double, 6, 3 > dStrainDStrainPlaneStrain()
    {
      Matrix< double, 6, 3 > T = Matrix< double, 6, 3 >::Zero();
      T( 0, 0 )                = 1;
      T( 1, 1 )                = 1;
      T( 3, 2 )                = 1;
      return T;
    }

  } // namespace PlaneStrain

  namespace PlaneStress {

    Marmot::Math::Tensors::Tensor322d compute_dStress_dDeformationGradient( const Marmot::Math::Tensors::Tensor633d& dS_dF_3D )
    {
      /*  dS^PS    dS^PS    dS    dF^Comp
       *  ----- == ----- * ---- * -------
       *  dF^PS    dS       dF    dF^PS
       *
       *  dF^Comp_ij = dF_ij^PS       if  ij =/= 33
       *             = dF_^Comp_33    else
       *
       *             with dS_33 = 0 =  dS33/dFij * dFij^PS + dS33/dF33 * dF^Comp_33
       *
       *             --> dF^Comp_33 = -dS33/dFij * dFij^PS * dF33/dS33
       *
       *             --> dF^Comp
       *                 -------
       *                  dF^PS
       * */

      // projection to plane
      Marmot::Math::Tensors::Tensor322d dS_dF = PlaneStrain::reduce3D_dStress_dDeformationGradient( dS_dF_3D );

      using namespace Marmot::Math::Tensors::IndexNotation;

      // clang-format off
                for ( int m = 0; m < 2; m ++ )
                    for ( int n = 0; n < 2; n ++ )
                        for ( int k = 0; k < 2; k ++ )
                            for ( int l = 0; l < 2; l ++ )
                                dS_dF(              toVoigt<2> (m,n), k, l )
                                    =  -  dS_dF_3D( toVoigt<3> (m,n), 2, 2 )
                                    * 1./ dS_dF_3D( toVoigt<3> (2,2), 2, 2 )
                                    *     dS_dF_3D( toVoigt<3> (2,2), k, l );
      // clang-format on
      return dS_dF;
    }

    Matrix< double, 3, 6 > dStressPlaneStressDStress()
    {
      Matrix< double, 3, 6 > T = Matrix< double, 3, 6 >::Zero();
      T( 0, 0 )                = 1;
      T( 1, 1 )                = 1;
      T( 2, 3 )                = 1;
      return T;
    }

    Marmot::Math::Matrix3d getPlaneStressTangent( const Marmot::Math::Matrix6d& C )
    {
      return dStressPlaneStressDStress() * C * dStrainDStrainPlaneStress( C );
    }

    Marmot::Math::Vector6d planeStressCompensationStrain( const Marmot::Math::Vector6d& strain, double nu )
    {
      const Marmot::Math::Vector6d& e                = strain;
      const double    strainCorrComp33 = -nu / ( 1 - nu ) * ( e( 0 ) + e( 1 ) ) - e( 2 );
      Marmot::Math::Vector6d        result;
      result << 0, 0, strainCorrComp33, 0, 0, 0;
      return result;
    }

    Marmot::Math::Matrix6d planeStressTangentTransformationMatrix( const Marmot::Math::Matrix6d& tangent )
    {
      /* Returns the transformation Matrix T which fullfills
       * planeStressIncrement = C : T * strainIncrement
       * for isotropic material behavior only!
       * */
      Marmot::Math::Matrix6d     T   = Marmot::Math::Matrix6d::Identity();
      const double t31 = -tangent( 2, 0 ) / tangent( 2, 2 );
      const double t32 = -tangent( 2, 1 ) / tangent( 2, 2 );
      T.row( 2 ).head( 3 ) << t31, t32, 0.0;
      return T;
    }

    Matrix< double, 6, 3 > dStrainDStrainPlaneStress( const Marmot::Math::Matrix6d& tangent )
    {
      Matrix< double, 6, 3 > T = Matrix< double, 6, 3 >::Zero();
      T( 0, 0 )                = 1;
      T( 1, 1 )                = 1;
      T( 3, 2 )                = 1;
      T( 2, 0 )                = -tangent( 2, 0 ) / tangent( 2, 2 );
      T( 2, 1 )                = -tangent( 2, 1 ) / tangent( 2, 2 );
      T( 2, 2 )                = -tangent( 2, 3 ) / tangent( 2, 2 );
      return T;
    }

  } // namespace PlaneStress
} // namespace Marmot::ContinuumMechanics
