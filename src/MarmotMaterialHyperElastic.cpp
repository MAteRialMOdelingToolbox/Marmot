#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotKinematics.h"
#include "Marmot/MarmotMaterialHyperElastic.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotVoigt.h"
#include <iostream>

using namespace Eigen;

void MarmotMaterialHyperElastic::computeStress( double*       Cauchy_,
                                                double*       dCauchy_d_F_np_,
                                                const double* F_n_,
                                                const double* F_np_,
                                                const double* timeOld_,
                                                const double  dT_,
                                                double&       pNewDT_ )
{
  using namespace Marmot;
  using namespace Marmot::ContinuumMechanics::TensorUtility::IndexNotation;

  Map< Vector6d >               Cauchy( Cauchy_ );
  Map< Matrix< double, 6, 9 > > dCauchy_d_F_np( dCauchy_d_F_np_ );
  const Map< const Matrix3d >   F_np( F_np_ );
  Vector6d                      E = ContinuumMechanics::Kinematics::Strain::GreenLagrange( F_np );

  Matrix6d dSdE;
  Vector6d S;

  computeStressPK2( S.data(), dSdE.data(), E.data(), timeOld_, dT_, pNewDT_ );

  double J = F_np.determinant();

  Matrix3d S_ = Marmot::ContinuumMechanics::VoigtNotation::voigtToStress( S );

  Cauchy = Marmot::ContinuumMechanics::VoigtNotation::stressToVoigt<double>( 1. / J * F_np * S_ * F_np.transpose() );

  TensorMap< Tensor< double, 3 > > dCauchydF( dCauchy_d_F_np_, 6, 3, 3 );
  dCauchydF.setZero();

  auto& F    = F_np;
  auto  FInv = F.inverse();

  EigenTensors::Tensor633d dSdF;
  dSdF.setZero();
  EigenTensors::Tensor633d dEdF = Marmot::ContinuumMechanics::Kinematics::Strain::dGreenLagrangedDeformationGradient(
    F_np );

  for ( int IJ = 0; IJ < 6; IJ++ )
    for ( int k = 0; k < 3; k++ )
      for ( int L = 0; L < 3; L++ )
        for ( int MN = 0; MN < 6; MN++ )
          dSdF( IJ, k, L ) += dSdE( IJ, MN ) * dEdF( MN, k, L );

  auto I = Matrix3d::Identity();

  for ( int ij = 0; ij < 6; ij++ ) {
    auto [i, j] = fromVoigt< 3 >( ij );
    for ( int k = 0; k < 3; k++ )
      for ( int L = 0; L < 3; L++ ) {
        dCauchydF( ij, k, L ) -= FInv( L, k ) * Cauchy( ij );
        for ( int N = 0; N < 3; N++ ) {
          dCauchydF( ij, k, L ) += 1. / J * +S_( L, N ) * ( F( j, N ) * I( i, k ) + F( i, N ) * I( j, k ) );
          for ( int M = 0; M < 3; M++ )
            dCauchydF( ij, k, L ) += 1. / J * F( i, M ) * dSdF( toVoigt< 3 >( M, N ), k, L ) * F( j, N );
        }
      }
  }
}

void MarmotMaterialHyperElastic::computePlaneStressPK2( double*       S2D,
                                                        double*       dSdE2D,
                                                        const double* E2D,
                                                        const double* timeOld,
                                                        const double  dT,
                                                        double&       pNewDT )

{
  throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << "not yet implemented" );
}

void MarmotMaterialHyperElastic::computeUniaxialStressPK2( double*       S1D,
                                                           double*       dSdE1D,
                                                           const double* E1D,
                                                           const double* timeOld,
                                                           const double  dT,
                                                           double&       pNewDT )
{
  throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << "not yet implemented" );
}
