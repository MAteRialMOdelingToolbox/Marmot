#include "Marmot/MarmotMaterialHyperElastic.h"
#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotKinematics.h"
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
  Vector6d                      E = ContinuumMechanics::Kinematics::strain::GreenLagrange( F_np );

  Matrix6d dSdE;
  Vector6d S;

  computeStressPK2( S.data(), dSdE.data(), E.data(), timeOld_, dT_, pNewDT_ );

  double J = F_np.determinant();

  Matrix3d S_ = Marmot::ContinuumMechanics::VoigtNotation::voigtToStress( S );

  Cauchy = Marmot::ContinuumMechanics::VoigtNotation::stressToVoigt( 1. / J * F_np * S_ * F_np.transpose() );

  TensorMap< Tensor< double, 3 > > dCauchydF( dCauchy_d_F_np_, 6, 3, 3 );
  dCauchydF.setZero();

  auto& F    = F_np;
  auto  FInv = F.inverse();

  EigenTensors::Tensor633d dSdF;
  dSdF.setZero();
  EigenTensors::Tensor633d dEdF = Marmot::ContinuumMechanics::Kinematics::strain::dGreenLagrangedDeformationGradient(
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

void MarmotMaterialHyperElastic::computePlaneStressPK2( double*       S,
                                                        double*       dSdE,
                                                        double*       E,
                                                        const double* timeOld,
                                                        const double  dT,
                                                        double&       pNewDT )

{
  // using namespace Marmot;

  // Map<Vector6d>  stress( stress_ );
  // Map<Matrix6d>  dStressDDStrain( dStressDDStrain_ );
  // Map<Vector6d>  dStrain( dStrain_ );
  // Map<VectorXd> stateVars( this->stateVars, this->nStateVars );

  // Vector6d  stressTemp;
  // VectorXd stateVarsOld = stateVars;
  // Vector6d  dStrainTemp  = dStrain;

  //// assumption of isochoric deformation for initial guess
  // dStrainTemp( 2 ) = ( -dStrain( 0 ) - dStrain( 1 ) );

  // int planeStressCount = 1;
  // while ( true ) {
  // stressTemp = stress;
  // stateVars  = stateVarsOld;

  // if ( pNewDT < 1.0 ) {
  // return;
  //}

  // double residual = stressTemp.array().abs()[2];

  // if ( residual < 1.e-10 || ( planeStressCount > 7 && residual < 1e-8 ) ) {
  // break;
  //}

  // double tangentCompliance = 1. / dStressDDStrain( 2, 2 );
  // if ( Math::isNaN( tangentCompliance ) || std::abs( tangentCompliance ) > 1e10 )
  // tangentCompliance = 1e10;

  // dStrainTemp[2] -= tangentCompliance * stressTemp[2];

  // planeStressCount += 1;
  // if ( planeStressCount > 13 ) {
  // pNewDT = 0.25;
  // MarmotJournal::warningToMSG( "PlaneStressWrapper requires cutback" );
  // return;
  //}
  //}

  // dStrain = dStrainTemp;
  // stress  = stressTemp;
}

void MarmotMaterialHyperElastic::computeUniaxialStressPK2( double*       S,
                                                           double*       dSdE,
                                                           double*       E,
                                                           const double* timeOld,
                                                           const double  dT,
                                                           double&       pNewDT )
{
  // using namespace Marmot;

  // Map<Vector6d>  stress( stress_ );
  // Map<Matrix6d>  dStressDDStrain( dStressDDStrain_ );
  // Map<Vector6d>  dStrain( dStrain_ );
  // Map<VectorXd> stateVars( this->stateVars, this->nStateVars );

  // Vector6d  stressTemp;
  // VectorXd stateVarsOld = stateVars;
  // Vector6d  dStrainTemp  = dStrain;

  // dStrainTemp( 2 ) = 0.0;
  // dStrainTemp( 1 ) = 0.0;

  // int count = 1;
  // while ( true ) {
  // stressTemp = stress;
  // stateVars  = stateVarsOld;

  // if ( pNewDT < 1.0 ) {
  // return;
  //}

  // double residual = stressTemp.array().abs()[1] + stressTemp.array().abs()[2];

  // if ( residual < 1.e-10 || ( count > 7 && residual < 1e-8 ) ) {
  // break;
  //}

  // dStrainTemp.segment<2>( 1 ) -= dStressDDStrain.block<2, 2>( 1, 1 ).colPivHouseholderQr().solve(
  // stressTemp.segment<2>( 1 ) );

  // count += 1;
  // if ( count > 13 ) {
  // pNewDT = 0.25;
  // MarmotJournal::warningToMSG( "UniaxialStressWrapper requires cutback" );
  // return;
  //}
  //}

  // dStrain = dStrainTemp;
  // stress  = stressTemp;
}
