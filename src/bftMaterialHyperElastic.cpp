#include "bftConstants.h"
#include "bftFunctions.h"
#include "bftKinematics.h"
#include "bftMaterialHyperElastic.h"
#include "bftTensor.h"
#include "bftVoigt.h"
#include <iostream>

using namespace Eigen;

void BftMaterialHyperElastic::computeStress( double*       Cauchy_,
                                             double*       dCauchy_d_F_np_,
                                             const double* F_n_,
                                             const double* F_np_,
                                             const double* timeOld_,
                                             const double  dT_,
                                             double&       pNewDT_ )
{
    using namespace bft;
    using namespace bft::TensorUtility::IndexNotation;

    Map<Vector6>              Cauchy( Cauchy_ );
    Map<Matrix<double, 6, 9>> dCauchy_d_F_np( dCauchy_d_F_np_ );
    const Map<const Matrix3d> F_np( F_np_ );
    Vector6                   E = kinematics::strain::GreenLagrange( F_np );

    Matrix6 dSdE;
    Vector6 S;

    computeStressPK2( S.data(), dSdE.data(), E.data(), timeOld_, dT_, pNewDT_ );

    double J = F_np.determinant();

    Matrix3d S_ = bft::Vgt::voigtToStress( S );

    Cauchy = bft::Vgt::stressToVoigt( 1. / J * F_np * S_ * F_np.transpose() );

    TensorMap<Tensor<double, 3>> dCauchydF( dCauchy_d_F_np_, 6, 3, 3 );
    dCauchydF.setZero();

    auto& F    = F_np;
    auto  FInv = F.inverse();

    Tensor633d dSdF;
    dSdF.setZero();
    Tensor633d dEdF = bft::kinematics::strain::dGreenLagrangedDeformationGradient( F_np );

    for ( int IJ = 0; IJ < 6; IJ++ )
        for ( int k = 0; k < 3; k++ )
            for ( int L = 0; L < 3; L++ )
                for ( int MN = 0; MN < 6; MN++ )
                    dSdF( IJ, k, L ) += dSdE( IJ, MN ) * dEdF( MN, k, L );

    auto I = Matrix3d::Identity();

    for ( int ij = 0; ij < 6; ij++ ) {
        auto [i, j] = fromVoigt<3>( ij );
        for ( int k = 0; k < 3; k++ )
            for ( int L = 0; L < 3; L++ ) {
                dCauchydF( ij, k, L ) -= FInv( L, k ) * Cauchy( ij );
                for ( int N = 0; N < 3; N++ ) {
                    dCauchydF( ij, k, L ) += 1. / J * +S_( L, N ) * ( F( j, N ) * I( i, k ) + F( i, N ) * I( j, k ) );
                    for ( int M = 0; M < 3; M++ )
                        dCauchydF( ij, k, L ) += 1. / J * F( i, M ) * dSdF( toVoigt<3>( M, N ), k, L ) * F( j, N );
                }
            }
    }
}

void BftMaterialHyperElastic::computePlaneStressPK2( double*       S,
                                                     double*       dSdE,
                                                     double*       E,
                                                     const double* timeOld,
                                                     const double  dT,
                                                     double&       pNewDT )

{
    // using namespace bft;

    // Map<Vector6>  stress( stress_ );
    // Map<Matrix6>  dStressDDStrain( dStressDDStrain_ );
    // Map<Vector6>  dStrain( dStrain_ );
    // Map<VectorXd> stateVars( this->stateVars, this->nStateVars );

    // Vector6  stressTemp;
    // VectorXd stateVarsOld = stateVars;
    // Vector6  dStrainTemp  = dStrain;

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
    // if ( isNaN( tangentCompliance ) || std::abs( tangentCompliance ) > 1e10 )
    // tangentCompliance = 1e10;

    // dStrainTemp[2] -= tangentCompliance * stressTemp[2];

    // planeStressCount += 1;
    // if ( planeStressCount > 13 ) {
    // pNewDT = 0.25;
    // BftJournal::warningToMSG( "PlaneStressWrapper requires cutback" );
    // return;
    //}
    //}

    // dStrain = dStrainTemp;
    // stress  = stressTemp;
}

void BftMaterialHyperElastic::computeUniaxialStressPK2( double*       S,
                                                        double*       dSdE,
                                                        double*       E,
                                                        const double* timeOld,
                                                        const double  dT,
                                                        double&       pNewDT )
{
    // using namespace bft;

    // Map<Vector6>  stress( stress_ );
    // Map<Matrix6>  dStressDDStrain( dStressDDStrain_ );
    // Map<Vector6>  dStrain( dStrain_ );
    // Map<VectorXd> stateVars( this->stateVars, this->nStateVars );

    // Vector6  stressTemp;
    // VectorXd stateVarsOld = stateVars;
    // Vector6  dStrainTemp  = dStrain;

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
    // BftJournal::warningToMSG( "UniaxialStressWrapper requires cutback" );
    // return;
    //}
    //}

    // dStrain = dStrainTemp;
    // stress  = stressTemp;
}
