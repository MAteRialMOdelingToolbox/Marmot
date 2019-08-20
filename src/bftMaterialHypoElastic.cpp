#include "bftFunctions.h"
#include "bftKinematics.h"
#include "bftMaterialHypoElastic.h"
#include "bftTensor.h"
#include "bftVoigt.h"
#include <iostream>

using namespace Eigen;

void BftMaterialHypoElastic::setCharacteristicElementLength( double length )

{
    characteristicElementLength = length;
}

void BftMaterialHypoElastic::computeStress( double*       stress_,
                                            double*       dStressDDDeformationGradient_,
                                            const double* FOld_,
                                            const double* FNew_,
                                            const double* timeOld,
                                            const double  dT,
                                            double&       pNewDT )
{
    // Standard implemenation of the Abaqus like Hughes-Winget algorithm
    // Approximation of the algorithmic tangent in order to
    // facilitate the dCauchy_dStrain tangent provided by
    // small strain material models
    //
    // strain increment is computed at midstep

    using namespace bft;
    using namespace bft::TensorUtility;
    using namespace kinematics::velocityGradient;

    const Map<const Matrix3d> FNew( FNew_ );
    const Map<const Matrix3d> FOld( FOld_ );
    bft::mVector6             stress( stress_ );

    Matrix3d FMidStep = 0.5 * ( FNew + FOld );

    Matrix3d l = ( FNew - FOld ) * FMidStep.inverse(); // actually l * dT

    Matrix3d dEps_  = 0.5 * ( l + l.transpose() ); // actually d * dT
    Matrix3d dOmega = 0.5 * ( l - l.transpose() ); // actually omega * dT

    Matrix3d dR = ( Matrix3d::Identity() - 0.5 * dOmega ).inverse() * ( Matrix3d::Identity() + 0.5 * dOmega );

    stress = Vgt::stressToVoigt( dR * Vgt::voigtToStress( stress ) * dR.transpose() );

    Matrix6 CJaumann;

    Vector6 dEps = Vgt::VoigtFromStrainMatrix( dEps_ );
    computeStress( stress.data(), CJaumann.data(), dEps.data(), timeOld, dT, pNewDT );

    Tensor633d                          dS_dl;
    Tensor633d                          dStressRotational_dl;
    Tensor633d                          dStressJaumann_dl;
    TensorMap<Eigen::Tensor<double, 3>> dS_dF( dStressDDDeformationGradient_, 6, 3, 3 );

    auto stressNew = Vgt::StressMatrixFromVoigt<3>( stress );

    dStressRotational_dl.setZero();
    for ( int ij = 0; ij < 6; ij++ ) {
        auto [i, j] = IndexNotation::fromVoigt<3>( ij );
        for ( int k = 0; k < 3; k++ )
            for ( int l = 0; l < 3; l++ )
                for ( int m = 0; m < 3; m++ )
                    dStressRotational_dl( ij, k, l ) += dOmega_dVelocityGradient( i, m, k, l ) * stressNew( m, j ) +
                                                        dOmega_dVelocityGradient( j, m, k, l ) * stressNew( i, m );
    }

    dStressJaumann_dl.setZero();
    for ( int ij = 0; ij < 6; ij++ )
        for ( int k = 0; k < 3; k++ )
            for ( int l = 0; l < 3; l++ )
                for ( int mn = 0; mn < 6; mn++ )
                    dStressJaumann_dl( ij, k, l ) += CJaumann( ij, mn ) * dStretchingRate_dVelocityGradient( mn, k, l );

    dS_dl     = dStressJaumann_dl + dStressRotational_dl;
    auto FInv = FNew.inverse();

    dS_dF.setZero();
    for ( int ij = 0; ij < 6; ij++ )
        for ( int k = 0; k < 3; k++ )
            for ( int l = 0; l < 3; l++ )
                for ( int m = 0; m < 3; m++ )
                    dS_dF( ij, k, l ) += dS_dl( ij, k, m ) * FInv( l, m );
}

void BftMaterialHypoElastic::computePlaneStress( double*       stress_,
                                                 double*       dStressDDStrain_,
                                                 double*       dStrain_,
                                                 const double* timeOld,
                                                 const double  dT,
                                                 double&       pNewDT )
{
    using namespace bft;

    Map<Vector6>  stress( stress_ );
    Map<Matrix6>  dStressDDStrain( dStressDDStrain_ );
    Map<Vector6>  dStrain( dStrain_ );
    Map<VectorXd> stateVars( this->stateVars, this->nStateVars );

    Vector6  stressTemp;
    VectorXd stateVarsOld = stateVars;
    Vector6  dStrainTemp  = dStrain;

    // assumption of isochoric deformation for initial guess
    dStrainTemp( 2 ) = ( -dStrain( 0 ) - dStrain( 1 ) );

    int planeStressCount = 1;
    while ( true ) {
        stressTemp = stress;
        stateVars  = stateVarsOld;

        computeStress( stressTemp.data(), dStressDDStrain.data(), dStrainTemp.data(), timeOld, dT, pNewDT );

        if ( pNewDT < 1.0 ) {
            return;
        }

        double residual = stressTemp.array().abs()[2];

        if ( residual < 1.e-10 || ( planeStressCount > 7 && residual < 1e-8 ) ) {
            break;
        }

        double tangentCompliance = 1. / dStressDDStrain( 2, 2 );
        if ( isNaN( tangentCompliance ) || std::abs( tangentCompliance ) > 1e10 )
            tangentCompliance = 1e10;

        dStrainTemp[2] -= tangentCompliance * stressTemp[2];

        planeStressCount += 1;
        if ( planeStressCount > 13 ) {
            pNewDT = 0.25;
            BftJournal::warningToMSG( "PlaneStressWrapper requires cutback" );
            return;
        }
    }

    dStrain = dStrainTemp;
    stress  = stressTemp;
}

void BftMaterialHypoElastic::computeUniaxialStress( double* stress_,
                                                    double* dStressDDStrain_,

                                                    double*       dStrain_,
                                                    const double* timeOld,
                                                    const double  dT,
                                                    double&       pNewDT )
{
    using namespace bft;

    Map<Vector6>  stress( stress_ );
    Map<Matrix6>  dStressDDStrain( dStressDDStrain_ );
    Map<Vector6>  dStrain( dStrain_ );
    Map<VectorXd> stateVars( this->stateVars, this->nStateVars );

    Vector6  stressTemp;
    VectorXd stateVarsOld = stateVars;
    Vector6  dStrainTemp  = dStrain;

    dStrainTemp( 2 ) = 0.0;
    dStrainTemp( 1 ) = 0.0;

    int count = 1;
    while ( true ) {
        stressTemp = stress;
        stateVars  = stateVarsOld;

        computeStress( stressTemp.data(), dStressDDStrain.data(), dStrainTemp.data(), timeOld, dT, pNewDT );

        if ( pNewDT < 1.0 ) {
            return;
        }

        double residual = stressTemp.array().abs()[1] + stressTemp.array().abs()[2];

        if ( residual < 1.e-10 || ( count > 7 && residual < 1e-8 ) ) {
            break;
        }

        dStrainTemp.segment<2>( 1 ) -= dStressDDStrain.block<2, 2>( 1, 1 ).colPivHouseholderQr().solve(
            stressTemp.segment<2>( 1 ) );

        count += 1;
        if ( count > 13 ) {
            pNewDT = 0.25;
            BftJournal::warningToMSG( "UniaxialStressWrapper requires cutback" );
            return;
        }
    }

    dStrain = dStrainTemp;
    stress  = stressTemp;
}
