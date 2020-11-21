#include "Marmot/HughesWinget.h"
#include "Marmot/MarmotKinematics.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotVoigt.h"

using namespace Eigen;

HughesWinget::HughesWinget( const Matrix3d& FOld, const Matrix3d& FNew, Formulation formulation )
    : theFormulation( formulation )
{

    Matrix3d FMidStep = 0.5 * ( FNew + FOld );

    l = ( FNew - FOld ) * FMidStep.inverse(); // actually l * dT

    Matrix3d dEps_ = 0.5 * ( l + l.transpose() );              // actually d * dT
    dOmega         = 0.5 * ( l - l.transpose() );              // actually omega * dT
    dEps           = Marmot::Vgt::VoigtFromStrainMatrix( dEps_ ); 
    dR             = ( Matrix3d::Identity() - 0.5 * dOmega ).inverse() * ( Matrix3d::Identity() + 0.5 * dOmega );
}

Marmot::Vector6 HughesWinget::getStrainIncrement()
{
    return dEps;
}

Matrix3d HughesWinget::getRotationIncrement()
{
    return dOmega;
}

Marmot::Vector6 HughesWinget::rotateTensor( const Marmot::Vector6& tensor )
{

    return Marmot::Vgt::stressToVoigt( dR * Marmot::Vgt::voigtToStress( tensor ) * dR.transpose() );
}

Marmot::Tensor633d HughesWinget::compute_dS_dF( const Marmot::Vector6& stress,
                                             const Matrix3d&     FInv,
                                             const Marmot::Matrix6& dChauchydEps )
{
    using namespace Marmot;
    using namespace Marmot::TensorUtility;
    using namespace Marmot::kinematics::velocityGradient;

    Tensor633d dS_dl;
    Tensor633d dS_dF;
    Tensor633d dStressRotational_dl;
    Tensor633d dStressJaumann_dl;
    auto       stressNew = Vgt::StressMatrixFromVoigt<3>( stress );

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
                    dStressJaumann_dl( ij, k, l ) += dChauchydEps( ij, mn ) *
                                                     dStretchingRate_dVelocityGradient( mn, k, l );

    dS_dl     = dStressJaumann_dl + dStressRotational_dl;

    dS_dF.setZero();
    for ( int ij = 0; ij < 6; ij++ )
        for ( int k = 0; k < 3; k++ )
            for ( int l = 0; l < 3; l++ )
                for ( int m = 0; m < 3; m++ )
                    dS_dF( ij, k, l ) += dS_dl( ij, k, m ) * FInv( l, m );

    return dS_dF;
}

Eigen::Matrix3d HughesWinget::compute_dScalar_dF(const Eigen::Matrix3d& FInv, const Marmot::Vector6& dScalarDEps){

    using namespace Marmot::kinematics::velocityGradient;
    Matrix3d dScalar_dl = Matrix3d::Zero();
    for ( int k = 0; k < 3; k++ )
        for ( int l = 0; l < 3; l++ )
            for ( int ij = 0; ij < 6; ij++ )
                dScalar_dl( k, l ) += dScalarDEps(ij) * dStretchingRate_dVelocityGradient(ij, k, l);
    
    return dScalar_dl * FInv;

}
