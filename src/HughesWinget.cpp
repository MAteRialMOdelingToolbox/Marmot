#include "Marmot/HughesWinget.h"
#include "Marmot/MarmotKinematics.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotVoigt.h"

namespace Marmot::NumericalAlgorithms {
  using namespace Eigen;
  HughesWinget::HughesWinget( const Matrix3d& FOld, const Matrix3d& FNew, Formulation formulation )
    : theFormulation( formulation )
  {

    Matrix3d FMidStep = 0.5 * ( FNew + FOld );

    l = ( FNew - FOld ) * FMidStep.inverse(); // actually l * dT

    Matrix3d dEps_ = 0.5 * ( l + l.transpose() ); // actually d * dT
    dOmega         = 0.5 * ( l - l.transpose() ); // actually omega * dT
    dEps           = Marmot::ContinuumMechanics::VoigtNotation::voigtFromStrainMatrix( dEps_ );
    dR             = ( Matrix3d::Identity() - 0.5 * dOmega ).inverse() * ( Matrix3d::Identity() + 0.5 * dOmega );
  }

  Marmot::Vector6d HughesWinget::getStrainIncrement() { return dEps; }

  Matrix3d HughesWinget::getRotationIncrement() { return dOmega; }

  Marmot::Vector6d HughesWinget::rotateTensor( const Marmot::Vector6d& tensor )
  {

    return Marmot::ContinuumMechanics::VoigtNotation::stressToVoigt(
      dR * Marmot::ContinuumMechanics::VoigtNotation::voigtToStress( tensor ) * dR.transpose() );
  }

  Marmot::EigenTensors::Tensor633d HughesWinget::compute_dS_dF( const Marmot::Vector6d& stress,
                                                                const Matrix3d&         FInv,
                                                                const Marmot::Matrix6d& dChauchydEps )
  {
    using namespace Marmot;
    using namespace Marmot::ContinuumMechanics::TensorUtility;
    using namespace Marmot::ContinuumMechanics::Kinematics::velocityGradient;

    EigenTensors::Tensor633d dS_dl;
    EigenTensors::Tensor633d dS_dF;
    EigenTensors::Tensor633d dStressRotational_dl;
    EigenTensors::Tensor633d dStressJaumann_dl;
    auto                     stressNew = ContinuumMechanics::VoigtNotation::stressMatrixFromVoigt< 3 >( stress );

    dStressRotational_dl.setZero();
    for ( int ij = 0; ij < 6; ij++ ) {
      auto [i, j] = IndexNotation::fromVoigt< 3 >( ij );
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
            dStressJaumann_dl( ij, k, l ) += dChauchydEps( ij, mn ) * dStretchingRate_dVelocityGradient( mn, k, l );

    dS_dl = dStressJaumann_dl + dStressRotational_dl;

    dS_dF.setZero();
    for ( int ij = 0; ij < 6; ij++ )
      for ( int k = 0; k < 3; k++ )
        for ( int l = 0; l < 3; l++ )
          for ( int m = 0; m < 3; m++ )
            dS_dF( ij, k, l ) += dS_dl( ij, k, m ) * FInv( l, m );

    return dS_dF;
  }

  Eigen::Matrix3d HughesWinget::compute_dScalar_dF( const Eigen::Matrix3d& FInv, const Marmot::Vector6d& dScalarDEps )
  {

    using namespace Marmot::ContinuumMechanics::Kinematics::velocityGradient;
    Matrix3d dScalar_dl = Matrix3d::Zero();
    for ( int k = 0; k < 3; k++ )
      for ( int l = 0; l < 3; l++ )
        for ( int ij = 0; ij < 6; ij++ )
          dScalar_dl( k, l ) += dScalarDEps( ij ) * dStretchingRate_dVelocityGradient( ij, k, l );

    return dScalar_dl * FInv;
  }
} // namespace Marmot::NumericalAlgorithms
