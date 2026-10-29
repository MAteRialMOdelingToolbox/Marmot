#include "Marmot/HughesWinget.h"
#include "Marmot/MarmotKinematics.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotVoigt.h"

namespace Marmot::ContinuumMechanics::Algorithms {
  using namespace Eigen;
  HughesWinget::HughesWinget( const Marmot::Math::Matrix3d& FOld,
                              const Marmot::Math::Matrix3d& FNew,
                              Formulation                   formulation )
    : theFormulation( formulation )
  {

    Marmot::Math::Matrix3d FMidStep = 0.5 * ( FNew + FOld );

    l = ( FNew - FOld ) * FMidStep.inverse();                   // actually l * dT

    Marmot::Math::Matrix3d dEps_ = 0.5 * ( l + l.transpose() ); // actually d * dT
    dOmega                       = 0.5 * ( l - l.transpose() ); // actually omega * dT
    dEps                         = Marmot::ContinuumMechanics::Voigt::voigtFromStrainMatrix( dEps_ );
    dR                           = ( Marmot::Math::Matrix3d::Identity() - 0.5 * dOmega ).inverse() *
         ( Marmot::Math::Matrix3d::Identity() + 0.5 * dOmega );
  }

  Marmot::Math::Vector6d HughesWinget::getStrainIncrement()
  {
    return dEps;
  }

  Marmot::Math::Matrix3d HughesWinget::getRotationIncrement()
  {
    return dR;
  }

  Marmot::Math::Vector6d HughesWinget::rotateTensor( const Marmot::Math::Vector6d& tensor )
  {

    return Marmot::ContinuumMechanics::Voigt::stressToVoigt< double >(
      dR * Marmot::ContinuumMechanics::Voigt::voigtToStress( tensor ) * dR.transpose() );
  }

  Marmot::Math::Tensors::Tensor633d HughesWinget::compute_dS_dF( const Marmot::Math::Vector6d& stress,
                                                                 const Marmot::Math::Matrix3d& FInv,
                                                                 const Marmot::Math::Matrix6d& dChauchydEps )
  {
    using namespace Marmot;
    using namespace Marmot::Math::Tensors;
    using namespace Marmot::ContinuumMechanics::Kinematics::VelocityGradient;

    Marmot::Math::Tensors::Tensor633d dS_dl;
    Marmot::Math::Tensors::Tensor633d dS_dF;
    Marmot::Math::Tensors::Tensor633d dStressRotational_dl;
    Marmot::Math::Tensors::Tensor633d dStressJaumann_dl;
    auto                              stressNew = ContinuumMechanics::Voigt::stressMatrixFromVoigt< 3 >( stress );

    dStressRotational_dl.setZero();
    for ( int ij = 0; ij < 6; ij++ ) {
      auto [i, j] = IndexNotation::fromVoigt< 3 >( ij );
      for ( int k = 0; k < 3; k++ )
        for ( int ll = 0; ll < 3; ll++ )
          for ( int m = 0; m < 3; m++ )
            dStressRotational_dl( ij, k, ll ) += dOmega_dVelocityGradient( i, m, k, ll ) * stressNew( m, j ) +
                                                 dOmega_dVelocityGradient( j, m, k, ll ) * stressNew( i, m );
    }

    dStressJaumann_dl.setZero();
    for ( int ij = 0; ij < 6; ij++ )
      for ( int k = 0; k < 3; k++ )
        for ( int ll = 0; ll < 3; ll++ )
          for ( int mn = 0; mn < 6; mn++ )
            dStressJaumann_dl( ij, k, ll ) += dChauchydEps( ij, mn ) * dStretchingRate_dVelocityGradient( mn, k, ll );

    dS_dl = dStressJaumann_dl + dStressRotational_dl;

    dS_dF.setZero();
    for ( int ij = 0; ij < 6; ij++ )
      for ( int k = 0; k < 3; k++ )
        for ( int ll = 0; ll < 3; ll++ )
          for ( int m = 0; m < 3; m++ )
            dS_dF( ij, k, ll ) += dS_dl( ij, k, m ) * FInv( ll, m );

    return dS_dF;
  }

  Eigen::Matrix3d HughesWinget::compute_dScalar_dF( const Eigen::Matrix3d&        FInv,
                                                    const Marmot::Math::Vector6d& dScalarDEps )
  {

    using namespace Marmot::ContinuumMechanics::Kinematics::VelocityGradient;
    Marmot::Math::Matrix3d dScalar_dl = Marmot::Math::Matrix3d::Zero();
    for ( int k = 0; k < 3; k++ )
      for ( int ll = 0; ll < 3; ll++ )
        for ( int ij = 0; ij < 6; ij++ )
          dScalar_dl( k, ll ) += dScalarDEps( ij ) * dStretchingRate_dVelocityGradient( ij, k, ll );

    return dScalar_dl * FInv;
  }
} // namespace Marmot::ContinuumMechanics::Algorithms
