#include "Marmot/HughesWinget.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotKinematics.h"
#include "Marmot/MarmotLowerDimensionalStress.h"
#include "Marmot/MarmotMaterialGradientEnhancedHypoElastic.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotVoigt.h"

using namespace Eigen;

void MarmotMaterialGradientEnhancedHypoElastic::computeStress( double*       stress_,
                                                               double&       K_local,
                                                               double&       nonLocalRadius,
                                                               double*       dStressDDDeformationGradient_,
                                                               double*       dK_localDDeformationGradient_,
                                                               double*       dStressDK,
                                                               const double* FOld_,
                                                               const double* FNew_,
                                                               const double  KOld,
                                                               const double  dK,
                                                               const double* timeOld,
                                                               const double  dT,
                                                               double&       pNewDT )
{
  // Standard implemenation of the Abaqus like Hughes-Winget algorithm
  // Approximation of the algorithmic tangent in order to
  // facilitate the dCauchy_dStrain tangent provided by
  // small strain material models

  using namespace Marmot;
  const Map< const Matrix3d > FNew( FNew_ );
  const Map< const Matrix3d > FOld( FOld_ );
  Marmot::mVector6d           stress( stress_ );

  Matrix6d CJaumann                = Matrix6d::Zero();
  Vector6d dK_LocalDStretchingRate = Vector6d::Zero();

  Marmot::NumericalAlgorithms::HughesWinget
    hughesWingetIntegrator( FOld, FNew, Marmot::NumericalAlgorithms::HughesWinget::Formulation::AbaqusLike );

  auto dEps = hughesWingetIntegrator.getStrainIncrement();
  stress    = hughesWingetIntegrator.rotateTensor( stress );

  computeStress( stress.data(),
                 K_local,
                 nonLocalRadius,
                 CJaumann.data(),
                 dK_LocalDStretchingRate.data(),
                 dStressDK,
                 dEps.data(),
                 KOld,
                 dK,
                 timeOld,
                 dT,
                 pNewDT );

  TensorMap< Eigen::Tensor< double, 3 > > dS_dF( dStressDDDeformationGradient_, 6, 3, 3 );
  Map< Matrix3d >                         dKLocal_dF( dK_localDDeformationGradient_ );

  Matrix3d FInv = FNew.inverse();
  dS_dF         = hughesWingetIntegrator.compute_dS_dF( stress, FInv, CJaumann );
  dKLocal_dF    = hughesWingetIntegrator.compute_dScalar_dF( FInv, dK_LocalDStretchingRate );
}

void MarmotMaterialGradientEnhancedHypoElastic::computePlaneStress( double*       stress2D_,
                                                                    double&       KLocal,
                                                                    double&       nonLocalRadius,
                                                                    double*       dStress_dStrain2D_,
                                                                    double*       dKLocal_dStrain2D_,
                                                                    double*       dStress_dK2D_,
                                                                    const double* dStrain2D_,
                                                                    double        KOld,
                                                                    double        dK,
                                                                    const double* timeOld,
                                                                    const double  dT,
                                                                    double&       pNewDT )
{
  using namespace Marmot;
  using namespace ContinuumMechanics::VoigtNotation;

  Map< const Matrix< double, 3, 1 > > dStrain2D( dStrain2D_ );
  Map< Matrix< double, 3, 1 > >       stress2D( stress2D_ );
  Map< Matrix< double, 3, 3 > >       dStress_dStrain2D( dStress_dStrain2D_ );
  Map< Matrix< double, 3, 1 > >       dStress_dK2D( dStress_dK2D_ );
  Map< Matrix< double, 3, 1 > >       dKLocal_dStrain2D( dKLocal_dStrain2D_ );
  Map< VectorXd >                     stateVars( this->stateVars, this->nStateVars );

  Matrix6d dStress_dStrain3D;
  Vector6d dKLocal_dStrain3D;
  Vector6d dStress_dK3D;

  Vector6d stressTemp3D;
  VectorXd stateVarsOld  = stateVars;
  Vector6d dStrain3DTemp = Marmot::ContinuumMechanics::VoigtNotation::make3DVoigt< VoigtSize::TwoD >( dStrain2D );

  // assumption of isochoric deformation for initial guess
  dStrain3DTemp( 2 ) = ( -dStrain2D( 0 ) - dStrain2D( 1 ) );

  int planeStressCount = 1;
  while ( true ) {
    stressTemp3D = Marmot::ContinuumMechanics::VoigtNotation::make3DVoigt< VoigtSize::TwoD >( stress2D );
    stateVars    = stateVarsOld;

    computeStress( stressTemp3D.data(),
                   KLocal,
                   nonLocalRadius,
                   dStress_dStrain3D.data(),
                   dKLocal_dStrain3D.data(),
                   dStress_dK3D.data(),
                   dStrain3DTemp.data(),
                   KOld,
                   dK,
                   timeOld,
                   dT,
                   pNewDT );

    if ( pNewDT < 1.0 ) {
      return;
    }

    double residual = stressTemp3D.array().abs()[2];

    if ( residual < 1.e-10 || ( planeStressCount > 7 && residual < 1e-5 ) ) {
      break;
    }

    double tangentCompliance = 1. / dStress_dStrain3D( 2, 2 );
    if ( Math::isNaN( tangentCompliance ) || std::abs( tangentCompliance ) > 1e10 )
      tangentCompliance = 1e10;

    dStrain3DTemp[2] -= tangentCompliance * stressTemp3D[2];

    planeStressCount += 1;
    if ( planeStressCount > 10 ) {
      pNewDT = 0.25;
      MarmotJournal::warningToMSG( "PlaneStressWrapper requires cutback" );
      return;
    }
  }

  stress2D = ContinuumMechanics::VoigtNotation::reduce3DVoigt< VoigtSize::TwoD >( stressTemp3D );

  dStress_dStrain2D = ContinuumMechanics::PlaneStress::getPlaneStressTangent( dStress_dStrain3D );
  dKLocal_dStrain2D = ContinuumMechanics::VoigtNotation::reduce3DVoigt< VoigtSize::TwoD >( dKLocal_dStrain3D );
  dStress_dK2D      = ContinuumMechanics::VoigtNotation::reduce3DVoigt< VoigtSize::TwoD >( dStress_dK3D );
}
