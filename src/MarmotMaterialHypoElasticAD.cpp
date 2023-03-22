#include "Marmot/MarmotMaterialHypoElasticAD.h"
#include "Marmot/HughesWinget.h"
#include "Marmot/MarmotAutomaticDifferentiation.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotKinematics.h"
#include "Marmot/MarmotLowerDimensionalStress.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include <iostream>

using namespace Eigen;
using namespace autodiff;

void MarmotMaterialHypoElasticAD::setCharacteristicElementLength( double length )
{
  characteristicElementLength = length;
}

void MarmotMaterialHypoElasticAD::computeStress( double*       stress_,
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

  using namespace Marmot;
  using namespace Marmot::ContinuumMechanics::TensorUtility;
  using namespace ContinuumMechanics::Kinematics::VelocityGradient;

  const Map< const Matrix3d > FOld( FOld_ );
  const Map< const Matrix3d > FNew( FNew_ );
  mVector6d                   stress( stress_ );
  Map< VectorXd >             stateVars( this->stateVars, this->nStateVars );

  Marmot::NumericalAlgorithms::HughesWinget
    hughesWingetIntegrator( FOld, FNew, Marmot::NumericalAlgorithms::HughesWinget::Formulation::AbaqusLike );

  auto dEps = hughesWingetIntegrator.getStrainIncrement();

  stress = hughesWingetIntegrator.rotateTensor( stress );

  // remember old state
  const VectorXd stateVarsOld = stateVars;
  const Vector6d stressOld    = stress;

  Matrix6d CJaumann;
  // ----------------------------------------
  // autodiff part
  // ----------------------------------------
  // compute stress and tangent with autodiff
  std::tie( stress, CJaumann ) = Marmot::AutomaticDifferentiation::dF_dX(
    [&]( const autodiff::VectorXdual& dE_ ) {
      // reset stateVars to old state
      stateVars = stateVarsOld;

      // set up dual vectors
      autodiff::VectorXdual s( stressOld );

      // compute stress
      computeStress( s.data(), dE_.data(), timeOld, dT, pNewDT );

      return s;
    },
    dEps );
  // ----------------------------------------

  TensorMap< Eigen::Tensor< double, 3 > > dS_dF( dStressDDDeformationGradient_, 6, 3, 3 );

  dS_dF = hughesWingetIntegrator.compute_dS_dF( stress, FNew.inverse(), CJaumann );
}

void MarmotMaterialHypoElasticAD::computePlaneStress( double*       stress2D_,
                                                      double*       dStress_dStrain2D_,
                                                      const double* dStrain2D_,
                                                      const double* timeOld,
                                                      const double  dT,
                                                      double&       pNewDT )
{
  using namespace Marmot;
  using namespace ContinuumMechanics::VoigtNotation;

  Map< const Matrix< double, 3, 1 > > dStrain2D( dStrain2D_ );
  Map< Matrix< double, 3, 1 > >       stress2D( stress2D_ );
  Map< Matrix< double, 3, 3 > >       dStress_dStrain2D( dStress_dStrain2D_ );
  Map< VectorXd >                     stateVars( this->stateVars, this->nStateVars );

  Matrix6d dStress_dStrain3D;

  Vector6d       stress3DTemp;
  const VectorXd stateVarsOld = stateVars;

  Vector6d dStrain3DTemp = Marmot::ContinuumMechanics::VoigtNotation::make3DVoigt< VoigtSize::TwoD >( dStrain2D );

  // assumption of isochoric deformation for initial guess
  dStrain3DTemp( 2 ) = ( -dStrain2D( 0 ) - dStrain2D( 1 ) );

  const Vector6d stress3DOld = Marmot::ContinuumMechanics::VoigtNotation::make3DVoigt< VoigtSize::TwoD >( stress2D );

  int planeStressCount = 1;
  while ( true ) {

    // ----------------------------------------
    // autodiff part
    // ----------------------------------------
    // compute stress and tangent with autodiff
    std::tie( stress3DTemp, dStress_dStrain3D ) = Marmot::AutomaticDifferentiation::dF_dX(
      [&]( const autodiff::VectorXdual& dE_ ) {
        // reset stateVars to old state
        stateVars = stateVarsOld;

        // std::cout << "stateVars = " << stateVars << std::endl;
        // set up dual vectors
        autodiff::VectorXdual s( stress3DOld );

        // compute stress
        computeStress( s.data(), dE_.data(), timeOld, dT, pNewDT );

        return s;
      },
      dStrain3DTemp );
    // ----------------------------------------

    double residual = stress3DTemp.array().abs()[2];

    if ( residual < 1.e-10 || ( planeStressCount > 7 && residual < 1e-8 ) ) {
      break;
    }

    double tangentCompliance = 1. / dStress_dStrain3D( 2, 2 );
    if ( Math::isNaN( tangentCompliance ) || std::abs( tangentCompliance ) > 1e10 )
      tangentCompliance = 1e10;

    dStrain3DTemp( 2 ) -= tangentCompliance * stress3DTemp( 2 );

    planeStressCount += 1;
    if ( planeStressCount > 13 ) {
      pNewDT = 0.25;
      MarmotJournal::warningToMSG( "PlaneStressWrapper requires cutback" );
      return;
    }
  }

  stress2D          = ContinuumMechanics::VoigtNotation::reduce3DVoigt< VoigtSize::TwoD >( stress3DTemp );
  dStress_dStrain2D = ContinuumMechanics::PlaneStress::getPlaneStressTangent( dStress_dStrain3D );
}

void MarmotMaterialHypoElasticAD::computeUniaxialStress( double* stress1D_,
                                                         double* dStress_dStrain1D_,

                                                         const double* dStrain1D_,
                                                         const double* timeOld,
                                                         const double  dT,
                                                         double&       pNewDT )
{
  using namespace Marmot;
  using namespace ContinuumMechanics::VoigtNotation;

  Map< const Matrix< double, 1, 1 > > dStrain1D( dStrain1D_ );
  Map< Matrix< double, 1, 1 > >       stress1D( stress1D_ );
  Map< VectorXd >                     stateVars( this->stateVars, this->nStateVars );

  Matrix6d dStress_dStrain3D;

  Vector6d dStrain3DTemp = Marmot::ContinuumMechanics::VoigtNotation::make3DVoigt< VoigtSize::OneD >( dStrain1D );
  Vector6d stress3DTemp;

  const VectorXd stateVarsOld = stateVars;
  const VectorXd stress3DOld  = Marmot::ContinuumMechanics::VoigtNotation::make3DVoigt< VoigtSize::OneD >( stress1D );
  ;

  int count = 1;
  while ( true ) {
    // ----------------------------------------
    // autodiff part
    // ----------------------------------------
    // compute stress and tangent with autodiff
    std::tie( stress3DTemp, dStress_dStrain3D ) = Marmot::AutomaticDifferentiation::dF_dX(
      [&]( const autodiff::VectorXdual& dE_ ) {
        // reset stateVars to old state
        stateVars = stateVarsOld;

        // set up dual vectors
        autodiff::VectorXdual s( stress3DOld );

        // compute stress
        computeStress( s.data(), dE_.data(), timeOld, dT, pNewDT );

        return s;
      },
      dStrain3DTemp );
    // ----------------------------------------

    if ( pNewDT < 1.0 ) {
      return;
    }

    if ( pNewDT < 1.0 ) {
      return;
    }

    const double residual = stress3DTemp.array().abs().segment( 1, 2 ).sum();

    if ( residual < 1.e-13 || ( count > 7 && residual < 1e-10 ) ) {
      break;
    }

    dStrain3DTemp.segment< 2 >( 1 ) -= dStress_dStrain3D.block< 2, 2 >( 1, 1 ).colPivHouseholderQr().solve(
      stress3DTemp.segment< 2 >( 1 ) );

    count += 1;
    if ( count > 13 ) {
      pNewDT = 0.25;
      MarmotJournal::warningToMSG( "UniaxialStressWrapper requires cutback" );
      return;
    }
  }

  stress1D              = ContinuumMechanics::VoigtNotation::reduce3DVoigt< VoigtSize::OneD >( stress3DTemp );
  dStress_dStrain1D_[0] = ContinuumMechanics::UniaxialStress::getUniaxialStressTangent( dStress_dStrain3D );
}
