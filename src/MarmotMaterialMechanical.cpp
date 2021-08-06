#include "Marmot/MarmotMaterialMechanical.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotLowerDimensionalStress.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotVoigt.h"
#include <iostream>

using namespace Eigen;

// void MarmotMaterialHypoElastic::computeStress( double*       stress_,
// double*       dStressDDStrain_,
// const double* FOld_,
// const double* FNew_,
// const double* timeOld,
// const double  dT,
// double&       pNewDT )
//{
// const Map<const Matrix3d> FNew( FNew_ );
// const Map<const Matrix3d> FOld( FOld_ );

////Marmot::Vector6d dEps = 1./2 * ( Marmot::ContinuumMechanics::VoigtNotation::StrainToVoigt( H + H.tranpose() ) );

////computeStress (stress_, dStressDDStrain_, dEps.data(), timeOld, dT, pNewDT);
//}

void MarmotMaterialMechanical::computePlaneStress( double*       stress2D_,
                                                   double*       dStress_dDeformationGradient2D_,
                                                   const double* FOld2D_,
                                                   const double* FNew2D_,
                                                   const double* timeOld,
                                                   const double  dT,
                                                   double&       pNewDT )
{

  using namespace Marmot;

  Map< Vector3d >               stress2D( stress2D_ );
  Eigen::TensorMap< Eigen::Tensor<double, 3 > > dStress_dDeformationGradient2D( dStress_dDeformationGradient2D_, 3, 2, 2);
  Map< const Matrix2d >         FNew2D( FNew2D_ );
  Map< const Matrix2d >         FOld2D( FOld2D_ );
  Map< VectorXd >               stateVars( this->stateVars, this->nStateVars );

  Vector6d stress3DTemp;
  VectorXd stateVarsOld = stateVars;

  Matrix3d FNew3D = Matrix3d::Identity();
  FNew3D.topLeftCorner(2,2) = FNew2D;

  Matrix3d FOld3D = Matrix3d::Identity();
  FOld3D .topLeftCorner(2,2) = FOld2D;

  // assumption of isochoric deformation for initial guess
  FNew3D( 2, 2 ) = 1. / ( FNew2D( 0, 0 ) * FNew2D( 1, 1 ) );

  EigenTensors::Tensor633d   dStress_dDeformationGradient3D;

  int planeStressCount = 1;
  while ( true ) {
    stress3DTemp = Marmot::ContinuumMechanics::VoigtNotation::make3DVoigt< Marmot::ContinuumMechanics::VoigtNotation::VoigtSize::TwoD >( stress2D );
    stateVars  = stateVarsOld;

    computeStress( stress3DTemp.data(),
                   dStress_dDeformationGradient3D.data(),
                   FOld3D.data(),
                   FNew3D.data(),
                   timeOld,
                   dT,
                   pNewDT );

    if ( pNewDT < 1.0 ) {
      return;
    }

    double residual = stress3DTemp.array().abs()[2];

    if ( residual < 1.e-10 || ( planeStressCount > 7 && residual < 1e-8 ) ) {
      break;
    }

    const double dS33_dF33 = dStress_dDeformationGradient2D( 2, 2, 2 );

    double tangentCompliance = 1. / dS33_dF33;
    if ( Math::isNaN( tangentCompliance ) || std::abs( tangentCompliance ) > 1e10 )
      tangentCompliance = 1e10;

    FNew3D( 2, 2 ) -= tangentCompliance * stress3DTemp( 2 );

    planeStressCount += 1;
    if ( planeStressCount > 13 ) {
      pNewDT = 0.25;
      MarmotJournal::warningToMSG( "PlaneStressWrapper requires cutback" );
      return;
    }
  }

  stress2D  = ContinuumMechanics::VoigtNotation::reduce3DVoigt< Marmot::ContinuumMechanics::VoigtNotation::VoigtSize::TwoD>( stress3DTemp );
  dStress_dDeformationGradient2D = ContinuumMechanics::PlaneStress::compute_dStress_dDeformationGradient( dStress_dDeformationGradient3D);
}

void MarmotMaterialMechanical::computeUniaxialStress( double*       stress_,
                                                      double*       dStressDDStrain_,
                                                      const double* FOld,
                                                      const double* FNew,
                                                      const double* timeOld,
                                                      const double  dT,
                                                      double&       pNewDT )
{

    throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << "not yet implemented");
}
