#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/HughesWinget.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotKinematics.h"
#include "Marmot/MarmotLowerDimensionalStress.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotVoigt.h"
#include <iostream>

using namespace Eigen;

void MarmotMaterialHypoElastic::setCharacteristicElementLength( double length )
{
  characteristicElementLength = length;
}

void MarmotMaterialHypoElastic::computePlaneStress( state2D&        state2D_,
                                                    double*         dStress_dStrain2D_,
                                                    const double*   dStrain2D_,
                                                    const timeInfo& timeInfo )
{
  using namespace Marmot;
  using namespace ContinuumMechanics::VoigtNotation;

  Map< const Matrix< double, 3, 1 > > dStrain2D( dStrain2D_ );
  Map< Matrix< double, 3, 1 > >       stress2D( state2D_.stress.data() );
  Map< Matrix< double, 3, 3 > >       dStress_dStrain2D( dStress_dStrain2D_ );
  Map< VectorXd >                     stateVars( state2D_.stateVars, this->getNumberOfRequiredStateVars() );

  Matrix6d dStress_dStrain3D;

  VectorXd stateVarsOld  = stateVars;
  Vector6d dStrain3DTemp = Marmot::ContinuumMechanics::VoigtNotation::make3DVoigt< VoigtSize::TwoD >( dStrain2D );

  // assumption of isochoric deformation for initial guess
  dStrain3DTemp( 2 ) = ( -dStrain2D( 0 ) - dStrain2D( 1 ) );

  state3D state;
  state.stress       = Marmot::ContinuumMechanics::VoigtNotation::make3DVoigt< VoigtSize::TwoD >( stress2D );
  state.stateVars    = stateVars.data();
  state.strainEnergy = state2D_.strainEnergy;

  int planeStressCount = 1;
  while ( true ) {

    stateVars = stateVarsOld;

    state.stress       = Marmot::ContinuumMechanics::VoigtNotation::make3DVoigt< VoigtSize::TwoD >( stress2D );
    state.strainEnergy = state2D_.strainEnergy;

    computeStress( state, dStress_dStrain3D.data(), dStrain3DTemp.data(), timeInfo );

    double residual = state.stress.array().abs()[2];

    if ( residual < 1.e-10 || ( planeStressCount > 7 && residual < 1e-8 ) ) {
      break;
    }

    double tangentCompliance = 1. / dStress_dStrain3D( 2, 2 );
    if ( Math::isNaN( tangentCompliance ) || std::abs( tangentCompliance ) > 1e10 )
      tangentCompliance = 1e10;

    dStrain3DTemp[2] -= tangentCompliance * state.stress[2];

    planeStressCount += 1;
    if ( planeStressCount > 13 ) {
      MarmotJournal::warningToMSG( "PlaneStressWrapper requires cutback" );
      throw std::runtime_error( "Plane stress iteration did not converge" );
      return;
    }
  }

  stress2D              = ContinuumMechanics::VoigtNotation::reduce3DVoigt< VoigtSize::TwoD >( state.stress );
  dStress_dStrain2D     = ContinuumMechanics::PlaneStress::getPlaneStressTangent( dStress_dStrain3D );
  state2D_.strainEnergy = state.strainEnergy;
}

void MarmotMaterialHypoElastic::computeUniaxialStress( state1D& state1D_,
                                                       double*  dStress_dStrain1D_,

                                                       const double*   dStrain1D_,
                                                       const timeInfo& timeInfo )
{
  using namespace Marmot;
  using namespace ContinuumMechanics::VoigtNotation;

  Map< const Matrix< double, 1, 1 > > dStrain1D( dStrain1D_ );
  Map< Matrix< double, 1, 1 > >       stress1D( &state1D_.stress );
  Map< VectorXd >                     stateVars( this->stateVars, this->nStateVars );

  Matrix6d dStress_dStrain3D;

  VectorXd stateVarsOld  = stateVars;
  Vector6d dStrain3DTemp = Marmot::ContinuumMechanics::VoigtNotation::make3DVoigt< VoigtSize::OneD >( dStrain1D );

  state3D state;
  state.stress       = Marmot::ContinuumMechanics::VoigtNotation::make3DVoigt< VoigtSize::OneD >( stress1D );
  state.stateVars    = stateVars.data();
  state.strainEnergy = state1D_.strainEnergy;
  int count          = 1;
  while ( true ) {
    stateVars          = stateVarsOld;
    state.stress       = Marmot::ContinuumMechanics::VoigtNotation::make3DVoigt< VoigtSize::OneD >( stress1D );
    state.strainEnergy = state1D_.strainEnergy;

    computeStress( state, dStress_dStrain3D.data(), dStrain3DTemp.data(), timeInfo );

    const double residual = state.stress.array().abs().segment( 1, 2 ).sum();

    if ( residual < 1.e-13 || ( count > 7 && residual < 1e-10 ) ) {
      break;
    }

    dStrain3DTemp.segment< 2 >( 1 ) -= dStress_dStrain3D.block< 2, 2 >( 1, 1 ).colPivHouseholderQr().solve(
      state.stress.segment< 2 >( 1 ) );

    count += 1;
    if ( count > 13 ) {
      MarmotJournal::warningToMSG( "UniaxialStressWrapper requires cutback" );
      throw std::runtime_error( "uniaxial stress iteration did not converge" );
      return;
    }
  }

  stress1D              = ContinuumMechanics::VoigtNotation::reduce3DVoigt< VoigtSize::OneD >( state.stress );
  dStress_dStrain1D_[0] = ContinuumMechanics::UniaxialStress::getUniaxialStressTangent( dStress_dStrain3D );
  state1D_.strainEnergy = state.strainEnergy;
}
