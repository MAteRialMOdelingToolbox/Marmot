#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotLowerDimensionalStress.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotVoigt.h"

namespace Marmot {

  using namespace Eigen;

  void MarmotMaterialHypoElastic::setCharacteristicElementLength( double length )
  {
    characteristicElementLength = length;
  }

  void MarmotMaterialHypoElastic::computePlaneStress( state2D&                state2D_,
                                                      Marmot::Matrix3d&       dStress_dStrain2D_,
                                                      const Marmot::Vector3d& dStrain2D_,
                                                      const timeInfo&         timeInfo ) const
  {
    using namespace Marmot;
    using namespace ContinuumMechanics::Voigt;

    Map< const Matrix< double, 3, 1 > > dStrain2D( dStrain2D_.data() );
    Map< Matrix< double, 3, 1 > >       stress2D( state2D_.stress.data() );
    Map< Matrix< double, 3, 3 > >       dStress_dStrain2D( dStress_dStrain2D_.data() );
    Map< VectorXd >                     stateVars( state2D_.stateVars, this->getNumberOfRequiredStateVars() );

    Matrix6d dStress_dStrain3D;

    VectorXd stateVarsOld  = stateVars;
    Vector6d dStrain3DTemp = Marmot::ContinuumMechanics::Voigt::make3DVoigt< VoigtSize::TwoD >( dStrain2D );

    // assumption of isochoric deformation for initial guess
    dStrain3DTemp( 2 ) = ( -dStrain2D( 0 ) - dStrain2D( 1 ) );

    state3D state( Marmot::ContinuumMechanics::Voigt::make3DVoigt< VoigtSize::TwoD >( stress2D ),
                   state2D_.elasticEnergyDensity,
                   state2D_.dissipation,
                   stateVars.data() );

    int planeStressCount = 1;
    while ( true ) {

      stateVars = stateVarsOld;

      state.stress               = Marmot::ContinuumMechanics::Voigt::make3DVoigt< VoigtSize::TwoD >( stress2D );
      state.elasticEnergyDensity = state2D_.elasticEnergyDensity;
      state.dissipation          = state2D_.dissipation;

      computeStress( state, dStress_dStrain3D, dStrain3DTemp, timeInfo );

      double residual = state.stress.array().abs()[2];

      if ( residual < 1.e-10 || ( planeStressCount > 7 && residual < 1e-8 ) ) {
        break;
      }

      double tangentCompliance = 1. / dStress_dStrain3D( 2, 2 );
      if ( Math::isNaN( tangentCompliance ) || std::abs( tangentCompliance ) > 1e10 )
        tangentCompliance = 1e10;

      dStrain3DTemp( 2 ) -= tangentCompliance * state.stress( 2 );

      planeStressCount += 1;
      if ( planeStressCount > 13 ) {
        MarmotJournal::warningToMSG( "PlaneStressWrapper requires cutback" );
        throw Marmot::StressUpdateFailed( "Plane stress iteration did not converge" );
      }
    }

    state2D_.stress               = ContinuumMechanics::Voigt::reduce3DVoigt< VoigtSize::TwoD >( state.stress );
    state2D_.elasticEnergyDensity = state.elasticEnergyDensity;
    state2D_.dissipation          = state.dissipation;
    state2D_.stateVars            = stateVars.data();

    dStress_dStrain2D_ = ContinuumMechanics::LowerOrder::PlaneStress::getPlaneStressTangent( dStress_dStrain3D );
  }

  void MarmotMaterialHypoElastic::computeUniaxialStress( state1D&        state1D_,
                                                         double&         dStress_dStrain1D_,
                                                         const double    dStrain1D_,
                                                         const timeInfo& timeInfo ) const
  {
    using namespace Marmot;
    using namespace ContinuumMechanics::Voigt;

    Map< const Matrix< double, 1, 1 > > dStrain1D( &dStrain1D_ );
    Map< Matrix< double, 1, 1 > >       stress1D( &state1D_.stress );
    Map< VectorXd >                     stateVars( state1D_.stateVars, stateLayout.totalSize() );

    Matrix6d dStress_dStrain3D;

    VectorXd stateVarsOld  = stateVars;
    Vector6d dStrain3DTemp = Marmot::ContinuumMechanics::Voigt::make3DVoigt< VoigtSize::OneD >( dStrain1D );

    state3D state( Marmot::ContinuumMechanics::Voigt::make3DVoigt< VoigtSize::OneD >( stress1D ),
                   state1D_.elasticEnergyDensity,
                   state1D_.dissipation,
                   stateVars.data() );

    int count = 1;
    while ( true ) {
      stateVars = stateVarsOld;

      state.stress               = Marmot::ContinuumMechanics::Voigt::make3DVoigt< VoigtSize::OneD >( stress1D );
      state.elasticEnergyDensity = state1D_.elasticEnergyDensity;
      state.dissipation          = state1D_.dissipation;

      computeStress( state, dStress_dStrain3D, dStrain3DTemp, timeInfo );

      const double residual = state.stress.array().abs().segment( 1, 2 ).sum();

      if ( residual < 1.e-13 || ( count > 7 && residual < 1e-10 ) ) {
        break;
      }

      dStrain3DTemp.segment< 2 >( 1 ) -= dStress_dStrain3D.block< 2, 2 >( 1, 1 ).colPivHouseholderQr().solve(
        state.stress.segment< 2 >( 1 ) );

      count += 1;
      if ( count > 13 ) {
        MarmotJournal::warningToMSG( "UniaxialStressWrapper requires cutback" );
        throw Marmot::StressUpdateFailed( "uniaxial stress iteration did not converge" );
      }
    }

    state1D_.stress               = ContinuumMechanics::Voigt::reduce3DVoigt< VoigtSize::OneD >( state.stress )[0];
    state1D_.elasticEnergyDensity = state.elasticEnergyDensity;
    state1D_.dissipation          = state.dissipation;
    state1D_.stateVars            = stateVars.data();

    dStress_dStrain1D_ = ContinuumMechanics::LowerOrder::UniaxialStress::getUniaxialStressTangent( dStress_dStrain3D );
  }
  double MarmotMaterialHypoElastic::getMaximumWaveSpeed( const state3D& state ) const
  {
    const int nStateVars = getNumberOfRequiredStateVars();

    std::vector< double > stateVarsCopy( nStateVars, 0.0 );
    if ( state.stateVars != nullptr && nStateVars > 0 ) {
      std::copy_n( state.stateVars, nStateVars, stateVarsCopy.begin() );
    }

    state3D stateCopy( state.stress, state.elasticEnergyDensity, state.dissipation, stateVarsCopy.data() );

    Marmot::Matrix6d dStress_dStrain = Marmot::Matrix6d::Zero();
    const auto       dStrain         = Marmot::Vector6d::Zero();
    const timeInfo   timeInfo{ 0.0, 1.0 };

    computeStress( stateCopy, dStress_dStrain, dStrain, timeInfo );

    const double maxStiffnessDiagonal = std::max( { dStress_dStrain( 0, 0 ),
                                                    dStress_dStrain( 1, 1 ),
                                                    dStress_dStrain( 2, 2 ),
                                                    dStress_dStrain( 3, 3 ),
                                                    dStress_dStrain( 4, 4 ),
                                                    dStress_dStrain( 5, 5 ) } );
    const double density              = getDensity( stateCopy.stateVars );

    return density > 0.0 ? std::sqrt( std::max( 0.0, maxStiffnessDiagonal ) / density ) : 0.0;
  }

} // namespace Marmot
