#include "Marmot/LinearViscoElasticWiechert.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotWiechert.h"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <stdexcept>

using namespace Marmot;

namespace Marmot::Materials {

  namespace {
    using VoigtVector6d       = Marmot::Vector6d;
    using VoigtMatrix6d       = Marmot::Matrix6d;
    using ConstVoigtVectorMap = Eigen::Map< const VoigtVector6d >;
    using TangentMatrixMap    = Eigen::Map< Eigen::Matrix< double, 6, 6, Eigen::RowMajor > >;
  } // namespace

  void LinearViscoElasticWiechert::initializeStateLayout()
  {
    const int nBranches = static_cast< int >( nMaxwell );

    stateLayout.add( "MaxwellStress", 6 * nBranches );
    stateLayout.finalize();
  }

  void LinearViscoElasticWiechert::computeStress( state3D&        state,
                                                  double*         C,
                                                  const double*   dStrain,
                                                  const timeInfo& timeInfo ) const
  {
    ConstVoigtVectorMap dE( dStrain );
    TangentMatrixMap    Cep( C );
    const VoigtVector6d dEVec = dE;

    Cep.setZero();

    if ( timeInfo.dT < 0.0 ) {
      throw std::runtime_error( "LinearViscoElasticWiechert: negative time increment." );
    }

    const double dt = std::max( timeInfo.dT, 0.0 );

    auto maxwellStressStateVars = stateLayout.getAs< Eigen::Map< Eigen::MatrixXd > >( state.stateVars,
                                                                                      "MaxwellStress",
                                                                                      6,
                                                                                      static_cast< int >( nMaxwell ) );

    const int nBranchesFromState = static_cast< int >( maxwellStressStateVars.cols() );
    const int nBranchesFromInput = static_cast< int >( nMaxwell );
    const int nBranches          = std::max( 0, std::min( nBranchesFromState, nBranchesFromInput ) );

    const VoigtMatrix6d Cel = Marmot::ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );

    const double equilibriumWeight = std::max( 0.0, 1.0 - branchElasticModuli.head( nBranches ).sum() );

    VoigtVector6d dSigma = equilibriumWeight * Cel * dEVec;

    Cep += equilibriumWeight * Cel;

    if ( nBranches > 0 ) {
      Eigen::Map< Wiechert::StateVarMatrix > activeMaxwellState( maxwellStressStateVars.data(), 6, nBranches );

      const Wiechert::StateVarMatrix qOld = activeMaxwellState;

      const Wiechert::Properties activeBranchElasticModuli   = branchElasticModuli.head( nBranches );
      const Wiechert::Properties activeBranchRelaxationTimes = branchRelaxationTimes.head( nBranches );

      Wiechert::updateStateVarMatrix( dt,
                                      activeBranchElasticModuli,
                                      activeBranchRelaxationTimes,
                                      activeMaxwellState,
                                      dEVec,
                                      Cel );

      dSigma += ( activeMaxwellState - qOld ).rowwise().sum();

      double        creepStiffness = 0.0;
      VoigtVector6d dStressDummy   = VoigtVector6d::Zero();

      Wiechert::evaluateWiechert( dt,
                                  activeBranchElasticModuli,
                                  activeBranchRelaxationTimes,
                                  Wiechert::StateVarMatrix::Zero( 6, nBranches ),
                                  creepStiffness,
                                  dStressDummy,
                                  1.0 );

      Cep += creepStiffness * Cel;
    }

    state.stress += dSigma;
    state.strainEnergyDensity += 0.5 * dE.dot( state.stress );
  }

  LinearViscoElasticWiechert::LinearViscoElasticWiechert( const double* materialProperties,
                                                          int           nMaterialProperties,
                                                          int           materialNumber )
    : MarmotMaterialHypoElastic( materialProperties, nMaterialProperties, materialNumber ),
      // clang-format off
      // elasticity parameters
      E( materialProperties[0] ),
      nu( materialProperties[1] ),
      m( materialProperties[2] ),
      n( materialProperties[3] ),
      nMaxwell( static_cast< size_t >( materialProperties[4] ) ),
      minTau( materialProperties[5] ),
      timeToDays( materialProperties[6] ),
      branchRelaxationTimes( static_cast< int >( nMaxwell ) ),
      branchElasticModuli( static_cast< int >( nMaxwell ) )
  // clang-format on
  {
    initializeStateLayout();

    constexpr double maxBranchFraction = 0.95;

    for ( int i = 0; i < static_cast< int >( nMaxwell ); ++i ) {
      const double tau     = minTau * std::pow( 10.0, static_cast< double >( i ) );
      const double tauDays = std::max( tau * timeToDays, 1e-30 );

      branchRelaxationTimes[i] = tauDays;
      branchElasticModuli[i]   = std::max( 0.0, m * std::pow( tauDays, -n ) );
    }

    const double weightSum = branchElasticModuli.sum();

    if ( weightSum > maxBranchFraction && weightSum > 0.0 ) {
      branchElasticModuli *= maxBranchFraction / weightSum;
    }
  }

} // namespace Marmot::Materials