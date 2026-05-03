#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <vector>
#include "Marmot/LinearViscoElasticWiechert.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotUtility.h"
#include "Marmot/MarmotViscoelasticity.h"
#include "Marmot/MarmotVoigt.h"
#include "Marmot/MarmotWiechert.h"

#include "Fastor/Fastor.h"
#include "Marmot/MarmotInterfaceMaterialHelperFunctions.h"
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <Fastor/expressions/linalg_ops/unary_norm_op.h>
#include <Fastor/tensor/TensorMap.h>

#include "autodiff/forward/real.hpp"
#include <iostream>
#include <map>
#include <string>

using namespace Marmot;
using namespace Eigen;

using Tensor1D = Fastor::Tensor< double, 3 >;
using Tensor2D = Fastor::Tensor< double, 3, 3 >;
using Tensor3D = Fastor::Tensor< double, 3, 3, 3 >;
using Tensor4D = Fastor::Tensor< double, 3, 3, 3, 3 >;

namespace Marmot::Materials {

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
    Eigen::Map< const Eigen::Matrix< double, 6, 1 > > dE( dStrain );
    Eigen::Map< Eigen::Matrix< double, 6, 6, Eigen::RowMajor > > Cep( C );

    Cep.setZero();

    if ( timeInfo.dT < 0.0 ) {
      throw std::runtime_error( "LinearViscoElasticWiechert: negative time increment." );
    }

    const double dt = std::max( timeInfo.dT, 0.0 );

    /*
     * v26.05-native Wiechert update.
     *
     * State convention:
     *   state.stress     = total Cauchy stress at n
     *   state.stateVars  = stress-like Maxwell branch history variables q_i
     *
     * Each Maxwell branch stores six Voigt components:
     *   q_i = state.stateVars[6*i : 6*i+6]
     *
     * Incremental update:
     *   q_i^{n+1} = a_i q_i^n + g_i b_i C_el Δε
     *   σ^{n+1}   = σ^n + C_inf Δε + Σ_i (q_i^{n+1} - q_i^n)
     *
     * where
     *   a_i = exp(-Δt / τ_i)
     *   b_i = (1 - a_i) / (Δt / τ_i)
     *
     * This is native Marmot v26.05 style: no old factory, no old state assignment,
     * no old six-argument material API.
     */

    const int nStateVars = getNumberOfRequiredStateVars();
    const int nBranchesFromState = nStateVars / 6;

    const int nBranchesFromInput = static_cast< int >( std::round( nMaxwell ) );
    const int nBranches = std::max( 0, std::min( nBranchesFromState, nBranchesFromInput ) );

    const double lambda = E * nu / ( ( 1.0 + nu ) * ( 1.0 - 2.0 * nu ) );
    const double mu     = E / ( 2.0 * ( 1.0 + nu ) );

    Eigen::Matrix< double, 6, 6, Eigen::RowMajor > Cel;
    Cel.setZero();

    Cel( 0, 0 ) = lambda + 2.0 * mu;
    Cel( 1, 1 ) = lambda + 2.0 * mu;
    Cel( 2, 2 ) = lambda + 2.0 * mu;

    Cel( 0, 1 ) = lambda;
    Cel( 0, 2 ) = lambda;
    Cel( 1, 0 ) = lambda;
    Cel( 1, 2 ) = lambda;
    Cel( 2, 0 ) = lambda;
    Cel( 2, 1 ) = lambda;

    Cel( 3, 3 ) = mu;
    Cel( 4, 4 ) = mu;
    Cel( 5, 5 ) = mu;

    /*
     * Power-law inspired branch placement.
     *
     * minTau and timeToDays are kept from the old material parameter list.
     * The branch relaxation times are logarithmically spaced:
     *   τ_i = minTau * 10^i
     *
     * The branch weights follow the old power-law parameters m and n in a
     * normalized positive distribution. The equilibrium fraction is whatever
     * remains after all branch weights are assigned.
     */
    std::vector< double > tau( nBranches, 0.0 );
    std::vector< double > weight( nBranches, 0.0 );

    double weightSum = 0.0;

    for ( int i = 0; i < nBranches; ++i ) {
      tau[i] = minTau * std::pow( 10.0, static_cast< double >( i ) );

      const double tauDays = std::max( tau[i] * timeToDays, 1e-30 );
      weight[i] = std::max( 0.0, m * std::pow( tauDays, -n ) );

      weightSum += weight[i];
    }

    const double maxBranchFraction = 0.95;

    if ( weightSum > maxBranchFraction && weightSum > 0.0 ) {
      for ( double& w : weight ) {
        w *= maxBranchFraction / weightSum;
      }
      weightSum = maxBranchFraction;
    }

    const double equilibriumWeight = std::max( 0.0, 1.0 - weightSum );

    Eigen::Matrix< double, 6, 1 > dSigma = equilibriumWeight * Cel * dE;

    Cep += equilibriumWeight * Cel;

    for ( int branch = 0; branch < nBranches; ++branch ) {
      double* qRaw = state.stateVars + 6 * branch;

      Eigen::Map< Eigen::Matrix< double, 6, 1 > > qOld( qRaw );
      Eigen::Matrix< double, 6, 1 > qPrevious = qOld;

      const double relaxationTime = std::max( tau[branch] * timeToDays, 1e-30 );
      const double x = dt / relaxationTime;

      const double a = std::exp( -x );

      double b = 1.0;
      if ( x > 1e-14 ) {
        b = ( 1.0 - a ) / x;
      }

      const Eigen::Matrix< double, 6, 1 > dqElastic = weight[branch] * b * Cel * dE;

      qOld = a * qOld + dqElastic;

      dSigma += qOld - qPrevious;

      Cep += weight[branch] * b * Cel;
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
      timeToDays( materialProperties[6] )
  // clang-format on
  {
    relaxationTimes = Marmot::Materials::Wiechert::initializeRelaxationTimes( nMaxwell, m );
    elasticModuli   = Marmot::Materials::Wiechert::initializeElasticModuli( nMaxwell, n );

    using namespace Marmot::ContinuumMechanics::Viscoelasticity;
    // elasticModuli =
    // Marmot::Materials::Wiechert::computeElasticModuli<powerLawApproximationOrder>(phi,
    // relaxationTimes);

    zerothWiechertStiffness = 0.0; // m_Ru*(1. - n_Ru )*pow( 2., n_Ru )*pow(minTau_Ru/sqrt(10.), n_Ru);
  }

  


  


  

} // namespace Marmot::Materials
