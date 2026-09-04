#include "Marmot/MarmotFiniteStrainViscoelasticity.h"
#include "Marmot/MarmotMath.h"
#include <algorithm>
#include <cstring>
#include <vector>

namespace Marmot::ContinuumMechanics::Viscoelasticity::FiniteStrain {

  MaxwellProperties createMaxwellProperties( int nMaxwell, const double* gammaTauPairVector )
  {
    if ( nMaxwell < 0 )
      throw std::invalid_argument( "Number of Maxwell elements cannot be negative." );
    if ( nMaxwell == 0 )
      return MaxwellProperties( 0, {}, 0.0, {} );
    if ( gammaTauPairVector == nullptr )
      throw std::invalid_argument( "gammaTauPairVector cannot be null when number of Maxwell elements is positive." );

    std::vector< double > gamma( nMaxwell );
    std::vector< double > tau( nMaxwell );
    double                sumGamma = 0.0;
    for ( int i = 0; i < nMaxwell; ++i ) {
      const double gammaValue = gammaTauPairVector[i * 2];
      const double tauValue   = gammaTauPairVector[i * 2 + 1];
      // A vanishing relaxation time is admissible and handled in the evaluation routines,
      // where the corresponding Maxwell element simply does not contribute.
      if ( tauValue < 0.0 )
        throw std::invalid_argument( "Relaxation times of Maxwell elements must not be negative." );

      gamma[i] = gammaValue;
      tau[i]   = tauValue;
      sumGamma += gammaValue;
    }
    return MaxwellProperties( nMaxwell, gamma, sumGamma, tau );
  }

  void evaluateGeneralizedMaxwellModel(
    TensorUtility::FastorTensors::StandardTensors::Tensor33d&           stress,
    TensorUtility::FastorTensors::StandardTensors::Tensor3333d&         tangent,
    const TensorUtility::FastorTensors::StandardTensors::Tensor333333d& dTangent_dDeformation,
    const TensorUtility::FastorTensors::StandardTensors::Tensor3333d&   initialCompliance,
    const TensorUtility::FastorTensors::StandardTensors::Tensor33d&     dStress,
    const double                                                        dT,
    const MaxwellProperties&                                            maxwellProperties,
    double*                                                             stateVars )
  {

    if ( maxwellProperties.nMaxwell == 0 )
      return;

    using namespace Fastor;
    using namespace Marmot::TensorUtility::FastorTensors::StandardTensors;
    using namespace Marmot::TensorUtility::FastorTensors::Indices;

    // copy of tangent to be incremented
    const Tensor3333d initialTangent = tangent;

    // scale equilibrium stress contribution
    stress  = stress * ( 1.0 - maxwellProperties.sumGamma );
    tangent = tangent * ( 1.0 - maxwellProperties.sumGamma );

    for ( int i = 0; i < maxwellProperties.nMaxwell; ++i ) {
      // get old  maxewell element stress from state variables
      const Tensor33d& Q_n = Tensor33d( stateVars + i * 9 );

      // get parameters of maxwell element
      const double& tau   = maxwellProperties.tau[i];
      const double& gamma = maxwellProperties.gamma[i];

      if ( std::abs( tau ) < 1e-12 ) {
        // Very small relaxation time effectively zero; skip contribution
        continue;
      }
      const double dT_tau    = std::max( dT / tau, 1e-15 );
      const double expFactor = Math::exp( -dT_tau );

      double alpha = expFactor;
      double beta  = gamma / dT_tau * ( 1.0 - expFactor );

      if ( dT_tau < 1e-6 ) {
        // use taylor expansion for small dt/tau
        alpha = 1.0 - dT_tau + 0.5 * dT_tau * dT_tau;
        beta  = gamma * ( 1.0 - 0.5 * dT_tau + 1.0 / 6.0 * dT_tau * dT_tau );
      }

      // compute new stress in maxwell element
      const Tensor33d Q_np = alpha * Q_n + beta * dStress;

      // add contribution to stress
      const Tensor33d H_np = einsum< ijkl, kl >( initialCompliance, Q_np );
      stress += einsum< ij, ijkl >( H_np, initialTangent );

      const Tensor3333d dH_np_dDeformation = einsum< ijmn, mnKL >( initialCompliance,
                                                                   evaluate( beta * initialTangent ) );

      tangent += einsum< ijkl, ijmn >( initialTangent, dH_np_dDeformation );
      tangent += einsum< ij, ijklmn >( H_np, dTangent_dDeformation );

      // update state variables
      memcpy( stateVars + i * 9, Q_np.data(), 9 * sizeof( double ) );
    }
  }

  void evaluateGeneralizedMaxwellModel( TensorUtility::FastorTensors::StandardTensors::Tensor33d&       stress,
                                        TensorUtility::FastorTensors::StandardTensors::Tensor3333d&     tangent,
                                        const TensorUtility::FastorTensors::StandardTensors::Tensor33d& dStress,
                                        const double                                                    dT,
                                        const MaxwellProperties& maxwellProperties,
                                        double*                  stateVars )
  {

    if ( maxwellProperties.nMaxwell == 0 )
      return;

    using namespace Fastor;
    using namespace Marmot::TensorUtility::FastorTensors::StandardTensors;
    using namespace Marmot::TensorUtility::FastorTensors::Indices;

    // copy of tangent to be incremented
    const Tensor3333d initialTangent = tangent;

    // scale equilibrium stress contribution
    stress  = stress * ( 1.0 - maxwellProperties.sumGamma );
    tangent = tangent * ( 1.0 - maxwellProperties.sumGamma );

    for ( int i = 0; i < maxwellProperties.nMaxwell; ++i ) {
      // get old  maxewell element stress from state variables
      const Tensor33d& Q_n = Tensor33d( stateVars + i * 9 );

      // get parameters of maxwell element
      const double& tau   = maxwellProperties.tau[i];
      const double& gamma = maxwellProperties.gamma[i];
      if ( std::abs( tau ) < 1e-12 ) {
        // Very small relaxation time effectively zero; skip contribution
        continue;
      }

      const double dT_tau    = std::max( dT / tau, 1e-15 );
      const double expFactor = Math::exp( -dT_tau );

      double alpha = expFactor;
      double beta  = gamma / dT_tau * ( 1.0 - expFactor );

      if ( dT_tau < 1e-6 ) {
        // use taylor expansion for small dt/tau
        alpha = 1.0 - dT_tau + 0.5 * dT_tau * dT_tau;
        beta  = gamma * ( 1.0 - 0.5 * dT_tau + 1.0 / 6.0 * dT_tau * dT_tau );
      }

      // compute new stress in maxwell element
      const Tensor33d Q_np = alpha * Q_n + beta * dStress;

      // add contribution to stress
      stress += Q_np;

      tangent += beta * initialTangent;

      // update state variables
      memcpy( stateVars + i * 9, Q_np.data(), 9 * sizeof( double ) );
    }
  }

} // namespace Marmot::ContinuumMechanics::Viscoelasticity::FiniteStrain
