#include "Marmot/Ogden.h"
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotEnergyDensityFunctions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotStressMeasures.h"
#include <Fastor/expressions/linalg_ops/unary_trans_op.h>
#include <Fastor/tensor/Tensor.h>
#include <Fastor/tensor_algebra/einsum_explicit.h>
#include <Fastor/tensor_algebra/indicial.h>
#include <autodiff/forward/dual/dual.hpp>
#include <stdexcept>

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Fastor;
  using namespace FastorIndices;
  using namespace FastorStandardTensors;

  namespace {
    int checkedNOgdenTerms( const double* materialProperties, int nMaterialProperties )
    {
      if ( nMaterialProperties < 1 )
        throw std::invalid_argument( "Ogden: expected the number of Ogden terms as the first material property" );
      const int n = static_cast< int >( materialProperties[0] );
      if ( n < 1 || nMaterialProperties < 1 + 2 * n )
        throw std::invalid_argument(
          "Ogden: not enough material properties provided for the requested number of Ogden terms" );
      return n;
    }
  } // namespace

  Ogden::Ogden( const double* materialProperties, int nMaterialProperties, int materialLabel )
    : MarmotMaterialFiniteStrainAD( materialProperties, nMaterialProperties, materialLabel ),
      nOgdenTerms( checkedNOgdenTerms( materialProperties, nMaterialProperties ) )
  {
    stateLayout.finalize();
  }

  double Ogden::getDensity( const double* stateVars ) const
  {
    const int densityIndex = 1 + 2 * nOgdenTerms;
    if ( nMaterialProperties < densityIndex + 1 )
      throw std::runtime_error( "Ogden: no density given!" );
    return materialProperties[densityIndex];
  }

  void Ogden::computeStressAD( ConstitutiveResponseAD< 3 >& response,
                               const DeformationAD< 3 >&    deformation,
                               const TimeIncrement&         timeIncrement ) const
  {
    using scalar  = autodiff::dual;
    const auto& F = deformation.F;

    using namespace ContinuumMechanics;

    // principal stretches lambda_i and eigenvectors Q of the right Cauchy-Green tensor C
    const auto [lambda, Q] = DeformationMeasures::principalStretches( F );
    const scalar J         = DeformationMeasures::volumeRatio( F );

    const double* mu    = &materialProperties[1];
    const double* alpha = &materialProperties[1 + nOgdenTerms];

    const auto [psi_,
                tauPrincipal] = EnergyDensityFunctions::Incompressible::OgdenPrincipalKirchhoffStress( lambda,
                                                                                                       J,
                                                                                                       mu,
                                                                                                       alpha,
                                                                                                       nOgdenTerms );

    const auto tau = StressMeasures::KirchhoffStressFromPrincipalKirchhoffStress( tauPrincipal, lambda, Q, F );

    response.tau                  = tau;
    response.elasticEnergyDensity = Math::makeReal( psi_ );
    response.dissipation          = 0.0;
  }
} // namespace Marmot::Materials
