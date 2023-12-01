#include "Marmot/CompressibleNeoHooke.h"
#include "Marmot/MarmotAutomaticDifferentiationForFastor.h"
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotEnergyDensityFunctions.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotMicromorphicTensorBasics.h"
#include "Marmot/MarmotNumericalDifferentiationForFastor.h"
#include "Marmot/MarmotStressMeasures.h"
#include "Marmot/MarmotVoigt.h"
#include <Fastor/expressions/linalg_ops/unary_trans_op.h>
#include <Fastor/tensor/Tensor.h>
#include <Fastor/tensor_algebra/einsum_explicit.h>
#include <Fastor/tensor_algebra/indicial.h>
#include <autodiff/forward/dual/dual.hpp>
#include <map>

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Fastor;
  using namespace FastorIndices;
  using namespace FastorStandardTensors;

  CompressibleNeoHooke::CompressibleNeoHooke( const double* materialProperties,
                                              int           nMaterialProperties,
                                              int           materialLabel )
    : MarmotMaterialFiniteStrain( materialProperties, nMaterialProperties, materialLabel )
  {
  }

  void CompressibleNeoHooke::computeStress( ConstitutiveResponse< 3 >& response,
                                            AlgorithmicModuli< 3 >&    tangents,
                                            const Deformation< 3 >&    deformation,
                                            const TimeIncrement&       timeIncrement )
  {
    const double& K = materialProperties[0];
    const double& G = materialProperties[1];

    const auto& F_ = deformation.F;
    double      psi_;
    Tensor33d   dPsi_dC;
    Tensor3333d d2Psi_dCdC;

    using namespace ContinuumMechanics;
    // compute Cauchy-Green deformation
    Tensor33d   C;
    Tensor3333d dC_dF;
    std::tie( C, dC_dF ) = DeformationMeasures::FirstOrderDerived::CauchyGreen( F_ );

    // compute energy density, first and second partial derivatives wrt Cauchy Green deformation
    std::tie( psi_, dPsi_dC, d2Psi_dCdC ) = EnergyDensityFunctions::SeconOrderDerived::PenceGouPotentialB( C, K, G );

    // compute Kirchhoff stress
    Tensor33d   PK2 = 2. * dPsi_dC;
    Tensor3333d dTau_dPK2, dTau_dF;
    std::tie( response.tau, dTau_dPK2, dTau_dF ) = StressMeasures::FirstOrderDerived::KirchhoffStressFromPK2( PK2, F_ );
    response.rho                                 = 1.0;
    response.elasticEnergyDensity                = psi_;

    // compute tangent operator
    tangents.dTau_dF = 2.0 * einsum< ijKL, KLMN >( einsum< ijKL, IJKL >( dTau_dPK2, d2Psi_dCdC ), dC_dF ) + dTau_dF;
  }

  StateView CompressibleNeoHooke::getStateView( const std::string& stateName )
  {
    static std::map< std::string, std::tuple< int, int > > stateMapping = {};

    const auto result = stateMapping.at( stateName );

    return { &this->stateVars[std::get< 0 >( result )], std::get< 1 >( result ) };
  }
} // namespace Marmot::Materials
