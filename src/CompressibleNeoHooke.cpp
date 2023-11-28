#include "Marmot/CompressibleNeoHooke.h"
#include "Marmot/MarmotAutomaticDifferentiationForFastor.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotMicromorphicTensorBasics.h"
#include "Marmot/MarmotNumericalDifferentiationForFastor.h"
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

    // compute Cauchy-Green deformation
    Tensor33d C = einsum< KI, KJ >( F_, F_ );

    // compute energy density, first and second partial derivatives wrt Cauchy Green deformation
    std::tie( psi_, dPsi_dC, d2Psi_dCdC ) = AutomaticDifferentiation::SecondOrder::d2f_dTensor_dTensor<
      3 >( [&]( const Fastor::Tensor< autodiff::dual2nd, 3, 3 >& Ce_ ) { return psi( Ce_, K, G ); }, C );

    // compute Kirchhoff stress
    response.S                    = einsum< iI, IJ, jJ, to_ij >( F_, 2. * dPsi_dC, F_ );
    response.rho                  = 1.0;
    response.elasticEnergyDensity = psi_;
    // compute tangent operator
    const auto&       I     = FastorStandardTensors::Spatial3D::I;
    const Tensor3333d dC_dF = einsum< LI, KJ, to_IJKL >( I, F_ ) + einsum< JL, KI, to_IJKL >( I, F_ );

    const Tensor3333d dS_dPK2 = einsum< iK, jL, to_ijKL >( F_, F_ );

    tangents.dS_dF = einsum< ijKL, KLMN >( einsum< ijKL, KLMN >( dS_dPK2, 2. * d2Psi_dCdC ), dC_dF );
  }

  StateView CompressibleNeoHooke::getStateView( const std::string& stateName )
  {
    static std::map< std::string, std::tuple< int, int > > stateMapping = {};

    const auto result = stateMapping.at( stateName );

    return { &this->stateVars[std::get< 0 >( result )], std::get< 1 >( result ) };
  }
} // namespace Marmot::Materials
