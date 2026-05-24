#include "Marmot/ADCompressibleNeoHooke.h"
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

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Fastor;
  using namespace FastorIndices;
  using namespace FastorStandardTensors;

  ADCompressibleNeoHooke::ADCompressibleNeoHooke( const double* materialProperties,
                                                  int           nMaterialProperties,
                                                  int           materialLabel )
    : MarmotMaterialFiniteStrainAD( materialProperties, nMaterialProperties, materialLabel )
  {
    stateLayout.finalize();
  }

  void ADCompressibleNeoHooke::computeStressAD( ConstitutiveResponseAD< 3 >& response,
                                                const DeformationAD< 3 >&    deformation,
                                                const TimeIncrement&         timeIncrement ) const
  {
    const double& K = materialProperties[0];
    const double& G = materialProperties[1];

    const auto& F_ = deformation.F;

    using namespace ContinuumMechanics;
    // compute Cauchy-Green deformation
    const auto C = DeformationMeasures::rightCauchyGreen( F_ );

    // compute energy density and first derivative w.r.t Cauchy Green deformation
    const auto [psi_, dPsi_dC] = EnergyDensityFunctions::FirstOrderDerived::PenceGouPotentialB( C, K, G );

    // compute Kirchhoff stress
    Tensor33t< autodiff::dual > PK2 = multiplyFastorTensorWithScalar( dPsi_dC, autodiff::dual( 2.0 ) );

    const auto tau                = StressMeasures::KirchhoffStressFromPK2( PK2, F_ );
    response.tau                  = tau;
    response.elasticEnergyDensity = Math::makeReal( psi_ );
  }
} // namespace Marmot::Materials
