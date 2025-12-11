#include "Marmot/CompressibleNeoHooke.h"
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotEnergyDensityFunctions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotStressMeasures.h"
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
    initializeStateLayout();
  }

  void CompressibleNeoHooke::computeStress( ConstitutiveResponse< 3 >& response,
                                            AlgorithmicModuli< 3 >&    tangents,
                                            const Deformation< 3 >&    deformation,
                                            const TimeIncrement&       timeIncrement ) const
  {
    const double& K = materialProperties[0];
    const double& G = materialProperties[1];

    const auto& F_ = deformation.F;

    using namespace ContinuumMechanics;
    // compute Cauchy-Green deformation
    const auto [C, dC_dF] = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( F_ );

    // compute energy density, first and second partial derivatives wrt Cauchy Green deformation
    const auto [psi_, dPsi_dC, d2Psi_dCdC] = EnergyDensityFunctions::SecondOrderDerived::PenceGouPotentialB( C, K, G );

    // compute Kirchhoff stress
    Tensor33d PK2 = 2. * dPsi_dC;

    const auto [tau, dTau_dPK2, dTau_dF] = StressMeasures::FirstOrderDerived::KirchhoffStressFromPK2( PK2, F_ );
    response.tau                         = tau;
    response.rho                         = 1.0;
    response.elasticEnergyDensity        = psi_;

    // compute tangent operator
    tangents.dTau_dF = 2.0 * einsum< ijKL, KLMN >( einsum< ijKL, IJKL >( dTau_dPK2, d2Psi_dCdC ), dC_dF ) + dTau_dF;
  }
} // namespace Marmot::Materials
