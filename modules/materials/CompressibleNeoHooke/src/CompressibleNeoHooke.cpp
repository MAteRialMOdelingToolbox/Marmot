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

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Fastor;
  using namespace Marmot::Math::Indices;
  using namespace Marmot::Math::FastorTensors;

  CompressibleNeoHooke::CompressibleNeoHooke( const double* materialProperties,
                                              int           nMaterialProperties,
                                              int           materialLabel )
    : Marmot::Materials::MaterialFiniteStrain( materialProperties, nMaterialProperties, materialLabel )
  {
    stateLayout.finalize();
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
    const auto [C, dC_dF] = Marmot::ContinuumMechanics::FiniteStrain::DeformationMeasures::FirstOrder::rightCauchyGreen(
      F_ );

    // compute energy density, first and second partial derivatives wrt Cauchy Green deformation
    const auto [psi_, dPsi_dC, d2Psi_dCdC] = Marmot::ContinuumMechanics::FiniteStrain::EnergyDensityFunctions::
      SecondOrder::PenceGouPotentialB( C, K, G );

    // compute Kirchhoff stress
    Tensor33d PK2 = 2. * dPsi_dC;

    const auto [tau, dTau_dPK2, dTau_dF] = Marmot::ContinuumMechanics::FiniteStrain::StressMeasures::FirstOrder::
      KirchhoffStressFromPK2( PK2, F_ );
    response.tau                  = tau;
    response.elasticEnergyDensity = psi_;
    response.dissipation          = 0.0;

    // compute tangent operator
    tangents.dTau_dF = 2.0 * einsum< ijKL, KLMN >( einsum< ijKL, IJKL >( dTau_dPK2, d2Psi_dCdC ), dC_dF ) + dTau_dF;
  }
} // namespace Marmot::Materials
