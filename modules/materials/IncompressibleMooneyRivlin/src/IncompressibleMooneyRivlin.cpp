#include "Marmot/IncompressibleMooneyRivlin.h"
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
  using namespace FastorIndices;
  using namespace FastorStandardTensors;

  IncompressibleMooneyRivlin::IncompressibleMooneyRivlin( const double* materialProperties,
                                                          int           nMaterialProperties,
                                                          int           materialLabel )
    : MarmotMaterialFiniteStrain( materialProperties, nMaterialProperties, materialLabel )
  {
    stateLayout.finalize();
  }

  void IncompressibleMooneyRivlin::computeStress( ConstitutiveResponse< 3 >& response,
                                                  AlgorithmicModuli< 3 >&    tangents,
                                                  const Deformation< 3 >&    deformation,
                                                  const TimeIncrement&       timeIncrement ) const
  {
    const double& C1 = materialProperties[0];
    const double& C2 = materialProperties[1];
    const auto&   F_ = deformation.F;

    using namespace ContinuumMechanics;

    // Unlike Ogden's general (non-integer) exponents, the classical Mooney-Rivlin potential is
    // quadratic in the isochoric stretches and can therefore be expressed directly in terms of the
    // invariants of C - no spectral decomposition of C is required.
    const auto [C, dC_dF] = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( F_ );

    const auto [psi_,
                dPsi_dC,
                d2Psi_dCdC] = EnergyDensityFunctions::Incompressible::SecondOrderDerived::MooneyRivlinPotential( C,
                                                                                                                 C1,
                                                                                                                 C2 );

    Tensor33d PK2 = 2. * dPsi_dC;

    const auto [tau, dTau_dPK2, dTau_dF] = StressMeasures::FirstOrderDerived::KirchhoffStressFromPK2( PK2, F_ );
    response.tau                         = tau;
    response.elasticEnergyDensity        = psi_;
    response.dissipation                 = 0.0;

    tangents.dTau_dF = 2.0 * einsum< ijKL, KLMN >( einsum< ijKL, IJKL >( dTau_dPK2, d2Psi_dCdC ), dC_dF ) + dTau_dF;
  }
} // namespace Marmot::Materials
