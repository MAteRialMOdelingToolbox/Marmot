#include "Marmot/CompressibleNeoHooke.h"
#include "Marmot/MarmotAutomaticDifferentiationForFastor.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotMicromorphicTensorBasics.h"
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

  CompressibleNeoHooke::CompressibleNeoHooke( const double* materialProperties,
                                              int           nMaterialProperties,
                                              int           materialLabel )
    : MarmotMaterialFiniteStrain( materialProperties, nMaterialProperties, materialLabel )
  {
  }

  void CompressibleNeoHooke::computeStress( ConstitutiveResponse< 3 >&       response,
                                            AlgorithmicModuli< 3 >&          tangents,
                                            const DeformationIncrement< 3 >& deformationIncrement,
                                            const TimeIncrement&             timeIncrement,
                                            double&                          pNewDT )
  {
    const double& K = materialProperties[0];
    const double& G = materialProperties[1];

    double                               psi_;
    Fastor::Tensor< double, 3, 3 >       dPsi_dC, PK2;
    Fastor::Tensor< double, 3, 3, 3, 3 > dS_dPK2, d2Psi_dCdC, dC_dF;
    Fastor::Tensor< double, 3, 3 >       C = einsum< KI, KJ >( deformationIncrement.F_np, deformationIncrement.F_np );

    const auto& I = FastorStandardTensors::Spatial3D::I;

    dC_dF = einsum< LI, KJ, to_IJKL >( I, deformationIncrement.F_np ) +
            einsum< JL, KI, to_IJKL >( I, deformationIncrement.F_np );
    dS_dPK2 = einsum< iK, jL, to_ijKL >( deformationIncrement.F_np, deformationIncrement.F_np );

    std::tie( psi_, dPsi_dC, d2Psi_dCdC ) = AutomaticDifferentiation::SecondOrder::d2f_dTensor_dTensor<
      3 >( [&]( const Fastor::Tensor< autodiff::dual2nd, 3, 3 >& Ce_ ) { return psi( Ce_, K, G ); }, C );

    PK2 = 2.0 * dPsi_dC;
    /* std::cout << "PK2 = " << PK2 << std::endl; */
    response.S = einsum< iI, IJ, jJ, to_ij >( deformationIncrement.F_np, PK2, deformationIncrement.F_np );

    std::cout << dC_dF << std::endl;

    tangents.dS_dF = dS_dPK2 * 2. * d2Psi_dCdC * dC_dF;
    std::exit( 1 );
  }

  StateView CompressibleNeoHooke::getStateView( const std::string& stateName )
  {
    static std::map< std::string, std::tuple< int, int > > stateMapping = {};

    const auto result = stateMapping.at( stateName );

    return { &this->stateVars[std::get< 0 >( result )], std::get< 1 >( result ) };
  }
} // namespace Marmot::Materials
