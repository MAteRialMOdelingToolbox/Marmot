#include "Marmot/MarmotFiniteStrainPlasticity.h"
#include "Marmot/MarmotTensorExponential.h"

namespace Marmot {
  namespace ContinuumMechanics::FiniteStrain::Plasticity {
    namespace FlowIntegration {

      namespace FirstOrderDerived {

        using namespace Fastor;
        using namespace FastorStandardTensors;
        using namespace FastorIndices;

        std::pair< Tensor33d, Tensor3333d > explicitIntegration( const Tensor33d& dGp )
        {
          const auto&              I       = Spatial3D::I;
          const Tensor33d          dFp     = permute< Index< 1, 0 > >( I + dGp );
          const static Tensor3333d dFp_dGp = einsum< in, Im, to_iImn >( I, I ); // TODO: Simply Transposition operator..
          return { dFp, dFp_dGp };
        }

        std::pair< Tensor33d, Tensor3333d > exponentialMap( const Tensor33d& dGp )
        {
          const auto [dFpT, dFpT_dGp] = TensorUtility::TensorExponential::FirstOrderDerived::
            computeTensorExponential( dGp, 15, 1e-14 );
          // clang-format off
            return { permute< Index< 1, 0 > >      ( dFpT     ),
                     permute< Index< 1, 0, 2, 3 > >( dFpT_dGp ) };
          // clang-format on
        }

      } // namespace FirstOrderDerived

    }   // namespace FlowIntegration

  }     // namespace ContinuumMechanics::FiniteStrain::Plasticity

} // namespace Marmot
