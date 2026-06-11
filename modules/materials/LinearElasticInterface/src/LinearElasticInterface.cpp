#include "Marmot/LinearElasticInterface.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotInterfaceMaterialHelperFunctions.h"

#include "Fastor/Fastor.h"

namespace Marmot::Materials {

  using namespace Marmot::FastorIndices;
  using namespace Marmot::FastorStandardTensors;

  LinearElasticInterface::LinearElasticInterface( const double* materialProperties,
                                                  int           nMaterialProperties,
                                                  int           materialNumber )
    : MarmotInterfaceMaterialHypoElastic::MarmotInterfaceMaterialHypoElastic( materialProperties,
                                                                              nMaterialProperties,
                                                                              materialNumber ),
      E_0( materialProperties[0] ),
      nu_0( materialProperties[1] ),
      h( materialProperties[2] )
  {
    assert( nMaterialProperties == 3 );
  }

  void LinearElasticInterface::computeStress( State&               state,
                                              Tangents&            tangents,
                                              const Deformation&   deformation,
                                              const TimeIncrement& timeIncrement )
  {
    using namespace Marmot::Materials::InterfaceMaterialHelperFunctions;

    // use Fastor because we really need to use the einsum
    auto&       force_ftensor            = state.force;
    auto&       surface_stress_ftensor   = state.surfaceStress;
    auto&       H_inv_ij_ftensor         = tangents.Q_ij;
    auto&       Z_ijkl_ftensor           = tangents.Z_ijkl;
    auto&       H_inv_nF_ijk_ftensor     = tangents.H_ijk;
    auto&       Yn_H_inv_Fn_ijkl_ftensor = tangents.Y_ijkl;
    auto        dU_ftensor               = deformation.dU;
    auto        dSurface_strain_ftensor  = deformation.dSurfaceStrain;
    const auto& normal_ftensor           = deformation.normal;

    auto [unitZ_ijkl,
          unitH_inv_ij,
          unitH_inv_Fn_ijk,
          unitYn_H_inv_Fn_ijkl] = calculateInterfaceMaterialParameters( normal_ftensor, nu_0 );

    Z_ijkl_ftensor           = h * E_0 * unitZ_ijkl;
    Yn_H_inv_Fn_ijkl_ftensor = h * E_0 * unitYn_H_inv_Fn_ijkl;
    H_inv_ij_ftensor         = 1. / h * E_0 * unitH_inv_ij;
    H_inv_nF_ijk_ftensor     = E_0 * unitH_inv_Fn_ijk;

    // handle zero strain increment
    if ( Fastor::norm( dU_ftensor ) < 1e-14 && Fastor::norm( dSurface_strain_ftensor ) < 1e-14 ) {
      return;
    }
    // Compute stress increment
    Tensor3d jumpU_ftensor = dU_ftensor( Fastor::seq( 0, 3 ), 0 ) - dU_ftensor( Fastor::seq( 3, Fastor::last ), 0 );

    Tensor91d average_dSurface_strain_ftensor = 1. / 2. *
                                                ( dSurface_strain_ftensor( Fastor::seq( 0, 9 ), 0 ) +
                                                  dSurface_strain_ftensor( Fastor::seq( 9, Fastor::last ), 0 ) );
    auto average_dSurface_strain_ftensor_reshape = Fastor::reshape< 3, 3 >( average_dSurface_strain_ftensor );
    force_ftensor += Fastor::einsum< ij, j, to_i >( H_inv_ij_ftensor, jumpU_ftensor );
    force_ftensor += Fastor::einsum< ijk, jk, to_i >( H_inv_nF_ijk_ftensor, average_dSurface_strain_ftensor_reshape );
    surface_stress_ftensor += Fastor::einsum< ijkl, kl, to_ij >( Z_ijkl_ftensor,
                                                                 average_dSurface_strain_ftensor_reshape );
    surface_stress_ftensor += Fastor::einsum< ijkl, kl, to_ij >( Yn_H_inv_Fn_ijkl_ftensor,
                                                                 average_dSurface_strain_ftensor_reshape );
    surface_stress_ftensor += Fastor::einsum< i, ijk, to_jk >( jumpU_ftensor, H_inv_nF_ijk_ftensor );
  };

} // namespace Marmot::Materials
