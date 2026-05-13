#include "Marmot/LinearElasticInterfaceBMGu.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotUtility.h"
#include "Marmot/MarmotVoigt.h"

#include "Fastor/Fastor.h"
#include "Marmot/MarmotInterfaceMaterialHelperFunctions.h"
#include <Fastor/tensor/TensorMap.h>

namespace Marmot::Materials {

  void LinearElasticInterfaceBMGu::initializeStateLayout() {}

  using namespace Marmot;
  using namespace Eigen;
  using namespace Fastor;

  using Tensor1D = Fastor::Tensor< double, 3 >;
  using Tensor2D = Fastor::Tensor< double, 3, 3 >;
  using Tensor3D = Fastor::Tensor< double, 3, 3, 3 >;
  using Tensor4D = Fastor::Tensor< double, 3, 3, 3, 3 >;

  LinearElasticInterfaceBMGu::LinearElasticInterfaceBMGu( const double* materialProperties,
                                                          int           nMaterialProperties,
                                                          int           materialNumber )
    : MarmotInterfaceMaterialHypoElastic::MarmotInterfaceMaterialHypoElastic( materialProperties,
                                                                              nMaterialProperties,
                                                                              materialNumber )
  {
    assert( nMaterialProperties == 7 );
  }

  void LinearElasticInterfaceBMGu::computeStress( double*       force,
                                                  double*       surface_stress,
                                                  double*       H_inv_ij,
                                                  double*       Z_ijkl,
                                                  double*       H_inv_nF_ijk,
                                                  double*       Yn_H_inv_Fn_ijkl,
                                                  const double* dU,
                                                  const double* dSurface_strain,
                                                  const double* normal,
                                                  const double* timeOld,
                                                  const double  dT,
                                                  double&       pNewDT )
  {
    using namespace Marmot::Materials::InterfaceMaterialHelperFunctions;

    // elasticity parameters
    const double& E_M  = this->materialProperties[0];
    const double& E_I  = this->materialProperties[2];
    const double& E_0  = this->materialProperties[4];
    const double& nu_0 = this->materialProperties[5];
    const double& h    = this->materialProperties[6];
    const double  Hbar = ( 2. / E_0 ) - ( 1. / E_M ) - ( 1. / E_I );
    const double  Zbar = ( E_M ) + ( E_I ) - ( 2. * E_0 );

    // map directly to force, surface stress, displacement, surface strain, normal and tangent stiffness
    // use Fastor because we really need to use the einsum
    auto force_ftensor            = Fastor::TensorMap< double, 3 >( force );
    auto surface_stress_ftensor   = Fastor::TensorMap< double, 3, 3 >( surface_stress );
    auto H_inv_ij_ftensor         = Fastor::TensorMap< double, 3, 3 >( H_inv_ij );
    auto Z_ijkl_ftensor           = Fastor::TensorMap< double, 3, 3, 3, 3 >( Z_ijkl );
    auto H_inv_nF_ijk_ftensor     = Fastor::TensorMap< double, 3, 3, 3 >( H_inv_nF_ijk );
    auto Yn_H_inv_Fn_ijkl_ftensor = Fastor::TensorMap< double, 3, 3, 3, 3 >( Yn_H_inv_Fn_ijkl );

    auto dU_ftensor_const              = Fastor::TensorMap< const double, 6, 1 >( dU );
    auto dSurface_strain_ftensor_const = Fastor::TensorMap< const double, 18, 1 >( dSurface_strain );
    auto normal_ftensor_const          = Fastor::TensorMap< const double, 3 >( normal );

    // Keep input temporaries for Fastor slicing/norm compatibility.
    Fastor::Tensor< double, 6, 1 >  dU_ftensor( dU_ftensor_const.data() );
    Fastor::Tensor< double, 18, 1 > dSurface_strain_ftensor( dSurface_strain_ftensor_const.data() );
    Fastor::Tensor< double, 3 >     normal_ftensor( normal_ftensor_const.data() );

    auto [unitZ_ijkl,
          unitH_inv_ij,
          unitH_inv_Fn_ijk,
          unitYn_H_inv_Fn_ijkl] = calculateInterfaceMaterialParameters( normal_ftensor, nu_0 );

    // Assign the material matrices to a larger structure. (Not necessary ...)
    Z_ijkl_ftensor           = -h / 2. * Zbar * unitZ_ijkl;
    Yn_H_inv_Fn_ijkl_ftensor = h * 0. * unitYn_H_inv_Fn_ijkl;
    H_inv_ij_ftensor         = 2. / h * 1. / Hbar * unitH_inv_ij;
    H_inv_nF_ijk_ftensor     = 0. * unitH_inv_Fn_ijk;

    // handle zero strain increment
    if ( Fastor::norm( dU_ftensor ) < 1e-14 && Fastor::norm( dSurface_strain_ftensor ) < 1e-14 ) {
      return;
    }
    // elastic step
    enum { i, j, k, l };

    Tensor1D jumpU_ftensor = dU_ftensor( Fastor::seq( 0, 3 ), 0 ) - dU_ftensor( Fastor::seq( 3, Fastor::last ), 0 );

    Fastor::Tensor< double, 9, 1 >
      average_dSurface_strain_ftensor = 1. / 2. *
                                        ( dSurface_strain_ftensor( Fastor::seq( 0, 9 ), 0 ) +
                                          dSurface_strain_ftensor( Fastor::seq( 9, Fastor::last ), 0 ) );
    auto average_dSurface_strain_ftensor_reshape = Fastor::reshape< 3, 3 >( average_dSurface_strain_ftensor );
    force_ftensor += Fastor::einsum< Fastor::Index< i, j >, Fastor::Index< j >, Fastor::OIndex< i > >( H_inv_ij_ftensor,
                                                                                                       jumpU_ftensor );
    force_ftensor -= Fastor::einsum< Fastor::Index< i, j, k >,
                                     Fastor::Index< j, k >,
                                     Fastor::OIndex< i > >( H_inv_nF_ijk_ftensor,
                                                            average_dSurface_strain_ftensor_reshape );
    surface_stress_ftensor -= Fastor::einsum< Fastor::Index< i, j, k, l >,
                                              Fastor::Index< k, l >,
                                              Fastor::OIndex< i, j > >( Z_ijkl_ftensor,
                                                                        average_dSurface_strain_ftensor_reshape );
    surface_stress_ftensor += Fastor::einsum< Fastor::Index< i, j, k, l >,
                                              Fastor::Index< k, l >,
                                              Fastor::OIndex< i, j > >( Yn_H_inv_Fn_ijkl_ftensor,
                                                                        average_dSurface_strain_ftensor_reshape );
    surface_stress_ftensor -= Fastor::einsum< Fastor::Index< i >,
                                              Fastor::Index< i, j, k >,
                                              Fastor::OIndex< j, k > >( jumpU_ftensor, H_inv_nF_ijk_ftensor );

    // std::cout << "jumpU_ftensor:\n" << jumpU_ftensor << std::endl;
    // std::cout << "force_ftensor:\n" << force_ftensor << std::endl;
  };

} // namespace Marmot::Materials
