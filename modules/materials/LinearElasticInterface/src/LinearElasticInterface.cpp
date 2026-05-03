#include "Marmot/LinearElasticInterface.h"
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

  void LinearElasticInterface::initializeStateLayout() {}

  using namespace Marmot;
  using namespace Eigen;
  using namespace Fastor;

  using Tensor1D = Fastor::Tensor< double, 3 >;
  using Tensor2D = Fastor::Tensor< double, 3, 3 >;
  using Tensor3D = Fastor::Tensor< double, 3, 3, 3 >;
  using Tensor4D = Fastor::Tensor< double, 3, 3, 3, 3 >;

  LinearElasticInterface::LinearElasticInterface( const double* materialProperties,
                                                  int           nMaterialProperties,
                                                  int           materialNumber )
    : MarmotMaterialHypoElasticInterface::MarmotMaterialHypoElasticInterface( materialProperties,
                                                                              nMaterialProperties,
                                                                              materialNumber )
  {
    assert( nMaterialProperties == 4 || nMaterialProperties == 6 || nMaterialProperties == 10 );
  }

  void LinearElasticInterface::computeStress( double*       force,
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
    const double& E_0  = this->materialProperties[0];
    const double& nu_0 = this->materialProperties[1];
    const double& h    = this->materialProperties[2];

    // map to force, surface stress, displacement, surface strain, normal and tangent stiffness
    // use Fastor because we really need to use the einsum
    Fastor::Tensor< double, 3 >          force_ftensor( force );
    Fastor::Tensor< double, 3, 3 >       surface_stress_ftensor( surface_stress );
    Fastor::Tensor< double, 3, 3 >       H_inv_ij_ftensor( H_inv_ij );
    Fastor::Tensor< double, 3, 3, 3, 3 > Z_ijkl_ftensor( Z_ijkl );
    Fastor::Tensor< double, 3, 3, 3 >    H_inv_nF_ijk_ftensor( H_inv_nF_ijk );
    Fastor::Tensor< double, 3, 3, 3, 3 > Yn_H_inv_Fn_ijkl_ftensor( Yn_H_inv_Fn_ijkl );

    auto dU_ftensor_const              = Fastor::TensorMap< const double, 6, 1 >( dU );
    auto dSurface_strain_ftensor_const = Fastor::TensorMap< const double, 18, 1 >( dSurface_strain );
    auto normal_ftensor_const          = Fastor::TensorMap< const double, 3 >( normal );

    Fastor::Tensor< double, 6, 1 >  dU_ftensor( dU_ftensor_const.data() );
    Fastor::Tensor< double, 18, 1 > dSurface_strain_ftensor( dSurface_strain_ftensor_const.data() );
    Fastor::Tensor< double, 3 >     normal_ftensor( normal_ftensor_const.data() );

    auto [unitZ_ijkl,
          unitH_inv_ij,
          unitH_inv_Fn_ijk,
          unitYn_H_inv_Fn_ijkl] = calculateInterfaceMaterialParameters( normal_ftensor, nu_0 );

    // Assign the material matrices to a larger structure. (Not necessary ...)
    Z_ijkl_ftensor           = h * E_0 * unitZ_ijkl;
    Yn_H_inv_Fn_ijkl_ftensor = h * E_0 * unitYn_H_inv_Fn_ijkl;
    H_inv_ij_ftensor         = 1. / h * E_0 * unitH_inv_ij;
    H_inv_nF_ijk_ftensor     = E_0 * unitH_inv_Fn_ijk;

    // handle zero strain increment
    if ( Fastor::norm( dU_ftensor ) < 1e-14 && Fastor::norm( dSurface_strain_ftensor ) < 1e-14 ) {
      std::copy( H_inv_ij_ftensor.data(), H_inv_ij_ftensor.data() + 9, H_inv_ij );
      std::copy( Z_ijkl_ftensor.data(), Z_ijkl_ftensor.data() + 81, Z_ijkl );
      std::copy( H_inv_nF_ijk_ftensor.data(), H_inv_nF_ijk_ftensor.data() + 27, H_inv_nF_ijk );
      std::copy( Yn_H_inv_Fn_ijkl_ftensor.data(), Yn_H_inv_Fn_ijkl_ftensor.data() + 81, Yn_H_inv_Fn_ijkl );
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
    force_ftensor += Fastor::einsum< Fastor::Index< i, j, k >,
                                     Fastor::Index< j, k >,
                                     Fastor::OIndex< i > >( H_inv_nF_ijk_ftensor,
                                                            average_dSurface_strain_ftensor_reshape );
    surface_stress_ftensor += Fastor::einsum< Fastor::Index< i, j, k, l >,
                                              Fastor::Index< k, l >,
                                              Fastor::OIndex< i, j > >( Z_ijkl_ftensor,
                                                                        average_dSurface_strain_ftensor_reshape );
    surface_stress_ftensor += Fastor::einsum< Fastor::Index< i, j, k, l >,
                                              Fastor::Index< k, l >,
                                              Fastor::OIndex< i, j > >( Yn_H_inv_Fn_ijkl_ftensor,
                                                                        average_dSurface_strain_ftensor_reshape );
    surface_stress_ftensor += Fastor::einsum< Fastor::Index< i >,
                                              Fastor::Index< i, j, k >,
                                              Fastor::OIndex< j, k > >( jumpU_ftensor, H_inv_nF_ijk_ftensor );

    std::copy( force_ftensor.data(), force_ftensor.data() + 3, force );
    std::copy( surface_stress_ftensor.data(), surface_stress_ftensor.data() + 3 * 3, surface_stress );
    std::copy( H_inv_ij_ftensor.data(), H_inv_ij_ftensor.data() + 9, H_inv_ij );
    std::copy( Z_ijkl_ftensor.data(), Z_ijkl_ftensor.data() + 81, Z_ijkl );
    std::copy( H_inv_nF_ijk_ftensor.data(), H_inv_nF_ijk_ftensor.data() + 27, H_inv_nF_ijk );
    std::copy( Yn_H_inv_Fn_ijkl_ftensor.data(), Yn_H_inv_Fn_ijkl_ftensor.data() + 81, Yn_H_inv_Fn_ijkl );
  };

} // namespace Marmot::Materials
