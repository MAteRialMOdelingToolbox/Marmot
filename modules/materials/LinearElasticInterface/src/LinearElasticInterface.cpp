#include "Marmot/LinearElasticInterface.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotUtility.h"
#include "Marmot/MarmotVoigt.h"

#include "Fastor/Fastor.h"
#include "Marmot/MarmotInterfaceMaterialHelperFunctions.h"
#include <Fastor/tensor/TensorMap.h>
#include <iostream>

namespace Marmot::Materials {

  void LinearElasticInterface::initializeStateLayout() {}

  using namespace Marmot;
  using namespace Eigen;
  using namespace Fastor;

  using Tensor1D = Marmot::FastorStandardTensors::Tensor3d;
  using Tensor2D = Marmot::FastorStandardTensors::Tensor33d;
  using Tensor3D = Marmot::FastorStandardTensors::Tensor333d;
  using Tensor4D = Marmot::FastorStandardTensors::Tensor3333d;

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

    std::cerr << "[LinearElasticInterface::computeStress] entered" << std::endl;
    std::cerr << "  this                    = " << this << std::endl;
    std::cerr << "  E_0                     = " << E_0 << std::endl;
    std::cerr << "  nu_0                    = " << nu_0 << std::endl;
    std::cerr << "  h                       = " << h << std::endl;
    std::cerr << "  state.force             = " << state.force << std::endl;
    std::cerr << "  state.surfaceStress     = " << state.surfaceStress << std::endl;
    std::cerr << "  tangents.Q_ij           = " << tangents.Q_ij << std::endl;
    std::cerr << "  tangents.Z_ijkl         = " << tangents.Z_ijkl << std::endl;
    std::cerr << "  tangents.H_ijk          = " << tangents.H_ijk << std::endl;
    std::cerr << "  tangents.Y_ijkl         = " << tangents.Y_ijkl << std::endl;
    std::cerr << "  deformation.dU          = " << deformation.dU << std::endl;
    std::cerr << "  deformation.dSurfaceStrain = " << deformation.dSurfaceStrain << std::endl;
    std::cerr << "  deformation.normal      = " << deformation.normal << std::endl;

    const double* timeOld = timeIncrement.timeOld;
    const double  dT      = timeIncrement.dT;
    double&       pNewDT  = timeIncrement.pNewDT;

    std::cerr << "  timeOld                 = " << timeOld << std::endl;
    std::cerr << "  dT                      = " << dT << std::endl;
    std::cerr << "  pNewDT                  = " << pNewDT << std::endl;

    (void)timeOld;
    (void)dT;
    (void)pNewDT;

    // map directly to force, surface stress, displacement, surface strain, normal and tangent stiffness
    // use Fastor because we really need to use the einsum
    std::cerr << "[LinearElasticInterface::computeStress] before TensorMap construction" << std::endl;
    auto force_ftensor                 = Fastor::TensorMap< double, 3 >( state.force );
    auto surface_stress_ftensor        = Fastor::TensorMap< double, 3, 3 >( state.surfaceStress );
    auto H_inv_ij_ftensor              = Fastor::TensorMap< double, 3, 3 >( tangents.Q_ij );
    auto Z_ijkl_ftensor                = Fastor::TensorMap< double, 3, 3, 3, 3 >( tangents.Z_ijkl );
    auto H_inv_nF_ijk_ftensor          = Fastor::TensorMap< double, 3, 3, 3 >( tangents.H_ijk );
    auto Yn_H_inv_Fn_ijkl_ftensor      = Fastor::TensorMap< double, 3, 3, 3, 3 >( tangents.Y_ijkl );
    auto dU_ftensor_const              = Fastor::TensorMap< const double, 6, 1 >( deformation.dU );
    auto dSurface_strain_ftensor_const = Fastor::TensorMap< const double, 18, 1 >( deformation.dSurfaceStrain );
    auto normal_ftensor_const          = Fastor::TensorMap< const double, 3 >( deformation.normal );

    std::cerr << "[LinearElasticInterface::computeStress] after TensorMap construction" << std::endl;

    Fastor::Tensor< double, 6, 1 >  dU_ftensor( dU_ftensor_const.data() );
    Fastor::Tensor< double, 18, 1 > dSurface_strain_ftensor( dSurface_strain_ftensor_const.data() );
    Fastor::Tensor< double, 3 >     normal_ftensor( normal_ftensor_const.data() );

    std::cerr << "[LinearElasticInterface::computeStress] copied input tensors" << std::endl;
    std::cerr << "  norm(dU)                 = " << Fastor::norm( dU_ftensor ) << std::endl;
    std::cerr << "  norm(dSurfaceStrain)     = " << Fastor::norm( dSurface_strain_ftensor ) << std::endl;
    std::cerr << "  norm(normal)             = " << Fastor::norm( normal_ftensor ) << std::endl;

    std::cerr << "[LinearElasticInterface::computeStress] before calculateInterfaceMaterialParameters" << std::endl;

    auto [unitZ_ijkl,
          unitH_inv_ij,
          unitH_inv_Fn_ijk,
          unitYn_H_inv_Fn_ijkl] = calculateInterfaceMaterialParameters( normal_ftensor, nu_0 );

    std::cerr << "[LinearElasticInterface::computeStress] after calculateInterfaceMaterialParameters" << std::endl;
    std::cerr << "[LinearElasticInterface::computeStress] before tangent assignment" << std::endl;

    Z_ijkl_ftensor           = h * E_0 * unitZ_ijkl;
    Yn_H_inv_Fn_ijkl_ftensor = h * E_0 * unitYn_H_inv_Fn_ijkl;
    H_inv_ij_ftensor         = 1. / h * E_0 * unitH_inv_ij;
    H_inv_nF_ijk_ftensor     = E_0 * unitH_inv_Fn_ijk;

    std::cerr << "[LinearElasticInterface::computeStress] after tangent assignment" << std::endl;

    // handle zero strain increment
    if ( Fastor::norm( dU_ftensor ) < 1e-14 && Fastor::norm( dSurface_strain_ftensor ) < 1e-14 ) {
      std::cerr << "[LinearElasticInterface::computeStress] zero increment, returning" << std::endl;
      return;
    }
    std::cerr << "[LinearElasticInterface::computeStress] nonzero increment, elastic step" << std::endl;
    // elastic step
    enum { i, j, k, l };

    std::cerr << "[LinearElasticInterface::computeStress] before jumpU computation" << std::endl;
    Tensor1D jumpU_ftensor = dU_ftensor( Fastor::seq( 0, 3 ), 0 ) - dU_ftensor( Fastor::seq( 3, Fastor::last ), 0 );
    std::cerr << "[LinearElasticInterface::computeStress] after jumpU computation, norm(jumpU) = "
              << Fastor::norm( jumpU_ftensor ) << std::endl;

    Fastor::Tensor< double, 9, 1 >
      average_dSurface_strain_ftensor = 1. / 2. *
                                        ( dSurface_strain_ftensor( Fastor::seq( 0, 9 ), 0 ) +
                                          dSurface_strain_ftensor( Fastor::seq( 9, Fastor::last ), 0 ) );
    auto average_dSurface_strain_ftensor_reshape = Fastor::reshape< 3, 3 >( average_dSurface_strain_ftensor );
    std::cerr << "[LinearElasticInterface::computeStress] after average surface strain, norm = "
              << Fastor::norm( average_dSurface_strain_ftensor_reshape ) << std::endl;
    std::cerr << "[LinearElasticInterface::computeStress] before force update 1" << std::endl;
    force_ftensor += Fastor::einsum< Fastor::Index< i, j >, Fastor::Index< j >, Fastor::OIndex< i > >( H_inv_ij_ftensor,
                                                                                                       jumpU_ftensor );
    std::cerr << "[LinearElasticInterface::computeStress] after force update 1" << std::endl;
    std::cerr << "[LinearElasticInterface::computeStress] before force update 2" << std::endl;
    force_ftensor += Fastor::einsum< Fastor::Index< i, j, k >,
                                     Fastor::Index< j, k >,
                                     Fastor::OIndex< i > >( H_inv_nF_ijk_ftensor,
                                                            average_dSurface_strain_ftensor_reshape );
    std::cerr << "[LinearElasticInterface::computeStress] after force update 2" << std::endl;
    std::cerr << "[LinearElasticInterface::computeStress] before surface stress update 1" << std::endl;
    surface_stress_ftensor += Fastor::einsum< Fastor::Index< i, j, k, l >,
                                              Fastor::Index< k, l >,
                                              Fastor::OIndex< i, j > >( Z_ijkl_ftensor,
                                                                        average_dSurface_strain_ftensor_reshape );
    std::cerr << "[LinearElasticInterface::computeStress] after surface stress update 1" << std::endl;
    std::cerr << "[LinearElasticInterface::computeStress] before surface stress update 2" << std::endl;
    surface_stress_ftensor += Fastor::einsum< Fastor::Index< i, j, k, l >,
                                              Fastor::Index< k, l >,
                                              Fastor::OIndex< i, j > >( Yn_H_inv_Fn_ijkl_ftensor,
                                                                        average_dSurface_strain_ftensor_reshape );
    std::cerr << "[LinearElasticInterface::computeStress] after surface stress update 2" << std::endl;
    std::cerr << "[LinearElasticInterface::computeStress] before surface stress update 3" << std::endl;
    surface_stress_ftensor += Fastor::einsum< Fastor::Index< i >,
                                              Fastor::Index< i, j, k >,
                                              Fastor::OIndex< j, k > >( jumpU_ftensor, H_inv_nF_ijk_ftensor );
    std::cerr << "[LinearElasticInterface::computeStress] after surface stress update 3" << std::endl;
    std::cerr << "[LinearElasticInterface::computeStress] leaving" << std::endl;
  };

} // namespace Marmot::Materials
