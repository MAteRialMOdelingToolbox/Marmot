#include "Marmot/LinearViscoElasticInterface.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialHypoElasticInterface.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotUtility.h"
#include "Marmot/MarmotViscoelasticity.h"
#include "Marmot/MarmotVoigt.h"
#include "Marmot/MarmotWiechertInterface.h"

#include "Fastor/Fastor.h"
#include "Marmot/MarmotInterfaceMaterialHelperFunctions.h"
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <Fastor/expressions/linalg_ops/unary_norm_op.h>
#include <Fastor/tensor/TensorMap.h>

#include "autodiff/forward/real.hpp"
#include <iostream>
#include <map>
#include <string>

using namespace Marmot;
using namespace Eigen;

using Tensor1D = Fastor::Tensor< double, 3 >;
using Tensor2D = Fastor::Tensor< double, 3, 3 >;
using Tensor3D = Fastor::Tensor< double, 3, 3, 3 >;
using Tensor4D = Fastor::Tensor< double, 3, 3, 3, 3 >;

namespace Marmot::Materials {

  void LinearViscoElasticInterface::initializeStateLayout()
  {
    // State variables are managed manually by LinearViscoElasticInterfaceStateVarManager.
  }



  LinearViscoElasticInterface::LinearViscoElasticInterface( const double* materialProperties,
                                                            int           nMaterialProperties,
                                                            int           materialNumber )
    : MarmotMaterialHypoElasticInterface( materialProperties, nMaterialProperties, materialNumber ),
      // clang-format off
      // elasticity parameters
      E_0( materialProperties[0] ),
      nu_0( materialProperties[1] ),
      h( materialProperties[2] ),
      m( materialProperties[3] ),
      n( materialProperties[4] ),
      nMaxwell( static_cast< size_t >( materialProperties[5] ) ),
      minTau( materialProperties[6] ),
      timeToDays( materialProperties[7] )
  // clang-format on
  {
    relaxationTimes = Marmot::Materials::WiechertInterface::initializeRelaxationTimes( nMaxwell, m );
    elasticModuli   = Marmot::Materials::WiechertInterface::initializeElasticModuli( nMaxwell, n );

    using namespace Marmot::ContinuumMechanics::Viscoelasticity;
    auto phi_ = [&]( autodiff::Real< powerLawApproximationOrder, double > tau ) {
      return ComplianceFunctions::powerLaw( tau, m, n );
    };

    // elasticModuli_Ru =
    // Marmot::Materials::WiechertInterface::computeElasticModuli_Ru<powerLawApproximationOrder>(phiRu_,
    // retardationTimes_Ru); elasticModuli_Rs =
    // Marmot::Materials::WiechertInterface::computeElasticModuli_Rs<powerLawApproximationOrder>(phiRs_,
    // retardationTimes_Rs);

    zerothWiechertStiffness = 0.0; // m_Ru*(1. - n_Ru )*pow( 2., n_Ru )*pow(minTau_Ru/sqrt(10.), n_Ru);
  }

  void LinearViscoElasticInterface::computeStress( double*       force,
                                                   double*       surfaceStress,
                                                   double*       H_inv_ij,
                                                   double*       Z_ijkl,
                                                   double*       H_inv_nF_ijk,
                                                   double*       Yn_H_inv_Fn_ijkl,
                                                   const double* dU,
                                                   const double* dSurfaceStrain,
                                                   const double* normal,
                                                   const double* timeOld,
                                                   const double  dT,
                                                   double&       pNewDT )
  {
    using namespace Marmot::Materials::InterfaceMaterialHelperFunctions;

    // map to force, surface stress, displacement, surface strain, normal and tangent stiffness
    // use Fastor because we really need to use the einsum

    Fastor::Tensor< double, 3 >          forceFtensor( force );
    Fastor::Tensor< double, 3, 3 >       surfaceStressFtensor( surfaceStress );
    Fastor::Tensor< double, 3, 3 >       H_inv_ij_Ftensor( H_inv_ij );
    Fastor::Tensor< double, 3, 3, 3, 3 > Z_ijkl_Ftensor( Z_ijkl );
    Fastor::Tensor< double, 3, 3, 3 >    H_inv_nF_ijk_Ftensor( H_inv_nF_ijk );
    Fastor::Tensor< double, 3, 3, 3, 3 > Yn_H_inv_Fn_ijkl_Ftensor( Yn_H_inv_Fn_ijkl );
    auto                                 dUFtensorConst = Fastor::TensorMap< const double, 6, 1 >( dU );
    auto dSurfaceStrainFtensorConst                     = Fastor::TensorMap< const double, 18, 1 >( dSurfaceStrain );
    auto normalFtensorConst                             = Fastor::TensorMap< const double, 3 >( normal );

    Fastor::Tensor< double, 6, 1 >  dUFtensor( dUFtensorConst.data() );
    Fastor::Tensor< double, 18, 1 > dSurfaceStrainFtensor( dSurfaceStrainFtensorConst.data() );
    Fastor::Tensor< double, 3 >     normalFtensor( normalFtensorConst.data() );

    auto [unitZ_ijkl,
          unitH_inv_ij,
          unitH_inv_nF_ijk,
          unitYn_H_inv_Fn_ijkl] = calculateInterfaceMaterialParameters( normalFtensor, nu_0 );

    // Convert unit tensors to Voigt full matrices using helper functions
    Eigen::Matrix< double, 3, 3 > unitH_inv_voigt_full      = convert2ndOrderTensorToMatrix_3x3( unitH_inv_ij );
    Eigen::Matrix< double, 9, 9 > unitZ_voigt_full          = convert4thOrderTensorToMatrix_9x9( unitZ_ijkl );
    Eigen::Matrix< double, 9, 9 > unitYn_H_inv_Fn_ijkl_full = convert4thOrderTensorToMatrix_9x9( unitYn_H_inv_Fn_ijkl );
    Eigen::Matrix< double, 3, 9 > unitH_inv_nF_ijk_full_3_9 = convert3rdOrderTensorToMatrix_3x9( unitH_inv_nF_ijk );
    Eigen::Matrix< double, 9, 3 > unitH_inv_nF_ijk_full_9_3 = convert3rdOrderTensorToMatrix_9x3( unitH_inv_nF_ijk );

    // Assign the material matrices to a larger structure. (Not necessary ...)
    Eigen::Matrix< double, 21, 21 > Cel = Eigen::Matrix< double, 21, 21 >::Zero();

    // handle zero strain increment
    if ( Fastor::norm( dUFtensor ) < 1e-14 && Fastor::norm( dSurfaceStrainFtensor ) < 1e-14 && timeOld == 0 ) {
      std::cout << "Zero strain increment in LinearViscoElasticInterface material.\n";
      Z_ijkl_Ftensor           = -h * E_0 * unitZ_ijkl;
      Yn_H_inv_Fn_ijkl_Ftensor = h * E_0 * unitYn_H_inv_Fn_ijkl;
      H_inv_ij_Ftensor         = 1. / h * E_0 * unitH_inv_ij;
      H_inv_nF_ijk_Ftensor     = E_0 * unitH_inv_nF_ijk;

      std::copy( H_inv_ij_Ftensor.data(), H_inv_ij_Ftensor.data() + 9, H_inv_ij );
      std::copy( Z_ijkl_Ftensor.data(), Z_ijkl_Ftensor.data() + 81, Z_ijkl );
      std::copy( H_inv_nF_ijk_Ftensor.data(), H_inv_nF_ijk_Ftensor.data() + 27, H_inv_nF_ijk );
      std::copy( Yn_H_inv_Fn_ijkl_Ftensor.data(), Yn_H_inv_Fn_ijkl_Ftensor.data() + 81, Yn_H_inv_Fn_ijkl );
      return;
    }

    // visco elastic step
    Eigen::Ref< WiechertInterface::mapStateVarMatrix_force_uu > creepStateVars_force_uu(
      stateVarManager->MaxwellStateVars_force_uu );
    Eigen::Ref< WiechertInterface::mapStateVarMatrix_force_us > creepStateVars_force_us(
      stateVarManager->MaxwellStateVars_force_us );
    Eigen::Ref< WiechertInterface::mapStateVarMatrix_surface_stress_Z > creepStateVars_surface_stress_Z(
      stateVarManager->MaxwellStateVars_surface_stress_Z );
    Eigen::Ref< WiechertInterface::mapStateVarMatrix_surface_stress_Y > creepStateVars_surface_stress_Y(
      stateVarManager->MaxwellStateVars_surface_stress_Y );
    Eigen::Ref< WiechertInterface::mapStateVarMatrix_surface_stress_us > creepStateVars_surface_stress_us(
      stateVarManager->MaxwellStateVars_surface_stress_us );

    const double dTimeDays = dT * timeToDays;

    Marmot::Vector3d creep_force_uu_Increment          = Marmot::Vector3d::Zero();
    Marmot::Vector3d creep_force_us_Increment          = Marmot::Vector3d::Zero();
    Marmot::Vector9d creep_surface_stress_Z_Increment  = Marmot::Vector9d::Zero();
    Marmot::Vector9d creep_surface_stress_Y_Increment  = Marmot::Vector9d::Zero();
    Marmot::Vector9d creep_surface_stress_us_Increment = Marmot::Vector9d::Zero();

    double creep_Stiffness = 0;

    WiechertInterface::evaluateWiechert( dTimeDays,
                                         elasticModuli,
                                         relaxationTimes,
                                         creepStateVars_force_uu,
                                         creepStateVars_force_us,
                                         creepStateVars_surface_stress_Z,
                                         creepStateVars_surface_stress_Y,
                                         creepStateVars_surface_stress_us,
                                         creep_Stiffness,
                                         creep_force_uu_Increment,
                                         creep_force_us_Increment,
                                         creep_surface_stress_Z_Increment,
                                         creep_surface_stress_Y_Increment,
                                         creep_surface_stress_us_Increment,
                                         1.0 );

    using namespace Marmot::ContinuumMechanics::Viscoelasticity;

    // Evaluate effective compliances due to the displacement jump and the surface stress

    double barE              = E_0 + zerothWiechertStiffness + creep_Stiffness;
    Z_ijkl_Ftensor           = -h * barE * unitZ_ijkl;
    Yn_H_inv_Fn_ijkl_Ftensor = h * barE * unitYn_H_inv_Fn_ijkl;
    H_inv_ij_Ftensor         = 1. / h * barE * unitH_inv_ij;
    H_inv_nF_ijk_Ftensor     = barE * unitH_inv_nF_ijk;

    enum { i, j, k, l };
    // Calculate jump increment
    Tensor1D jumpUFtensor = dUFtensor( Fastor::seq( 0, 3 ), 0 ) - dUFtensor( Fastor::seq( 3, Fastor::last ), 0 );

    // Calculate average surface strain increment
    Fastor::Tensor< double, 9, 1 > averageDsurfaceStrainFtensor = 1. / 2. *
                                                                  ( dSurfaceStrainFtensor( Fastor::seq( 0, 9 ), 0 ) +
                                                                    dSurfaceStrainFtensor( Fastor::seq( 9,
                                                                                                        Fastor::last ),
                                                                                           0 ) );
    auto averageDsurfaceStrainFtensorReshape = Fastor::reshape< 3, 3 >( averageDsurfaceStrainFtensor );

    // Wrap Eigen vectors in Fastor TensorMaps (no copy)
    Fastor::TensorMap< double, 3 > creep_force_uu_IncrementFastor( creep_force_uu_Increment.data() );
    Fastor::TensorMap< double, 3 > creep_force_us_IncrementFastor( creep_force_us_Increment.data() );

    // Convert Eigen 9-element vectors to Fastor 3x3 tensors
    // Note: Fastor::reshape only works on Tensor, not TensorMap
    // Step 1: Wrap Eigen data in TensorMaps (no copy, just wraps the pointer)
    Fastor::TensorMap< double, 9 > creep_surface_stress_Z_Increment_map( creep_surface_stress_Z_Increment.data() );
    Fastor::TensorMap< double, 9 > creep_surface_stress_Y_Increment_map( creep_surface_stress_Y_Increment.data() );
    Fastor::TensorMap< double, 9 > creep_surface_stress_us_Increment_map( creep_surface_stress_us_Increment.data() );
    // Step 2: Copy TensorMap data into actual Tensors (required for reshape to work)
    Fastor::Tensor< double, 9 > creep_surface_stress_Z_Increment_tensor  = creep_surface_stress_Z_Increment_map;
    Fastor::Tensor< double, 9 > creep_surface_stress_Y_Increment_tensor  = creep_surface_stress_Y_Increment_map;
    Fastor::Tensor< double, 9 > creep_surface_stress_us_Increment_tensor = creep_surface_stress_us_Increment_map;
    // Step 3: Reshape 9-element tensors to 3×3 matrices
    auto creep_surface_stress_Z_IncrementFastor  = Fastor::reshape< 3, 3 >( creep_surface_stress_Z_Increment_tensor );
    auto creep_surface_stress_Y_IncrementFastor  = Fastor::reshape< 3, 3 >( creep_surface_stress_Y_Increment_tensor );
    auto creep_surface_stress_us_IncrementFastor = Fastor::reshape< 3, 3 >( creep_surface_stress_us_Increment_tensor );

    // std::cout<<"creep_force_uu_IncrementFastor:\n"<<creep_force_uu_IncrementFastor<<'\n';
    // std::cout<<"creep_force_us_IncrementFastor:\n"<<creep_force_us_IncrementFastor<<'\n';
    // std::cout<<"creep_surface_stress_Z_IncrementFastor:\n"<<creep_surface_stress_Z_IncrementFastor<<'\n';
    // std::cout<<"creep_surface_stress_Y_IncrementFastor:\n"<<creep_surface_stress_Y_IncrementFastor<<'\n';
    // std::cout<<"creep_surface_stress_us_IncrementFastor:\n"<<creep_surface_stress_us_IncrementFastor<<'\n';

    // std::cout<<"creep_Ru_increment_fastor:\n"<<creep_Ru_increment_fastor<<'\n';
    Tensor1D
      dForce_uu = Fastor::einsum< Fastor::Index< i, j >, Fastor::Index< j >, Fastor::OIndex< i > >( H_inv_ij_Ftensor,
                                                                                                    jumpUFtensor ) -
                  1. / h * creep_force_uu_IncrementFastor;

    Tensor1D dForce_us = Fastor::einsum< Fastor::Index< i, j, k >,
                                         Fastor::Index< j, k >,
                                         Fastor::OIndex< i > >( H_inv_nF_ijk_Ftensor,
                                                                averageDsurfaceStrainFtensorReshape ) -
                         creep_force_us_IncrementFastor;

    // std::cout<<"creep_Rs_increment_fastor:\n"<<creep_Rs_increment_fastor<<'\n';
    Tensor2D dSurfaceStress_Z_ij = Fastor::einsum< Fastor::Index< i, j, k, l >,
                                                   Fastor::Index< k, l >,
                                                   Fastor::OIndex< i, j > >( Z_ijkl_Ftensor,
                                                                             averageDsurfaceStrainFtensorReshape ) +
                                   h * creep_surface_stress_Z_IncrementFastor;

    Tensor2D dSurfaceStress_Y_ij = Fastor::einsum< Fastor::Index< i, j, k, l >,
                                                   Fastor::Index< k, l >,
                                                   Fastor::OIndex< i, j > >( Yn_H_inv_Fn_ijkl_Ftensor,
                                                                             averageDsurfaceStrainFtensorReshape ) -
                                   h * creep_surface_stress_Y_IncrementFastor;

    Tensor2D dSurfaceStress_us_ij = Fastor::einsum< Fastor::Index< i >,
                                                    Fastor::Index< i, j, k >,
                                                    Fastor::OIndex< j, k > >( jumpUFtensor, H_inv_nF_ijk_Ftensor ) -
                                    creep_surface_stress_us_IncrementFastor;

    forceFtensor += dForce_uu;
    forceFtensor -= dForce_us;

    surfaceStressFtensor -= dSurfaceStress_Z_ij;
    surfaceStressFtensor += dSurfaceStress_Y_ij;
    surfaceStressFtensor -= dSurfaceStress_us_ij;

    std::copy( forceFtensor.data(), forceFtensor.data() + 3, force );
    std::copy( surfaceStressFtensor.data(), surfaceStressFtensor.data() + 3 * 3, surfaceStress );

    std::copy( H_inv_ij_Ftensor.data(), H_inv_ij_Ftensor.data() + 9, H_inv_ij );
    std::copy( Z_ijkl_Ftensor.data(), Z_ijkl_Ftensor.data() + 81, Z_ijkl );
    std::copy( H_inv_nF_ijk_Ftensor.data(), H_inv_nF_ijk_Ftensor.data() + 27, H_inv_nF_ijk );
    std::copy( Yn_H_inv_Fn_ijkl_Ftensor.data(), Yn_H_inv_Fn_ijkl_Ftensor.data() + 81, Yn_H_inv_Fn_ijkl );

    // Use already available functionality convert Fastor tensors to Eigen matricfes/vectors
    // Tranform to Eigen matrices to work with the internal machinery of KelvinChainInterface ...
    Eigen::Map< Eigen::Matrix< double, 3, Eigen::RowMajor > > jumpUVoigtFull( jumpUFtensor.data() );
    Eigen::Map< Eigen::Matrix< double, 9, Eigen::RowMajor > > averageDsurfaceStrainFtensorReshapeVoigtFull(
      averageDsurfaceStrainFtensorReshape.data() );

    WiechertInterface::updateStateVarMatrix_force_uu( dTimeDays,
                                                      elasticModuli,
                                                      relaxationTimes,
                                                      creepStateVars_force_uu,
                                                      jumpUVoigtFull,
                                                      unitH_inv_voigt_full );

    WiechertInterface::updateStateVarMatrix_force_us( dTimeDays,
                                                      elasticModuli,
                                                      relaxationTimes,
                                                      creepStateVars_force_us,
                                                      averageDsurfaceStrainFtensorReshapeVoigtFull,
                                                      unitH_inv_nF_ijk_full_3_9 );

    WiechertInterface::updateStateVarMatrix_surface_stress_Z( dTimeDays,
                                                              elasticModuli,
                                                              relaxationTimes,
                                                              creepStateVars_surface_stress_Z,
                                                              averageDsurfaceStrainFtensorReshapeVoigtFull,
                                                              unitZ_voigt_full );

    WiechertInterface::updateStateVarMatrix_surface_stress_Y( dTimeDays,
                                                              elasticModuli,
                                                              relaxationTimes,
                                                              creepStateVars_surface_stress_Y,
                                                              averageDsurfaceStrainFtensorReshapeVoigtFull,
                                                              unitYn_H_inv_Fn_ijkl_full );

    WiechertInterface::updateStateVarMatrix_surface_stress_us( dTimeDays,
                                                               elasticModuli,
                                                               relaxationTimes,
                                                               creepStateVars_surface_stress_us,
                                                               jumpUVoigtFull,
                                                               unitH_inv_nF_ijk_full_3_9 );

    // std::cout<<"Updated state variables successfully.\n";
    // std::cout<<"creepStateVars_force_uu:\n"<<creepStateVars_force_uu<<'\n';
    // std::cout<<"creepStateVars_force_us:\n"<<creepStateVars_force_us<<'\n';
    // std::cout<<"creepStateVars_surface_stress_Z:\n"<<creepStateVars_surface_stress_Z<<'\n';
    // std::cout<<"creepStateVars_surface_stress_Y:\n"<<creepStateVars_surface_stress_Y<<'\n';
    // std::cout<<"creepStateVars_surface_stress_us:\n"<<creepStateVars_surface_stress_us<<'\n';
    return;
  };

  void LinearViscoElasticInterface::assignStateVars( double* stateVars_, int nStateVars )
  {
    if ( nStateVars < getNumberOfRequiredStateVars() )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": Not sufficient stateVars!" );

    this->stateVarManager = std::make_unique< LinearViscoElasticInterfaceStateVarManager >( stateVars_, nMaxwell );

    MarmotMaterialHypoElasticInterface::assignStateVars( stateVars_, nStateVars );
  }

  StateView LinearViscoElasticInterface::getStateView( const std::string& stateName )
  {
    return stateVarManager->getStateView( stateName );
  }

  int LinearViscoElasticInterface::getNumberOfRequiredStateVars() const
  {
    return LinearViscoElasticInterfaceStateVarManager::layout.nRequiredStateVars + 2 * 3 * nMaxwell + 3 * 9 * nMaxwell;
  }
} // namespace Marmot::Materials
