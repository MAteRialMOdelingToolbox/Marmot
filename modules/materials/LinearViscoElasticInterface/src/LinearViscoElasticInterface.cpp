#include "Marmot/LinearViscoElasticInterface.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotInterfaceMaterialHelperFunctions.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotWiechertInterface.h"

#include "Fastor/Fastor.h"
#include <Eigen/Dense>
#include <stdexcept>

using namespace Marmot::FastorStandardTensors;
using namespace Marmot::FastorIndices;

namespace Marmot::Materials {

  LinearViscoElasticInterface::LinearViscoElasticInterface( const double* materialProperties,
                                                            int           nMaterialProperties,
                                                            int           materialNumber )
    : MarmotInterfaceMaterialHypoElastic( materialProperties, nMaterialProperties, materialNumber ),
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
    if ( nMaxwell != 1 ) {
      throw std::invalid_argument( "LinearViscoElasticInterface currently supports exactly one Maxwell element." );
    }

    stateLayout.add( "MaxwellStateVars_force_uu", 3 * nMaxwell );
    stateLayout.add( "MaxwellStateVars_force_us", 3 * nMaxwell );
    stateLayout.add( "MaxwellStateVars_surface_stress_Z", 9 * nMaxwell );
    stateLayout.add( "MaxwellStateVars_surface_stress_Y", 9 * nMaxwell );
    stateLayout.add( "MaxwellStateVars_surface_stress_us", 9 * nMaxwell );
    stateLayout.finalize();

    relaxationTimes.resize( nMaxwell );
    relaxationTimes.setConstant( m );
    elasticModuli.resize( nMaxwell );
    elasticModuli.setConstant( n );

    // TODO(v26.05): Support a generalized Maxwell chain by generating distinct branch moduli and relaxation times
    // from the power-law parameters and minTau. Repeating the same (n, m) pair would only duplicate one branch.

    zerothWiechertStiffness = 0.0; // m_Ru*(1. - n_Ru )*pow( 2., n_Ru )*pow(minTau_Ru/sqrt(10.), n_Ru);
  }

  void LinearViscoElasticInterface::computeStress( State&               state,
                                                   Tangents&            tangents,
                                                   const Deformation&   deformation,
                                                   const TimeIncrement& timeIncrement )
  {
    using namespace Marmot::Materials::InterfaceMaterialHelperFunctions;

    double*      stateVars = state.stateVars;
    const double timeOld   = timeIncrement.timeOld;
    const double dT        = timeIncrement.dT;

    if ( stateVars == nullptr && getNumberOfRequiredStateVars() > 0 ) {
      throw std::runtime_error( "LinearViscoElasticInterface: state variables not provided." );
    }

    // use Fastor because we really need to use the einsum

    auto&       forceFtensor             = state.force;
    auto&       surfaceStressFtensor     = state.surfaceStress;
    auto&       H_inv_ij_Ftensor         = tangents.Q_ij;
    auto&       Z_ijkl_Ftensor           = tangents.Z_ijkl;
    auto&       H_inv_nF_ijk_Ftensor     = tangents.H_ijk;
    auto&       Yn_H_inv_Fn_ijkl_Ftensor = tangents.Y_ijkl;
    auto        dUFtensor                = deformation.dU;
    auto        dSurfaceStrainFtensor    = deformation.dSurfaceStrain;
    const auto& normalFtensor            = deformation.normal;

    auto [unitZ_ijkl,
          unitH_inv_ij,
          unitH_inv_nF_ijk,
          unitYn_H_inv_Fn_ijkl] = calculateInterfaceMaterialParameters( normalFtensor, nu_0 );

    // Convert unit tensors to Voigt full matrices using helper functions
    Eigen::Matrix< double, 3, 3 >
      unitH_inv_voigt_full = Marmot::ContinuumMechanics::TensorUtility::convert2ndOrderTensorToMatrix_3x3(
        unitH_inv_ij );
    Eigen::Matrix< double, 9, 9 >
      unitZ_voigt_full = Marmot::ContinuumMechanics::TensorUtility::convert4thOrderTensorToMatrix_9x9( unitZ_ijkl );
    Eigen::Matrix< double, 9, 9 >
      unitYn_H_inv_Fn_ijkl_full = Marmot::ContinuumMechanics::TensorUtility::convert4thOrderTensorToMatrix_9x9(
        unitYn_H_inv_Fn_ijkl );
    Eigen::Matrix< double, 3, 9 >
      unitH_inv_nF_ijk_full_3_9 = Marmot::ContinuumMechanics::TensorUtility::convert3rdOrderTensorToMatrix_3x9(
        unitH_inv_nF_ijk );
    // Assign the material matrices to a larger structure. (Not necessary ...)
    // handle zero strain increment
    if ( Fastor::norm( dUFtensor ) < 1e-14 && Fastor::norm( dSurfaceStrainFtensor ) < 1e-14 && timeOld == 0.0 ) {
      Z_ijkl_Ftensor           = -h * E_0 * unitZ_ijkl;
      Yn_H_inv_Fn_ijkl_Ftensor = h * E_0 * unitYn_H_inv_Fn_ijkl;
      H_inv_ij_Ftensor         = 1. / h * E_0 * unitH_inv_ij;
      H_inv_nF_ijk_Ftensor     = E_0 * unitH_inv_nF_ijk;

      return;
    }

    // visco elastic step
    auto creepStateVars_force_uu_map         = stateLayout.getAs< Eigen::Map< Eigen::MatrixXd > >( stateVars,
                                                                                           "MaxwellStateVars_force_uu",
                                                                                           3,
                                                                                           nMaxwell );
    auto creepStateVars_force_us_map         = stateLayout.getAs< Eigen::Map< Eigen::MatrixXd > >( stateVars,
                                                                                           "MaxwellStateVars_force_us",
                                                                                           3,
                                                                                           nMaxwell );
    auto creepStateVars_surface_stress_Z_map = stateLayout.getAs<
      Eigen::Map< Eigen::MatrixXd > >( stateVars, "MaxwellStateVars_surface_stress_Z", 9, nMaxwell );
    auto creepStateVars_surface_stress_Y_map = stateLayout.getAs<
      Eigen::Map< Eigen::MatrixXd > >( stateVars, "MaxwellStateVars_surface_stress_Y", 9, nMaxwell );
    auto creepStateVars_surface_stress_us_map = stateLayout.getAs<
      Eigen::Map< Eigen::MatrixXd > >( stateVars, "MaxwellStateVars_surface_stress_us", 9, nMaxwell );

    Eigen::Ref< WiechertInterface::StateVarMatrix_force_uu > creepStateVars_force_uu( creepStateVars_force_uu_map );
    Eigen::Ref< WiechertInterface::StateVarMatrix_force_us > creepStateVars_force_us( creepStateVars_force_us_map );
    Eigen::Ref< WiechertInterface::StateVarMatrix_surface_stress_Z > creepStateVars_surface_stress_Z(
      creepStateVars_surface_stress_Z_map );
    Eigen::Ref< WiechertInterface::StateVarMatrix_surface_stress_Y > creepStateVars_surface_stress_Y(
      creepStateVars_surface_stress_Y_map );
    Eigen::Ref< WiechertInterface::StateVarMatrix_surface_stress_us > creepStateVars_surface_stress_us(
      creepStateVars_surface_stress_us_map );

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

    // Evaluate effective compliances due to the displacement jump and the surface stress

    double barE              = E_0 + zerothWiechertStiffness + creep_Stiffness;
    Z_ijkl_Ftensor           = -h * barE * unitZ_ijkl;
    Yn_H_inv_Fn_ijkl_Ftensor = h * barE * unitYn_H_inv_Fn_ijkl;
    H_inv_ij_Ftensor         = 1. / h * barE * unitH_inv_ij;
    H_inv_nF_ijk_Ftensor     = barE * unitH_inv_nF_ijk;

    // Calculate jump increment
    Tensor3d jumpUFtensor = dUFtensor( Fastor::seq( 0, 3 ), 0 ) - dUFtensor( Fastor::seq( 3, Fastor::last ), 0 );

    // Calculate average surface strain increment
    Tensor91d averageDsurfaceStrainFtensor = 1. / 2. *
                                             ( dSurfaceStrainFtensor( Fastor::seq( 0, 9 ), 0 ) +
                                               dSurfaceStrainFtensor( Fastor::seq( 9, Fastor::last ), 0 ) );
    auto averageDsurfaceStrainFtensorReshape = Fastor::reshape< 3, 3 >( averageDsurfaceStrainFtensor );

    // Wrap Eigen vectors in Fastor TensorMaps (no copy)
    Tensor3d creep_force_uu_IncrementFastor( creep_force_uu_Increment.data() );
    Tensor3d creep_force_us_IncrementFastor( creep_force_us_Increment.data() );

    // Convert Eigen 9-element vectors to Fastor 3x3 tensors
    // Note: Fastor::reshape only works on Tensor, not TensorMap
    // Step 1: Wrap Eigen data in TensorMaps (no copy, just wraps the pointer)
    TensorMap9d creep_surface_stress_Z_Increment_map( creep_surface_stress_Z_Increment.data() );
    TensorMap9d creep_surface_stress_Y_Increment_map( creep_surface_stress_Y_Increment.data() );
    TensorMap9d creep_surface_stress_us_Increment_map( creep_surface_stress_us_Increment.data() );
    // Step 2: Copy TensorMap data into actual Tensors (required for reshape to work)
    Tensor9d creep_surface_stress_Z_Increment_tensor  = creep_surface_stress_Z_Increment_map;
    Tensor9d creep_surface_stress_Y_Increment_tensor  = creep_surface_stress_Y_Increment_map;
    Tensor9d creep_surface_stress_us_Increment_tensor = creep_surface_stress_us_Increment_map;
    // Step 3: Reshape 9-element tensors to 3×3 matrices
    auto creep_surface_stress_Z_IncrementFastor  = Fastor::reshape< 3, 3 >( creep_surface_stress_Z_Increment_tensor );
    auto creep_surface_stress_Y_IncrementFastor  = Fastor::reshape< 3, 3 >( creep_surface_stress_Y_Increment_tensor );
    auto creep_surface_stress_us_IncrementFastor = Fastor::reshape< 3, 3 >( creep_surface_stress_us_Increment_tensor );

    Tensor3d dForce_uu = Fastor::einsum< ij, j, to_i >( H_inv_ij_Ftensor, jumpUFtensor ) -
                         1. / h * creep_force_uu_IncrementFastor;

    Tensor3d dForce_us = Fastor::einsum< ijk, jk, to_i >( H_inv_nF_ijk_Ftensor, averageDsurfaceStrainFtensorReshape ) -
                         creep_force_us_IncrementFastor;

    // std::cout<<"creep_Rs_increment_fastor:\n"<<creep_Rs_increment_fastor<<'\n';
    Tensor33d dSurfaceStress_Z_ij = Fastor::einsum< ijkl, kl, to_ij >( Z_ijkl_Ftensor,
                                                                       averageDsurfaceStrainFtensorReshape ) +
                                    h * creep_surface_stress_Z_IncrementFastor;

    Tensor33d dSurfaceStress_Y_ij = Fastor::einsum< ijkl, kl, to_ij >( Yn_H_inv_Fn_ijkl_Ftensor,
                                                                       averageDsurfaceStrainFtensorReshape ) -
                                    h * creep_surface_stress_Y_IncrementFastor;

    Tensor33d dSurfaceStress_us_ij = Fastor::einsum< i, ijk, to_jk >( jumpUFtensor, H_inv_nF_ijk_Ftensor ) -
                                     creep_surface_stress_us_IncrementFastor;

    forceFtensor += dForce_uu;
    forceFtensor -= dForce_us;

    surfaceStressFtensor -= dSurfaceStress_Z_ij;
    surfaceStressFtensor += dSurfaceStress_Y_ij;
    surfaceStressFtensor -= dSurfaceStress_us_ij;

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

    return;
  };

  int LinearViscoElasticInterface::getNumberOfRequiredStateVars() const
  {
    return stateLayout.totalSize();
  }
} // namespace Marmot::Materials
