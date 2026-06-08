#include "Marmot/VonMisesInterface.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotUtility.h"
#include "Marmot/MarmotVoigt.h"
#include "Marmot/VonMises.h"

#include "Fastor/Fastor.h"
#include "Marmot/MarmotInterfaceMaterialHelperFunctions.h"
#include "autodiff/forward/real.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <map>
#include <string>

using namespace Marmot;
using namespace Eigen;

namespace Marmot::Materials {

  VonMisesInterface::VonMisesInterface( const double* materialProperties, int nMaterialProperties, int materialNumber )
    : MarmotInterfaceMaterialHypoElastic( materialProperties, nMaterialProperties, materialNumber ),
      // clang-format off
      // elasticity parameters
      E_0( materialProperties[0] ),
      nu_0( materialProperties[1] ),
      h( materialProperties[2] ),
      // plasticity parameters
      yieldStress( materialProperties[3] ),
      HLin( materialProperties[4] ),
      deltaYieldStress( materialProperties[5] ),
      delta( materialProperties[6] ),
      // Re-map properties for VonMisesModel: [E, nu, yieldStress, HLin, deltaYieldStress, delta]
      // (skip h at index 2)
      vonMisesProps{ materialProperties[0], materialProperties[1],
                     materialProperties[3], materialProperties[4],
                     materialProperties[5], materialProperties[6] }
  // clang-format on
  {
    vonMisesModel = std::make_unique< VonMisesModel >( vonMisesProps.data(), 6, materialNumber );
    initializeStateLayout();
  }

  VonMisesInterface::~VonMisesInterface() = default;

  void VonMisesInterface::computeStress( State&               state,
                                         Tangents&            tangents,
                                         const Deformation&   deformation,
                                         const TimeIncrement& timeIncrement )
  {
    using namespace Marmot::Materials::InterfaceMaterialHelperFunctions;
    enum { i, j, k, l };
    auto scaled_forceFtensor         = Fastor::TensorMap< double, 3 >( state.force );
    auto scaled_averageStressFtensor = Fastor::TensorMap< double, 3, 3 >( state.surfaceStress );
    auto Q_ij_Ftensor_scaled         = Fastor::TensorMap< double, 3, 3 >( tangents.Q_ij );
    auto Z_ijkl_Ftensor_scaled       = Fastor::TensorMap< double, 3, 3, 3, 3 >( tangents.Z_ijkl );
    auto H_ijk_Ftensor_scaled        = Fastor::TensorMap< double, 3, 3, 3 >( tangents.H_ijk );
    auto Y_ijkl_Ftensor_scaled       = Fastor::TensorMap< double, 3, 3, 3, 3 >( tangents.Y_ijkl );

    double*       stateVars = state.stateVars;
    const double* timeOld   = timeIncrement.timeOld;
    const double  dT        = timeIncrement.dT;

    // map to force, surface stress, displacement, surface strain, normal and tangent stiffness
    // use Fastor because we really need to use the einsum

    Fastor::Tensor< double, 6, 1 >  dUFtensor( deformation.dU );
    Fastor::Tensor< double, 18, 1 > dSurfaceDispGradientFtensor( deformation.dSurfaceStrain );
    Fastor::Tensor< double, 3 >     normalFtensor( deformation.normal );

    // Evaluate average stress on the layer using Von Mises yield criterion.
    // vonMisesModel.computeStress updates averageStress and writes the new
    // elastoplastic tangent into C_ep (state var), which persists across increments.

    // Displacement jump: top(0:3) - bottom(3:6)
    Fastor::Tensor< double, 3 > dJumpU = dUFtensor( Fastor::seq( 0, 3 ), 0 ) -
                                         dUFtensor( Fastor::seq( 3, Fastor::last ), 0 );

    // Average surface strain: 0.5*(top(0:9) + bottom(9:18)), reshaped to 3x3
    Fastor::Tensor< double, 9, 1 >
         dSurfaceDispGradientAvgFlat = 0.5 * ( dSurfaceDispGradientFtensor( Fastor::seq( 0, 9 ), 0 ) +
                                            dSurfaceDispGradientFtensor( Fastor::seq( 9, Fastor::last ), 0 ) );
    auto dSurfaceDispGradientAvg     = Fastor::Tensor< double, 3, 3 >(
      Fastor::reshape< 3, 3 >( dSurfaceDispGradientAvgFlat ) );
    // Average displacement gradient: jump contribution (normal-to-layer) + surface strain
    Fastor::Tensor< double, 3, 3 >
      dU_kl_Jump = ( 1. / h ) *
                   Fastor::einsum< Fastor::Index< i >, Fastor::Index< j >, Fastor::OIndex< i, j > >( dJumpU,
                                                                                                     normalFtensor );

    // Symmetrize and include surface strain to obtain the full average strain increment
    Fastor::Tensor< double, 3, 3 > dDispGradAvg = ( dU_kl_Jump + dSurfaceDispGradientAvg );
    Fastor::Tensor< double, 3, 3 > dStrainAvg   = 0.5 * ( dDispGradAvg + Fastor::transpose( dDispGradAvg ) );

    // Convert 3x3 strain tensor to Voigt 6-vector (with factor 2 on shear components)
    Eigen::Map< const Eigen::Matrix< double, 3, 3, Eigen::RowMajor > > dStrainAvgEigen( dStrainAvg.data() );
    const Marmot::Vector6d dStrainAvgVoigt = Marmot::ContinuumMechanics::VoigtNotation::strainToVoigt(
      dStrainAvgEigen );

    // VonMisesModel updates stress in-place (incremental hypoelastic-plastic);
    // read the current 3x3 averageStress, symmetrize, convert to Voigt 6-vector
    Eigen::Map< const Eigen::Matrix< double, 3, 3, Eigen::RowMajor > > scaled_averageStressCurrent(
      state.surfaceStress );
    const Eigen::Matrix< double, 3, 3 > scaled_averageStressSym = 0.5 * ( scaled_averageStressCurrent +
                                                                          scaled_averageStressCurrent.transpose() );
    Marmot::Vector6d                    averageStressVoigt      = 1. / h *
                                          Marmot::ContinuumMechanics::VoigtNotation::stressToVoigt(
                                            scaled_averageStressSym );

    if ( stateVars == nullptr && getNumberOfRequiredStateVars() > 0 ) {
      throw std::runtime_error( MakeString() << __PRETTY_FUNCTION__ << ": state vars not provided." );
    }

    Marmot::Matrix6d C_ep  = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E_0, nu_0 );
    double&          kappa = stateLayout.getAs< double& >( stateVars, "kappa" );

    MarmotMaterialHypoElastic::state3D  vonMisesState{ averageStressVoigt, 0.0, &kappa };
    MarmotMaterialHypoElastic::timeInfo vonMisesTimeInfo{ timeOld[0], dT };

    vonMisesModel->computeStress( vonMisesState, C_ep, dStrainAvgVoigt, vonMisesTimeInfo );
    averageStressVoigt = vonMisesState.stress;

    auto [Z_ijkl_ep, Q_ij_ep, H_ijk_ep, Y_ijkl_ep] = calculateInterfaceMaterialParameters( normalFtensor, C_ep );

    Q_ij_Ftensor_scaled   = ( 1.0 / h ) * Q_ij_ep;
    Z_ijkl_Ftensor_scaled = (h)*Z_ijkl_ep;
    H_ijk_Ftensor_scaled  = H_ijk_ep;
    Y_ijkl_Ftensor_scaled = h * Y_ijkl_ep;

    // Expand the updated Voigt 6-vector back to full 3x3 and write into averageStress for the FE model
    // Write updated Voigt stress directly into state.surfaceStress.
    // Use RowMajor because Fastor::TensorMap reads the same memory layout.
    const Eigen::Matrix< double, 3, 3, Eigen::RowMajor >
      scaled_averageStressFull = h * Marmot::ContinuumMechanics::VoigtNotation::voigtToStress( averageStressVoigt );

    auto scaled_averageStressFullFtensor = Fastor::TensorMap< const double, 3, 3 >( scaled_averageStressFull.data() );

    scaled_averageStressFtensor = scaled_averageStressFullFtensor;

    scaled_forceFtensor = ( 1.0 / h ) *
                          Fastor::einsum< Fastor::Index< i, j >,
                                          Fastor::Index< j >,
                                          Fastor::OIndex< i > >( scaled_averageStressFtensor, normalFtensor );

    return;
  };
  double VonMisesInterface::getDensity()
  {
    if ( this->nMaterialProperties < 8 )
      throw std::runtime_error( MakeString() << __PRETTY_FUNCTION__ << ": No density given! nMaterialProperties < 8" );
    return this->materialProperties[7];
  }

} // namespace Marmot::Materials
