#include "Marmot/VonMisesInterface.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotUtility.h"
#include "Marmot/MarmotViscoelasticity.h"
#include "Marmot/MarmotVoigt.h"
#include "Marmot/MarmotWiechertInterface.h"
#include "Marmot/VonMises.h"

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

  void VonMisesInterface::initializeStateLayout()
  {
    stateLayout.add( "kappa", 1 );
    stateLayout.finalize();
  }

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
                     materialProperties[5], materialProperties[6] },
      // Instantiate VonMisesModel once using the re-mapped properties stored in vonMisesProps
      vonMisesModel( vonMisesProps.data(), 6, materialNumber )
  // clang-format on
  {
    initializeStateLayout();
  }
  void VonMisesInterface::computeStress( double*       scaled_force,
                                         double*       scaled_averageStress,
                                         double*       Q_ij,
                                         double*       Z_ijkl,
                                         double*       H_ijk,
                                         double*       Y_ijkl,
                                         const double* dU,
                                         const double* dSurfaceDispGradient,
                                         const double* normal,
                                         const double* timeOld,
                                         const double  dT,
                                         double&       pNewDT )
  {
    using namespace Marmot::Materials::InterfaceMaterialHelperFunctions;
    enum { i, j, k, l };

    // map to force, surface stress, displacement, surface strain, normal and tangent stiffness
    // use Fastor because we really need to use the einsum

    Fastor::Tensor< double, 6, 1 >  dUFtensor( dU );
    Fastor::Tensor< double, 18, 1 > dSurfaceDispGradientFtensor( dSurfaceDispGradient );
    Fastor::Tensor< double, 3 >     normalFtensor( normal );

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
      scaled_averageStress );
    const Eigen::Matrix< double, 3, 3 > scaled_averageStressSym = 0.5 * ( scaled_averageStressCurrent +
                                                                          scaled_averageStressCurrent.transpose() );
    Marmot::Vector6d                    averageStressVoigt      = 1. / h *
                                          Marmot::ContinuumMechanics::VoigtNotation::stressToVoigt(
                                            scaled_averageStressSym );

    if ( stateVars == nullptr ) {
      throw std::runtime_error( MakeString() << __PRETTY_FUNCTION__ << ": state vars not assigned." );
    }

    Eigen::Matrix< double, 6, 6, Eigen::RowMajor >
            C_ep  = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E_0, nu_0 );
    double& kappa = stateLayout.getAs< double& >( stateVars, "kappa" );

    MarmotMaterialHypoElastic::state3D  vonMisesState{ averageStressVoigt, 0.0, &kappa };
    MarmotMaterialHypoElastic::timeInfo vonMisesTimeInfo{ timeOld[0], dT };

    vonMisesModel.computeStress( vonMisesState, C_ep.data(), dStrainAvgVoigt.data(), vonMisesTimeInfo );
    averageStressVoigt = vonMisesState.stress;

    auto [Z_ijkl_ep, Q_ij_ep, H_ijk_ep, Y_ijkl_ep] = calculateInterfaceMaterialParameters( normalFtensor, C_ep );

    Fastor::Tensor< double, 3, 3 >       Q_ij_Ftensor_scaled   = ( 1.0 / h ) * Q_ij_ep;
    Fastor::Tensor< double, 3, 3, 3, 3 > Z_ijkl_Ftensor_scaled = (h)*Z_ijkl_ep;
    Fastor::Tensor< double, 3, 3, 3 >    H_ijk_Ftensor_scaled  = H_ijk_ep;
    Fastor::Tensor< double, 3, 3, 3, 3 > Y_ijkl_Ftensor_scaled = h * Y_ijkl_ep;

    // Expand the updated Voigt 6-vector back to full 3x3 and write into averageStress for the FE model
    // voigtToStress returns column-major Eigen matrix; convert to row-major before copying
    // so that Fastor (row-major) reads the layout correctly
    const Eigen::Matrix< double, 3, 3, Eigen::RowMajor >
      scaled_averageStressFull = Marmot::ContinuumMechanics::VoigtNotation::voigtToStress( averageStressVoigt );

    std::copy( scaled_averageStressFull.data(), scaled_averageStressFull.data() + 9, scaled_averageStress );
    Fastor::TensorMap< const double, 3, 3 > scaled_averageStressFtensor( scaled_averageStress );
    Fastor::TensorMap< double, 3 >          scaled_forceFtensor( scaled_force );

    scaled_forceFtensor = Fastor::einsum< Fastor::Index< i, j >,
                                          Fastor::Index< j >,
                                          Fastor::OIndex< i > >( scaled_averageStressFtensor, normalFtensor );

    std::copy( Q_ij_Ftensor_scaled.data(), Q_ij_Ftensor_scaled.data() + 9, Q_ij );
    std::copy( Z_ijkl_Ftensor_scaled.data(), Z_ijkl_Ftensor_scaled.data() + 81, Z_ijkl );
    std::copy( H_ijk_Ftensor_scaled.data(), H_ijk_Ftensor_scaled.data() + 27, H_ijk );
    std::copy( Y_ijkl_Ftensor_scaled.data(), Y_ijkl_Ftensor_scaled.data() + 81, Y_ijkl );

    return;
  };
  void VonMisesInterface::assignStateVars( double* stateVars, int nStateVars )
  {
    if ( nStateVars < getNumberOfRequiredStateVars() )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": Not sufficient stateVars!" );

    MarmotInterfaceMaterialHypoElastic::assignStateVars( stateVars, nStateVars );
    this->stateVars = stateVars;
  }

  StateView VonMisesInterface::getStateView( const std::string& stateName )
  {
    if ( stateVars == nullptr ) {
      throw std::runtime_error( MakeString() << __PRETTY_FUNCTION__ << ": state vars not assigned." );
    }
    return MarmotInterfaceMaterialHypoElastic::getStateView( stateName, stateVars );
  }

  double VonMisesInterface::getDensity()
  {
    if ( this->nMaterialProperties < 8 )
      throw std::runtime_error( MakeString() << __PRETTY_FUNCTION__ << ": No density given! nMaterialProperties < 8" );
    return this->materialProperties[7];
  }

} // namespace Marmot::Materials
