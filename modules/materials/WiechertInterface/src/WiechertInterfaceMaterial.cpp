#include "Marmot/WiechertInterfaceMaterial.h"

#include "Marmot/LinearViscoElasticWiechert.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotInterfaceMaterialHelperFunctions.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

#include "Fastor/Fastor.h"
#include <Eigen/Dense>

#include <stdexcept>

using namespace Marmot;
using namespace Marmot::FastorIndices;
using namespace Marmot::FastorStandardTensors;

namespace Marmot::Materials {

  WiechertInterfaceMaterial::WiechertInterfaceMaterial( const double* materialProperties,
                                                        int           nMaterialProperties,
                                                        int           materialNumber )
    : MarmotInterfaceMaterialHypoElastic( materialProperties, nMaterialProperties, materialNumber ),
      h( materialProperties[2] )
  {
    if ( nMaterialProperties < 8 ) {
      throw std::invalid_argument( "WiechertInterfaceMaterial requires at least 8 material properties." );
    }

    bulkMaterialProperties.reserve( nMaterialProperties - 1 );
    bulkMaterialProperties.push_back( materialProperties[0] );
    bulkMaterialProperties.push_back( materialProperties[1] );
    bulkMaterialProperties.insert( bulkMaterialProperties.end(),
                                   materialProperties + 3,
                                   materialProperties + nMaterialProperties );

    bulkMaterial = std::make_unique< LinearViscoElasticWiechert >( bulkMaterialProperties.data(),
                                                                   static_cast< int >( bulkMaterialProperties.size() ),
                                                                   materialNumber );

    stateLayout.add( "maxwellStateVars", bulkMaterial->getNumberOfRequiredStateVars() );
    stateLayout.finalize();
  }

  WiechertInterfaceMaterial::~WiechertInterfaceMaterial() = default;

  void WiechertInterfaceMaterial::computeStress( State&               state,
                                                 Tangents&            tangents,
                                                 const Deformation&   deformation,
                                                 const TimeIncrement& timeIncrement )
  {
    using namespace Marmot::Materials::InterfaceMaterialHelperFunctions;

    auto&       force          = state.force;
    auto&       surfaceStress  = state.surfaceStress;
    auto&       Q_ij           = tangents.Q_ij;
    auto&       Z_ijkl         = tangents.Z_ijkl;
    auto&       H_ijk          = tangents.H_ijk;
    auto&       Y_ijkl         = tangents.Y_ijkl;
    const auto& normal         = deformation.normal;
    auto        dU             = deformation.dU;
    auto        dSurfaceStrain = deformation.dSurfaceStrain;

    const Tensor3d jumpU = dU( Fastor::seq( 0, 3 ), 0 ) - dU( Fastor::seq( 3, Fastor::last ), 0 );

    const Tensor91d averageSurfaceGradientFlat = 0.5 * ( dSurfaceStrain( Fastor::seq( 0, 9 ), 0 ) +
                                                         dSurfaceStrain( Fastor::seq( 9, Fastor::last ), 0 ) );
    const Tensor33d averageSurfaceGradient( Fastor::reshape< 3, 3 >( averageSurfaceGradientFlat ) );
    const Tensor33d displacementGradient = ( 1. / h ) * Fastor::einsum< i, j, to_ij >( jumpU, normal ) +
                                           averageSurfaceGradient;
    const Tensor33d strainIncrement = 0.5 * ( displacementGradient + Fastor::transpose( displacementGradient ) );

    const Eigen::Map< const Eigen::Matrix< double, 3, 3, Eigen::RowMajor > > strainIncrementEigen(
      strainIncrement.data() );
    const Vector6d strainIncrementVoigt = ContinuumMechanics::VoigtNotation::strainToVoigt( strainIncrementEigen );

    const Eigen::Map< const Eigen::Matrix< double, 3, 3, Eigen::RowMajor > > scaledStressCurrent(
      state.surfaceStress.data() );
    const Eigen::Matrix3d scaledStressSym = 0.5 * ( scaledStressCurrent + scaledStressCurrent.transpose() );
    const Vector6d stressVoigt = ( 1. / h ) * ContinuumMechanics::VoigtNotation::stressToVoigt( scaledStressSym );

    Matrix6d                           tangent          = Matrix6d::Zero();
    double*                            maxwellStateVars = stateLayout.getPtr( state.stateVars, "maxwellStateVars" );
    MarmotMaterialHypoElastic::state3D bulkState{ stressVoigt, 0.0, 0.0, maxwellStateVars };
    const MarmotMaterialHypoElastic::timeInfo bulkTime{ timeIncrement.timeOld, timeIncrement.dT };

    bulkMaterial->computeStress( bulkState, tangent, strainIncrementVoigt, bulkTime );

    auto [Z, Q, H, Y] = calculateInterfaceMaterialParameters( normal, tangent );
    Q_ij              = ( 1. / h ) * Q;
    Z_ijkl            = h * Z;
    H_ijk             = H;
    Y_ijkl            = h * Y;

    const Eigen::Matrix< double, 3, 3, Eigen::RowMajor >
      scaledStress = h * ContinuumMechanics::VoigtNotation::voigtToStress( bulkState.stress );
    surfaceStress  = Tensor33d( scaledStress.data() );
    force          = ( 1. / h ) * Fastor::einsum< ij, j, to_i >( surfaceStress, normal );
  }

  void WiechertInterfaceMaterial::initializeYourself( double* stateVars, int )
  {
    bulkMaterial->initializeYourself( stateLayout.getPtr( stateVars, "maxwellStateVars" ),
                                      bulkMaterial->getNumberOfRequiredStateVars() );
  }

  double WiechertInterfaceMaterial::getDensity()
  {
    return bulkMaterial->getDensity( nullptr );
  }

} // namespace Marmot::Materials
