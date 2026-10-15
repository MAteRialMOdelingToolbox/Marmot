#include "Marmot/MarmotInterfaceMaterialHypoElastic.h"

#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotInterfaceMaterialHelperFunctions.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotUtils.h"
#include "Marmot/MarmotVoigt.h"

#include "Fastor/Fastor.h"
#include <Eigen/Dense>

#include <stdexcept>
#include <string>

using namespace Marmot;
using namespace Marmot::Math::Indices;
using namespace Marmot::Math::FastorTensors;

Marmot::Materials::InterfaceMaterialHypoElastic::InterfaceMaterialHypoElastic( const std::string& materialName,
                                                                               const double*      matProperties_,
                                                                               int                nMaterialProperties_,
                                                                               int                materialNumber_ )
  : materialProperties( matProperties_ ),
    nMaterialProperties( nMaterialProperties_ ),
    h( Marmot::Utilities::checkedMaterialProperty( matProperties_, nMaterialProperties_, 2 ) ),
    materialNumber( materialNumber_ )
{
  baseMaterialProperties.reserve( nMaterialProperties - 1 );
  baseMaterialProperties.push_back( materialProperties[0] );
  baseMaterialProperties.push_back( materialProperties[1] );
  baseMaterialProperties.insert( baseMaterialProperties.end(),
                                 materialProperties + 3,
                                 materialProperties + nMaterialProperties );

  baseMaterial = std::unique_ptr< Marmot::Materials::MaterialHypoElastic >(
    Marmot::Registry::MaterialHypoElasticFactory::createMaterial( materialName,
                                                                  baseMaterialProperties.data(),
                                                                  static_cast< int >( baseMaterialProperties.size() ),
                                                                  materialNumber ) );
  if ( !baseMaterial ) {
    throw std::invalid_argument( "Unknown base material for Marmot::Materials::InterfaceMaterialHypoElastic: " +
                                 materialName );
  }

  stateLayout.add( "baseMaterialStateVars", baseMaterial->getNumberOfRequiredStateVars() );
  stateLayout.finalize();
}

void Marmot::Materials::InterfaceMaterialHypoElastic::setCharacteristicElementLength( double length )
{
  characteristicElementLength = length;
  baseMaterial->setCharacteristicElementLength( length );
}

void Marmot::Materials::InterfaceMaterialHypoElastic::computeStress( State&               state,
                                                                     Tangents&            tangents,
                                                                     const Deformation&   deformation,
                                                                     const TimeIncrement& timeIncrement )
{
  if ( !baseMaterial ) {
    throw std::logic_error( "Marmot::Materials::InterfaceMaterialHypoElastic has no base material." );
  }

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

  const Tensor3d jumpU = dU( Fastor::seq( 0, 3 ) ) - dU( Fastor::seq( 3, Fastor::last ) );

  const Tensor9d  averageSurfaceGradientFlat = 0.5 * ( dSurfaceStrain( Fastor::seq( 0, 9 ) ) +
                                                      dSurfaceStrain( Fastor::seq( 9, Fastor::last ) ) );
  const Tensor33d averageSurfaceGradient( Fastor::reshape< 3, 3 >( averageSurfaceGradientFlat ) );
  const Tensor33d displacementGradient = ( 1. / h ) * Fastor::einsum< i, j, to_ij >( jumpU, normal ) +
                                         averageSurfaceGradient;
  const Tensor33d strainIncrement = 0.5 * ( displacementGradient + Fastor::transpose( displacementGradient ) );

  const Marmot::Math::Vector6d strainIncrementVoigt = ContinuumMechanics::VoigtNotation::strainToVoigt(
    Marmot::Math::mapEigenToFastor( strainIncrement ) );

  const auto                   scaledStressCurrent = Marmot::Math::mapEigenToFastor( state.surfaceStress );
  const Eigen::Matrix3d        scaledStressSym     = 0.5 * ( scaledStressCurrent + scaledStressCurrent.transpose() );
  const Marmot::Math::Vector6d stressVoigt         = ( 1. / h ) *
                                             ContinuumMechanics::VoigtNotation::stressToVoigt( scaledStressSym );

  Marmot::Math::Matrix6d tangent               = Marmot::Math::Matrix6d::Zero();
  double*                baseMaterialStateVars = stateLayout.getPtr( state.stateVars, "baseMaterialStateVars" );
  Marmot::Materials::MaterialHypoElastic::state3D        baseState{ stressVoigt, 0.0, 0.0, baseMaterialStateVars };
  const Marmot::Materials::MaterialHypoElastic::timeInfo timeInfo{ timeIncrement.timeOld + timeIncrement.dT,
                                                                   timeIncrement.dT };

  baseMaterial->computeStress( baseState, tangent, strainIncrementVoigt, timeInfo );

  auto [Z, Q, H, Y] = calculateInterfaceMaterialParameters( normal, tangent );
  Q_ij              = ( 1. / h ) * Q;
  Z_ijkl            = h * Z;
  H_ijk             = H;
  Y_ijkl            = h * Y;

  Marmot::Math::mapEigenToFastor( surfaceStress ) = h * ContinuumMechanics::VoigtNotation::voigtToStress(
                                                          baseState.stress );
  force = ( 1. / h ) * Fastor::einsum< ij, j, to_i >( surfaceStress, normal );
}

void Marmot::Materials::InterfaceMaterialHypoElastic::initializeYourself( double* stateVars, int )
{
  baseMaterial->initializeYourself( stateLayout.getPtr( stateVars, "baseMaterialStateVars" ),
                                    baseMaterial->getNumberOfRequiredStateVars() );
}

double Marmot::Materials::InterfaceMaterialHypoElastic::getDensity()
{
  return baseMaterial->getDensity( nullptr );
}
