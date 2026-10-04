#include "Marmot/MarmotInterfaceMaterialHypoElastic.h"

#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotInterfaceMaterialHelperFunctions.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

#include "Fastor/Fastor.h"
#include <Eigen/Dense>

#include <stdexcept>
#include <string>

using namespace Marmot;
using namespace Marmot::FastorIndices;
using namespace Marmot::FastorStandardTensors;

MarmotInterfaceMaterialHypoElastic::MarmotInterfaceMaterialHypoElastic( const std::string& materialName,
                                                                        const double*      matProperties_,
                                                                        int                nMaterialProperties_,
                                                                        int                materialNumber_ )
  : materialProperties( matProperties_ ), nMaterialProperties( nMaterialProperties_ ), materialNumber( materialNumber_ )
{
  if ( nMaterialProperties < 3 ) {
    throw std::invalid_argument(
      "MarmotInterfaceMaterialHypoElastic requires at least E, nu, and interface thickness h." );
  }

  h = materialProperties[2];
  baseMaterialProperties.reserve( nMaterialProperties - 1 );
  baseMaterialProperties.push_back( materialProperties[0] );
  baseMaterialProperties.push_back( materialProperties[1] );
  baseMaterialProperties.insert( baseMaterialProperties.end(),
                                 materialProperties + 3,
                                 materialProperties + nMaterialProperties );

  baseMaterial = std::unique_ptr< MarmotMaterialHypoElastic >(
    MarmotLibrary::MarmotMaterialHypoElasticFactory::createMaterial( materialName,
                                                                     baseMaterialProperties.data(),
                                                                     static_cast< int >(
                                                                       baseMaterialProperties.size() ),
                                                                     materialNumber ) );
  if ( !baseMaterial ) {
    throw std::invalid_argument( "Unknown base material for MarmotInterfaceMaterialHypoElastic: " + materialName );
  }

  stateLayout.add( "baseMaterialStateVars", baseMaterial->getNumberOfRequiredStateVars() );
  stateLayout.finalize();
}

void MarmotInterfaceMaterialHypoElastic::setCharacteristicElementLength( double length )
{
  characteristicElementLength = length;
  baseMaterial->setCharacteristicElementLength( length );
}

void MarmotInterfaceMaterialHypoElastic::computeStress( State&               state,
                                                        Tangents&            tangents,
                                                        const Deformation&   deformation,
                                                        const TimeIncrement& timeIncrement )
{
  if ( !baseMaterial ) {
    throw std::logic_error( "MarmotInterfaceMaterialHypoElastic has no base material." );
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

  const Eigen::Map< const Eigen::Matrix< double, 3, 3, Eigen::RowMajor > > strainIncrementEigen(
    strainIncrement.data() );
  const Vector6d strainIncrementVoigt = ContinuumMechanics::VoigtNotation::strainToVoigt( strainIncrementEigen );

  const Eigen::Map< const Eigen::Matrix< double, 3, 3, Eigen::RowMajor > > scaledStressCurrent(
    state.surfaceStress.data() );
  const Eigen::Matrix3d scaledStressSym = 0.5 * ( scaledStressCurrent + scaledStressCurrent.transpose() );
  const Vector6d        stressVoigt = ( 1. / h ) * ContinuumMechanics::VoigtNotation::stressToVoigt( scaledStressSym );

  Matrix6d tangent               = Matrix6d::Zero();
  double*  baseMaterialStateVars = stateLayout.getPtr( state.stateVars, "baseMaterialStateVars" );
  MarmotMaterialHypoElastic::state3D        baseState{ stressVoigt, 0.0, 0.0, baseMaterialStateVars };
  const MarmotMaterialHypoElastic::timeInfo timeInfo{ timeIncrement.timeOld + timeIncrement.dT, timeIncrement.dT };

  baseMaterial->computeStress( baseState, tangent, strainIncrementVoigt, timeInfo );

  auto [Z, Q, H, Y] = calculateInterfaceMaterialParameters( normal, tangent );
  Q_ij              = ( 1. / h ) * Q;
  Z_ijkl            = h * Z;
  H_ijk             = H;
  Y_ijkl            = h * Y;

  const Eigen::Matrix< double, 3, 3, Eigen::RowMajor > scaledStress = h *
                                                                      ContinuumMechanics::VoigtNotation::voigtToStress(
                                                                        baseState.stress );
  surfaceStress = Tensor33d( scaledStress.data() );
  force         = ( 1. / h ) * Fastor::einsum< ij, j, to_i >( surfaceStress, normal );
}

void MarmotInterfaceMaterialHypoElastic::initializeYourself( double* stateVars, int )
{
  baseMaterial->initializeYourself( stateLayout.getPtr( stateVars, "baseMaterialStateVars" ),
                                    baseMaterial->getNumberOfRequiredStateVars() );
}

double MarmotInterfaceMaterialHypoElastic::getDensity()
{
  return baseMaterial->getDensity( nullptr );
}
