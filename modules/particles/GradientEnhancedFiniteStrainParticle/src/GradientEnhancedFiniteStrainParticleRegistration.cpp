#include "Marmot/GradientEnhancedFiniteStrainParticle.h"
#include "Marmot/MarmotMaterialPoint.h"
#include "Marmot/MarmotParticleLibrary.h"

namespace Marmot::Meshfree {

  using namespace MarmotLibrary;

  const static bool GradientEnhancedFiniteStrainParticle_PlaneStrain_isRegistered = MarmotLibrary::MarmotParticleFactory::
    registerParticle( "GradientEnhancedFiniteStrain/PlaneStrain/Point",
                      []( int                                                  cellID,
                          const double*                                        nodeCoordinates,
                          int                                                  sizeNodeCoordinates,
                          double                                               volume,
                          // MarmotMaterialPoint&                                 mp,
                          const std::string&                                 materialName,
                          const double*                                 materialProperties,
                          int                                         sizeMaterialProperties,
                          const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                        -> Marmot::Meshfree::MarmotParticle* {
                        return new GradientEnhancedFiniteStrainParticle< 2 >( cellID,
                                                                            nodeCoordinates,
                                                                            sizeNodeCoordinates,
                                                                            volume,
                                                                            // mp,
                                                                            materialName,
                                                                            materialProperties,
                                                                            sizeMaterialProperties,
                                                                            approximation );
                      } );

  const static bool GradientEnhancedFiniteStrainParticle_3D_isRegistered = MarmotLibrary::MarmotParticleFactory::
    registerParticle( "GradientEnhancedFiniteStrain/3D/Point",
                      []( int                                                  cellID,
                          const double*                                        nodeCoordinates,
                          int                                                  sizeNodeCoordinates,
                          double                                               volume,
                          const std::string&                                   materialName,
                          const double*                                        materialProperties,
                          int                                                  sizeMaterialProperties,
                          const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                        -> Marmot::Meshfree::MarmotParticle* {
                        return new GradientEnhancedFiniteStrainParticle< 3 >( cellID,
                                                                            nodeCoordinates,
                                                                            sizeNodeCoordinates,
                                                                            volume,
                                                                            materialName,
                                                                            materialProperties,
                                                                            sizeMaterialProperties,
                                                                            approximation );
                      } );

} // namespace Marmot::Meshfree
