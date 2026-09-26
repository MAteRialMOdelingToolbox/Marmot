#include "Marmot/DisplacementParticle.h"
#include "Marmot/MarmotMaterialPoint.h"
#include "Marmot/MarmotParticleLibrary.h"

namespace Marmot::Meshfree {

  using namespace MarmotLibrary;

  const static bool DisplacementParticle_PlaneStrain_isRegistered = MarmotLibrary::MarmotParticleFactory::
    registerParticle( "Displacement/PlaneStrain/Point",
                      []( int           cellID,
                          const double* nodeCoordinates,
                          int           sizeNodeCoordinates,
                          double        volume,
                          // MarmotMaterialPoint&                                 mp,
                          const std::string&                                   materialName,
                          const double*                                        materialProperties,
                          int                                                  sizeMaterialProperties,
                          const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                        -> Marmot::Meshfree::MarmotParticle* {
                        return new DisplacementParticle< 2 >( cellID,
                                                              nodeCoordinates,
                                                              sizeNodeCoordinates,
                                                              volume,
                                                              // mp,
                                                              materialName,
                                                              materialProperties,
                                                              sizeMaterialProperties,
                                                              approximation );
                      } );

} // namespace Marmot::Meshfree
