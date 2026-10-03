#include "Marmot/DisplacementMaterialPoint.h"
#include "Marmot/MarmotMPMLibrary.h"

namespace Marmot::MaterialPoints::Registration {

  const static bool DisplacementPlaneStrainMaterialPoint_isRegistered = MarmotLibrary::MarmotMaterialPointFactory::
    registerMaterialPoint( "Displacement/PlaneStrain",
                           []( int           materialPointNumber,
                               const double* vertexCoordinates,
                               int           nVertexCoordinates,
                               double        volume ) -> MarmotMaterialPoint* {
                             return new DisplacementMaterialPoint2D( materialPointNumber,
                                                                     vertexCoordinates,
                                                                     nVertexCoordinates,
                                                                     volume );
                           } );

  const static bool Displacement3DMaterialPoint_isRegistered = MarmotLibrary::MarmotMaterialPointFactory::
    registerMaterialPoint( "Displacement/3D",
                           []( int           materialPointNumber,
                               const double* vertexCoordinates,
                               int           nVertexCoordinates,
                               double        volume ) -> MarmotMaterialPoint* {
                             return new DisplacementMaterialPoint3D( materialPointNumber,
                                                                     vertexCoordinates,
                                                                     nVertexCoordinates,
                                                                     volume );
                           } );

} // namespace Marmot::MaterialPoints::Registration
