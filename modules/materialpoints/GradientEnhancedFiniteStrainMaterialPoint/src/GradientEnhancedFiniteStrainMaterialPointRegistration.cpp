#include "Marmot/GradientEnhancedFiniteStrainMaterialPoint.h"
#include "Marmot/MarmotMPMLibrary.h"

namespace Marmot::MaterialPoints::Registration {

  const static bool
    GradientEnhancedFiniteStrainMaterialPoint2D_isRegistered = MarmotLibrary::MarmotMaterialPointFactory::
      registerMaterialPoint( "GradientEnhancedFiniteStrain/PlaneStrain",
                             []( int mpNumber, const double* vertexCoordinates, int nVertexCoordinates, double volume )
                               -> MarmotMaterialPoint* {
                               return new GradientEnhancedFiniteStrainMaterialPoint2D( mpNumber,
                                                                                       vertexCoordinates,
                                                                                       nVertexCoordinates,
                                                                                       volume );
                             } );

  const static bool
    GradientEnhancedFiniteStrainMaterialPoint3D_isRegistered = MarmotLibrary::MarmotMaterialPointFactory::
      registerMaterialPoint( "GradientEnhancedFiniteStrain/3D",
                             []( int mpNumber, const double* vertexCoordinates, int nVertexCoordinates, double volume )
                               -> MarmotMaterialPoint* {
                               return new GradientEnhancedFiniteStrainMaterialPoint3D( mpNumber,
                                                                                       vertexCoordinates,
                                                                                       nVertexCoordinates,
                                                                                       volume );
                             } );

} // namespace Marmot::MaterialPoints::Registration
