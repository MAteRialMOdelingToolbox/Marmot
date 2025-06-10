#include "Marmot/MarmotGeostaticStress.h"
#include "Marmot/MarmotMath.h"
#include <Eigen/Dense>

namespace Marmot::GeostaticStress {

  std::tuple< double, double, double > getGeostaticStressFromLinearDistribution(
    const double* geostaticStressDefinition,
    double        coordinate_y )
  {

    const double sigY1 = geostaticStressDefinition[0];
    const double sigY2 = geostaticStressDefinition[2];
    const double y1    = geostaticStressDefinition[1];
    const double y2    = geostaticStressDefinition[3];

    const double S22Geostatic = Marmot::Math::linearInterpolation( coordinate_y, y1, y2, sigY1, sigY2 );
    const double S11Geostatic = geostaticStressDefinition[4] * S22Geostatic;
    const double S33Geostatic = geostaticStressDefinition[5] * S22Geostatic;

    return { S11Geostatic, S22Geostatic, S33Geostatic };
  }
} // namespace Marmot::GeostaticStress
