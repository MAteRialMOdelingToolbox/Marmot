#pragma once

#include "Marmot/MarmotJournal.h"
#include <Eigen/Core>
#include <functional>

namespace Marmot::GeostaticStress {

  std::tuple< double, double, double > getGeostaticStressFromLinearDistribution( const double* geostaticStressDefintion,
                                                                                 double        coordinate_y );

}
