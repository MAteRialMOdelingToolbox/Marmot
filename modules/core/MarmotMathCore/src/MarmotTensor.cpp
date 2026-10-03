#include "Marmot/MarmotTensor.h"

namespace Marmot {
  namespace ContinuumMechanics::TensorUtility {
    Eigen::Matrix3d dyadicProduct( const Eigen::Vector3d& vector1, const Eigen::Vector3d& vector2 )
    {
      Eigen::Matrix3d dyade;

      for ( int i = 0; i < vector1.rows(); i++ )
        for ( int j = 0; j < vector1.rows(); j++ )
          dyade( i, j ) = vector1( i ) * vector2( j );

      return dyade;
    }
  } // namespace ContinuumMechanics::TensorUtility
} // namespace Marmot
