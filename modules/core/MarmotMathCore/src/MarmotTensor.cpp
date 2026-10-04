#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotFastorTensorBasics.h"

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

    Matrix99d convert4thOrderTensorToMatrix_9x9( const Fastor::Tensor< double, 3, 3, 3, 3 >& tensor )
    {
      const FastorStandardTensors::Tensor99d tensorAsFastorMatrix = Fastor::reshape< 9, 9 >( tensor );

      // Eigen maps Fastor's row-major storage as column-major, so transpose back to preserve (ij)(kl) ordering.
      return Eigen::Map< const Matrix99d >( tensorAsFastorMatrix.data() ).transpose();
    }

    Matrix93d convert3rdOrderTensorToMatrix_9x3( const Fastor::Tensor< double, 3, 3, 3 >& tensor )
    {
      const FastorStandardTensors::Tensor93d tensorAsFastorMatrix = Fastor::reshape< 9, 3 >( tensor );

      return Eigen::Map< const Matrix39d >( tensorAsFastorMatrix.data() ).transpose();
    }

    Matrix39d convert3rdOrderTensorToMatrix_3x9( const Fastor::Tensor< double, 3, 3, 3 >& tensor )
    {
      const FastorStandardTensors::Tensor39d tensorAsFastorMatrix = Fastor::reshape< 3, 9 >( tensor );

      return Eigen::Map< const Matrix93d >( tensorAsFastorMatrix.data() ).transpose();
    }

    Matrix3d convert2ndOrderTensorToMatrix_3x3( const Fastor::Tensor< double, 3, 3 >& tensor )
    {
      const FastorStandardTensors::Tensor33d tensorAsFastorMatrix = Fastor::reshape< 3, 3 >( tensor );

      return Eigen::Map< const Matrix3d >( tensorAsFastorMatrix.data() ).transpose();
    }

  } // namespace ContinuumMechanics::TensorUtility
} // namespace Marmot
