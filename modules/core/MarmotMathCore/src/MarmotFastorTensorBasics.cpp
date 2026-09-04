#include "Marmot/MarmotFastorTensorBasics.h"
namespace Marmot::TensorUtility::FastorTensors {
  StandardTensors::Tensor3333d invertMinorSymmetricFourthOrderTensor( const StandardTensors::Tensor3333d& C )
  {
    const double& sqrt2 = Constants::sqrt2;

    Fastor::Tensor< double, 6, 6 > C66;
    Fastor::Tensor< double, 9, 9 > C99 = Fastor::reshape< 9, 9 >( C );
    C99                                = 0.5 * ( C99 + transpose( C99 ) ); // make sure it's symmetric

    // Mapping for Mandel 6D basis: 11, 12, 13, 22, 23, 33
    int rowColMap[6] = { 0, 1, 2, 4, 5, 8 };

    // Mandel scaling: 1.0 for normal strains, sqrt(2) for shear strains
    double mandelScale[6] = { 1.0, sqrt2, sqrt2, 1.0, sqrt2, 1.0 };

    // Map to 6x6 Mandel stiffness symmetrically
    for ( size_t i = 0; i < 6; ++i ) {
      for ( size_t j = 0; j < 6; ++j ) {
        C66( i, j ) = C99( rowColMap[i], rowColMap[j] ) * mandelScale[i] * mandelScale[j];
      }
    }

    // invert reduced stiffness
    Fastor::Tensor< double, 6, 6 > C66inv = inverse< Fastor::InvCompType::BlockLU >( C66 );

    Fastor::Tensor< double, 9, 9 > C99inv( 0.0 );

    // 9x9 layout: 11, 12, 13, 21, 22, 23, 31, 32, 33
    int rowColMapFull[9] = { 0, 1, 2, 1, 3, 4, 2, 4, 5 };

    // Inverse Mandel scaling for mapping 6x6 back to 9x9: 1.0 for normal, 1/sqrt(2) for shear
    double invScaleFull[9] =
      { 1.0, 1.0 / sqrt2, 1.0 / sqrt2, 1.0 / sqrt2, 1.0, 1.0 / sqrt2, 1.0 / sqrt2, 1.0 / sqrt2, 1.0 };

    // Map back to 9x9 compliance symmetrically
    for ( size_t i = 0; i < 9; ++i ) {
      for ( size_t j = 0; j < 9; ++j ) {
        C99inv( i, j ) = C66inv( rowColMapFull[i], rowColMapFull[j] ) * invScaleFull[i] * invScaleFull[j];
      }
    }

    // reshape back to 4th order
    StandardTensors::Tensor3333d Cinv = Fastor::reshape< 3, 3, 3, 3 >( C99inv );
    return Cinv;
  }
} // namespace Marmot::TensorUtility::FastorTensors
