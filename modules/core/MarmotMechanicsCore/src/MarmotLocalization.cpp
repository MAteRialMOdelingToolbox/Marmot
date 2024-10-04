#include "Marmot/MarmotLocalization.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTensor.h"
#include <cmath>
#include <iostream>

namespace Marmot {
  namespace ContinuumMechanics::LocalizationAnalysis {

    Marmot::Matrix3d computeAcousticTensor( const Marmot::Matrix6d& materialTangent,
                                            const Marmot::Vector3d& normalVector )
    {
      Marmot::Matrix3d acousticTensor;
      int              indexRow, indexCol;

      acousticTensor.setZero();

      for ( int i = 0; i < 3; i++ ) {
        for ( int j = 0; j < 3; j++ ) {
          indexRow = Marmot::ContinuumMechanics::TensorUtility::IndexNotation::toVoigt< 3 >( i, j );
          for ( int k = 0; k < 3; k++ ) {
            for ( int l = 0; l < 3; l++ ) {
              indexCol = Marmot::ContinuumMechanics::TensorUtility::IndexNotation::toVoigt< 3 >( k, l );

              acousticTensor( j, k ) += normalVector( i ) * materialTangent( indexRow, indexCol ) * normalVector( l );
            }
          }
        }
      }

      return acousticTensor;
    }

    bool localizationChecker( const Marmot::Matrix3d& acousticTensor )
    {
      double Tol = 1e-12;
      return acousticTensor.determinant() <= Tol ? true : false;
    }

    Marmot::Vector3d computeNormalVector( double alpha, double beta )
    {
      double           alphaRad = Marmot::Math::degToRad( alpha );
      double           betaRad  = Marmot::Math::degToRad( beta );
      Marmot::Vector3d nVec;

      nVec << cos( betaRad ) * cos( alphaRad ), cos( betaRad ) * sin( alphaRad ), sin( betaRad );
      return nVec;
    }

    double minimumDeterminantAcousticTensor( const Marmot::Matrix6d& materialTangent )
    {
      Marmot::Matrix3d acousticTensor;
      Marmot::Vector3d nVec;
      double           detQ;
      double           detQNew;

      detQ = 1e300;
      for ( double alpha = 0.; alpha <= 360.; alpha += 10. ) {
        for ( double beta = 0.; beta <= 90.; beta += 10. ) {
          nVec           = computeNormalVector( alpha, beta );
          acousticTensor = computeAcousticTensor( materialTangent, nVec );
          detQNew        = acousticTensor.determinant();

          if ( detQNew <= detQ )
            detQ = detQNew;
        }
      }

      return detQ;
    }
  } // namespace ContinuumMechanics::LocalizationAnalysis
} // namespace Marmot
