#pragma once
#include "bftFunctions.h"
#include "bftMath.h"
#include "bftTypedefs.h"

/* Substepper for elastoplastic materials, implicit return mapping version
 * Matthias Neuner (2016)
 * */

namespace bft {

    template <int nSizeMatTangent>
    class PerezFougetSubstepper {

      public:
        typedef Eigen::Matrix<double, nSizeMatTangent, nSizeMatTangent> TangentSizedMatrix;
        PerezFougetSubstepper( double         initialStepSize,
                               double         minimumStepSize,
                               double         scaleUpFactor,
                               double         scaleDownFactor,
                               int            nPassesToIncrease,
                               const Matrix6& Cel );
        bool   isFinished();
        double getNextSubstep();
        bool   decreaseSubstepSize();

        void    finishElasticSubstep();
        void    finishSubstep( const TangentSizedMatrix& dXdY );
        Matrix6 consistentStiffness();

      private:
        const double initialStepSize, minimumStepSize, scaleUpFactor, scaleDownFactor;
        const int    nPassesToIncrease;

        double currentProgress;
        double currentSubstepSize;
        int    passedSubsteps;

        const Matrix6& Cel;

        TangentSizedMatrix consistentTangent;
    };
} // namespace bft

namespace bft {
    template <int n>
    PerezFougetSubstepper<n>::PerezFougetSubstepper( double         initialStepSize,
                                                     double         minimumStepSize,
                                                     double         scaleUpFactor,
                                                     double         scaleDownFactor,
                                                     int            nPassesToIncrease,
                                                     const Matrix6& Cel )
        :

          initialStepSize( initialStepSize ),
          minimumStepSize( minimumStepSize ),
          scaleUpFactor( scaleUpFactor ),
          scaleDownFactor( scaleDownFactor ),
          nPassesToIncrease( nPassesToIncrease ),
          currentProgress( 0.0 ),
          currentSubstepSize( initialStepSize ),
          passedSubsteps( 0 ),
          Cel( Cel )
    {
        consistentTangent = TangentSizedMatrix::Zero();
    }

    template <int n>
    bool PerezFougetSubstepper<n>::isFinished()
    {
        // this is due to numerical accuracy ...
        return ( 1.0 - currentProgress ) <= 2e-16;
    }

    template <int n>
    double PerezFougetSubstepper<n>::getNextSubstep()
    {
        if ( passedSubsteps >= nPassesToIncrease )
            currentSubstepSize *= scaleUpFactor;

        const double remainingProgress = 1.0 - currentProgress;
        if ( remainingProgress < currentSubstepSize )
            currentSubstepSize = remainingProgress;

        passedSubsteps++;
        currentProgress += currentSubstepSize;

        return currentSubstepSize;
    }

    template <int n>
    bool PerezFougetSubstepper<n>::decreaseSubstepSize()
    {
        currentProgress -= currentSubstepSize;
        passedSubsteps = 0;

        currentSubstepSize *= scaleDownFactor;

        if ( currentSubstepSize < minimumStepSize )
            return BftJournal::warningToMSG( "UMAT: Substepper: Minimal stepzsize reached" );
        else
            return BftJournal::notificationToMSG( "UMAT: Substepper: Decreasing stepsize" );
    }

    template <int n>
    void PerezFougetSubstepper<n>::finishElasticSubstep()
    {
        consistentTangent += currentSubstepSize * TangentSizedMatrix::Identity();
    }

    template <int n>
    void PerezFougetSubstepper<n>::finishSubstep( const TangentSizedMatrix& dXdY )
    {
        finishElasticSubstep();
        consistentTangent.applyOnTheLeft( dXdY );

        for ( int i = 0; i < consistentTangent.rows(); i++ )
            for ( int j = 0; j < consistentTangent.cols(); j++ )
                if ( std::abs( consistentTangent( i, j ) ) < 1e-12 )
                    consistentTangent( i, j ) = 0.0;
    }

    template <int n>
    Matrix6 PerezFougetSubstepper<n>::consistentStiffness()
    {
        return consistentTangent.topLeftCorner( 6, 6 ) * Cel;
    }
} // namespace bft
