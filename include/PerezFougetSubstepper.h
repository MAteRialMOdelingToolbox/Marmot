#pragma once
#include "bftFunctions.h"
#include "bftTypedefs.h"

/* Deprecated, do not use anymore (!)
 * Matthias Neuner (2015)
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

        void    extendConsistentTangent();
        void    extendConsistentTangent( const TangentSizedMatrix& Tangent );
        Matrix6 consistentStiffness();

      private:
        const double initialStepSize, minimumStepSize, scaleUpFactor, scaleDownFactor;
        const int    nPassesToIncrease;

        double         currentProgress;
        double         currentSubstepSize;
        int            passedSubsteps;
        const Matrix6& Cel;

        TangentSizedMatrix elasticTangent;
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
        elasticTangent                       = TangentSizedMatrix::Identity();
        elasticTangent.topLeftCorner( 6, 6 ) = Cel;
        consistentTangent                    = TangentSizedMatrix::Zero();
    }

    template <int n>
    bool PerezFougetSubstepper<n>::isFinished()
    {
        return currentProgress >= 1.0;
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
    void PerezFougetSubstepper<n>::extendConsistentTangent()
    {
        consistentTangent += currentSubstepSize * elasticTangent;
    }

    template <int n>
    void PerezFougetSubstepper<n>::extendConsistentTangent( const TangentSizedMatrix& matTangent )
    {
        extendConsistentTangent();
        consistentTangent.applyOnTheLeft( matTangent );
    }

    template <int n>
    Matrix6 PerezFougetSubstepper<n>::consistentStiffness()
    {
        return consistentTangent.topLeftCorner( 6, 6 );
    }
} // namespace bft
