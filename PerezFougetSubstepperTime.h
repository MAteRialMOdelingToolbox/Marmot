#pragma once
#include "MarmotTypedefs.h"

/*
11-22-2015 matthias neuner:
Modified Version of the Perez-Fouget Substepper,
to account for time-variant elastic Stiffness Tensor Cel(t_n+1)
NO changes in algorithmic formulation!

modifications:
(*) elastic substep consistent tangent matrix update needs current Cel(t) for update
(*) No need for Cel for object construction anymore
*/

namespace Marmot {

    template <int nSizeTangent>
    class PerezFougetSubstepperTime {

      public:
        typedef Eigen::Matrix<double, nSizeTangent, nSizeTangent> TangentSizedMatrix;

        PerezFougetSubstepperTime( double initialStepSize,
                                   double minimumStepSize,
                                   double scaleUpFactor,
                                   double scaleDownFactor,
                                   int    nPassesToIncrease );
        bool   isFinished();
        double getNextSubstep();
        double getFinishedProgress();
        bool   decreaseSubstepSize();

        void    extendConsistentTangent( const Matrix6& CelT );
        void    extendConsistentTangent( const Matrix6& CelT, const TangentSizedMatrix& matTangent );
        Matrix6 consistentStiffness();

      private:
        const double initialStepSize, minimumStepSize, scaleUpFactor, scaleDownFactor;
        const int    nPassesToIncrease;

        double currentProgress;
        double currentSubstepSize;
        int    passedSubsteps;

        TangentSizedMatrix elasticTangent;
        TangentSizedMatrix consistentTangent;
    };
} // namespace Marmot

#include "MarmotFunctions.h"

namespace Marmot {
    template <int s>
    PerezFougetSubstepperTime<s>::PerezFougetSubstepperTime( double initialStepSize,
                                                             double minimumStepSize,
                                                             double scaleUpFactor,
                                                             double scaleDownFactor,
                                                             int    nPassesToIncrease )
        : initialStepSize( initialStepSize ),
          minimumStepSize( minimumStepSize ),
          scaleUpFactor( scaleUpFactor ),
          scaleDownFactor( scaleDownFactor ),
          nPassesToIncrease( nPassesToIncrease ),
          currentProgress( 0.0 ),
          currentSubstepSize( initialStepSize ),
          passedSubsteps( 0 )

    {
        elasticTangent    = TangentSizedMatrix::Identity();
        consistentTangent = TangentSizedMatrix::Zero();
    }

    template <int s>
    bool PerezFougetSubstepperTime<s>::isFinished()
    {
        return currentProgress >= 1.0;
    }

    template <int s>
    double PerezFougetSubstepperTime<s>::getNextSubstep()
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

    template <int s>
    double PerezFougetSubstepperTime<s>::getFinishedProgress()
    {
        return currentProgress - currentSubstepSize;
    }

    template <int s>
    bool PerezFougetSubstepperTime<s>::decreaseSubstepSize()
    {
        currentProgress -= currentSubstepSize;
        passedSubsteps = 0;

        currentSubstepSize *= scaleDownFactor;

        if ( currentSubstepSize < minimumStepSize )
            return MarmotJournal::warningToMSG( "UMAT: Substepper: Minimal stepzsize reached" );
        else
            return MarmotJournal::notificationToMSG( "UMAT: Substepper: Decreasing stepsize" );
    }

    template <int s>
    void PerezFougetSubstepperTime<s>::extendConsistentTangent( const Matrix6& CelT )
    {

        elasticTangent.topLeftCorner( 6, 6 ) = CelT;
        consistentTangent += currentSubstepSize * elasticTangent;
    }

    template <int s>
    void PerezFougetSubstepperTime<s>::extendConsistentTangent( const Matrix6&            CelT,
                                                                const TangentSizedMatrix& matTangent )
    {
        extendConsistentTangent( CelT );
        consistentTangent.applyOnTheLeft( matTangent );
    }

    template <int s>
    Matrix6 PerezFougetSubstepperTime<s>::consistentStiffness()
    {
        return consistentTangent.topLeftCorner( 6, 6 );
    }
} // namespace Marmot
