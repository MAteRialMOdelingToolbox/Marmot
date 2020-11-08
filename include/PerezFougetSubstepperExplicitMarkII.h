#pragma once
#include "MarmotFunctions.h"
#include "MarmotTypedefs.h"

/* Substepper for semi-explicit elastoplastic materials
 * Matthias Neuner (2016)
 * */

namespace marmot {

    template <int nSizeMatTangent>
    class PerezFougetSubstepper {

      public:
        typedef Eigen::Matrix<double, nSizeMatTangent, nSizeMatTangent> TangentSizedMatrix;
        typedef Eigen::Matrix<double, nSizeMatTangent, 6>               MatrixStateStrain;
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
        void    finishSubstep( const TangentSizedMatrix& dXdY, const TangentSizedMatrix& dYdXOld );
        Matrix6 consistentStiffness();

      private:
        const double initialStepSize, minimumStepSize, scaleUpFactor, scaleDownFactor;
        const int    nPassesToIncrease;

        double currentProgress;
        double currentSubstepSize;
        int    passedSubsteps;

        const Matrix6&     Cel;
        MatrixStateStrain  I76;
        TangentSizedMatrix I77;

        MatrixStateStrain consistentTangent;
    };
} // namespace marmot

namespace marmot {
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
        consistentTangent = MatrixStateStrain::Zero();
        I76.setZero();
        I76.topLeftCorner( 6, 6 ) = Matrix6::Identity();
        I77                       = TangentSizedMatrix::Identity();
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
            return MarmotJournal::warningToMSG( "UMAT: Substepper: Minimal stepzsize reached" );
        else
            return MarmotJournal::notificationToMSG( "UMAT: Substepper: Decreasing stepsize" );
    }

    template <int n>
    void PerezFougetSubstepper<n>::finishElasticSubstep()
    {
        consistentTangent += currentSubstepSize * I76 * Cel;
    }

    template <int n>
    void PerezFougetSubstepper<n>::finishSubstep( const TangentSizedMatrix& dXdY, const TangentSizedMatrix& dYdXOld )
    {
        consistentTangent.applyOnTheLeft( I77 - dYdXOld );
        consistentTangent += currentSubstepSize * I76 * Cel;
        consistentTangent.applyOnTheLeft( dXdY );
    }

    template <int n>
    Matrix6 PerezFougetSubstepper<n>::consistentStiffness()
    {
        return consistentTangent.topLeftCorner( 6, 6 );
    }
} // namespace marmot
