#pragma once 
#include "bftTypedefs.h"
//#include "bftVoigt.h"
#include "bftFunctions.h"

namespace bft{

    template <size_t materialTangentSize, size_t nIntegrationDependentStateVars>
        class NeunerAdaptiveSubstepper
        {
            public:
                typedef Matrix<double, materialTangentSize, materialTangentSize> TangentSizedMatrix;
                typedef Matrix<double, nIntegrationDependentStateVars, 1> IntegrationStateVector;
                enum ResultAccuracyModes{SingleStep, DoubleStep, RichardsonExtrapolation};

                NeunerAdaptiveSubstepper(double initialStepSize, double minimumStepSize, double scaleUpFactor, double scaleDownFactor, 
                                            int nPassesToIncrease, const Matrix6& Cel, double errorTolerance, ResultAccuracyModes resultMode,
                                                    const Vector6& stressOld,
                                                    const IntegrationStateVector& stateVarsOld);
                bool isFinished();
                double getNextSubstep();
                bool finishSubstep(const Vector6& resultStress, const TangentSizedMatrix& Tangent, const IntegrationStateVector& stateVars); 
                bool decreaseSubstepSize();
                
                void getAcceptedOldStressAndStateVars(Vector6& stress, IntegrationStateVector& stateVars);

                bft::Vector6 getResultingStress();
                bft::Matrix6 getResultingAlgorithmicTangent();
                IntegrationStateVector getResultingStateVars();


            private:
                enum SubsteppingState{FullStep, FirstHalfStep, SecondHalfStep};

                const double    initialStepSize,
                                minimumStepSize,
                                scaleUpFactor,
                                scaleDownFactor;
                const int       nPassesToIncrease;

                double currentProgress;
                double currentSubstepSize;
                int passedSubsteps;
                const Matrix6& Cel;

                const ResultAccuracyModes resultAccuracy;

                Vector6                 stressProgressHalf, stressProgressFull;
                IntegrationStateVector  stateProgressHalf,  stateProgressFull;
                TangentSizedMatrix      consistentTangentProgressHalf, consistentTangentProgressFull;

                TangentSizedMatrix elasticTangent;

                SubsteppingState currentState;

        };

}

namespace bft{
    template <size_t n, size_t nState>
    NeunerAdaptiveSubstepper<n, nState>::NeunerAdaptiveSubstepper(double initialStepSize, 
                                                    double minimumStepSize, 
                                                    double scaleUpFactor, 
                                                    double scaleDownFactor, 
                                                    int nPassesToIncrease,
                                                    const Matrix6& Cel, 
                                                    double errorTolerance,
                                                    ResultAccuracyModes resultMode,
                                                    const Vector6& stressOld,
                                                    const IntegrationStateVector& stateVarsOld):
                                                    initialStepSize(initialStepSize),
                                                    minimumStepSize(minimumStepSize),
                                                    scaleUpFactor(scaleUpFactor),
                                                    scaleDownFactor(scaleDownFactor),
                                                    nPassesToIncrease(nPassesToIncrease),
                                                    currentProgress(0.0),
                                                    currentSubstepSize(initialStepSize),
                                                    passedSubsteps(0),
                                                    Cel(Cel),
                                                    resultAccuracy(resultAccuracy)
    {
            elasticTangent = TangentSizedMatrix::Identity();
            elasticTangent.topLeftCorner(6,6) = Cel;

            consistentTangentProgressHalf= TangentSizedMatrix::Zero();
            consistentTangentProgressFull= TangentSizedMatrix::Zero();

            stressProgressHalf = stressOld;
            stressProgressFull = stressOld;

            stateProgressHalf = stateVarsOld;
            stateProgressFull = stateVarsOld;

            currentState = FullStep;
    }

    template<size_t n, size_t nState>
    bool NeunerAdaptiveSubstepper<n, nState>::isFinished()
    {
            return currentProgress >= 1.0;
    }

    template<size_t n, size_t nState>
    double NeunerAdaptiveSubstepper<n, nState>::getNextSubstep()
    {
        switch(currentState){
            case FullStep:{

                if(passedSubsteps >= nPassesToIncrease)
                        currentSubstepSize *= scaleUpFactor;

                const double remainingProgress = 1.0 - currentProgress;
                if( remainingProgress < currentSubstepSize)
                        currentSubstepSize = remainingProgress;            

                return currentSubstepSize;
                break;}
            case FirstHalfStep:{
                return currentSubstepSize/2;
                break;}
            case SecondHalfStep:{
                passedSubsteps++;
                return currentSubstepSize/2;
                break;}
        }
                //currentProgress += currentSubstepSize;
    }

    template<size_t n, size_t nState>
    bool NeunerAdaptiveSubstepper<n, nState>::decreaseSubstepSize()
    {
        switch(currentState)
        {
            case FullStep: {
                     
                    break; }
            case FirstHalfStep:

                           currentState = FullStep;
                           break;
            case SecondHalfStep:
                           currentState = FirstHalfStep;
                           break;
        }
            passedSubsteps = 0;
            currentSubstepSize *= scaleDownFactor;

            if(currentSubstepSize < minimumStepSize)
                    return warningToMSG("UMAT: Substepper: Minimal stepzsize reached");
            else
                    return notificationToMSG("UMAT: Substepper: Decreasing stepsize");
    }

    template<size_t n, size_t nState>
    bool NeunerAdaptiveSubstepper<n, nState>::finishSubstep()
    {
        switch(currentState){
            case FullStep:{
                          currentState = FirstHalfStep;
                          break;
                          }
            case FirstHalfStep:{
                           currentState = SecondHalfStep;
                           break;}
            case SecondHalfStep:{
                            //error Estimation
                            //
                            //Error large than tolerance?
                            //
                            //YES: decrease SubstepSize
                            //
                            //ELSE:

                            currentProgress += currentSubstepSize;
                            currentState = FullStep;
                            break;}
        }

    }

    //template<size_t n, size_t nState>
    //void NeunerAdaptiveSubstepper<n, nState>::extendConsistentTangent()
    //{
            //consistentTangent += currentSubstepSize * elasticTangent;
    //}

    template<int n>
    void PerezFougetSubstepper<n>::extendConsistentTangent(const TangentSizedMatrix& matTangent)
    {
            extendConsistentTangent();
            consistentTangent.applyOnTheLeft(matTangent);
    }

    template <int n>
    Matrix6 PerezFougetSubstepper<n>::consistentStiffness()
    {
            return consistentTangent.topLeftCorner(6,6);
    }
}
