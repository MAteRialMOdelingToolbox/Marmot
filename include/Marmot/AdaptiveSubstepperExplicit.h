/* ---------------------------------------------------------------------
 *                                       _   
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_ 
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_ 
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 * 
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck, 
 * 2020 - today
 * 
 * festigkeitslehre@uibk.ac.at
 * 
 * Matthias Neuner matthias.neuner@uibk.ac.at
 * 
 * This file is part of the MAteRialMOdellingToolbox (marmot).
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of marmot.
 * ---------------------------------------------------------------------
 */

#pragma once
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTypedefs.h"

/* Adaptive Substepper, employing an error estimation and
 * Richardson Extrapolation for an (semi)-Explicit Return  Mapping algorithm
 *
 * Matthias Neuner (2016), developed within the DK-CIM collaboration
 * */

namespace Marmot::NumericalAlgorithms {

    template <size_t materialTangentSize, size_t nIntegrationDependentStateVars>
    class AdaptiveSubstepperExplicit {
      public:
        typedef Eigen::Matrix<double, materialTangentSize, materialTangentSize> TangentSizedMatrix;
        typedef Eigen::Matrix<double, materialTangentSize, 6>                   MatrixStateStrain;
        typedef Eigen::Matrix<double, nIntegrationDependentStateVars, 1>        IntegrationStateVector;

        AdaptiveSubstepperExplicit( double         initialStepSize,
                                    double         minimumStepSize,
                                    double         maxScaleUpFactor,
                                    double         scaleDownFactor,
                                    double         integrationErrorTolerance,
                                    int            nPassesToIncrease,
                                    const Matrix6d& Cel );
        void   setConvergedProgress( const Marmot::Vector6d& stressOld, const IntegrationStateVector& stateVarsOld );
        bool   isFinished();
        double getNextSubstep();
        int    getNumberOfSubsteps();
        int    getNumberDiscardedSubstepsDueToError();
        bool   finishSubstep( const Marmot::Vector6d&           resultStress,
                              const TangentSizedMatrix&     dXdY,
                              const TangentSizedMatrix&     dYdXOld,
                              const IntegrationStateVector& stateVars );

        void finishElasticSubstep( const Marmot::Vector6d& resultStress );
        bool discardSubstep();
        bool repeatSubstep( double decrementationFactor );

        void    getConvergedProgress( Marmot::Vector6d& stress, IntegrationStateVector& stateVars );
        Matrix6d getCurrentTangentOperator();
        void    getResults( Marmot::Vector6d& stress, Matrix6d& consistentTangent, IntegrationStateVector& stateVars );

      private:
        const double initialStepSize, minimumStepSize, maxScaleUpFactor, scaleDownFactor, integrationErrorTolerance;
        const int    nPassesToIncrease;
        const bool   ignoreErrorToleranceOnMinimumStepSize;

        const Matrix6d& Cel;
        double         currentProgress;
        double         currentSubstepSize;
        int            passedSubsteps;
        int            substepIndex;
        int            discardedDueToError;

        // internal storages for the progress of the total increment
        Marmot::Vector6d           stressProgress;
        IntegrationStateVector stateProgress;
        MatrixStateStrain      consistentTangentProgress;

        // temporal storages, which are used until a cycle full/half/half has finished successfully
        Marmot::Vector6d           stressProgressHalfTemp, stressProgressFullTemp;
        IntegrationStateVector stateProgressHalfTemp, stateProgressFullTemp;
        MatrixStateStrain      consistentTangentProgressHalfTemp, consistentTangentProgressFullTemp;

        const TangentSizedMatrix& IXX();
        const MatrixStateStrain&  IX6();

        enum SubsteppingState { FullStep, FirstHalfStep, SecondHalfStep };
        SubsteppingState currentState;

        bool acceptSubstepWithFullStepOnly();
        bool splitCurrentSubstep();
    };

} // namespace Marmot::NumericalAlgorithms

namespace Marmot::NumericalAlgorithms {

    template <size_t n, size_t nState>
    AdaptiveSubstepperExplicit<n, nState>::AdaptiveSubstepperExplicit( double         initialStepSize,
                                                                       double         minimumStepSize,
                                                                       double         maxScaleUpFactor,
                                                                       double         scaleDownFactor,
                                                                       double         integrationErrorTolerance,
                                                                       int            nPassesToIncrease,
                                                                       const Matrix6d& Cel )
        : initialStepSize( initialStepSize ),
          minimumStepSize( minimumStepSize ),
          maxScaleUpFactor( maxScaleUpFactor ),
          scaleDownFactor( scaleDownFactor ),
          integrationErrorTolerance( integrationErrorTolerance ),
          nPassesToIncrease( nPassesToIncrease ),
          ignoreErrorToleranceOnMinimumStepSize( true ),
          Cel( Cel ),
          currentProgress( 0.0 ),
          currentSubstepSize( initialStepSize ),
          passedSubsteps( 0 ),
          substepIndex( -1 ),
          discardedDueToError( 0 )
    {
        consistentTangentProgress         = MatrixStateStrain::Zero();
        consistentTangentProgressFullTemp = MatrixStateStrain::Zero();
        consistentTangentProgressHalfTemp = MatrixStateStrain::Zero();

        stressProgressHalfTemp.setZero();
        stressProgressFullTemp.setZero();
        stateProgressHalfTemp.setZero();
        stateProgressFullTemp.setZero();
        consistentTangentProgressHalfTemp.setZero();
        consistentTangentProgressFullTemp.setZero();

        currentState = FullStep;
    }

    template <size_t n, size_t nState>
    const typename AdaptiveSubstepperExplicit<n, nState>::MatrixStateStrain& AdaptiveSubstepperExplicit<n,
                                                                                                        nState>::IX6()
    {
        static MatrixStateStrain IX6 = [] {
            MatrixStateStrain tmp     = MatrixStateStrain::Zero();
            tmp.topLeftCorner( 6, 6 ) = Matrix6d::Identity();
            return tmp;
        }();
        return IX6;
    }

    template <size_t n, size_t nState>
    const typename AdaptiveSubstepperExplicit<n, nState>::TangentSizedMatrix& AdaptiveSubstepperExplicit<n,
                                                                                                         nState>::IXX()
    {
        static TangentSizedMatrix IXX = [] {
            TangentSizedMatrix tmp = TangentSizedMatrix::Identity();
            return tmp;
        }();
        return IXX;
    }

    template <size_t n, size_t nState>
    void AdaptiveSubstepperExplicit<n, nState>::setConvergedProgress( const Marmot::Vector6d&           stressOld,
                                                                      const IntegrationStateVector& stateVarsOld )
    {
        stressProgress = stressOld;
        stateProgress  = stateVarsOld;
    }

    template <size_t n, size_t nState>
    bool AdaptiveSubstepperExplicit<n, nState>::isFinished()
    {
        return ( ( 1.0 - currentProgress ) <= 2e-16 && currentState == FullStep );
    }

    template <size_t n, size_t nState>
    double AdaptiveSubstepperExplicit<n, nState>::getNextSubstep()
    {
        switch ( currentState ) {
        case FullStep: {
            const double remainingProgress = 1.0 - currentProgress;
            if ( remainingProgress < currentSubstepSize )
                currentSubstepSize = remainingProgress;
            substepIndex++;
            return currentSubstepSize;
            break;
        }
        case FirstHalfStep: {
            return 0.5 * currentSubstepSize;
            break;
        }
        case SecondHalfStep: {
            return 0.5 * currentSubstepSize;
            break;
        }
        }
        return 0;
    }
    template <size_t n, size_t nState>
    void AdaptiveSubstepperExplicit<n, nState>::getConvergedProgress( Marmot::Vector6d&           stress,
                                                                      IntegrationStateVector& stateVars )
    {
        switch ( currentState ) {
        case FullStep: {
            stress    = stressProgress;
            stateVars = stateProgress;
            break;
        }
        case FirstHalfStep: {
            stress    = stressProgress;
            stateVars = stateProgress;
            break;
        }
        case SecondHalfStep: {
            stress    = stressProgressHalfTemp;
            stateVars = stateProgressHalfTemp;
            break;
        }
        }
    }

    template <size_t n, size_t nState>
    bool AdaptiveSubstepperExplicit<n, nState>::discardSubstep()
    {
        passedSubsteps = 0;
        switch ( currentState ) {
        case FullStep: {
            currentSubstepSize *= scaleDownFactor; // we use the scale factor only here
            break;
        }
        // these cases should actually never happen, as the full step has already converged!
        case FirstHalfStep:
            MarmotJournal::warningToMSG( "UMAT: warning, 1th half sub step has not converged after already "
                          "converged full step" );
            return acceptSubstepWithFullStepOnly();

        case SecondHalfStep:
            MarmotJournal::warningToMSG( "UMAT: warning, 2th half sub step has not converged after already "
                          "converged full step" );
            return acceptSubstepWithFullStepOnly();
        }

        currentState = FullStep;

        if ( currentSubstepSize < minimumStepSize )
            return MarmotJournal::warningToMSG( "UMAT: Substepper: Minimal stepzsize reached" );
        else
            return true;
    }

    template <size_t n, size_t nState>
    bool AdaptiveSubstepperExplicit<n, nState>::repeatSubstep( double factorNew )
    {
        currentState   = FullStep;
        passedSubsteps = 0;

        currentSubstepSize *= factorNew; // we use the scale factor only here

        if ( currentSubstepSize < minimumStepSize )
            return MarmotJournal::warningToMSG( "UMAT: Substepper: Minimal stepzsize reached" );
        else
            return true;
    }

    template <size_t n, size_t nState>
    bool AdaptiveSubstepperExplicit<n, nState>::finishSubstep( const Marmot::Vector6d&           resultStress,
                                                               const TangentSizedMatrix&     dXdY,
                                                               const TangentSizedMatrix&     dYdXOld,
                                                               const IntegrationStateVector& stateVars )
    {
        if ( currentState == FullStep ) {
            stressProgressFullTemp            = resultStress;
            stateProgressFullTemp             = stateVars;
            consistentTangentProgressFullTemp = consistentTangentProgress;
            consistentTangentProgressFullTemp.applyOnTheLeft( IXX() - dYdXOld );
            consistentTangentProgressFullTemp += currentSubstepSize * IX6() * Cel;
            consistentTangentProgressFullTemp.applyOnTheLeft( dXdY );
            currentState = FirstHalfStep;
            return true;
        }
        else if ( currentState == FirstHalfStep ) {

            stressProgressHalfTemp            = resultStress;
            stateProgressHalfTemp             = stateVars;
            consistentTangentProgressHalfTemp = consistentTangentProgress;
            consistentTangentProgressHalfTemp.applyOnTheLeft( IXX() - dYdXOld );
            consistentTangentProgressHalfTemp += 0.5 * currentSubstepSize * IX6() * Cel;
            consistentTangentProgressHalfTemp.applyOnTheLeft( dXdY );
            currentState = SecondHalfStep;
            return true;
        }

        else if ( currentState == SecondHalfStep ) {
            // error Estimation

            currentState             = FullStep;
            const double error       = ( resultStress - stressProgressFullTemp ).norm();
            const double errorRatio  = error / integrationErrorTolerance;
            double       scaleFactor = 1.0;
            if ( errorRatio > 1e-10 )
                scaleFactor = 0.9 * std::sqrt( 1. / errorRatio );

            // saturations
            if ( scaleFactor < 0.1 )
                scaleFactor = 0.1;
            if ( scaleFactor * currentSubstepSize < minimumStepSize )
                scaleFactor = minimumStepSize / currentSubstepSize;
            if ( scaleFactor > maxScaleUpFactor )
                scaleFactor = maxScaleUpFactor;
            if ( scaleFactor > 10 )
                scaleFactor = 10;

            // Error large than tolerance?
            if ( error > integrationErrorTolerance ) {
                discardedDueToError++;
                passedSubsteps = 0;
                if ( errorRatio < 2 ) {
                    return splitCurrentSubstep();
                }
                else {
                    return repeatSubstep( scaleFactor );
                }
                // return splitCurrentSubstep(); }
            }
            else {

                stressProgressHalfTemp = resultStress;
                stateProgressHalfTemp  = stateVars;
                consistentTangentProgressHalfTemp.applyOnTheLeft( IXX() - dYdXOld );
                consistentTangentProgressHalfTemp += 0.5 * currentSubstepSize * IX6() * Cel;
                consistentTangentProgressHalfTemp.applyOnTheLeft( dXdY );

                consistentTangentProgress = 2 * consistentTangentProgressHalfTemp - consistentTangentProgressFullTemp;
                stressProgress            = 2 * stressProgressHalfTemp - stressProgressFullTemp;
                stateProgress             = 2 * stateProgressHalfTemp - stateProgressFullTemp;
                currentProgress += currentSubstepSize;

                passedSubsteps++;
                currentSubstepSize *= scaleFactor;

                return true;
            }
        }
        return false;
    }

    template <size_t n, size_t nState>
    void AdaptiveSubstepperExplicit<n, nState>::finishElasticSubstep( const Marmot::Vector6d& newStress )
    {
        switch ( currentState ) {
        case FullStep: {
            // this means, that the complete current cycle is already successfull,
            // as the two half steps must also be elastic!
            consistentTangentProgress += currentSubstepSize * IX6() * Cel;
            stressProgress = newStress;
            // no need for two half steps if full step was already elastic
            currentProgress += currentSubstepSize;
            currentState = FullStep;
            passedSubsteps++;
            break;
        }

        case FirstHalfStep: {
            consistentTangentProgressHalfTemp = consistentTangentProgress;
            consistentTangentProgressHalfTemp += 0.5 * currentSubstepSize * IX6() * Cel;
            stressProgressHalfTemp = newStress;
            currentState           = SecondHalfStep;
            break;
        }
        case SecondHalfStep: {
            acceptSubstepWithFullStepOnly();
            break;
        }
        }
    }

    template <size_t n, size_t nState>
    void AdaptiveSubstepperExplicit<n, nState>::getResults( Marmot::Vector6d&           stress,
                                                            Matrix6d&                consistentTangentOperator,
                                                            IntegrationStateVector& stateVars )
    {
        stress                    = stressProgress;
        stateVars                 = stateProgress;
        consistentTangentOperator = consistentTangentProgress.topLeftCorner( 6, 6 );
    }

    template <size_t n, size_t nState>
    bool AdaptiveSubstepperExplicit<n, nState>::acceptSubstepWithFullStepOnly()
    {
        consistentTangentProgress = consistentTangentProgressFullTemp;
        stressProgress            = stressProgressFullTemp;
        stateProgress             = stateProgressFullTemp;

        currentProgress += currentSubstepSize;
        currentState = FullStep;

        return true;
    }

    template <size_t n, size_t nState>
    bool AdaptiveSubstepperExplicit<n, nState>::splitCurrentSubstep()
    {
        if ( currentSubstepSize < 2 * minimumStepSize ) {
            if ( ignoreErrorToleranceOnMinimumStepSize )
                return acceptSubstepWithFullStepOnly();
            else
                return false;
        }

        consistentTangentProgressFullTemp = consistentTangentProgressHalfTemp;
        stressProgressFullTemp            = stressProgressHalfTemp;
        stateProgressFullTemp             = stateProgressHalfTemp;
        currentSubstepSize *= 0.5;
        currentState = FirstHalfStep;

        return true;
    }
    template <size_t n, size_t nState>
    int AdaptiveSubstepperExplicit<n, nState>::getNumberOfSubsteps()
    {
        return substepIndex;
    }
    template <size_t n, size_t nState>
    int AdaptiveSubstepperExplicit<n, nState>::getNumberDiscardedSubstepsDueToError()
    {
        return discardedDueToError;
    }
} // namespace Marmot::NumericalAlgorithms
