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
#include "Marmot/MarmotTypedefs.h"

/**
 * @file PerezFougetSubstepperTime.h
 * @brief Time-variant variant of the Pérez–Fouget sub-stepper.
 *
 * A modified version of the Pérez–Fouget sub-stepper that accounts for a
 * time-variant elastic stiffness tensor \f$ \mathbb{C}^\mathrm{el}(t_{n+1}) \f$.
 * The algorithmic formulation is unchanged; the only difference is that the
 * current elastic tangent must be supplied when extending the consistent tangent.
 */

namespace Marmot::NumericalAlgorithms {
  /**
   * @brief Pérez–Fouget sub-stepper with time-dependent elastic stiffness.
   *
   * Compared to @ref PerezFougetSubstepper, this variant does not require the
   * elastic stiffness @c Cel at construction time; instead it is passed
   * explicitly to @c extendConsistentTangent at each sub-step.
   *
   * @tparam sizeMaterialState  Total size of the material state vector.
   */
  template < int sizeMaterialState >
  class PerezFougetSubstepperTime {

  public:
    /// Matrix to carry the Jacobian of a material state
    typedef Eigen::Matrix< double, sizeMaterialState, sizeMaterialState > TangentSizedMatrix;

    /**
     * @brief Construct a PerezFougetSubstepperTime.
     * @param initialStepSize   Initial sub-step size as a fraction of the total increment (0, 1].
     * @param minimumStepSize   Minimum allowable sub-step size; triggers a warning if reached.
     * @param scaleUpFactor     Factor by which the sub-step size is increased after @p nPassesToIncrease successes.
     * @param scaleDownFactor   Factor applied to the sub-step size when convergence fails.
     * @param nPassesToIncrease Number of successful sub-steps before the step size may be increased.
     */
    PerezFougetSubstepperTime( double initialStepSize,
                               double minimumStepSize,
                               double scaleUpFactor,
                               double scaleDownFactor,
                               int    nPassesToIncrease );
    /// Check if the subincrementation has finished
    bool isFinished();
    /// get the next subincrement size
    double getNextSubstep();
    /// get the total finished progress of the subincrementation process
    double getFinishedProgress();
    /// decrease the next subincrement
    bool decreaseSubstepSize();

    /**
     * @brief Accumulate the consistent tangent for an elastic sub-step.
     * @param CelT  Elastic stiffness matrix (6×6) at the current (end-of-sub-step) time.
     */
    void extendConsistentTangent( const Matrix6d& CelT );

    /**
     * @brief Accumulate the consistent tangent for an elastoplastic sub-step.
     * @param CelT        Elastic stiffness matrix (6×6) at the current time.
     * @param matTangent  Inverse algorithmic tangent (`TangentSizedMatrix`); left-applied to the current accumulation.
     */
    void extendConsistentTangent( const Matrix6d& CelT, const TangentSizedMatrix& matTangent );

    /**
     * @brief Return the assembled consistent algorithmic stiffness after all sub-steps are finished.
     * @return  The 6×6 consistent tangent operator.
     */
    Matrix6d consistentStiffness();

  private:
    const double initialStepSize, minimumStepSize, scaleUpFactor, scaleDownFactor;
    const int    nPassesToIncrease;

    double currentProgress;
    double currentSubstepSize;
    int    passedSubsteps;

    TangentSizedMatrix elasticTangent;
    TangentSizedMatrix consistentTangent;
  };
} // namespace Marmot::NumericalAlgorithms

#include "Marmot/MarmotJournal.h"

namespace Marmot::NumericalAlgorithms {
  template < int s >
  PerezFougetSubstepperTime< s >::PerezFougetSubstepperTime( double initialStepSize,
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

  template < int s >
  bool PerezFougetSubstepperTime< s >::isFinished()
  {
    return currentProgress >= 1.0;
  }

  template < int s >
  double PerezFougetSubstepperTime< s >::getNextSubstep()
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

  template < int s >
  double PerezFougetSubstepperTime< s >::getFinishedProgress()
  {
    return currentProgress - currentSubstepSize;
  }

  template < int s >
  bool PerezFougetSubstepperTime< s >::decreaseSubstepSize()
  {
    currentProgress -= currentSubstepSize;
    passedSubsteps = 0;

    currentSubstepSize *= scaleDownFactor;

    if ( currentSubstepSize < minimumStepSize )
      return MarmotJournal::warningToMSG( "UMAT: Substepper: Minimal stepzsize reached" );
    else
      return MarmotJournal::notificationToMSG( "UMAT: Substepper: Decreasing stepsize" );
  }

  template < int s >
  void PerezFougetSubstepperTime< s >::extendConsistentTangent( const Matrix6d& CelT )
  {

    elasticTangent.topLeftCorner( 6, 6 ) = CelT;
    consistentTangent += currentSubstepSize * elasticTangent;
  }

  template < int s >
  void PerezFougetSubstepperTime< s >::extendConsistentTangent( const Matrix6d&           CelT,
                                                                const TangentSizedMatrix& matTangent )
  {
    extendConsistentTangent( CelT );
    consistentTangent.applyOnTheLeft( matTangent );
  }

  template < int s >
  Matrix6d PerezFougetSubstepperTime< s >::consistentStiffness()
  {
    return consistentTangent.topLeftCorner( 6, 6 );
  }
} // namespace Marmot::NumericalAlgorithms
