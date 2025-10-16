#include "Marmot/NewtonConvergenceChecker.h"

using namespace Eigen;

namespace Marmot::NumericalAlgorithms {

  NewtonConvergenceChecker::NewtonConvergenceChecker( const VectorXd& residualScaleVector,
                                                      int             nMaxNewtonCycles,
                                                      int             nMaxNewtonCyclesAlt,
                                                      double          newtonTol,
                                                      double          newtonRTol,
                                                      double          newtonTolAlt,
                                                      double          newtonRTolAlt )
    : residualScaleVector( residualScaleVector ),
      nMaxNewtonCycles( nMaxNewtonCycles ),
      nMaxNewtonCyclesAlt( nMaxNewtonCyclesAlt ),
      newtonTol( newtonTol ),
      newtonRTol( newtonRTol ),
      newtonTolAlt( newtonTolAlt ),
      newtonRTolAlt( newtonRTolAlt )
  {
  }

  double NewtonConvergenceChecker::relativeNorm( const VectorXd& increment, const VectorXd& reference )
  {
    double incNorm = increment.norm();
    double refNorm = reference.norm();

    if ( incNorm < 1e-14 )
      return incNorm;

    if ( refNorm < 1e-12 )
      // for a too small reference norm, a reasonable relative norm cannot be computed
      return 0.0;

    return increment.norm() / refNorm;
  }

  double NewtonConvergenceChecker::residualNorm( const VectorXd& residual )
  {
    return ( residualScaleVector.array() * residual.array() ).matrix().norm();
  }

  bool NewtonConvergenceChecker::iterationFinished( const VectorXd& residual,
                                                    const VectorXd& X,
                                                    const VectorXd& dX,
                                                    int             numberOfIterations )
  {
    if ( isConverged( residual, X, dX, numberOfIterations ) || numberOfIterations > nMaxNewtonCyclesAlt )
      return true;
    else
      return false;
  }

  bool NewtonConvergenceChecker::isConverged( const VectorXd& residual,
                                              const VectorXd& X,
                                              const VectorXd& dX,
                                              int             numberOfIterations )
  {
    const double resNorm = residualNorm( residual );
    const double relNorm = relativeNorm( dX, X );

    if ( numberOfIterations <= nMaxNewtonCycles ) {
      if ( resNorm <= newtonTol && relNorm <= newtonRTol )
        return true;
    }
    else if ( numberOfIterations <= nMaxNewtonCyclesAlt + 1 ) {
      if ( resNorm <= newtonTolAlt && relNorm <= newtonRTolAlt )
        return true;
    }
    return false;
  }
} // namespace Marmot::NumericalAlgorithms
