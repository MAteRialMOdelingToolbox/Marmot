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
#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotTypedefs.h"
#include <algorithm>
#include <tuple>

namespace Marmot::NumericalAlgorithms::TimeIntegration {

  /** @brief Explicit Euler time integrator
   * @param yN Current value
   * @param dt Time step size
   * @param fRate Function that computes the rate of change
   * @param fRateArgs Additional arguments for the rate function
   * @return Updated value after time step
   *
   * This function computes one single time step using the explicit Euler method:
   * \f[ \boldsymbol{y}_{n+1} = \boldsymbol{y}_n + f(\boldsymbol{y}_n) * \Delta t \f]
   * where \f$ \boldsymbol{y}_n \f$ is the current value, \f$ \Delta t \f$ is the time step size,
   * and \f$ f(\boldsymbol{y}_n) \f$ is the rate of change.
   *
   */
  template < typename functionType, typename yType, typename... Args >
  yType explicitEuler( yType yN, const double dt, functionType fRate, Args&&... fRateArgs )
  {
    return yN + fRate( yN, fRateArgs... ) * dt;
  }

  /**
   * @brief Implicit Euler time integrator for a linear vector-valued rate equation
   * @param yN Current value
   * @param dt Time step size
   * @param fRate Function that computes the rate of change
   * @param fRateArgs Additional arguments for the rate function
   * @return Updated value after time step
   *
   * This function computes one single time step using the implicit Euler method for a linear vector-valued rate
   * equation: \f[ \boldsymbol{y}_{n+1} = \boldsymbol{y}_n + f(\boldsymbol{y}_{n+1}) \Delta t \f] where \f$
   * \boldsymbol{y}_n \f$ is the current value, \f$ \Delta t \f$ is the time step size, and \f$ f(\boldsymbol{y}) \f$
   * is the linear rate of change evaluated at the next time step.
   *
   *
   * @todo: Is this really semi-implicit? It looks like a fully implicit Euler step for linear systems.
   * @todo: Use external central difference function?
   *
   * */
  template < int ySize, typename functionType, typename... Args >
  Eigen::Matrix< double, ySize, 1 > semiImplicitEuler( Eigen::Matrix< double, ySize, 1 > yN,
                                                       const double                      dt,
                                                       functionType                      fRate,
                                                       Args&&... fRateArgs )
  {
    Eigen::Matrix< double, ySize, ySize > fS = Eigen::Matrix< double, ySize, ySize >::Zero();
    Eigen::Matrix< double, ySize, ySize > Iy = Eigen::Matrix< double, ySize, ySize >::Identity();

    Eigen::VectorXd leftX( ySize );
    Eigen::VectorXd rightX( ySize );

    for ( size_t i = 0; i < ySize; i++ ) {
      double volatile h = std::max( 1.0, std::abs( yN( i ) ) ) * Constants::cubicRootEps();
      leftX             = yN;
      leftX( i ) -= h;
      rightX = yN;
      rightX( i ) += h;
      fS.col( i ) = 1. / ( 2. * h ) * ( fRate( rightX, fRateArgs... ) - fRate( leftX, fRateArgs... ) );
    }

    const auto A = Iy - dt * fS;
    const auto r = fRate( yN, fRateArgs... );

    if constexpr ( ySize == 1 ) {
      Eigen::Matrix< double, 1, 1 > out;
      out( 0 ) = yN( 0 ) + dt * r( 0 ) / A( 0, 0 );
      return out;
    }

    return yN + A.colPivHouseholderQr().solve( dt * r );
  }

  /**
   * @brief Explicit Euler integration with error estimation based on Richardson extrapolation
   * @tparam functionType Type of the rate function
   * @tparam yType Type of the current value
   * @tparam Args Additional argument types for the rate function
   * @param yN Current value
   * @param dt Time step size
   * @param fRate Function that computes the rate of change
   * @param fRateArgs Additional arguments for the rate function
   * @return Updated value after time step
   *
   * This function computes one single time step using the explicit Euler method with Richardson extrapolation for
   * error estimation: \f[ \boldsymbol{y}_{n+1} = 2 \left( \boldsymbol{y}_n + f\left(\boldsymbol{y}_n
   * + \frac{f(\boldsymbol{y}_n) \Delta t}{2}\right) \frac{\Delta t}{2} \right) - \left( \boldsymbol{y}_n +
   * f(\boldsymbol{y}_n) \Delta t \right) \f] where \f$ \boldsymbol{y}_n \f$ is the current value, \f$ \Delta t \f$ is
   * the time step size, and \f$ f(\boldsymbol{y}) \f$ is the rate of change.
   *
   */
  template < typename functionType, typename yType, typename... Args >
  yType explicitEulerRichardson( yType yN, const double dt, functionType fRate, Args&&... fRateArgs )
  {
    yType fN = fRate( yN, fRateArgs... );
    yType u  = yN + fN * dt;
    yType v  = yN + fN * dt / 2.;
    yType w  = v + fRate( v, fRateArgs... ) * dt / 2.;

    return 2. * w - u;
  }

  /**
   * @brief Explicit Euler integration based on Richardson extrapolation with error estimation and time step
   * estimation
   * @tparam ySize Size of the state vector
   * @tparam functionType Type of the rate function
   * @tparam Args Additional argument types for the rate function
   * @param yN Current value
   * @param dt Current time step size
   * @param TOL Desired tolerance for the error estimation
   * @param fRate Function that computes the rate of change
   * @param fRateArgs Additional arguments for the rate function
   * @return Tuple containing the updated value after time step and the new time step size
   *
   * This function computes one single time step using the explicit Euler method with Richardson extrapolation for
   * error estimation and adaptive time stepping: \f[ \boldsymbol{y}_{n+1} = 2 \left( \boldsymbol{y}_n +
   * f\left(\boldsymbol{y}_n
   * + \frac{f(\boldsymbol{y}_n) \Delta t}{2}\right) \frac{\Delta t}{2} \right) - \left( \boldsymbol{y}_n +
   * f(\boldsymbol{y}_n) \Delta t \right) \f] where \f$ \boldsymbol{y}_n \f$ is the current value, \f$ \Delta t \f$ is
   * the current time step size, and \f$ f(\boldsymbol{y}) \f$ is the rate of change.
   *
   * The function also estimates the error of the time step and adjusts the time step size for the next iteration
   * based on the desired tolerance @p TOL. The new time step size is computed as: \f[ \Delta t_{\text{new}} =
   * \Delta t \cdot \min\left(2, \max\left(0.2, 0.9 \sqrt{\frac{TOL}{EST}}\right)\right) \f] where \f$ EST \f$ is the
   * estimated error.
   *
   * @todo: reuse explicitEulerRichardson function?
   */
  template < int ySize, typename functionType, typename... Args >
  std::tuple< Eigen::Matrix< double, ySize, 1 >, double > explicitEulerRichardsonWithErrorEstimator(
    Eigen::Matrix< double, ySize, 1 > yN,
    const double                      dt,
    const double                      TOL,
    functionType                      fRate,
    Args&&... fRateArgs )
  {

    typedef Eigen::Matrix< double, ySize, 1 > ySized;
    ySized                                    fN   = fRate( yN, fRateArgs... );
    ySized                                    u    = yN + fN * dt;
    ySized                                    v    = yN + fN * dt / 2.;
    ySized                                    w    = v + fRate( v, fRateArgs... ) * dt / 2.;
    ySized                                    yNew = 2. * w - u;

    // error estimator
    const double AERR    = 1.0;
    const double aI      = AERR / TOL;
    const double rI      = 1.0;
    double       scaling = 0;
    ySized       ESTVec  = ySized::Zero();

    for ( int i = 0; i < ySize; i++ ) {
      scaling     = aI + rI * std::max( abs( yNew( i ) ), abs( yN( i ) ) );
      ESTVec( i ) = abs( w( i ) - u( i ) ) / abs( scaling );
    }

    const double EST    = ESTVec.maxCoeff();
    const double tauNew = dt * std::min( 2., std::max( 0.2, 0.9 * std::sqrt( TOL / EST ) ) );

    return std::make_tuple( yNew, tauNew );
  }

} // namespace Marmot::NumericalAlgorithms::TimeIntegration
