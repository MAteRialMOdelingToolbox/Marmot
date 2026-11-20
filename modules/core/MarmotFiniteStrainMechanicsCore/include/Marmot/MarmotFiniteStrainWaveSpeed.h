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
#include <Fastor/tensor/Tensor.h>
#include <algorithm>
#include <cmath>

namespace Marmot::ContinuumMechanics::FiniteStrain {

  /**
   * @brief Largest longitudinal wave speed of a finite-strain material, from a central-difference estimate of the
   * diagonal of its stiffness.
   *
   * @tparam KirchhoffStressAt Callable `Fastor::Tensor<double, 3, 3>( const Fastor::Tensor<double, 3, 3>& F )`
   *                           returning the Kirchhoff stress of the material at @p F, evaluated on a fresh copy of
   *                           the state variables, so that the estimate leaves the state untouched.
   * @param[in] F                 Deformation gradient the stiffness is taken at.
   * @param[in] density           Mass density in the reference configuration.
   * @param[in] kirchhoffStressAt The stress evaluation.
   * @return @f$ \sqrt{\max_i C_{ii} / \rho_0} @f$, with @f$ C_{ii} @f$ the change of @f$ \tau_{ii} @f$ per unit
   *         change of the logarithmic strain @f$ \ln F_{ii} @f$ (stretches @f$ 1 \pm 10^{-6} @f$); zero for a
   *         non-positive density.
   *
   * @details The default of both finite-strain material interfaces
   * (Marmot::Materials::MaterialFiniteStrain::getMaximumWaveSpeed and
   * Marmot::Materials::MaterialGradientEnhancedFiniteStrain::getMaximumWaveSpeed), so that the two critical time steps
   * of the explicit solvers -- the one of the finite elements and the one of the particles -- read the material the
   * same way.
   */
  template < typename KirchhoffStressAt >
  double maximumWaveSpeedByCentralDifferences( const Fastor::Tensor< double, 3, 3 >& F,
                                               double                                density,
                                               KirchhoffStressAt&&                   kirchhoffStressAt )
  {
    constexpr double dEps = 1e-6;

    double maxStiffnessDiagonal = 0.0;

    for ( int i = 0; i < 3; ++i ) {
      auto FPlus  = F;
      auto FMinus = F;
      FPlus( i, i ) *= ( 1.0 + dEps );
      FMinus( i, i ) *= ( 1.0 - dEps );

      const Fastor::Tensor< double, 3, 3 > tauPlus  = kirchhoffStressAt( FPlus );
      const Fastor::Tensor< double, 3, 3 > tauMinus = kirchhoffStressAt( FMinus );

      const double dLogStrain = std::log( 1.0 + dEps ) - std::log( 1.0 - dEps );
      const double Cii        = ( tauPlus( i, i ) - tauMinus( i, i ) ) / dLogStrain;
      maxStiffnessDiagonal    = std::max( maxStiffnessDiagonal, Cii );
    }

    return density > 0.0 ? std::sqrt( std::max( 0.0, maxStiffnessDiagonal ) / density ) : 0.0;
  }

} // namespace Marmot::ContinuumMechanics::FiniteStrain
