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
#include "Marmot/YieldSurfaceCombinationManager.h"
#include <cmath>
#include <functional>
#include <string>

/**
 * @file RankineSoftening.h
 * @brief Building blocks shared by the multi-surface Rankine softening plasticity models.
 */

namespace Marmot::Materials::RankineSoftening {

  using YieldSurfaceCombinationManager = Marmot::ContinuumMechanics::Algorithms::YieldSurfaceCombinationManager< 3 >;
  using YieldSurfFlagArr               = YieldSurfaceCombinationManager::YieldSurfFlagArr;

  /**
   * @brief Isotropic hardening/softening law of VonMisesModel,
   * \f$ f_\mathrm{y}(\kappa) = f_\mathrm{y}^0 + H_\mathrm{iso,lin}\kappa
   *     + \Delta f_\mathrm{y}^{0,\infty}\left(1-\exp(-\delta\kappa)\right) \f$.
   */
  struct Hardening {
    double fy0;
    double HLin;
    double deltaYieldStress;
    double delta;

    template < typename T >
    T fy( const T& kappa ) const
    {
      using std::exp;
      return fy0 + HLin * kappa + deltaYieldStress * ( 1. - exp( -delta * kappa ) );
    }

    template < typename T >
    T dfy_dKappa( const T& kappa ) const
    {
      using std::exp;
      return HLin + deltaYieldStress * delta * exp( -delta * kappa );
    }
  };

  /**
   * @brief Active-set strategy for the three Rankine surfaces: starting from @p activeSurfaces, all combinations are
   * tried until @p returnMappingAttempt reports an admissible solution.
   * @throws Marmot::Utilities::StressUpdateFailed if no combination succeeds.
   */
  void searchActiveSurfaces( YieldSurfFlagArr                                        activeSurfaces,
                             const std::function< bool( const YieldSurfFlagArr& ) >& returnMappingAttempt,
                             const std::string&                                      caller );

  /**
   * @brief Orthonormal tensor basis of the principal frame, as flattened (column-major) 3x3 tensors.
   */
  struct SpectralBasis {
    /// eigenprojections \f$ \mathbf{E}_A = \boldsymbol{n}_A\otimes\boldsymbol{n}_A \f$
    Eigen::Matrix< double, 9, 3 > E;
    /// \f$ \mathbf{N}_{AB} = \sqrt{2}\,\mathrm{sym}(\boldsymbol{n}_A\otimes\boldsymbol{n}_B) \f$ for (A,B) = (1,2),
    /// (1,3), (2,3)
    Eigen::Matrix< double, 9, 3 > N;
  };

  SpectralBasis spectralBasis( const Eigen::Matrix3d& eigenVectors );

  /**
   * @brief Coefficients \f$ c_{AB} = (\sigma_A-\sigma_B)/(\sigma^\mathrm{trial}_A-\sigma^\mathrm{trial}_B) \f$ of the
   * principal-frame rotation for (A,B) = (1,2), (1,3), (2,3), with the limit for repeated trial eigenvalues.
   *
   * Together with spectralBasis, the derivative of the coaxial update \f$ \boldsymbol{\sigma} = \sum_A \sigma_A
   * \mathbf{E}_A \f$ wrt the symmetric trial stress reads \f$ \sum_{A,B} D_{AB}\,\mathbf{E}_A\otimes\mathbf{E}_B +
   * \sum_{A<B} c_{AB}\,\mathbf{N}_{AB}\otimes\mathbf{N}_{AB} \f$.
   *
   * @param sigTrial principal trial stresses
   * @param sig principal stresses after the return mapping
   * @param D \f$ \partial\sigma_A / \partial\sigma^\mathrm{trial}_B \f$
   * @param tol absolute tolerance for treating trial eigenvalues as repeated
   */
  Eigen::Vector3d spinCoefficients( const Eigen::Vector3d& sigTrial,
                                    const Eigen::Vector3d& sig,
                                    const Eigen::Matrix3d& D,
                                    double                 tol );

} // namespace Marmot::Materials::RankineSoftening
