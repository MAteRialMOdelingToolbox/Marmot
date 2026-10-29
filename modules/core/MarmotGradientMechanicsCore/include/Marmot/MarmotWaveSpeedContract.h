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
#include "Marmot/MarmotMaterialGeneralGradientEnhancedHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <vector>

namespace Marmot::Testing {

  /**
   * @brief Check that an override of getMaximumWaveSpeed() keeps the contract of the base class' default.
   * @details A material overrides getMaximumWaveSpeed() with a closed form to avoid the default's full constitutive
   * evaluation. The artificial bulk viscosity's damage degradation divides the wave speed at the current non-local
   * field by the one at a zero field, so an override that ignores the field returns the same speed for both: the
   * ratio is 1 and the degradation silently does nothing. This function compares the override to the default
   * (called non-virtually on the same object) at every given field, and checks that neither touches the state.
   * @tparam nNonlocalVariables Number of non-local variables of the material.
   * @param[in] material The material whose override is checked.
   * @param[in] currentResponse A response, including its state variables, at which to evaluate both.
   * @param[in] fields The non-local fields to evaluate at; they should include fields at which the material is
   * damaged, at least one of them non-zero.
   * @param[in] relativeTolerance 0, the default, demands the bits of the default.
   * @throws std::runtime_error if the override differs from the default or changes the state variables.
   */
  template < int nNonlocalVariables >
  void checkMaximumWaveSpeedOverrideMatchesDefault(
    const MarmotMaterialGeneralGradientEnhancedHypoElastic< nNonlocalVariables >&                    material,
    const typename MarmotMaterialGeneralGradientEnhancedHypoElastic< nNonlocalVariables >::response& currentResponse,
    const std::vector< Eigen::Vector< double, nNonlocalVariables > >&                                fields,
    double relativeTolerance = 0.0 )
  {
    const int                   nStateVars = material.getNumberOfRequiredStateVars();
    const std::vector< double > stateVarsBefore( currentResponse.stateVars, currentResponse.stateVars + nStateVars );

    for ( const auto& K : fields ) {
      const double override_ = material.getMaximumWaveSpeed( currentResponse, K );
      const double default_  = material.MarmotMaterialGeneralGradientEnhancedHypoElastic<
        nNonlocalVariables >::getMaximumWaveSpeed( currentResponse, K );

      throwExceptionOnFailure( std::abs( override_ - default_ ) <= relativeTolerance * std::abs( default_ ),
                               MakeString() << __PRETTY_FUNCTION__ << ": getMaximumWaveSpeed() returns " << override_
                                            << " at the non-local field " << K.transpose() << ", the default returns "
                                            << default_ );
    }

    throwExceptionOnFailure( std::equal( stateVarsBefore.begin(), stateVarsBefore.end(), currentResponse.stateVars ),
                             MakeString()
                               << __PRETTY_FUNCTION__ << ": getMaximumWaveSpeed() changed the state variables" );
  }

} // namespace Marmot::Testing
