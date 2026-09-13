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

#include "Marmot/AT2PhaseField.h"
#include "Marmot/MarmotMaterialGeneralGradientEnhancedHypoElasticFactory.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool AT2PhaseFieldIsRegistered = MarmotMaterialGeneralGradientEnhancedHypoElasticFactory<
      1 >::registerMaterial< AT2PhaseField >( "AT2PHASEFIELD" );

    // NOT registered with GradientEnhancedHughesWingetWrapper, deliberately. That wrapper hands the
    // wrapped model the forward-rotated stress of the last increment, and this model ignores its
    // `res.stress` argument entirely -- it returns g(phi) * C : eps from its own stored `strain` state.
    // Wrapped, its Kirchhoff stress does not rotate at all under a rigid rotation. See the note in
    // MarmotMaterialGradientEnhancedHughesWinget.h.

  } // namespace Registration

} // namespace Marmot::Materials
