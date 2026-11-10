/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck
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

#include "Marmot/CompressibleFiniteStrainFractionalSchapery.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
#include "Marmot/MarmotMaterialFiniteStrainSubstepped.h"

namespace {
  using namespace Marmot::Materials;

  const static bool CompressibleFiniteStrainFractionalSchaperyRegistered = Marmot::Registry::
    MaterialFiniteStrainFactory::registerMaterial< CompressibleFiniteStrainFractionalSchapery >(
      "COMPRESSIBLEFINITESTRAINFRACTIONALSCHAPERY" );

  const static bool CompressibleFiniteStrainFractionalSchaperySubsteppedRegistered = Marmot::Registry::
    MaterialFiniteStrainFactory::registerMaterial<
      Marmot::Materials::MaterialFiniteStrainSubstepped< CompressibleFiniteStrainFractionalSchapery > >(
      "COMPRESSIBLEFINITESTRAINFRACTIONALSCHAPERY_SUBSTEPPED" );
} // namespace
