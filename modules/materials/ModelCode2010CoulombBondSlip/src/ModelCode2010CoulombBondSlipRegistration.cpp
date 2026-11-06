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
#include "Marmot/MarmotBondSlipLawFactory.h"
#include "Marmot/ModelCode2010CoulombBondSlip.h"

namespace Marmot::Materials::Registration {
  using namespace Marmot::Registry;

  const static bool ModelCode2010CoulombBondSlipIsRegistered = Marmot::Registry::BondSlipLawFactory::
    registerBondSlipLaw< ModelCode2010CoulombBondSlip >( "MODELCODE2010COULOMBBONDSLIP" );
} // namespace Marmot::Materials::Registration
