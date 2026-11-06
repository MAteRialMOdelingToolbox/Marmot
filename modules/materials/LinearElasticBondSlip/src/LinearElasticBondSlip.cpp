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
#include "Marmot/LinearElasticBondSlip.h"
#include "Marmot/MarmotJournal.h"
#include <stdexcept>

namespace Marmot::Materials {

  LinearElasticBondSlip::LinearElasticBondSlip( const double* materialProperties,
                                                int           nMaterialProperties,
                                                int           materialNumber )
    : Marmot::Materials::BondSlipLaw( materialProperties, nMaterialProperties, materialNumber ),
      Kt( materialProperties[0] ),
      Kn( materialProperties[1] )
  {
    if ( nMaterialProperties < 2 )
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__ << ": requires 2 properties, K_t and K_n" );
    stateLayout.finalize();
  }

  void LinearElasticBondSlip::computeBondStress( State&                        state,
                                                 Marmot::Math::Matrix3d&       dBondStress_dSlip,
                                                 const Marmot::Math::Vector3d& slip,
                                                 const Marmot::Math::Vector3d&,
                                                 const TimeIncrement& ) const
  {
    dBondStress_dSlip = Marmot::Math::Vector3d( Kt, Kn, Kn ).asDiagonal();
    state.bondStress  = dBondStress_dSlip * slip;

    state.elasticEnergyDensity = 0.5 * slip.dot( state.bondStress );
  }
} // namespace Marmot::Materials
