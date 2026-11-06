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
#include "Marmot/EmbeddedBondElement.h"
#include "Marmot/MarmotElementFactory.h"
#include <string>

namespace Marmot::Discretizations::Elements::Registration {

  using namespace Marmot::Registry;

  template < int nDim, int nRebarNodes, int nHostNodes >
  bool registerBond( const std::string& name )
  {
    return Marmot::Registry::ElementFactory::
      registerElement( name, []( int elementID ) -> Marmot::Discretizations::Elements::Element* {
        return new EmbeddedBondElement< nDim, nRebarNodes, nHostNodes >( elementID );
      } );
  }

  // EB<nDim>D<number of rebar nodes><host shape>
  const static bool EB2D2Q4_isRegistered  = registerBond< 2, 2, 4 >( "EB2D2Q4" );
  const static bool EB2D2Q8_isRegistered  = registerBond< 2, 2, 8 >( "EB2D2Q8" );
  const static bool EB2D3Q4_isRegistered  = registerBond< 2, 3, 4 >( "EB2D3Q4" );
  const static bool EB2D3Q8_isRegistered  = registerBond< 2, 3, 8 >( "EB2D3Q8" );
  const static bool EB3D2T4_isRegistered  = registerBond< 3, 2, 4 >( "EB3D2T4" );
  const static bool EB3D2T10_isRegistered = registerBond< 3, 2, 10 >( "EB3D2T10" );
  const static bool EB3D2H8_isRegistered  = registerBond< 3, 2, 8 >( "EB3D2H8" );
  const static bool EB3D2H20_isRegistered = registerBond< 3, 2, 20 >( "EB3D2H20" );
  const static bool EB3D3T4_isRegistered  = registerBond< 3, 3, 4 >( "EB3D3T4" );
  const static bool EB3D3T10_isRegistered = registerBond< 3, 3, 10 >( "EB3D3T10" );
  const static bool EB3D3H8_isRegistered  = registerBond< 3, 3, 8 >( "EB3D3H8" );
  const static bool EB3D3H20_isRegistered = registerBond< 3, 3, 20 >( "EB3D3H20" );

} // namespace Marmot::Discretizations::Elements::Registration
