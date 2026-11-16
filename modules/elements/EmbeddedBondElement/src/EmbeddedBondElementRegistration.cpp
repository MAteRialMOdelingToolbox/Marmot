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
#include "Marmot/EmbeddedLargeSlipBondElement.h"
#include "Marmot/MarmotElementFactory.h"
#include <string>
#include <utility>

namespace Marmot::Discretizations::Elements::Registration {

  using namespace Marmot::Registry;

  template < int nDim, int nBarNodes, int nHostNodes >
  bool registerBond( const std::string& name )
  {
    return Marmot::Registry::ElementFactory::
      registerElement( name, []( int elementID ) -> Marmot::Discretizations::Elements::Element* {
        return new EmbeddedBondElement< nDim, nBarNodes, nHostNodes >( elementID );
      } );
  }

  // EB<nDim>D<number of bar nodes><host shape>
  const static bool EB2D2T3_isRegistered  = registerBond< 2, 2, 3 >( "EB2D2T3" );
  const static bool EB2D3T3_isRegistered  = registerBond< 2, 3, 3 >( "EB2D3T3" );
  const static bool EB2D2Q4_isRegistered  = registerBond< 2, 2, 4 >( "EB2D2Q4" );
  const static bool EB2D2Q8_isRegistered  = registerBond< 2, 2, 8 >( "EB2D2Q8" );
  const static bool EB2D3Q4_isRegistered  = registerBond< 2, 3, 4 >( "EB2D3Q4" );
  const static bool EB2D3Q8_isRegistered  = registerBond< 2, 3, 8 >( "EB2D3Q8" );
  const static bool EB3D2T4_isRegistered  = registerBond< 3, 2, 4 >( "EB3D2T4" );
  const static bool EB3D2H8_isRegistered  = registerBond< 3, 2, 8 >( "EB3D2H8" );
  const static bool EB3D2H20_isRegistered = registerBond< 3, 2, 20 >( "EB3D2H20" );
  const static bool EB3D3T4_isRegistered  = registerBond< 3, 3, 4 >( "EB3D3T4" );
  const static bool EB3D3H8_isRegistered  = registerBond< 3, 3, 8 >( "EB3D3H8" );
  const static bool EB3D3H20_isRegistered = registerBond< 3, 3, 20 >( "EB3D3H20" );

  // large slip: EBLS<nDim>D<number of bar nodes><host shape>W<number of bar elements in the window>.
  // The window is a runtime size of the element (one class per bar and host shape); one name per window size is
  // registered only because the element factory takes a plain function pointer of the element label, which cannot
  // carry the size. The largest window is the one documented cap: a larger one is reported by the factory as an
  // unknown element.
  constexpr int maxWindow = 48;

  template < int nDim, int nBarNodes, int nHostNodes, int nBarElements >
  Marmot::Discretizations::Elements::Element* makeLargeSlipBond( int elementID )
  {
    return new EmbeddedLargeSlipBondElement< nDim, nBarNodes, nHostNodes >( elementID, nBarElements );
  }

  template < int nDim, int nBarNodes, int nHostNodes, int... windows >
  bool registerLargeSlipBonds( const std::string& prefix, std::integer_sequence< int, windows... > )
  {
    return ( Marmot::Registry::ElementFactory::
               registerElement( prefix + "W" + std::to_string( windows + 1 ),
                                &makeLargeSlipBond< nDim, nBarNodes, nHostNodes, windows + 1 > ) &&
             ... );
  }

  template < int nDim, int nBarNodes, int nHostNodes >
  bool registerLargeSlipBonds( const std::string& prefix )
  {
    return registerLargeSlipBonds< nDim, nBarNodes, nHostNodes >( prefix,
                                                                  std::make_integer_sequence< int, maxWindow >{} );
  }

  const static bool EBLS2D2T3_isRegistered  = registerLargeSlipBonds< 2, 2, 3 >( "EBLS2D2T3" );
  const static bool EBLS2D3T3_isRegistered  = registerLargeSlipBonds< 2, 3, 3 >( "EBLS2D3T3" );
  const static bool EBLS2D2Q4_isRegistered  = registerLargeSlipBonds< 2, 2, 4 >( "EBLS2D2Q4" );
  const static bool EBLS2D2Q8_isRegistered  = registerLargeSlipBonds< 2, 2, 8 >( "EBLS2D2Q8" );
  const static bool EBLS2D3Q4_isRegistered  = registerLargeSlipBonds< 2, 3, 4 >( "EBLS2D3Q4" );
  const static bool EBLS2D3Q8_isRegistered  = registerLargeSlipBonds< 2, 3, 8 >( "EBLS2D3Q8" );
  const static bool EBLS3D2T4_isRegistered  = registerLargeSlipBonds< 3, 2, 4 >( "EBLS3D2T4" );
  const static bool EBLS3D2H8_isRegistered  = registerLargeSlipBonds< 3, 2, 8 >( "EBLS3D2H8" );
  const static bool EBLS3D2H20_isRegistered = registerLargeSlipBonds< 3, 2, 20 >( "EBLS3D2H20" );
  const static bool EBLS3D3T4_isRegistered  = registerLargeSlipBonds< 3, 3, 4 >( "EBLS3D3T4" );
  const static bool EBLS3D3H8_isRegistered  = registerLargeSlipBonds< 3, 3, 8 >( "EBLS3D3H8" );
  const static bool EBLS3D3H20_isRegistered = registerLargeSlipBonds< 3, 3, 20 >( "EBLS3D3H20" );

} // namespace Marmot::Discretizations::Elements::Registration
