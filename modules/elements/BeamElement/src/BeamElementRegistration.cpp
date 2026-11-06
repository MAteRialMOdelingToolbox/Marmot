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
#include "Marmot/BeamElement.h"
#include "Marmot/CorotationalBeamElement.h"
#include "Marmot/CorotationalBeamElement3D.h"
#include "Marmot/MarmotElementFactory.h"

namespace Marmot::Discretizations::Elements::Registration {

  template < int nDim, int nNodes >
  bool registerBeam( const std::string& name )
  {
    return Marmot::Registry::ElementFactory::registerElement( name,
                                                              []( int elementID )
                                                                -> Marmot::Discretizations::Elements::Element* {
                                                                return new BeamElement< nDim, nNodes >( elementID );
                                                              } );
  }

  // Euler-Bernoulli beams B<nDim>D<nNodes>: 2 nodes (cubic Hermite), 3 nodes (quintic Hermite; end, end, mid)
  const static bool B2D2_isRegistered = registerBeam< 2, 2 >( "B2D2" );
  const static bool B2D3_isRegistered = registerBeam< 2, 3 >( "B2D3" );
  const static bool B3D2_isRegistered = registerBeam< 3, 2 >( "B3D2" );
  const static bool B3D3_isRegistered = registerBeam< 3, 3 >( "B3D3" );

  // co-rotational (finite rotation) Euler-Bernoulli beam, 2D, 2 nodes
  const static bool B2D2CR_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "B2D2CR", []( int elementID ) -> Marmot::Discretizations::Elements::Element* {
      return new CorotationalBeamElement2D( elementID );
    } );

  // co-rotational (finite rotation) Euler-Bernoulli beam, 3D, 2 nodes
  const static bool B3D2CR_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "B3D2CR", []( int elementID ) -> Marmot::Discretizations::Elements::Element* {
      return new CorotationalBeamElement3D( elementID );
    } );

} // namespace Marmot::Discretizations::Elements::Registration
