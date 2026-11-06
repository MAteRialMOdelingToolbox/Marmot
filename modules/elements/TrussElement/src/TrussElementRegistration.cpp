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
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/TrussElement.h"

namespace Marmot::Discretizations::Elements::Registration {

  template < int nDim, int nNodes, typename TrussElement< nDim, nNodes >::Kinematics kinematics >
  bool registerTruss( const std::string& name )
  {
    return Marmot::Registry::ElementFactory::registerElement( name,
                                                              []( int elementID )
                                                                -> Marmot::Discretizations::Elements::Element* {
                                                                return new TrussElement< nDim, nNodes >( elementID,
                                                                                                         kinematics );
                                                              } );
  }

  // small strain, hypoelastic materials
  const static bool T2D2_isRegistered = registerTruss< 2, 2, TrussElement< 2, 2 >::SmallStrain >( "T2D2" );
  const static bool T2D3_isRegistered = registerTruss< 2, 3, TrussElement< 2, 3 >::SmallStrain >( "T2D3" );
  const static bool T3D2_isRegistered = registerTruss< 3, 2, TrussElement< 3, 2 >::SmallStrain >( "T3D2" );
  const static bool T3D3_isRegistered = registerTruss< 3, 3, TrussElement< 3, 3 >::SmallStrain >( "T3D3" );

  // finite strain, finite strain materials
  const static bool T2D2CR_isRegistered = registerTruss< 2, 2, TrussElement< 2, 2 >::FiniteStrain >( "T2D2CR" );
  const static bool T2D3CR_isRegistered = registerTruss< 2, 3, TrussElement< 2, 3 >::FiniteStrain >( "T2D3CR" );
  const static bool T3D2CR_isRegistered = registerTruss< 3, 2, TrussElement< 3, 2 >::FiniteStrain >( "T3D2CR" );
  const static bool T3D3CR_isRegistered = registerTruss< 3, 3, TrussElement< 3, 3 >::FiniteStrain >( "T3D3CR" );

} // namespace Marmot::Discretizations::Elements::Registration
