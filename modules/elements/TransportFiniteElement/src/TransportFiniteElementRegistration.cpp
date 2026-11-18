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
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/TransportFiniteElement.h"

/* Naming: DC<section><nodes><fields>, after the diffusive heat transfer elements DC2D4, DC3D8, DCAX4 of Abaqus, with a
 * suffix for the nodal fields: T = temperature. The integration is full: the capacity N^T N is integrated exactly on
 * an undistorted element (degree 2 p).
 */
namespace {
  using namespace Marmot::Discretizations::Elements;

  template < class T,
             Marmot::Discretizations::Elements::Quadrature::QuadratureRule quadratureRule,
             typename T::SectionType                                       sectionType >
  Marmot::Registry::ElementFactory::elementFactoryFunction makeHeatFactoryFunction()
  {
    return []( int elementID ) -> Marmot::Discretizations::Elements::Element* {
      return new T( elementID, quadratureRule, sectionType, { "temperature" } );
    };
  }

  using Q = Marmot::Discretizations::Elements::Quadrature::QuadratureRule;

  const static bool DC2D3T_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "DC2D3T",
                     makeHeatFactoryFunction< TransportFiniteElement< 2, 3, 1 >,
                                              Quadrature::gaussExactToDegree( 2 ),
                                              TransportFiniteElement< 2, 3, 1 >::Plane >() );

  const static bool DC2D4T_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "DC2D4T",
                     makeHeatFactoryFunction< TransportFiniteElement< 2, 4, 1 >,
                                              Quadrature::gaussExactToDegree( 2 ),
                                              TransportFiniteElement< 2, 4, 1 >::Plane >() );

  const static bool DC2D8T_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "DC2D8T",
                     makeHeatFactoryFunction< TransportFiniteElement< 2, 8, 1 >,
                                              Quadrature::gaussExactToDegree( 4 ),
                                              TransportFiniteElement< 2, 8, 1 >::Plane >() );

  const static bool DCAX4T_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "DCAX4T",
                     makeHeatFactoryFunction< TransportFiniteElement< 2, 4, 1 >,
                                              Quadrature::gaussExactToDegree( 2 ),
                                              TransportFiniteElement< 2, 4, 1 >::Axisymmetric >() );

  const static bool DCAX8T_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "DCAX8T",
                     makeHeatFactoryFunction< TransportFiniteElement< 2, 8, 1 >,
                                              Quadrature::gaussExactToDegree( 4 ),
                                              TransportFiniteElement< 2, 8, 1 >::Axisymmetric >() );

  const static bool DC3D4T_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "DC3D4T",
                     makeHeatFactoryFunction< TransportFiniteElement< 3, 4, 1 >,
                                              Quadrature::gaussExactToDegree( 2 ),
                                              TransportFiniteElement< 3, 4, 1 >::Solid >() );

  const static bool DC3D8T_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "DC3D8T",
                     makeHeatFactoryFunction< TransportFiniteElement< 3, 8, 1 >,
                                              Quadrature::gaussExactToDegree( 2 ),
                                              TransportFiniteElement< 3, 8, 1 >::Solid >() );

  const static bool DC3D20T_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "DC3D20T",
                     makeHeatFactoryFunction< TransportFiniteElement< 3, 20, 1 >,
                                              Quadrature::gaussExactToDegree( 4 ),
                                              TransportFiniteElement< 3, 20, 1 >::Solid >() );

} // namespace
