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
 * suffix for the nodal fields: T = temperature, H = (pore) relative humidity, TH = both (coupled heat and moisture
 * transport, fields in this order). The integration is full: the capacity
 * N^T N is integrated exactly on an undistorted element (degree 2 p).
 */
namespace {
  using namespace Marmot::Discretizations::Elements;

  /// the nodal field of the single field transport elements, by suffix
  struct Temperature {
    static constexpr const char* name = "temperature";
  };
  struct RelativeHumidity {
    static constexpr const char* name = "relative humidity";
  };

  /// the nodal fields of a transport element: the names, in the order of the dofs at a node
  template < class... Field >
  struct Fields {
    static constexpr int                      nFields = sizeof...( Field );
    static std::array< std::string, nFields > names() { return { Field::name... }; }
  };

  template < class T,
             class Field,
             Marmot::Discretizations::Elements::Quadrature::QuadratureRule quadratureRule,
             typename T::SectionType                                       sectionType >
  Marmot::Registry::ElementFactory::elementFactoryFunction makeTransportFactoryFunction()
  {
    return []( int elementID ) -> Marmot::Discretizations::Elements::Element* {
      return new T( elementID, quadratureRule, sectionType, Field::names() );
    };
  }

  /// register the eight element types DC2D3, DC2D4, DC2D8, DCAX4, DCAX8, DC3D4, DC3D8, DC3D20 + suffix
  template < class Field >
  bool registerElements( const std::string& suffix )
  {
    using Marmot::Registry::ElementFactory;
    constexpr int nF = Field::nFields;
    using P3         = TransportFiniteElement< 2, 3, nF >;
    using P4         = TransportFiniteElement< 2, 4, nF >;
    using P8         = TransportFiniteElement< 2, 8, nF >;
    using S4         = TransportFiniteElement< 3, 4, nF >;
    using S8         = TransportFiniteElement< 3, 8, nF >;
    using S20        = TransportFiniteElement< 3, 20, nF >;

    constexpr auto q2 = Quadrature::gaussExactToDegree( 2 );
    constexpr auto q4 = Quadrature::gaussExactToDegree( 4 );

    bool ok = true;
    ok &= ElementFactory::registerElement( "DC2D3" + suffix,
                                           makeTransportFactoryFunction< P3, Field, q2, P3::Plane >() );
    ok &= ElementFactory::registerElement( "DC2D4" + suffix,
                                           makeTransportFactoryFunction< P4, Field, q2, P4::Plane >() );
    ok &= ElementFactory::registerElement( "DC2D8" + suffix,
                                           makeTransportFactoryFunction< P8, Field, q4, P8::Plane >() );
    ok &= ElementFactory::registerElement( "DCAX4" + suffix,
                                           makeTransportFactoryFunction< P4, Field, q2, P4::Axisymmetric >() );
    ok &= ElementFactory::registerElement( "DCAX8" + suffix,
                                           makeTransportFactoryFunction< P8, Field, q4, P8::Axisymmetric >() );
    ok &= ElementFactory::registerElement( "DC3D4" + suffix,
                                           makeTransportFactoryFunction< S4, Field, q2, S4::Solid >() );
    ok &= ElementFactory::registerElement( "DC3D8" + suffix,
                                           makeTransportFactoryFunction< S8, Field, q2, S8::Solid >() );
    ok &= ElementFactory::registerElement( "DC3D20" + suffix,
                                           makeTransportFactoryFunction< S20, Field, q4, S20::Solid >() );
    return ok;
  }

  // heat conduction: DC2D3T ... DC3D20T
  const static bool heatElementsAreRegistered = registerElements< Fields< Temperature > >( "T" );
  // moisture diffusion: DC2D3H ... DC3D20H
  const static bool moistureElementsAreRegistered = registerElements< Fields< RelativeHumidity > >( "H" );
  // coupled heat and moisture transport: DC2D3TH ... DC3D20TH, fields temperature, relative humidity
  const static bool thermoHygroElementsAreRegistered = registerElements< Fields< Temperature, RelativeHumidity > >(
    "TH" );

} // namespace
