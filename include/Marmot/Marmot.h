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
 * Matthias Neuner matthias.neuner@uibk.ac.at
 * Magdalena Schreter magdalena.schreter@uibk.ac.at
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
#pragma once
#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotMaterial.h"
#include <functional>
#include <iostream>
#include <string>
#include <unordered_map>

namespace MarmotLibrary {

  // MaterialFactory
  //
  // - Allows materials to register themselve with their name and ID
  // - Allows the user to create instances of materials

  class MarmotMaterialFactory {
  public:
    using materialFactoryFunction = MarmotMaterial* (*)( const double* materialProperties,
                                                         int           nMaterialProperties,
                                                         int           materialNumber );
    MarmotMaterialFactory()       = delete;

    static int getMaterialCodeFromName( const std::string& materialName );

    static MarmotMaterial* createMaterial( int           materialCode,
                                           const double* materialProperties,
                                           int           nMaterialProperties,
                                           int           materialNumber );

    static bool registerMaterial( int                     materialCode,
                                  const std::string&      materialName,
                                  materialFactoryFunction factoryFunction );

  private:
    static std::unordered_map< std::string, int >             materialNameToCodeAssociation;
    static std::unordered_map< int, materialFactoryFunction > materialFactoryFunctionByCode;
  };

  // ElementFactory
  //
  // - Allows elements to register themselve with their name and ID
  // - Allows the user to create instances of elements

  class MarmotElementFactory {
  public:
    using elementFactoryFunction = MarmotElement* (*)( int elementNumber );
    MarmotElementFactory()       = delete;

    static int getElementCodeFromName( const std::string& elementName );

    static MarmotElement* createElement( int elementCode, int elementNumber );

    static bool registerElement( const std::string&     elementName,
                                 int                    elementCode,
                                 elementFactoryFunction factoryFunction );

  private:
    static std::unordered_map< std::string, int >            elementNameToCodeAssociation;
    static std::unordered_map< int, elementFactoryFunction > elementFactoryFunctionByCode;
  };

} // namespace MarmotLibrary
