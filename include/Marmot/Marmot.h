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

/**
 * @class MarmotMaterialFactory
 * @brief Factory class for creating material instances.
 *
 * This class provides a mechanism to register materials by their code and name,
 * and to create material instances based on their properties.
 * It allows for dynamic material creation without hardcoding specific material types.
 */
  class MarmotMaterialFactory {
  public:
    using materialFactoryFunction = MarmotMaterial* (*)( const double* materialProperties,
                                                         int           nMaterialProperties,
                                                         int           materialNumber );
    MarmotMaterialFactory()       = delete;

    /**
     * @brief Get the unique material code from its name.
     * @param materialName Name of the material.
     * @return Unique code associated with the material name, or -1 if not found.
     */
    static int getMaterialCodeFromName( const std::string& materialName );

    /**
     * @brief Create a material instance based on its code and properties.
     * @param materialCode Unique code for the material.
     * @param materialProperties Array of material properties.
     * @param nMaterialProperties Number of properties in the array.
     * @param materialNumber Unique identifier for the material instance.
     * @return Pointer to the created MarmotMaterial instance, or nullptr if creation failed.
     */
    static MarmotMaterial* createMaterial( int           materialCode,
                                           const double* materialProperties,
                                           int           nMaterialProperties,
                                           int           materialNumber );

    /**
     * @brief Register a material with its code and factory function.
     * @param materialCode Unique code for the material.
     * @param materialName Name of the material.
     * @param factoryFunction Function to create material instances.
     * @return True if registration was successful, false if the code already exists.
     */
    static bool registerMaterial( int                     materialCode,
                                  const std::string&      materialName,
                                  materialFactoryFunction factoryFunction );

  private:
    static std::unordered_map< std::string, int >             materialNameToCodeAssociation;
    static std::unordered_map< int, materialFactoryFunction > materialFactoryFunctionByCode;
  };

/**
 * @class MarmotElementFactory
 * @brief Factory class for creating element instances.
 * This class provides a mechanism to register elements by their code and name,
 * and to create element instances based on their properties.
 */

  class MarmotElementFactory {
  public:
    using elementFactoryFunction = MarmotElement* (*)( int elementNumber );
    MarmotElementFactory()       = delete;

    /**
     * @brief Get the unique element code from its name.
     * @param elementName Name of the element.
     * @return Unique code associated with the element name, or throws an exception if not found.
     */
    static int getElementCodeFromName( const std::string& elementName );

    /**
     * @brief Create an element instance based on its code and number.
     * @param elementCode Unique code for the element.
     * @param elementNumber Unique identifier for the element instance.
     * @return Pointer to the created MarmotElement instance, or nullptr if creation failed.
     */
    static MarmotElement* createElement( int elementCode, int elementNumber );

    /**
     * @brief Register an element with its code and factory function.
     * @param elementName Name of the element.
     * @param elementCode Unique code for the element.
     * @param factoryFunction Function to create element instances.
     * @return True if registration was successful, false if the code already exists.
     */
    static bool registerElement( const std::string&     elementName,
                                 int                    elementCode,
                                 elementFactoryFunction factoryFunction );

  private:
    static std::unordered_map< std::string, int >            elementNameToCodeAssociation;
    static std::unordered_map< int, elementFactoryFunction > elementFactoryFunctionByCode;
  };

} // namespace MarmotLibrary
