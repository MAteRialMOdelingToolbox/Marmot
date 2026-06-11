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
#pragma once
#include "Marmot/MarmotTypedefs.h"
#include <map>
#include <vector>

namespace Marmot {

  /**
   * @brief Utilities for constructing DOF-layout structures in finite element assemblies.
   */
  namespace FiniteElement {

    /**
     * @brief Builds a node-field layout from a map of field names to their DOF dimensions.
     *
     * @param fieldSizes Map from field name to a pair (nodesPerElement, dofsPerNode).
     * @return A 2-D vector where each row corresponds to a node and each entry is a field name.
     */
    const std::vector< std::vector< std::string > > makeNodeFieldLayout(
      const std::map< std::string, std::pair< int, int > >& fieldSizes );

    /**
     * @brief Computes the DOF permutation pattern for a blocked (field-major) layout.
     *
     * @param nodeFields     Node-field layout as returned by @ref makeNodeFieldLayout.
     * @param fieldSizes     Map from field name to a pair (nodesPerElement, dofsPerNode).
     * @return Integer permutation vector that reorders DOFs from node-major to field-major order.
     */
    std::vector< int > makeBlockedLayoutPermutationPattern(
      const std::vector< std::vector< std::string > >&      nodeFields,
      const std::map< std::string, std::pair< int, int > >& fieldSizes );

  } // namespace FiniteElement
} // namespace Marmot
