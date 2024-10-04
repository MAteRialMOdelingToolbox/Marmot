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

namespace Marmot {
  namespace NumericalAlgorithms {
    /** Manager for yield surface combinations for multisurface plasticity:
     * Try different yield surface combinations and track already used combinations
     * */
    template < int nYieldSurfaces >
    class YieldSurfaceCombinationManager {
      /// Column in the YieldSurfFlagArr for marking a combination as used
      const int idxUsedFlag;

    public:
      /// An array which contains every possible (reasonable) combination of yield surfaces
      Eigen::Array< bool, ( 1 << nYieldSurfaces ) - 1, ( nYieldSurfaces + 1 ) > yieldSurfaceCombinations;
      /// An array to carry active/nonactive states of yield surfaces
      typedef Eigen::Array< bool, 1, nYieldSurfaces > YieldSurfFlagArr;
      /// An array to carry values of yield functions
      typedef Eigen::Array< double, 1, nYieldSurfaces > YieldSurfResArr;

      YieldSurfaceCombinationManager();
      /// get another unused combination of yield surfaces
      bool getAnotherYieldFlagCombination( YieldSurfFlagArr& activeSurfaces );
      /// set the current combination as used
      void markYieldFlagCombinationAsUsed( const YieldSurfFlagArr& activeSurfaces );
      /// reset all yieldsurfaces as unused
      void resetUsedYieldFlagCombinations();
    };
  } // namespace NumericalAlgorithms
} // namespace Marmot

namespace Marmot {
  namespace NumericalAlgorithms {
    template < int n >
    YieldSurfaceCombinationManager< n >::YieldSurfaceCombinationManager() : idxUsedFlag( n )
    {
      yieldSurfaceCombinations.setZero();
      const int numRows = ( 1 << n ) - 1;
      const int numCols = n;

      for ( int row = 0; row < numRows; row++ )
        for ( int col = 0; col < numCols; col++ )
          yieldSurfaceCombinations( row, col ) = ( row + 1 ) & 1 << ( col );
      return;
    }

    template < int n >
    bool YieldSurfaceCombinationManager< n >::getAnotherYieldFlagCombination( YieldSurfFlagArr& activeSurfaces )
    {
      for ( int i = 0; i < yieldSurfaceCombinations.rows(); i++ ) {
        bool             alreadyUsed          = yieldSurfaceCombinations.row( i )( idxUsedFlag );
        YieldSurfFlagArr combinationCandidate = yieldSurfaceCombinations.row( i ).head( n );
        if ( !alreadyUsed && ( combinationCandidate == true ).any() ) {
          activeSurfaces = combinationCandidate;
          return true;
        }
      }
      return false;
    }

    template < int n >
    void YieldSurfaceCombinationManager< n >::resetUsedYieldFlagCombinations()
    {
      yieldSurfaceCombinations.col( idxUsedFlag ).setConstant( false );
    }

    template < int n >
    void YieldSurfaceCombinationManager< n >::markYieldFlagCombinationAsUsed( const YieldSurfFlagArr& activeSurfaces )
    {
      for ( int i = 0; i < yieldSurfaceCombinations.rows(); i++ )
        if ( ( yieldSurfaceCombinations.row( i ).head( n ) == activeSurfaces ).all() )
          yieldSurfaceCombinations.row( i )( idxUsedFlag ) = true;
    }
  } // namespace NumericalAlgorithms
} // namespace Marmot
