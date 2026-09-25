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
#include "Marmot/MarmotMaterialPoint.h"
#include "Marmot/MarmotJournal.h"
#include <string>
#include <unordered_map>
#include <vector>
#include <stdexcept>

class MarmotCell {

public:
  virtual ~MarmotCell() = default;

  virtual const std::vector< std::vector< std::string > >& getNodeFields() const = 0;

  virtual const std::vector< int >& getDofIndicesPermutationPattern() const = 0;

  virtual int getNNodes() const = 0;

  virtual int getNDofPerCell() const = 0;

  virtual std::string getCellShape() const = 0;

  virtual bool isCoordinateInCell( const double* coordinates ) const = 0;

  virtual void getBoundingBox( double* boundingBoxMin, double* boundingBoxMax ) const = 0;

  virtual void assignMaterialPoints( const std::vector< MarmotMaterialPoint* >& materialPoints ) = 0;

  virtual void computeMaterialPointKernels( const double* Q,
                                            double*       fInt,
                                            double*       dFInt_dQ,
                                            double        timeNew,
                                            double        dT ) const = 0;
  
  virtual void computeLumpedInertia( double* I ) {
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << "not yet implemented" );
  }
  
  virtual void computeConsistentInertia( double* I ) {
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << "not yet implemented" );
  }

  virtual void computeBodyLoad( int type, const double* load, double* fExt, double* dExt_dQ, double timeNew, double dT )
    const = 0;

  virtual void computeDistributedLoad( int           type,
                                       int           surfaceID,
                                       int           materialPointNumber,
                                       const double* load,
                                       double*       fExt,
                                       double*       dExt_dQ,
                                       double        timeNew,
                                       double        dT ) const = 0;

  virtual void getInterpolationVector( double* vec, const double* coordinates ) const = 0;

  virtual const std::unordered_map< std::string, int >& getSupportedBodyLoadTypes() const = 0;

  virtual const std::unordered_map< std::string, int >& getSupportedDistributedLoadTypes() const = 0;

  virtual void interpolateFieldsToMaterialPoints( const double* Q ) const = 0;
};
