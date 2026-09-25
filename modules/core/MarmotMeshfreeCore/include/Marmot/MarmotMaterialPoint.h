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
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotUtils.h"
#include <string>

class MarmotMaterialPoint {

public:
  virtual ~MarmotMaterialPoint(){};

  virtual void assignStateVars( double* stateVars, int nStateVars ) = 0;

  virtual void assignMaterial( const MarmotMaterialSection& material ) = 0;

  virtual void initializeYourself() = 0;

  virtual void prepareYourself( double timeNew, double dT ) = 0;

  virtual void computeYourself( double timeNew, double dT ) = 0;

  virtual void acceptStateAndPosition(){};

  virtual StateView getStateView( const std::string& stateName ) const = 0;

  virtual int getNumberOfRequiredStateVars() const = 0;

  virtual void getCoordinatesAtCenter( double* coordinates ) const = 0;

  virtual int getMaterialPointNumber() const = 0;

  virtual int getDimension() const = 0;

  virtual int getNumberOfVertices() const = 0;

  virtual std::string getMaterialPointShape() const = 0;

  virtual void getVertexCoordinates( double* coordinates ) const = 0;

  virtual void getCenterDisplacement( double* displacement ) const = 0;

  virtual double getVolumeUndeformed() const = 0;

  virtual double getDensityUndeformed() const = 0;

  virtual void setInitialCondition( const std::string& conditionName, const double* value ) = 0;
};
