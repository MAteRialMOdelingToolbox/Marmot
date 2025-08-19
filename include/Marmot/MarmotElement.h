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
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotUtils.h"
#include <stdexcept>
#include <string>
#include <vector>

class MarmotElement {

public:
  enum StateTypes {
    Sigma11,
    Sigma22,
    Sigma33,
    HydrostaticStress,
    GeostaticStress,
    MarmotMaterialStateVars,
    MarmotMaterialInitialization,
    HasEigenDeformation,
  };

  enum DistributedLoadTypes {
    Pressure,
    SurfaceTorsion,
    SurfaceTraction,
  };

  virtual ~MarmotElement();

  virtual int getNumberOfRequiredStateVars() = 0;

  virtual std::vector< std::vector< std::string > > getNodeFields() = 0;

  virtual std::vector< int > getDofIndicesPermutationPattern() = 0;

  virtual int getNNodes() = 0;

  virtual int getNSpatialDimensions() = 0;

  virtual int getNDofPerElement() = 0;

  virtual std::string getElementShape() = 0;

  virtual void assignStateVars( double* stateVars, int nStateVars ) = 0;

  virtual void assignProperty( const ElementProperties& property );

  virtual void assignProperty( const MarmotMaterialSection& property );

  virtual void assignNodeCoordinates( const double* coordinates ) = 0;

  virtual void initializeYourself() = 0;

  virtual void setInitialConditions( StateTypes state, const double* values ) = 0;

  virtual void computeYourself( const double* QTotal,
                                const double* dQ,
                                double*       Pint,
                                double*       K,
                                const double* time,
                                double        dT,
                                double&       pNewdT ) = 0;

  virtual void computeDistributedLoad( DistributedLoadTypes loadType,
                                       double*              Pext,
                                       double*              K,
                                       int                  elementFace,
                                       const double*        load,
                                       const double*        QTotal,
                                       const double*        time,
                                       double               dT ) = 0;

  virtual void computeBodyForce( double*       Pext,
                                 double*       K,
                                 const double* load,
                                 const double* QTotal,
                                 const double* time,
                                 double        dT ) = 0;

  virtual void computeLumpedInertia( double* I )
  {
    throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << "not yet implemented" );
  };

  virtual void computeConsistentInertia( double* I )
  {
    throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << "not yet implemented" );
  };

  virtual StateView getStateView( const std::string& stateName, int quadraturePoint ) = 0;

  virtual std::vector< double > getCoordinatesAtCenter() = 0;

  virtual std::vector< std::vector< double > > getCoordinatesAtQuadraturePoints() = 0;

  virtual int getNumberOfQuadraturePoints() = 0;
};
