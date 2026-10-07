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

#include "Marmot/MarmotMeshfreeKernelFunction.h"
#include <cmath>

namespace Marmot::Meshfree {

  class MarmotMeshfreeKernelFunctionBSpline2ndOrderBoxed : public MarmotMeshfreeKernelFunction {

  private:
    double*      _centerCoord;
    const double _supportRadius;
    const int    _dim;

  public:
    MarmotMeshfreeKernelFunctionBSpline2ndOrderBoxed( double* centerCoord, int dim, double supportRadius );

    virtual ~MarmotMeshfreeKernelFunctionBSpline2ndOrderBoxed() = default;

    double computeKernelFunction( const double* coord ) const override;

    void computeKernelFunctionGradient( const double* coord, double* grad ) const override;

    double computeBSpline2ndOrder( double coord_minus_center ) const;

    double computeBSpline2ndOrderGradient( double coord_minus_center ) const;

    const double* getCenterCoordinates() const override;

    void moveTo( const double* coord ) override;

    bool isInSupport( const double* coord ) const override;

    void getBoundingBox( double* min, double* max ) const override;
  };

} // namespace Marmot::Meshfree
