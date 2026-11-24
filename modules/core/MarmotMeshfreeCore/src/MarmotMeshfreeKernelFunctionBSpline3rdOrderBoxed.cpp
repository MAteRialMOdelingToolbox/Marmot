/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck
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
#include "Marmot/MarmotMeshfreeKernelFunctionBSpline3rdOrderBoxed.h"
#include "Marmot/MarmotJournal.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <stdexcept>

namespace Marmot::Discretizations::Meshfree {

  KernelFunctionBSpline3rdOrderBoxed::KernelFunctionBSpline3rdOrderBoxed( double* centerCoord,
                                                                          int     dim,
                                                                          double  supportRadius )
    : _centerCoord( centerCoord ), _supportRadius( supportRadius ), _dim( dim )
  {
    if ( dim < 1 || dim > maxDim )
      throw std::invalid_argument( Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << ": dimension " << dim
                                                                   << " is not supported (1 to " << maxDim << ")" );
  }

  double KernelFunctionBSpline3rdOrderBoxed::computeKernelFunction( const double* coord ) const
  {
    double phi = 1.0;
    for ( int i = 0; i < _dim; i++ )
      phi *= computeBSpline3rdOrder( coord[i] - _centerCoord[i] );
    return phi;
  }

  double KernelFunctionBSpline3rdOrderBoxed::computeKernelFunctionAndGradient( const double* coord, double* grad ) const
  {
    // the one-dimensional factors w(x_i - x_A,i) and their derivatives, once per direction
    std::array< double, maxDim > w, dw_dr;
    for ( int i = 0; i < _dim; i++ )
      w[i] = computeBSpline3rdOrder( coord[i] - _centerCoord[i], &dw_dr[i] );

    double phi = 1.0;
    for ( int i = 0; i < _dim; i++ )
      phi *= w[i];

    // product rule: dphi/dx_i = w'_i prod_(j != i) w_j
    for ( int i = 0; i < _dim; i++ ) {
      grad[i] = 1.0;
      for ( int j = 0; j < _dim; j++ )
        grad[i] *= i == j ? dw_dr[j] : w[j];
    }
    return phi;
  }

  double KernelFunctionBSpline3rdOrderBoxed::computeBSpline3rdOrder( double coord_minus_center, double* dw_dr ) const
  {
    const double z = std::abs( coord_minus_center ) / _supportRadius;
    double       w, dw_dz;
    if ( z <= 1. / 2 ) {
      w     = 2.0 / 3.0 - 4.0 * z * z + 4.0 * z * z * z;
      dw_dz = -8.0 * z + 12.0 * z * z;
    }
    else if ( z <= 1 ) {
      w     = 4.0 / 3.0 - 4.0 * z + 4.0 * z * z - 4.0 / 3.0 * z * z * z;
      dw_dz = -4.0 + 8.0 * z - 4.0 * z * z;
    }
    else {
      w     = 0;
      dw_dz = 0;
    }
    if ( dw_dr )
      *dw_dr = dw_dz * ( coord_minus_center > 0 ? 1.0 / _supportRadius : -1.0 / _supportRadius );
    return w;
  }

  const double* KernelFunctionBSpline3rdOrderBoxed::getCenterCoordinates() const
  {
    return _centerCoord;
  }

  void KernelFunctionBSpline3rdOrderBoxed::moveTo( const double* coordinate )
  {
    for ( int i = 0; i < _dim; i++ ) {
      _centerCoord[i] = coordinate[i];
    }
  }

  bool KernelFunctionBSpline3rdOrderBoxed::isInSupport( const double* coord ) const
  {
    return computeKernelFunction( coord ) > 0;
  }

  void KernelFunctionBSpline3rdOrderBoxed::getBoundingBox( double* min, double* max ) const
  {
    for ( int i = 0; i < _dim; i++ ) {
      min[i] = _centerCoord[i] - _supportRadius;
      max[i] = _centerCoord[i] + _supportRadius;
    }
  }

}; // namespace Marmot::Discretizations::Meshfree
