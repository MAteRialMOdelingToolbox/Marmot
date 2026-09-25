#include "Marmot/MarmotMeshfreeKernelFunctionBSpline3rdOrderBoxed.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <cmath>

namespace Marmot::Meshfree {

  MarmotMeshfreeKernelFunctionBSpline3rdOrderBoxed::MarmotMeshfreeKernelFunctionBSpline3rdOrderBoxed(
    double* centerCoord,
    int     dim,
    double  supportRadius )
    : _centerCoord( centerCoord ), _supportRadius( supportRadius ), _dim( dim )
  {
  }

  double MarmotMeshfreeKernelFunctionBSpline3rdOrderBoxed::computeKernelFunction( const double* coord ) const
  {

    double res = 1.0;

    for ( int i = 0; i < _dim; i++ ) {
      res *= computeBSpline3rdOrder( coord[i] - _centerCoord[i] );
    }

    return res;
  }

  void MarmotMeshfreeKernelFunctionBSpline3rdOrderBoxed::computeKernelFunctionGradient( const double* coord,
                                                                                        double*       grad ) const
  {
    for ( int i = 0; i < _dim; i++ ) {
      grad[i] = 0;
    }

    for ( int i = 0; i < _dim; i++ ) {
      double res = 1.0;
      for ( int j = 0; j < _dim; j++ ) {
        if ( i == j )
          res *= computeBSpline3rdOrderGradient( coord[j] - _centerCoord[j] );
        else
          res *= computeBSpline3rdOrder( coord[j] - _centerCoord[j] );
      }
      grad[i] = res;
    }
  }

  double MarmotMeshfreeKernelFunctionBSpline3rdOrderBoxed::computeBSpline3rdOrder( double coord_minus_center ) const
  {
    const double z = std::abs( coord_minus_center ) / _supportRadius;

    if ( z <= 1. / 2 )
      return 2.0 / 3.0 - 4.0 * z * z + 4.0 * z * z * z;
    if ( z <= 1 )
      return 4.0 / 3.0 - 4.0 * z + 4.0 * z * z - 4.0 / 3.0 * z * z * z;
    return 0;
  }

  double MarmotMeshfreeKernelFunctionBSpline3rdOrderBoxed::computeBSpline3rdOrderGradient(
    double coord_minus_center ) const
  {
    const double z         = std::abs( coord_minus_center ) / _supportRadius;
    const double dz_dcoord = coord_minus_center > 0 ? 1.0 / _supportRadius : -1.0 / _supportRadius;

    if ( z <= 1. / 2 )
      return ( -8.0 * z + 12.0 * z * z ) * dz_dcoord;
    if ( z <= 1 )
      return ( -4.0 + 8.0 * z - 4.0 * z * z ) * dz_dcoord;
    return 0;
  }

  const double* MarmotMeshfreeKernelFunctionBSpline3rdOrderBoxed::getCenterCoordinates() const
  {
    return _centerCoord;
  }

  void MarmotMeshfreeKernelFunctionBSpline3rdOrderBoxed::moveTo( const double* coordinate )
  {
    for ( int i = 0; i < _dim; i++ ) {
      _centerCoord[i] = coordinate[i];
    }
  }

  bool MarmotMeshfreeKernelFunctionBSpline3rdOrderBoxed::isInSupport( const double* coord ) const
  {
    return computeKernelFunction( coord ) > 0;
  }

  void MarmotMeshfreeKernelFunctionBSpline3rdOrderBoxed::getBoundingBox( double* min, double* max ) const
  {
    for ( int i = 0; i < _dim; i++ ) {
      min[i] = _centerCoord[i] - _supportRadius;
      max[i] = _centerCoord[i] + _supportRadius;
    }
  }

}; // namespace Marmot::Meshfree
