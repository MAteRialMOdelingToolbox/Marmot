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
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotJournal.h"
#include <Eigen/Dense>
#include <array>
#include <stdexcept>
#include <vector>

/**
 * @brief Kinematics shared by the bond elements of embedded bars (EmbeddedBondElement, EmbeddedLargeSlipBondElement):
 * the integration rule along the bar, the local frame of the bar axis and its transport to the current axis.
 */
namespace Marmot::Discretizations::Elements::EmbeddedBondKinematics {

  /**
   * @brief The Gauss-Lobatto rule along a bar part, from the quadrature catalogue: its end points sit at the ends of
   * the part, which avoids the traction oscillations of Gauss rules for stiff bond.
   * @param n number of points, 2 to 5 (degree of exactness 2n - 3)
   * @return points and weights on [-1, 1]
   * @throws std::invalid_argument for other numbers of points
   */
  inline const std::vector< Quadrature::QuadraturePointInfo >& gaussLobatto( int n )
  {
    if ( n < 2 || n > 5 )
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << "bond elements: 2 to 5 integration points are supported, " << n
                                   << " were requested" );
    return Quadrature::getQuadraturePoints( ReferenceDomain::Line, Quadrature::gaussLobattoExactToDegree( 2 * n - 3 ) );
  }

  /**
   * @brief Local frame of the unit tangent t, columns [t, n1, (n2)].
   *
   * In 2D, n1 is t rotated by +90 degrees. In 3D, the normals are an arbitrary right-handed pair; the global axis least
   * aligned with t gives a well conditioned first normal. They may differ between integration points, so the normal
   * components of slip and bond stress must not be compared between points, and bond-slip laws must be isotropic in
   * the normal plane.
   */
  template < int nDim >
  Eigen::Matrix< double, nDim, nDim > localFrame( const Eigen::Matrix< double, nDim, 1 >& t )
  {
    Eigen::Matrix< double, nDim, nDim > R;
    if constexpr ( nDim == 2 ) {
      R.col( 0 ) = t;
      R.col( 1 ) = Eigen::Vector2d( -t( 1 ), t( 0 ) );
    }
    else {
      int axis;
      t.cwiseAbs().minCoeff( &axis );
      const Eigen::Vector3d n1 = t.cross( Eigen::Vector3d::Unit( axis ) ).normalized();
      R.col( 0 )               = t;
      R.col( 1 )               = n1;
      R.col( 2 )               = t.cross( n1 );
    }
    return R;
  }

  /// A frame [t, n1, (n2)] following the current bar axis, and the derivatives of its normals w.r.t. the tangent
  template < int nDim >
  struct TransportedFrame {
    Eigen::Matrix< double, nDim, nDim >                         R;                ///< columns [t, n1, (n2)]
    std::array< Eigen::Matrix< double, nDim, nDim >, nDim - 1 > dNormal_dTangent; ///< d n_i / d t
  };

  /**
   * @brief The reference frame carried to the current unit tangent t by the smallest rotation taking the reference
   * tangent to t.
   *
   * In 2D, \f$\boldsymbol{n}_1\f$ is t rotated by +90 degrees. In 3D, the rotation (Rodrigues, about
   * \f$\boldsymbol{v} = \boldsymbol{t}_0 \times \boldsymbol{t}\f$ with \f$c = \boldsymbol{t}_0 \cdot \boldsymbol{t}\f$)
   * \f[ \boldsymbol{Q} = c\,\boldsymbol{I} + [\boldsymbol{v}]_\times + \frac{\boldsymbol{v} \otimes \boldsymbol{v}}{1 +
   * c}, \qquad \boldsymbol{n}_i = \boldsymbol{Q}\, \boldsymbol{n}_i^0 \f]
   * depends on t only: it follows the bending of the bar, but not a rotation of the bar about its own axis. The
   * normal components are therefore objective for laws isotropic in the normal plane (all laws here); a stored
   * normal history (e.g., the plastic normal slip) is not rotated with a twist of the bar.
   *
   * @param R0 the reference frame [t0, n1^0, (n2^0)]
   * @param t  the current unit tangent
   * @throws std::runtime_error if t is opposite to t0 (the bar turned around)
   */
  template < int nDim >
  TransportedFrame< nDim > transportFrame( const Eigen::Matrix< double, nDim, nDim >& R0,
                                           const Eigen::Matrix< double, nDim, 1 >&    t )
  {
    TransportedFrame< nDim > frame;
    frame.R.col( 0 ) = t;
    if constexpr ( nDim == 2 ) {
      frame.R.col( 1 ) = Eigen::Vector2d( -t( 1 ), t( 0 ) );
      frame.dNormal_dTangent[0] << 0, -1, 1, 0;
    }
    else {
      const Eigen::Vector3d t0 = R0.col( 0 );
      const Eigen::Vector3d v  = t0.cross( t );
      const double          c  = t0.dot( t );
      if ( c <= -1 + 1e-8 )
        throw std::runtime_error( "bond element: the bar axis turned around against its reference direction" );
      auto skew = []( const Eigen::Vector3d& a ) {
        Eigen::Matrix3d S;
        S << 0, -a( 2 ), a( 1 ), a( 2 ), 0, -a( 0 ), -a( 1 ), a( 0 ), 0;
        return S;
      };
      const Eigen::Matrix3d Q = c * Eigen::Matrix3d::Identity() + skew( v ) + v * v.transpose() / ( 1 + c );
      // dc/dt = t0^T, dv/dt = [t0]x
      const Eigen::Matrix3d dv_dt = skew( t0 );
      for ( int i = 0; i < 2; i++ ) {
        const Eigen::Vector3d n0  = R0.col( 1 + i );
        frame.R.col( 1 + i )      = Q * n0;
        const double vn           = v.dot( n0 );
        frame.dNormal_dTangent[i] = n0 * t0.transpose() - skew( n0 ) * dv_dt +
                                    ( vn * dv_dt + v * ( n0.transpose() * dv_dt ) ) / ( 1 + c ) -
                                    v * vn * t0.transpose() / ( ( 1 + c ) * ( 1 + c ) );
      }
    }
    return frame;
  }

} // namespace Marmot::Discretizations::Elements::EmbeddedBondKinematics
