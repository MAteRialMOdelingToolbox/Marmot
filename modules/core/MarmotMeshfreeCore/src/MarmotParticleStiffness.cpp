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

#include "Marmot/MarmotParticleStiffness.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace Marmot::Discretizations::Meshfree {

  template < int nDim >
  void AbsoluteStiffnessBound< nDim >::visitShapeFunctions( int           nDimOfParticle,
                                                            int           nDofPerNode,
                                                            int           nNodes,
                                                            const double* N,
                                                            const double* dN_dY,
                                                            const double* d2N_dYdY )
  {
    if ( nDimOfParticle != nDim )
      throw std::invalid_argument( "AbsoluteStiffnessBound: a particle of dimension " +
                                   std::to_string( nDimOfParticle ) + " for a bound of dimension " +
                                   std::to_string( nDim ) );
    _nDofPerNode         = nDofPerNode;
    _nDof                = nDofPerNode * nNodes;
    _withSecondGradients = d2N_dYdY != nullptr;
    _moments.assign( _nColumns * nDofPerNode, Moments() );

    for ( int c = 0; c < _nColumns; c++ )
      for ( int k = 0; k < nDofPerNode; k++ ) {
        Moments& m = _moments[c * nDofPerNode + k];
        for ( int B = 0; B < nNodes; B++ ) {
          const double w_B = _columnWeights[c * _nDof + nDofPerNode * B + k];
          if ( w_B == 0.0 )
            continue;
          m.active = true;
          m.sumOfValues += w_B * std::abs( N[B] );
          const Eigen::Matrix< double, nDim, 1 > g_B = Eigen::Map< const Eigen::Matrix< double, nDim, 1 > >( dN_dY +
                                                                                                             nDim * B );
          addMoments( g_B, w_B, m.sumOfGradientNorms, m.gradientMoments );
          if ( _withSecondGradients )
            addMoments( symmetricComponents( d2N_dYdY + nDim * nDim * B ),
                        w_B,
                        m.sumOfSecondGradientNorms,
                        m.secondGradientMoments );
        }
      }
  }

  template < int nDim >
  void AbsoluteStiffnessBound< nDim >::visitStiffnessOfNode( int           A,
                                                             const double* df_dIncrement,
                                                             const double* df_dIncrementGradient,
                                                             const double* df_dIncrementSecondGradient )
  {
    // the dofs per node of the particles of Marmot: the displacement, and with one or six nonlocal fields
    switch ( _nDofPerNode ) {
    case nDim:
      return addBoundOfRowsOfNode< nDim >( A, df_dIncrement, df_dIncrementGradient, df_dIncrementSecondGradient );
    case nDim + 1:
      return addBoundOfRowsOfNode< nDim + 1 >( A, df_dIncrement, df_dIncrementGradient, df_dIncrementSecondGradient );
    case nDim + 6:
      return addBoundOfRowsOfNode< nDim + 6 >( A, df_dIncrement, df_dIncrementGradient, df_dIncrementSecondGradient );
    default: return addBoundOfRowsOfNode< 0 >( A, df_dIncrement, df_dIncrementGradient, df_dIncrementSecondGradient );
    }
  }

  template < int nDim >
  template < int nDofPerNodeAtCompileTime >
  inline void AbsoluteStiffnessBound< nDim >::addBoundOfRowsOfNode( int           A,
                                                                    const double* df_dIncrement,
                                                                    const double* df_dIncrementGradient,
                                                                    const double* df_dIncrementSecondGradient )
  {
    constexpr int nAtCompileTime = nDofPerNodeAtCompileTime > 0 ? nDofPerNodeAtCompileTime : Eigen::Dynamic;
    using Block                  = Eigen::Map< const Eigen::Matrix< double, nAtCompileTime, nAtCompileTime > >;
    const int n                  = nDofPerNodeAtCompileTime > 0 ? nDofPerNodeAtCompileTime : _nDofPerNode;
    // block l of a derivative, in the layout of StiffnessVisitor
    const auto block = [n]( const double* derivative, int l ) { return Block( derivative + n * n * l, n, n ); };

    for ( int i = 0; i < n; i++ )
      for ( int k = 0; k < n; k++ ) {
        const double                     a = block( df_dIncrement, 0 )( i, k );
        Eigen::Matrix< double, nDim, 1 > b;
        for ( int l = 0; l < nDim; l++ )
          b( l ) = block( df_dIncrementGradient, l )( i, k );
        Eigen::Matrix< double, nSym, 1 > h = Eigen::Matrix< double, nSym, 1 >::Zero();
        if ( _withSecondGradients )
          for ( int l = 0, s = 0; l < nDim; l++ )
            for ( int m = l; m < nDim; m++, s++ )
              h( s ) = block( df_dIncrementSecondGradient, l * nDim + m )( i, k ) +
                       ( m != l ? block( df_dIncrementSecondGradient, m * nDim + l )( i, k ) : 0.0 );

        // Most factors of a node's rows are structurally zero (e.g. a displacement row has no gradient factor for a
        // nonlocal column, and only the stabilization has second-gradient factors): a zero factor adds nothing.
        const bool withValueTerm          = a != 0.0;
        const bool withGradientTerm       = ( b.array() != 0.0 ).any();
        const bool withSecondGradientTerm = ( h.array() != 0.0 ).any();
        if ( !withValueTerm && !withGradientTerm && !withSecondGradientTerm )
          continue;

        for ( int c = 0; c < _nColumns; c++ ) {
          const Moments& m = _moments[c * n + k];
          if ( !m.active )
            continue;
          double sum = std::abs( a ) * m.sumOfValues;
          if ( withGradientTerm )
            sum += std::sqrt( std::max( 0.0, m.sumOfGradientNorms * b.dot( m.gradientMoments * b ) ) );
          if ( withSecondGradientTerm )
            sum += std::sqrt( std::max( 0.0, m.sumOfSecondGradientNorms * h.dot( m.secondGradientMoments * h ) ) );
          _bound[c * _nDof + n * A + i] += sum;
        }
      }
  }

  template < int nDim >
  template < int n >
  inline void AbsoluteStiffnessBound< nDim >::addMoments( const Eigen::Matrix< double, n, 1 >& g,
                                                          double                               w,
                                                          double&                              sumOfNorms,
                                                          Eigen::Matrix< double, n, n >&       moments )
  {
    const double norm = g.norm();
    if ( norm == 0.0 )
      return;
    sumOfNorms += w * norm;
    moments += ( w / norm ) * g * g.transpose();
  }

  template < int nDim >
  inline auto AbsoluteStiffnessBound< nDim >::symmetricComponents( const double* d2N_dYdY )
    -> Eigen::Matrix< double, nSym, 1 >
  {
    Eigen::Matrix< double, nSym, 1 > components;
    for ( int l = 0, s = 0; l < nDim; l++ )
      for ( int m = l; m < nDim; m++, s++ )
        components( s ) = d2N_dYdY[l * nDim + m];
    return components;
  }

  template class AbsoluteStiffnessBound< 2 >;
  template class AbsoluteStiffnessBound< 3 >;

  void addBoundOfAbsoluteStiffnessTimes( Particle&     particle,
                                         const double* columnWeights,
                                         int           nColumns,
                                         double*       bound,
                                         double        timeNew,
                                         double        dT )
  {
    switch ( particle.getDimension() ) {
    case 2: {
      AbsoluteStiffnessBound< 2 > absoluteStiffness( columnWeights, nColumns, bound );
      particle.describeStiffness( absoluteStiffness, timeNew, dT );
      break;
    }
    case 3: {
      AbsoluteStiffnessBound< 3 > absoluteStiffness( columnWeights, nColumns, bound );
      particle.describeStiffness( absoluteStiffness, timeNew, dT );
      break;
    }
    default:
      throw std::invalid_argument( "addBoundOfAbsoluteStiffnessTimes: a particle of dimension " +
                                   std::to_string( particle.getDimension() ) );
    }
  }

} // namespace Marmot::Discretizations::Meshfree
