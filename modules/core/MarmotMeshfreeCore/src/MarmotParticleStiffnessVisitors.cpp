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

#include "Marmot/MarmotParticleStiffnessVisitors.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

namespace Marmot::Discretizations::Meshfree {

  namespace {

    /**
     * @brief Call @p function with the number of dofs per node at compile time, as a std::integral_constant.
     * @tparam n The number tried first; the next ones up to maxDofsPerNode follow.
     * @tparam Function A callable taking a std::integral_constant< int, nDofPerNode >.
     * @param[in] nDofPerNode The number of dofs per node.
     * @param[in] visitor The name of the visitor, for the message.
     * @param[in] function The callable.
     * @throws std::invalid_argument for more than maxDofsPerNode dofs per node.
     */
    template < int n = 1, typename Function >
    void withDofsPerNodeAtCompileTime( int nDofPerNode, const char* visitor, Function&& function )
    {
      if constexpr ( n > maxDofsPerNode )
        throw std::invalid_argument( std::string( visitor ) + ": a particle with " + std::to_string( nDofPerNode ) +
                                     " dofs per node; the stiffness visitors are compiled for 1 to " +
                                     std::to_string( maxDofsPerNode ) );
      else if ( nDofPerNode == n )
        function( std::integral_constant< int, n >() );
      else
        withDofsPerNodeAtCompileTime< n + 1 >( nDofPerNode, visitor, std::forward< Function >( function ) );
    }

    /**
     * @brief Let a particle describe its stiffness to a visitor of its dimension.
     * @tparam Visitor The visitor template, e.g. AbsoluteStiffnessBound.
     * @tparam Arguments The types of the arguments of its constructor.
     * @param[in,out] particle The particle.
     * @param[in] timeNew The time of the zero increment.
     * @param[in] dT The time increment of the zero increment.
     * @param[in] arguments The arguments of the constructor of the visitor.
     * @throws std::invalid_argument for a particle of a dimension other than 2 or 3.
     */
    template < template < int > class Visitor, typename... Arguments >
    void describeStiffnessTo( Particle& particle, double timeNew, double dT, Arguments... arguments )
    {
      switch ( particle.getDimension() ) {
      case 2: {
        Visitor< 2 > visitor( arguments... );
        particle.describeStiffness( visitor, timeNew, dT );
        break;
      }
      case 3: {
        Visitor< 3 > visitor( arguments... );
        particle.describeStiffness( visitor, timeNew, dT );
        break;
      }
      default:
        throw std::invalid_argument( "the stiffness visitors: a particle of dimension " +
                                     std::to_string( particle.getDimension() ) );
      }
    }

  } // namespace

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
    withDofsPerNodeAtCompileTime( nDofPerNode, "AbsoluteStiffnessBound", [this]( auto nAtCompileTime ) {
      _addBoundOfRowsOfNode = &AbsoluteStiffnessBound::template addBoundOfRowsOfNode<
        decltype( nAtCompileTime )::value >;
    } );
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
    ( this->*_addBoundOfRowsOfNode )( A, df_dIncrement, df_dIncrementGradient, df_dIncrementSecondGradient );
  }

  template < int nDim >
  template < int nDofPerNode >
  void AbsoluteStiffnessBound< nDim >::addBoundOfRowsOfNode( int           A,
                                                             const double* df_dIncrement,
                                                             const double* df_dIncrementGradient,
                                                             const double* df_dIncrementSecondGradient )
  {
    constexpr int n = nDofPerNode;
    using Block     = Eigen::Map< const Eigen::Matrix< double, n, n > >;
    // block l of a derivative, in the layout of StiffnessVisitor
    const auto block = []( const double* derivative, int l ) { return Block( derivative + n * n * l ); };

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
  void AbsoluteStiffnessBound< nDim >::addMoments( const Eigen::Matrix< double, n, 1 >& g,
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
  auto AbsoluteStiffnessBound< nDim >::symmetricComponents( const double* d2N_dYdY ) -> Eigen::Matrix< double, nSym, 1 >
  {
    Eigen::Matrix< double, nSym, 1 > components;
    for ( int l = 0, s = 0; l < nDim; l++ )
      for ( int m = l; m < nDim; m++, s++ )
        components( s ) = d2N_dYdY[l * nDim + m];
    return components;
  }

  template < int nDim >
  void StiffnessAssembler< nDim >::visitShapeFunctions( int           nDimOfParticle,
                                                        int           nDofPerNode,
                                                        int           nNodes,
                                                        const double* N,
                                                        const double* dN_dY,
                                                        const double* d2N_dYdY )
  {
    if ( nDimOfParticle != nDim )
      throw std::invalid_argument( "StiffnessAssembler: a particle of dimension " + std::to_string( nDimOfParticle ) +
                                   " for an assembler of dimension " + std::to_string( nDim ) );
    withDofsPerNodeAtCompileTime( nDofPerNode, "StiffnessAssembler", [this]( auto nAtCompileTime ) {
      _addRowsOfNodeToK = &StiffnessAssembler::template addRowsOfNodeToK< decltype( nAtCompileTime )::value >;
    } );
    _nDof     = nDofPerNode * nNodes;
    _N        = Eigen::Map< const Eigen::MatrixXd >( N, 1, nNodes );
    _dN_dY    = Eigen::Map< const Eigen::MatrixXd >( dN_dY, nDim, nNodes );
    _d2N_dYdY = d2N_dYdY ? Eigen::MatrixXd( Eigen::Map< const Eigen::MatrixXd >( d2N_dYdY, nDim * nDim, nNodes ) )
                         : Eigen::MatrixXd();
  }

  template < int nDim >
  void StiffnessAssembler< nDim >::visitStiffnessOfNode( int           A,
                                                         const double* df_dIncrement,
                                                         const double* df_dIncrementGradient,
                                                         const double* df_dIncrementSecondGradient )
  {
    ( this->*_addRowsOfNodeToK )( A, df_dIncrement, df_dIncrementGradient, df_dIncrementSecondGradient );
  }

  template < int nDim >
  template < int nDofPerNode >
  void StiffnessAssembler< nDim >::addRowsOfNodeToK( int           A,
                                                     const double* df_dIncrement,
                                                     const double* df_dIncrementGradient,
                                                     const double* df_dIncrementSecondGradient )
  {
    using OfNode = StiffnessOfNode< nDim, nDofPerNode >;
    OfNode ofNode;
    ofNode.df_dIncrement         = Eigen::Map< const typename OfNode::Block >( df_dIncrement );
    ofNode.df_dIncrementGradient = Eigen::Map< const decltype( OfNode::df_dIncrementGradient ) >(
      df_dIncrementGradient );
    if ( df_dIncrementSecondGradient )
      ofNode.df_dIncrementSecondGradient = Eigen::Map< const decltype( OfNode::df_dIncrementSecondGradient ) >(
        df_dIncrementSecondGradient );
    else
      ofNode.df_dIncrementSecondGradient.setZero();
    addRowsOfNode( ofNode, A, _N, _dN_dY, _d2N_dYdY, Eigen::Map< Eigen::MatrixXd >( _K, _nDof, _nDof ) );
  }

  template class AbsoluteStiffnessBound< 2 >;
  template class AbsoluteStiffnessBound< 3 >;
  template class StiffnessAssembler< 2 >;
  template class StiffnessAssembler< 3 >;

  void addBoundOfAbsoluteStiffnessTimes( Particle&     particle,
                                         const double* columnWeights,
                                         int           nColumns,
                                         double*       bound,
                                         double        timeNew,
                                         double        dT )
  {
    describeStiffnessTo< AbsoluteStiffnessBound >( particle, timeNew, dT, columnWeights, nColumns, bound );
  }

  void assembleStiffness( Particle& particle, double* K, double timeNew, double dT )
  {
    describeStiffnessTo< StiffnessAssembler >( particle, timeNew, dT, K );
  }

} // namespace Marmot::Discretizations::Meshfree
