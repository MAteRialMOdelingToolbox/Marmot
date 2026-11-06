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
#include "Marmot/MarmotBondSlipLaw.h"
#include "Marmot/MarmotBondSlipLawFactory.h"
#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryElement.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace Marmot::Discretizations::Elements {

  /**
   * @class Marmot::Discretizations::Elements::EmbeddedLargeSlipBondElement
   * @brief Bond-slip of an embedded reinforcement bar with large slip: the bar may slide out of (or into) its host.
   *
   * @tparam nDim        number of spatial dimensions, 2 or 3
   * @tparam nHostNodes  number of nodes of the host element (Quad4, Quad8, Tetra4, Tetra10, Hexa8, Hexa20)
   * @tparam nBarNodes   number of nodes of the bar elements, 2 or 3 (end, end, mid)
   *
   * Unlike the small-slip EmbeddedBondElement, the bond is integrated over the **channel**: the part of the original
   * bar path inside the host element, fixed in the host. The bond history (the state of the bond-slip law) belongs to
   * the channel points, i.e., to the concrete. At every evaluation, the partner of a channel point on the bar is
   * searched -- the bar point currently next to it (closest point projection onto the current bar axis) -- among a
   * window of consecutive bar elements. The tangential slip is the reference arc length the bar has slid past the
   * channel point,
   * \f[ s_t = S_c - S^*, \f]
   * with the reference arc lengths \f$S_c\f$ of the channel point and \f$S^*\f$ of its current bar partner, and the
   * normal slip is the relative position of the bar partner w.r.t. the channel point, normal to the axis.
   *
   * Each channel point represents a tributary part of the channel (the Gauss-Lobatto weights tile the channel). Its
   * bond is scaled by \f$f = 3 x^2 - 2 x^3\f$ of the fraction \f$x \in [0, 1]\f$ of the tributary part still covered
   * by the bar, with the positions of the bar ends found from the slip of the point. As a bar end slides along the
   * channel, the bond fades smoothly (\f$C^1\f$, a kink would make Newton's method cycle); a channel point the bar has
   * left completely (\f$f = 0\f$) has no bond and its history is frozen.
   * Parts of the bar outside of all channels -- e.g., pulled out of the concrete -- are unbonded. The bonded length
   * therefore shrinks, exactly, and the free length grows as the bar is pulled out. For small slips (\f$f = 1\f$),
   * the element is identical to the small-slip EmbeddedBondElement.
   *
   * Nodes: the host element nodes, followed by the unique nodes of the window of nBarElements consecutive bar
   * elements in chain order (bar2: e0, e1, ..., eN; bar3: e0, m0, e1, m1, ..., eN).
   *
   * Element properties: [perimeter, c_s, c_e, start is a bar end, end is a bar end, (number of integration points)],
   * with the channel \f$[c_s, c_e]\f$ given in the chain parameter \f$c \in [0, n_{\mathrm{bar\,elements}}]\f$ of
   * the window (\f$c = k + (\xi + 1)/2\f$ in bar element \f$k\f$). If the partner of a channel point lies beyond an
   * end of the window that is not an end of the bar, the window is too small for the slip and an exception is
   * thrown.
   *
   * The geometry of the channel and of the frame is the reference one; the tangent is computed by central
   * differences of the internal forces.
   *
   * Material: a bond-slip law registered in the Marmot::Registry::BondSlipLawFactory.
   */
  template < int nDim, int nHostNodes, int nBarNodes >
  class EmbeddedLargeSlipBondElement : public Marmot::Discretizations::Elements::Element {

    static_assert( nDim == 2 || nDim == 3, "EmbeddedLargeSlipBondElement: nDim must be 2 or 3" );
    static_assert( nBarNodes == 2 || nBarNodes == 3, "EmbeddedLargeSlipBondElement: nBarNodes must be 2 or 3" );

  public:
    using HostGeometry = Marmot::Discretizations::Elements::GeometryElement< nDim, nHostNodes >;
    using VectorDim    = Eigen::Matrix< double, nDim, 1 >;
    using MatrixDim    = Eigen::Matrix< double, nDim, nDim >;
    using BarN         = Eigen::Matrix< double, 1, nBarNodes >;

    const int elLabel;
    const int nBarElements;
    const int nChainNodes;
    const int nNodes;
    const int nDofs;

    /// Element properties
    Eigen::Map< const Eigen::VectorXd > elementProperties;
    /// The geometry of the host element
    HostGeometry hostGeometry;
    /// Reference coordinates of the chain nodes, [dim, chain node]
    Eigen::Map< const Eigen::MatrixXd > chainCoordinates;

    /// The state of a channel point
    class QPStateVarManager : public Marmot::Utilities::StateVarVectorManager {
      /// \hideinitializer
      inline const static auto layout = makeLayout( {
        { .name = "slip", .length = 3 },
        { .name = "bond stress", .length = 3 },
        { .name = "elastic energy", .length = 1 },
        { .name = "dissipation", .length = 1 },
        { .name = "partner", .length = 1 },
        { .name = "active", .length = 1 },
        { .name = "covered fraction", .length = 1 },
        { .name = "begin of material state", .length = 0 },
      } );

    public:
      Eigen::Map< Eigen::Vector3d > slip;
      Eigen::Map< Eigen::Vector3d > bondStress;
      double&                       elasticEnergy;
      double&                       dissipation;
      double&                       partner;
      double&                       active;
      double&                       coveredFraction;
      Eigen::Map< Eigen::VectorXd > materialStateVars;

      static int getNumberOfRequiredStateVarsQuadraturePointOnly() { return layout.nRequiredStateVars; };

      QPStateVarManager( double* theStateVarVector, int nStateVars )
        : Marmot::Utilities::StateVarVectorManager( theStateVarVector, layout ),
          slip( &find( "slip" ) ),
          bondStress( &find( "bond stress" ) ),
          elasticEnergy( find( "elastic energy" ) ),
          dissipation( find( "dissipation" ) ),
          partner( find( "partner" ) ),
          active( find( "active" ) ),
          coveredFraction( find( "covered fraction" ) ),
          materialStateVars( &find( "begin of material state" ),
                             nStateVars - getNumberOfRequiredStateVarsQuadraturePointOnly() ) {};
    };

    /// A channel point
    struct QuadraturePoint {
      double                                            c   = 0; ///< chain parameter of the channel point
      double                                            S   = 0; ///< reference arc length along the chain
      double                                            pdS = 0; ///< perimeter x reference length of the point
      double                                            Sa  = 0; ///< start of the tributary part (reference arc length)
      double                                            Sb  = 0; ///< end of the tributary part (reference arc length)
      VectorDim                                         X;       ///< reference position
      VectorDim                                         hostXi;  ///< parent coordinates in the host element
      Eigen::Matrix< double, 1, nHostNodes >            hostN;   ///< host shape functions
      MatrixDim                                         R;       ///< local frame [t, n1, (n2)] as columns
      std::unique_ptr< Marmot::Materials::BondSlipLaw > bondSlipLaw;
    };

    std::vector< QuadraturePoint > qps;

  private:
    double*               stateVars_  = nullptr;
    int                   nStateVars_ = 0;
    std::vector< double > elementLengths_, cumulativeLengths_;
    double                characteristicLength_ = 1.0;

  public:
    EmbeddedLargeSlipBondElement( int elementID, int nBarElements_ )
      : elLabel( elementID ),
        nBarElements( nBarElements_ ),
        nChainNodes( nBarNodes == 2 ? nBarElements_ + 1 : 2 * nBarElements_ + 1 ),
        nNodes( nHostNodes + nChainNodes ),
        nDofs( nDim * ( nHostNodes + nChainNodes ) ),
        elementProperties( nullptr, 0 ),
        chainCoordinates( nullptr, nDim, nBarNodes == 2 ? nBarElements_ + 1 : 2 * nBarElements_ + 1 )
    {
      qps.resize( nBarNodes == 2 ? 2 : 3 );
    }

    /// the chain node indices of the nodes of bar element k, in the node order of the bar shape (end, end, mid)
    std::array< int, nBarNodes > barElementNodes( int k ) const
    {
      if constexpr ( nBarNodes == 2 )
        return { k, k + 1 };
      else
        return { 2 * k, 2 * k + 2, 2 * k + 1 };
    }

    static BarN barN( double xi )
    {
      if constexpr ( nBarNodes == 2 )
        return Marmot::Discretizations::Elements::Spatial1D::Bar2::N( xi );
      else
        return Marmot::Discretizations::Elements::Spatial1D::Bar3::N( xi );
    }

    static BarN barDNdXi( double xi )
    {
      if constexpr ( nBarNodes == 2 )
        return Marmot::Discretizations::Elements::Spatial1D::Bar2::dNdXi( xi );
      else
        return Marmot::Discretizations::Elements::Spatial1D::Bar3::dNdXi( xi );
    }

    /// nodal coordinates of bar element k, [dim, node], from the chain coordinates x (reference or current)
    Eigen::Matrix< double, nDim, nBarNodes > barElementCoordinates( const Eigen::MatrixXd& x, int k ) const
    {
      Eigen::Matrix< double, nDim, nBarNodes > xe;
      const auto                               idx = barElementNodes( k );
      for ( int i = 0; i < nBarNodes; i++ )
        xe.col( i ) = x.col( idx[i] );
      return xe;
    }

    /// reference arc length from the chain start to the point xi of bar element k
    double arcLength( int k, double xi ) const
    {
      const auto   Xe   = barElementCoordinates( chainCoordinates, k );
      const auto   rule = Marmot::Discretizations::Elements::Quadrature::Spatial1D::gaussPointList3;
      double       S    = cumulativeLengths_[k];
      const double a = -1.0, b = xi;
      for ( const auto& gp : rule ) {
        const double z = a + 0.5 * ( b - a ) * ( gp.xi( 0 ) + 1 );
        S += ( Xe * barDNdXi( z ).transpose() ).norm() * 0.5 * ( b - a ) * gp.weight;
      }
      return S;
    }

    /// reference arc length of the point xi of bar element k, extrapolated linearly beyond the element ends
    double arcLengthExtrapolated( int k, double xi ) const
    {
      const double xiClamped = std::clamp( xi, -1.0, 1.0 );
      double       S         = arcLength( k, xiClamped );
      if ( xi != xiClamped ) {
        const auto Xe = barElementCoordinates( chainCoordinates, k );
        S += ( xi - xiClamped ) * ( Xe * barDNdXi( xiClamped ).transpose() ).norm();
      }
      return S;
    }

    static std::vector< std::pair< double, double > > gaussLobatto( int n )
    {
      switch ( n ) {
      case 2: return { { -1., 1. }, { 1., 1. } };
      case 3: return { { -1., 1. / 3 }, { 0., 4. / 3 }, { 1., 1. / 3 } };
      case 4: {
        const double a = std::sqrt( 1. / 5 );
        return { { -1., 1. / 6 }, { -a, 5. / 6 }, { a, 5. / 6 }, { 1., 1. / 6 } };
      }
      case 5: {
        const double a = std::sqrt( 3. / 7 );
        return { { -1., 0.1 }, { -a, 49. / 90 }, { 0., 32. / 45 }, { a, 49. / 90 }, { 1., 0.1 } };
      }
      default:
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << "EmbeddedLargeSlipBondElement: 2 to 5 integration points are "
                                        "supported, "
                                     << n << " were requested" );
      }
    }

    static MatrixDim localFrame( const VectorDim& t )
    {
      MatrixDim R;
      if constexpr ( nDim == 2 ) {
        R.col( 0 ) = t;
        R.col( 1 ) = VectorDim( -t( 1 ), t( 0 ) );
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

    VectorDim findHostParentCoordinates( const VectorDim& X ) const
    {
      const bool tet = hostGeometry.shape == Marmot::Discretizations::Elements::ElementShapes::Tetra4 ||
                       hostGeometry.shape == Marmot::Discretizations::Elements::ElementShapes::Tetra10;
      VectorDim  xi  = tet ? VectorDim::Constant( 0.25 ) : VectorDim::Zero();
      for ( int iteration = 0;; iteration++ ) {
        const VectorDim residual = hostGeometry.NB( hostGeometry.N( xi ) ) * hostGeometry.coordinates - X;
        const VectorDim dXi      = -hostGeometry.Jacobian( hostGeometry.dNdXi( xi ) ).inverse() * residual;
        xi += dXi;
        if ( dXi.norm() < 1e-13 )
          break;
        if ( iteration > 50 )
          throw std::invalid_argument( Marmot::Utilities::MakeString()
                                       << "EmbeddedLargeSlipBondElement " << elLabel
                                       << ": inverse mapping into the host element did not converge" );
      }
      constexpr double tol = 1e-6;
      const bool inside = tet ? ( xi.minCoeff() >= -tol && xi.sum() <= 1 + tol ) : xi.cwiseAbs().maxCoeff() <= 1 + tol;
      if ( !inside )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << "EmbeddedLargeSlipBondElement " << elLabel
                                     << ": a channel point lies outside of the host element, parent "
                                        "coordinates "
                                     << xi.transpose() );
      return xi;
    }

    /// the partner of a point on the current bar axis
    struct Partner {
      int    element     = 0;
      double xi          = 0;     ///< clamped to the element
      double xiUnclamped = 0;     ///< extrapolated beyond the chain ends
      bool   beforeStart = false; ///< the projection lies before the first chain node
      bool   afterEnd    = false; ///< the projection lies behind the last chain node
    };

    Partner findPartner( const Eigen::MatrixXd& x, const VectorDim& point ) const
    {
      Partner best;
      double  bestDistance = std::numeric_limits< double >::infinity();
      double  xiFirst = 0, xiLast = 0;
      for ( int k = 0; k < nBarElements; k++ ) {
        const auto xe = barElementCoordinates( x, k );
        double     xi = 0.0;
        for ( int it = 0; it < 30; it++ ) {
          const VectorDim d   = xe * barN( xi ).transpose() - point;
          const VectorDim dx  = xe * barDNdXi( xi ).transpose();
          VectorDim       ddx = VectorDim::Zero();
          if constexpr ( nBarNodes == 3 )
            ddx = xe * Eigen::Vector3d( 1., 1., -2. );
          const double f   = d.dot( dx );
          const double df  = dx.dot( dx ) + d.dot( ddx );
          const double dxi = -f / df;
          xi += dxi;
          if ( std::abs( dxi ) < 1e-14 )
            break;
        }
        if ( k == 0 )
          xiFirst = xi;
        if ( k == nBarElements - 1 )
          xiLast = xi;
        const double    xiClamped = std::clamp( xi, -1.0, 1.0 );
        const VectorDim d         = xe * barN( xiClamped ).transpose() - point;
        const double    distance  = d.norm() + ( std::abs( xi - xiClamped ) > 1e-12 ? 1e-12 : 0.0 );
        if ( distance < bestDistance ) {
          bestDistance = distance;
          best         = { k, xiClamped, xi, false, false };
        }
      }
      best.beforeStart = best.element == 0 && xiFirst < -1.0 - 1e-12;
      best.afterEnd    = best.element == nBarElements - 1 && xiLast > 1.0 + 1e-12;
      return best;
    }

    int getNumberOfRequiredStateVars() override
    {
      return static_cast< int >( qps.size() ) * ( QPStateVarManager::getNumberOfRequiredStateVarsQuadraturePointOnly() +
                                                  qps[0].bondSlipLaw->getNumberOfRequiredStateVars() );
    }

    std::vector< std::vector< std::string > > getNodeFields() override
    {
      return std::vector< std::vector< std::string > >( nNodes, { "displacement" } );
    }

    std::vector< int > getDofIndicesPermutationPattern() override
    {
      std::vector< int > pattern( nDofs );
      for ( int i = 0; i < nDofs; i++ )
        pattern[i] = i;
      return pattern;
    }

    int getNNodes() override { return nNodes; }

    int getNSpatialDimensions() override { return nDim; }

    int getNDofPerElement() override { return nDofs; }

    /// visualized as its host element (the leading nodes), on which the channel lies
    std::string getElementShape() override { return hostGeometry.getElementShape(); }

    void assignStateVars( double* stateVars, int nStateVars ) override
    {
      stateVars_  = stateVars;
      nStateVars_ = nStateVars;
    }

    void assignProperty( const Marmot::Discretizations::ElementProperties& property ) override
    {
      if ( property.nElementProperties < 5 )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__
                                     << ": requires the properties [perimeter, c_s, c_e, start is a bar "
                                        "end, end is a bar end, (number of integration points)]" );
      new ( &elementProperties )
        Eigen::Map< const Eigen::VectorXd >( property.elementProperties, property.nElementProperties );
      const int nQps = property.nElementProperties > 5 ? static_cast< int >( std::lround( elementProperties[5] ) )
                                                       : ( nBarNodes == 2 ? 2 : 3 );
      gaussLobatto( nQps );
      qps.clear();
      qps.resize( nQps );
    }

    void assignProperty( const Marmot::Discretizations::MaterialSection& section ) override
    {
      for ( auto& qp : qps )
        qp.bondSlipLaw = std::unique_ptr< Marmot::Materials::BondSlipLaw >(
          Marmot::Registry::BondSlipLawFactory::createBondSlipLaw( section.materialName,
                                                                   section.materialProperties,
                                                                   section.nMaterialProperties,
                                                                   elLabel ) );
    }

    void assignNodeCoordinates( const double* coordinates ) override
    {
      hostGeometry.assignNodeCoordinates( coordinates );
      new ( &chainCoordinates )
        Eigen::Map< const Eigen::MatrixXd >( coordinates + nHostNodes * nDim, nDim, nChainNodes );
    }

    void initializeYourself() override
    {
      // reference lengths of the bar elements
      elementLengths_.assign( nBarElements, 0.0 );
      cumulativeLengths_.assign( nBarElements + 1, 0.0 );
      for ( int k = 0; k < nBarElements; k++ ) {
        cumulativeLengths_[k + 1] = 0.0; // arcLength reads cumulativeLengths_[k] only
        elementLengths_[k]        = arcLength( k, 1.0 ) - cumulativeLengths_[k];
        cumulativeLengths_[k + 1] = cumulativeLengths_[k] + elementLengths_[k];
      }
      characteristicLength_ = std::max( 1e-12, cumulativeLengths_.back() / nBarElements );

      const double perimeter = elementProperties[0];
      const double cStart    = elementProperties[1];
      const double cEnd      = elementProperties[2];
      const auto   rule      = gaussLobatto( static_cast< int >( qps.size() ) );

      // the reference arc lengths of the channel ends, to tile the channel into the tributary parts of the points
      auto arcLengthAt = [&]( double c ) {
        const int k = std::clamp( static_cast< int >( std::floor( c ) ), 0, nBarElements - 1 );
        return arcLength( k, 2.0 * ( c - k ) - 1.0 );
      };
      const double SStart = arcLengthAt( std::min( cStart, cEnd ) ), SEnd = arcLengthAt( std::max( cStart, cEnd ) );
      double       cumulativeWeight = 0.0;

      for ( size_t q = 0; q < qps.size(); q++ ) {
        auto& qp = qps[q];
        qp.c     = cStart + 0.5 * ( cEnd - cStart ) * ( rule[q].first + 1 );
        qp.Sa    = SStart + ( SEnd - SStart ) * cumulativeWeight / 2.0;
        cumulativeWeight += rule[q].second;
        qp.Sb = SStart + ( SEnd - SStart ) * cumulativeWeight / 2.0;
        // the channel lies in one bar element; taken from its midpoint, since a point at the end of the channel (an
        // integer c) would otherwise be attributed to the neighboring bar element, with its (possibly different) length
        const int    k = std::clamp( static_cast< int >( std::floor( 0.5 * ( cStart + cEnd ) ) ), 0, nBarElements - 1 );
        const double xi = 2.0 * ( qp.c - k ) - 1.0;

        const auto      Xe = barElementCoordinates( chainCoordinates, k );
        const VectorDim G  = Xe * barDNdXi( xi ).transpose();
        qp.X               = Xe * barN( xi ).transpose();
        qp.S               = arcLength( k, xi );
        qp.pdS             = perimeter * G.norm() * std::abs( cEnd - cStart ) * rule[q].second; // dS = |G| 2 dc
        qp.R               = localFrame( G / G.norm() );
        qp.hostXi          = findHostParentCoordinates( qp.X );
        qp.hostN           = hostGeometry.N( qp.hostXi );
      }
    }

    void setInitialConditions( StateTypes state, const double* ) override
    {
      if ( state != Marmot::Discretizations::Elements::Element::MarmotMaterialInitialization )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__ << ": invalid initial condition" );
      const int nQp = nStateVars_ / static_cast< int >( qps.size() );
      for ( size_t i = 0; i < qps.size(); i++ ) {
        QPStateVarManager sv( stateVars_ + i * nQp, nQp );
        sv.partner         = qps[i].c;
        sv.active          = 1.0;
        sv.coveredFraction = 1.0;
        qps[i].bondSlipLaw->initializeYourself( sv.materialStateVars.data(),
                                                static_cast< int >( sv.materialStateVars.size() ) );
      }
    }

    /**
     * @brief The internal forces for the total displacements Q, from the states at the beginning of the increment in
     * `states`, which are updated in place.
     */
    void evaluate( const Eigen::VectorXd& Q, Eigen::VectorXd& P, double* states, double time, double dT ) const
    {
      P.setZero( nDofs );
      const int nQp = nStateVars_ / static_cast< int >( qps.size() );

      Eigen::Map< const Eigen::Matrix< double, nDim, Eigen::Dynamic > > Uhost( Q.data(), nDim, nHostNodes );
      Eigen::Map< const Eigen::MatrixXd > Uchain( Q.data() + nHostNodes * nDim, nDim, nChainNodes );
      const Eigen::MatrixXd               x = chainCoordinates + Uchain;

      const double startIsBarEnd = elementProperties[3];
      const double endIsBarEnd   = elementProperties[4];

      for ( size_t i = 0; i < qps.size(); i++ ) {
        const auto&       qp = qps[i];
        QPStateVarManager sv( states + i * nQp, nQp );

        const VectorDim xc      = qp.X + Uhost * qp.hostN.transpose();
        const Partner   partner = findPartner( x, xc );

        if ( ( partner.beforeStart && startIsBarEnd < 0.5 ) || ( partner.afterEnd && endIsBarEnd < 0.5 ) )
          throw std::runtime_error( Marmot::Utilities::MakeString()
                                    << "EmbeddedLargeSlipBondElement " << elLabel
                                    << ": the slip exceeds the window of bar elements; increase maxSlip" );

        const auto      xe       = barElementCoordinates( x, partner.element );
        const VectorDim xb       = xe * barN( partner.xi ).transpose();
        const VectorDim relative = qp.R.transpose() * ( xb - xc );

        // tangential slip: the arc length the bar slid past the channel point (extrapolated beyond the bar ends)
        Eigen::Vector3d slip = Eigen::Vector3d::Zero();
        slip( 0 )            = qp.S - arcLengthExtrapolated( partner.element, partner.xiUnclamped );
        for ( int d = 1; d < nDim; d++ )
          slip( d ) = relative( d );

        // the fraction of the tributary part still covered by the bar, whose ends moved by the slip
        double coveredStart = qp.Sa, coveredEnd = qp.Sb;
        if ( startIsBarEnd > 0.5 )
          coveredStart = std::max( coveredStart, cumulativeLengths_.front() + slip( 0 ) );
        if ( endIsBarEnd > 0.5 )
          coveredEnd = std::min( coveredEnd, cumulativeLengths_.back() + slip( 0 ) );
        // C1 smooth fade: a kink in the bond force when a bar end leaves a tributary part makes Newton cycle
        const double covered  = std::clamp( ( coveredEnd - coveredStart ) / ( qp.Sb - qp.Sa ), 0.0, 1.0 );
        const double fraction = covered * covered * ( 3.0 - 2.0 * covered );

        sv.partner         = partner.element + 0.5 * ( partner.xi + 1 );
        sv.coveredFraction = covered;
        if ( fraction <= 0.0 ) {
          // the bar has left the channel point: no bond, the history is kept
          sv.active        = 0.0;
          sv.bondStress    = Eigen::Vector3d::Zero();
          sv.elasticEnergy = 0.0;
          continue;
        }

        Marmot::Materials::BondSlipLaw::State state{ sv.bondStress,
                                                     sv.elasticEnergy / qp.pdS,
                                                     sv.dissipation / qp.pdS,
                                                     sv.materialStateVars.data() };
        Eigen::Matrix3d                       C;
        qp.bondSlipLaw->computeBondStress( state, C, slip, slip - sv.slip, { time, dT } );

        sv.slip          = slip;
        sv.bondStress    = state.bondStress;
        sv.elasticEnergy = state.elasticEnergyDensity * qp.pdS * fraction;
        sv.dissipation   = state.dissipation * qp.pdS;
        sv.active        = 1.0;

        const VectorDim T   = qp.R * state.bondStress.template head< nDim >() * fraction;
        const BarN      M   = barN( partner.xi );
        const auto      idx = barElementNodes( partner.element );
        for ( int a = 0; a < nHostNodes; a++ )
          P.template segment< nDim >( a * nDim ) -= qp.hostN( a ) * T * qp.pdS;
        for ( int b = 0; b < nBarNodes; b++ )
          P.template segment< nDim >( ( nHostNodes + idx[b] ) * nDim ) += M( b ) * T * qp.pdS;
      }
    }

    void computeKernels( const double* QTotal_, const double*, double* Pe_, double* Ke_, double time, double dT )
      override
    {
      const Eigen::Map< const Eigen::VectorXd > QTotal( QTotal_, nDofs );
      Eigen::Map< Eigen::VectorXd >             Pe( Pe_, nDofs );
      Eigen::Map< Eigen::MatrixXd >             Ke( Ke_, nDofs, nDofs );

      // the tangent by central differences, each evaluation from the states at the beginning of the increment
      const std::vector< double > statesOld( stateVars_, stateVars_ + nStateVars_ );
      std::vector< double >       scratch( nStateVars_ );
      Eigen::VectorXd             Q = QTotal, Pplus, Pminus;
      const double                h = 1e-7 * characteristicLength_;
      for ( int j = 0; j < nDofs; j++ ) {
        Q( j )  = QTotal( j ) + h;
        scratch = statesOld;
        evaluate( Q, Pplus, scratch.data(), time, dT );
        Q( j )  = QTotal( j ) - h;
        scratch = statesOld;
        evaluate( Q, Pminus, scratch.data(), time, dT );
        Q( j ) = QTotal( j );
        Ke.col( j ) += ( Pplus - Pminus ) / ( 2 * h );
      }

      Eigen::VectorXd P;
      evaluate( QTotal, P, stateVars_, time, dT );
      Pe += P;
    }

    void computeKernelsExplicit( const double* QTotal_, const double*, double* Pe_, double time, double dT ) override
    {
      const Eigen::Map< const Eigen::VectorXd > QTotal( QTotal_, nDofs );
      Eigen::VectorXd                           P;
      evaluate( QTotal, P, stateVars_, time, dT );
      Eigen::Map< Eigen::VectorXd >( Pe_, nDofs ) += P;
    }

    void computeDistributedLoad( DistributedLoadTypes,
                                 double*,
                                 double*,
                                 int,
                                 const double*,
                                 const double*,
                                 double,
                                 double ) override
    {
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__ << ": distributed loads are not supported by bond elements" );
    }

    void computeBodyForce( double*, double*, const double*, const double*, double, double ) override {}

    void computeLumpedInertia( double* ) override {}

    void computeConsistentInertia( double* ) override {}

    void computeInternalEnergy( double& internalEnergy ) override
    {
      const int nQp = nStateVars_ / static_cast< int >( qps.size() );
      for ( size_t i = 0; i < qps.size(); i++ )
        internalEnergy += QPStateVarManager( stateVars_ + i * nQp, nQp ).elasticEnergy;
    }

    Marmot::Utilities::StateView getStateView( const std::string& stateName, int qpNumber ) override
    {
      const int         nQp  = nStateVars_ / static_cast< int >( qps.size() );
      double*           base = stateVars_ + qpNumber * nQp;
      QPStateVarManager sv( base, nQp );
      if ( sv.contains( stateName ) )
        return sv.getStateView( stateName );
      return qps[qpNumber].bondSlipLaw->getStateView( stateName, sv.materialStateVars.data() );
    }

    std::vector< double > getCoordinatesAtCenter() override
    {
      VectorDim X = VectorDim::Zero();
      for ( const auto& qp : qps )
        X += qp.X / qps.size();
      return std::vector< double >( X.data(), X.data() + nDim );
    }

    std::vector< std::vector< double > > getCoordinatesAtQuadraturePoints() override
    {
      std::vector< std::vector< double > > result;
      for ( const auto& qp : qps )
        result.emplace_back( qp.X.data(), qp.X.data() + nDim );
      return result;
    }

    int getNumberOfQuadraturePoints() override { return static_cast< int >( qps.size() ); }
  };

} // namespace Marmot::Discretizations::Elements
