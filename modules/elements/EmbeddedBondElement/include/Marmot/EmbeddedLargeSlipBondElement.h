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
#include "Marmot/EmbeddedBondKinematics.h"
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
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace Marmot::Discretizations::Elements {

  /**
   * @class Marmot::Discretizations::Elements::EmbeddedLargeSlipBondElement
   * @brief Bond-slip of an embedded bar with large slip: the bar may slide out of (or into) its host element.
   *
   * @tparam nDim        number of spatial dimensions, 2 or 3
   * @tparam nBarNodes   number of nodes of the bar elements, 2 or 3 (end, end, mid)
   * @tparam nHostNodes  number of nodes of the host element (Tri3, Quad4, Quad8, Tetra4, Hexa8, Hexa20)
   *
   * **Glossary.**
   *  - *bar*: the 1D (truss or beam) elements of one unbranched bar; \f$\eta \in [-1, 1]\f$ is the parameter of a bar
   *    element, \f$S\f$ the reference arc length along the bar.
   *  - *window*: the nBarElements consecutive bar elements this element may couple to; \f$c \in [0,
   *    n_{\mathrm{bar\,elements}}]\f$ is the window parameter, \f$c = k + (\eta + 1)/2\f$ in bar element \f$k\f$ of the
   *    window.
   *  - *channel*: the part \f$[c_s, c_e]\f$ of the original bar path inside the host element, fixed in the host. It
   *    owns the bond history: the bond-slip law states belong to the concrete, not to the bar.
   *  - *channel point*: an integration point of the channel (Gauss-Lobatto); its *tributary part* is the part of the
   *    channel its weight represents, \f$[S_a, S_b]\f$ in reference arc length.
   *  - *partner*: the point of the current bar next to a channel point (closest point projection).
   *  - *covered fraction* \f$x\f$: the part of the tributary part still covered by the bar; *fade*
   *    \f$f(x) = 3 x^2 - 2 x^3\f$.
   *
   * **Contact mechanics.** This is a point-to-segment contact formulation (Wriggers, Computational Contact Mechanics,
   * 2006, ch. 4; Laursen, Computational Contact and Impact Mechanics, 2002): channel point = slave point, window =
   * master segments, partner = closest point projection \f$\bar\eta\f$, normal slip = normal gap, tangential slip =
   * slip in the convected coordinate of the master.
   *
   * **Kinematics.** A channel point at \f$\boldsymbol{x}_c = \boldsymbol{X}_c + \sum_a N_a \boldsymbol{u}_a\f$ (host)
   * has its partner at \f$\bar\eta\f$, the solution of the closest point problem
   * \f[ g(\bar\eta) = \big(\boldsymbol{x}(\bar\eta) - \boldsymbol{x}_c\big) \cdot \boldsymbol{x}_{,\eta}(\bar\eta) = 0
   * \f]
   * on the current bar \f$\boldsymbol{x}(\eta) = \sum_b M_b(\eta)\, \boldsymbol{x}_b\f$, extended linearly beyond the
   * window ends. With the gap \f$\boldsymbol{d} = \boldsymbol{x}(\bar\eta) - \boldsymbol{x}_c\f$ (bar minus host, as in
   * the small-slip element) and the current frame \f$[\boldsymbol{t}, \boldsymbol{n}_1, \boldsymbol{n}_2]\f$ at the
   * partner (\f$\boldsymbol{t} = \boldsymbol{x}_{,\eta}/|\boldsymbol{x}_{,\eta}|\f$, the normals carried along from the
   * reference frame of the channel point, EmbeddedBondKinematics::transportFrame), the slip is
   * \f[ s_t = S_c - S(\bar\eta), \qquad s_{n,i} = \boldsymbol{n}_i \cdot \boldsymbol{d}, \f]
   * the reference arc length the bar has slid past the channel point, and the normal gap. Both are invariant under a
   * rigid motion of host and bar (objective); the normal frame follows the bending of the bar but not a twist about its
   * axis, which the laws, isotropic in the normal plane, do not need.
   *
   * **Bond and fade.** The tributary part of a channel point is covered by the bar between its ends, which have moved
   * by \f$s_t\f$: \f$x = \big(\min(S_b, S_{\mathrm{end}} + s_t) - \max(S_a, S_{\mathrm{start}} + s_t)\big) / (S_b -
   * S_a)\f$, clamped to [0, 1] (only at bar ends). The bond is scaled by \f$f(x)\f$: \f$C^1\f$, as a kink in the bond
   * force when a bar end leaves a tributary part makes Newton's method cycle. A channel point the bar has left
   * completely (\f$f = 0\f$) has no bond and its history is frozen until the bar comes back. Parts of the bar outside
   * of all channels -- e.g., pulled out of the concrete -- are unbonded, so the bonded length shrinks with the
   * pull-out, up to the smoothing of \f$f\f$ over the tributary part at the bar end (a bar slid by 1 into a channel of
   * length 3 with 3 points is bonded over \f$0.5 f(1) + 2 f(0.75) = 2.1875\f$ instead of 2). For small displacements
   * and a covered channel, the element agrees with the small-slip EmbeddedBondElement to first order.
   *
   * **Residual and tangent.** With \f$p\,\Delta S\f$ the perimeter times the tributary length,
   * \f[ \boldsymbol{f} = \sum_{\mathrm{points}} f\, p\, \Delta S\; \boldsymbol{B}^T \boldsymbol{\tau}, \qquad
   *     \boldsymbol{B} = \frac{\partial \boldsymbol{s}}{\partial \boldsymbol{u}}, \f]
   * over the host dofs and the dofs of the partner's bar element (the work conjugate of the slip). The linearization
   * of the closest point gives
   * \f[ \delta\bar\eta = -\frac{(\delta\boldsymbol{x}_b - \delta\boldsymbol{x}_c) \cdot \boldsymbol{x}_{,\eta} +
   *     \boldsymbol{d} \cdot \delta\boldsymbol{x}_{,\eta}}{\boldsymbol{x}_{,\eta} \cdot \boldsymbol{x}_{,\eta} +
   *     \boldsymbol{d} \cdot \boldsymbol{x}_{,\eta\eta}}, \qquad \delta s_t = -S_{,\eta}\, \delta\bar\eta, \qquad
   *     \delta s_{n,i} = \boldsymbol{n}_i \cdot \delta\boldsymbol{d} + \boldsymbol{d} \cdot \delta\boldsymbol{n}_i \f]
   * (\f$\delta\boldsymbol{x}_b, \delta\boldsymbol{x}_{,\eta}\f$ at fixed \f$\eta\f$), coded analytically in
   * localKinematics(). The tangent
   * \f[ \boldsymbol{K} = p\,\Delta S \Big( f\, \boldsymbol{B}^T \frac{\partial\boldsymbol{\tau}}{\partial
   *     \boldsymbol{s}} \boldsymbol{B} + \frac{\mathrm{d}f}{\mathrm{d}x} \frac{\mathrm{d}x}{\mathrm{d}s_t}
   *     \boldsymbol{B}^T \boldsymbol{\tau} \otimes \boldsymbol{B}_t + f\, \frac{\partial (\boldsymbol{B}^T
   *     \boldsymbol{\tau})}{\partial \boldsymbol{u}} \Big|_{\boldsymbol{\tau}} \Big) \f]
   * takes the material part from the algorithmic tangent of the bond-slip law, the fade part analytically, and the
   * geometric part (the second derivatives of the closest point projection and of the frame, times the bond stress)
   * by central differences of the analytic \f$\boldsymbol{B}\f$ over the host and partner dofs only (geometry, no law
   * evaluation): geometricStiffness().
   *
   * **Bars.** The bar elements should form a straight or smoothly curved bar: at a kink, the partner of a point on the
   * convex side of the kink is the kink node itself (its projection lies on neither element), and the slip jumps by
   * about \f$2 g \tan(\alpha/2)\f$ for a gap g and a kink angle \f$\alpha\f$ as the point crosses the kink.
   *
   * **Nodes**: the host element nodes, followed by the unique nodes of the window in bar order (bar2: e0, e1, ...,
   * eN; bar3: e0, m0, e1, m1, ..., eN). Unlike the small-slip element (bar nodes first), the host comes first, so that
   * the element is drawn as its host element (the leading nodes), on which the channel lies.
   *
   * **Element properties**: [perimeter, c_s, c_e, start is a bar end, end is a bar end, (number of channel points
   * 2..5)], with \f$0 \leq c_s < c_e \leq n_{\mathrm{bar\,elements}}\f$ within one bar element, and the flags 0 or 1.
   * If the partner of a channel point lies beyond an end of the window that is not an end of the bar, the window is
   * too small for the slip and an exception is thrown.
   *
   * Material: a bond-slip law registered in the Marmot::Registry::BondSlipLawFactory.
   */
  template < int nDim, int nBarNodes, int nHostNodes >
  class EmbeddedLargeSlipBondElement : public Marmot::Discretizations::Elements::Element {

    static_assert( nDim == 2 || nDim == 3, "EmbeddedLargeSlipBondElement: nDim must be 2 or 3" );
    static_assert( nBarNodes == 2 || nBarNodes == 3, "EmbeddedLargeSlipBondElement: nBarNodes must be 2 or 3" );

  public:
    using HostGeometry   = Marmot::Discretizations::Elements::GeometryElement< nDim, nHostNodes >; ///< host shape
    using VectorDim      = Eigen::Matrix< double, nDim, 1 >;                                       ///< spatial vector
    using MatrixDim      = Eigen::Matrix< double, nDim, nDim >;                                    ///< spatial tensor
    using BarN           = Eigen::Matrix< double, 1, nBarNodes >;    ///< shape functions of a bar element
    using BarCoordinates = Eigen::Matrix< double, nDim, nBarNodes >; ///< nodal coordinates of a bar element

    /// the dofs a channel point couples: the host element, then the partner's bar element
    static constexpr int nLocalDofs = ( nHostNodes + nBarNodes ) * nDim;
    using LocalVector               = Eigen::Matrix< double, nLocalDofs, 1 >;          ///< over the local dofs
    using LocalMatrix               = Eigen::Matrix< double, nLocalDofs, nLocalDofs >; ///< over the local dofs
    using SlipJacobian              = Eigen::Matrix< double, nDim, nLocalDofs >;       ///< d slip / d local dofs

    const int elLabel;                                                                 ///< element label
    const int nBarElements; ///< number of bar elements in the window
    const int nWindowNodes; ///< number of unique nodes of the window
    const int nNodes;       ///< host nodes + window nodes
    const int nDofs;        ///< number of element dofs

    /// The element properties as given
    Eigen::Map< const Eigen::VectorXd > elementProperties;
    double                              perimeter     = 0;     ///< perimeter of the bar
    double                              cStart        = 0;     ///< start of the channel, window parameter
    double                              cEnd          = 0;     ///< end of the channel, window parameter
    bool                                startIsBarEnd = false; ///< the window starts at an end of the bar
    bool                                endIsBarEnd   = false; ///< the window ends at an end of the bar

    /// The geometry of the host element
    HostGeometry hostGeometry;
    /// Reference coordinates of the window nodes, [dim, window node]
    Eigen::Map< const Eigen::MatrixXd > windowCoordinates;

    /// The state of a channel point
    class QPStateVarManager : public Marmot::Utilities::StateVarVectorManager {
      /// \hideinitializer
      inline const static auto layout = makeLayout( {
        { .name = "slip", .length = 3 },
        { .name = "bond stress", .length = 3 },
        { .name = "elastic energy density", .length = 1 },
        { .name = "dissipation density", .length = 1 },
        { .name = "dissipation", .length = 1 },
        { .name = "partner", .length = 1 },
        { .name = "active", .length = 1 },
        { .name = "covered fraction", .length = 1 },
        { .name = "begin of material state", .length = 0 },
      } );

    public:
      Eigen::Map< Eigen::Vector3d > slip;                 ///< [s_t, s_n1, s_n2]
      Eigen::Map< Eigen::Vector3d > bondStress;           ///< bond stress of the law (not faded)
      double&                       elasticEnergyDensity; ///< of the law, per bond surface (not faded)
      double&                       dissipationDensity;   ///< of the law, per bond surface (not faded)
      double&                       dissipation;          ///< dissipated energy of the point, faded increments
      double&                       partner;              ///< window parameter of the partner
      double&                       active;               ///< 1 if the bar covers the point (f > 0), else 0
      double&                       coveredFraction;      ///< x
      Eigen::Map< Eigen::VectorXd > materialStateVars;    ///< the state of the bond-slip law

      /// number of state variables without the bond-slip law
      static int getNumberOfRequiredStateVarsQuadraturePointOnly() { return layout.nRequiredStateVars; };

      /// a view of the state of a channel point
      QPStateVarManager( double* theStateVarVector, int nStateVars )
        : Marmot::Utilities::StateVarVectorManager( theStateVarVector, layout ),
          slip( &find( "slip" ) ),
          bondStress( &find( "bond stress" ) ),
          elasticEnergyDensity( find( "elastic energy density" ) ),
          dissipationDensity( find( "dissipation density" ) ),
          dissipation( find( "dissipation" ) ),
          partner( find( "partner" ) ),
          active( find( "active" ) ),
          coveredFraction( find( "covered fraction" ) ),
          materialStateVars( &find( "begin of material state" ),
                             nStateVars - getNumberOfRequiredStateVarsQuadraturePointOnly() ){};
    };

    /// A channel point
    struct ChannelPoint {
      double                                            c                    = 0; ///< window parameter
      double                                            S                    = 0; ///< reference arc length
      double                                            perimeterTimesLength = 0; ///< p x tributary length
      double                                            tributaryStart       = 0; ///< reference arc length
      double                                            tributaryEnd         = 0; ///< reference arc length
      VectorDim                                         X;                        ///< reference position
      VectorDim                                         hostXi;      ///< parent coordinates in the host element
      Eigen::Matrix< double, 1, nHostNodes >            hostN;       ///< host shape functions
      MatrixDim                                         R;           ///< reference frame [t, n1, (n2)] as columns
      std::unique_ptr< Marmot::Materials::BondSlipLaw > bondSlipLaw; ///< the bond-slip law of the point
    };

    /// The channel points
    std::vector< ChannelPoint > channelPoints;

    /// The tolerance of "inside the host element" in parent coordinates (as the small-slip element)
    static constexpr double insideTolerance = 1e-6;

  private:
    double*               stateVars_  = nullptr;
    int                   nStateVars_ = 0;
    std::vector< double > cumulativeLengths_; ///< reference arc length at the start of each bar element, and the end
    double                characteristicLength_ = 1.0;

  public:
    /**
     * @brief Construct the element.
     * @param elementID     element label
     * @param nBarElements_ number of bar elements in the window: a runtime size, so that the registrations of all
     * window sizes share one class
     */
    EmbeddedLargeSlipBondElement( int elementID, int nBarElements_ )
      : elLabel( elementID ),
        nBarElements( nBarElements_ ),
        nWindowNodes( nBarNodes == 2 ? nBarElements_ + 1 : 2 * nBarElements_ + 1 ),
        nNodes( nHostNodes + nWindowNodes ),
        nDofs( nDim * ( nHostNodes + nWindowNodes ) ),
        elementProperties( nullptr, 0 ),
        windowCoordinates( nullptr, nDim, nBarNodes == 2 ? nBarElements_ + 1 : 2 * nBarElements_ + 1 )
    {
      channelPoints.resize( defaultNumberOfChannelPoints() );
    }

    /// 2 channel points for linear bars, 3 for quadratic bars
    static constexpr int defaultNumberOfChannelPoints() { return nBarNodes == 2 ? 2 : 3; }

    /// the window node indices of the nodes of bar element k, in the node order of the bar shape (end, end, mid)
    std::array< int, nBarNodes > barElementNodes( int k ) const
    {
      if constexpr ( nBarNodes == 2 )
        return { k, k + 1 };
      else
        return { 2 * k, 2 * k + 2, 2 * k + 1 };
    }

    /// shape functions of a bar element
    static BarN barN( double eta )
    {
      if constexpr ( nBarNodes == 2 )
        return Spatial1D::Bar2::N( eta );
      else
        return Spatial1D::Bar3::N( eta );
    }

    /// first derivatives of the shape functions of a bar element
    static BarN barDNdEta( double eta )
    {
      if constexpr ( nBarNodes == 2 )
        return Spatial1D::Bar2::dNdXi( eta );
      else
        return Spatial1D::Bar3::dNdXi( eta );
    }

    /// second derivatives of the shape functions of a bar element
    static BarN barD2NdEta2( double eta )
    {
      if constexpr ( nBarNodes == 2 )
        return Spatial1D::Bar2::d2NdXi2( eta );
      else
        return Spatial1D::Bar3::d2NdXi2( eta );
    }

    /// nodal coordinates of bar element k, [dim, node], from the window coordinates x (reference or current)
    BarCoordinates barElementCoordinates( const Eigen::MatrixXd& x, int k ) const
    {
      BarCoordinates xe;
      const auto     idx = barElementNodes( k );
      for ( int i = 0; i < nBarNodes; i++ )
        xe.col( i ) = x.col( idx[i] );
      return xe;
    }

    /// reference arc length from the window start to the point eta of bar element k, eta in [-1, 1]
    double arcLength( int k, double eta ) const
    {
      const auto& rule = Quadrature::getQuadraturePoints( ReferenceDomain::Line, Quadrature::gaussExactToDegree( 5 ) );
      const BarCoordinates Xe = barElementCoordinates( windowCoordinates, k );
      double               S  = cumulativeLengths_[k];
      for ( const auto& gp : rule ) {
        const double z = -1.0 + 0.5 * ( eta + 1.0 ) * ( gp.xi( 0 ) + 1 );
        S += ( Xe * barDNdEta( z ).transpose() ).norm() * 0.5 * ( eta + 1.0 ) * gp.weight;
      }
      return S;
    }

    /// A point of bar element k of the window: shape functions and reference arc length
    struct BarPoint {
      BarN   N, dN_dEta, d2N_dEta2; ///< shape functions and derivatives
      double S       = 0;           ///< reference arc length
      double dS_dEta = 0;           ///< derivative of the reference arc length
    };

    /// the point eta of bar element k; beyond the first or last node of the window, on its linear extension
    BarPoint barPoint( int k, double eta ) const
    {
      BarPoint             p;
      const BarCoordinates Xe     = barElementCoordinates( windowCoordinates, k );
      const bool           before = k == 0 && eta < -1.0;
      const bool           after  = k == nBarElements - 1 && eta > 1.0;
      if ( before || after ) {
        const double end = before ? -1.0 : 1.0;
        p.dN_dEta        = barDNdEta( end );
        p.N              = barN( end ) + ( eta - end ) * p.dN_dEta;
        p.d2N_dEta2      = BarN::Zero();
        p.dS_dEta        = ( Xe * p.dN_dEta.transpose() ).norm();
        p.S              = ( before ? cumulativeLengths_[k] : cumulativeLengths_[k + 1] ) + ( eta - end ) * p.dS_dEta;
      }
      else {
        p.N         = barN( eta );
        p.dN_dEta   = barDNdEta( eta );
        p.d2N_dEta2 = barD2NdEta2( eta );
        p.dS_dEta   = ( Xe * p.dN_dEta.transpose() ).norm();
        p.S         = arcLength( k, eta );
      }
      return p;
    }

    /**
     * @brief The closest point on bar element k (extended linearly at the window ends) to the point xc: Newton on the
     * orthogonality \f$g(\eta) = (\boldsymbol{x}(\eta) - \boldsymbol{x}_c) \cdot \boldsymbol{x}_{,\eta}\f$, \f$g' =
     * \boldsymbol{x}_{,\eta} \cdot \boldsymbol{x}_{,\eta} + (\boldsymbol{x} - \boldsymbol{x}_c) \cdot
     * \boldsymbol{x}_{,\eta\eta}\f$.
     * @param k        bar element of the window
     * @param xe       current nodal coordinates of the bar element
     * @param xc       the point
     * @param etaStart start of the iteration
     * @return the parameter eta of the closest point, possibly outside of [-1, 1]
     * @throws std::runtime_error if g' <= 0 (not a minimum: the bar is curved more strongly than its distance to the
     * point) or without convergence
     */
    double closestPoint( int k, const BarCoordinates& xeGlobal, const VectorDim& xcGlobal, double etaStart ) const
    {
      // relative to the first node: the round-off of the gap is then about eps (|gap| + L), not eps |x|
      const VectorDim      origin = xeGlobal.col( 0 );
      const BarCoordinates xe     = xeGlobal.colwise() - origin;
      const VectorDim      xc     = xcGlobal - origin;
      double               eta    = etaStart;
      for ( int iteration = 0; iteration < 50; iteration++ ) {
        const BarPoint  p                   = barPoint( k, eta );
        const VectorDim gap                 = xe * p.N.transpose() - xc;
        const VectorDim dx_dEta             = xe * p.dN_dEta.transpose();
        const VectorDim d2x_dEta2           = xe * p.d2N_dEta2.transpose();
        const double    orthogonality       = gap.dot( dx_dEta );
        const double    dOrthogonality_dEta = dx_dEta.dot( dx_dEta ) + gap.dot( d2x_dEta2 );
        if ( !( dOrthogonality_dEta > 0 ) )
          throw std::runtime_error( Marmot::Utilities::MakeString()
                                    << "EmbeddedLargeSlipBondElement " << elLabel
                                    << ": no closest point on the bar (it is curved more strongly than its distance to "
                                       "the channel point)" );
        const double step = -orthogonality / dOrthogonality_dEta;
        eta += step;
        // relative to |eta|: beyond the element, x(eta) is the difference of large terms of the extrapolated shape
        // functions
        if ( std::abs( step ) < 1e-12 * std::max( 1.0, std::abs( eta ) ) )
          return eta;
      }
      throw std::runtime_error( Marmot::Utilities::MakeString()
                                << "EmbeddedLargeSlipBondElement " << elLabel
                                << ": the closest point projection onto the bar did not converge" );
    }

    /// Tolerance in the bar parameter for a partner to belong to a bar element (at a node shared with the next one)
    static constexpr double etaTolerance = 1e-9;

    /// the partner of a channel point on the current bar
    struct Partner {
      int    element = 0;     ///< bar element of the window
      double eta     = 0;     ///< parameter in the bar element; outside of [-1, 1] only at the window ends
      bool   atKink  = false; ///< the partner is a node at a kink (eta fixed at -1 or 1)
    };

    /**
     * @brief The partner of the current channel point position xc: the closest point on the window, searched from the
     * bar element of the previous partner towards the side the projection points to.
     * @param x            current window coordinates
     * @param xc           current position of the channel point
     * @param startElement bar element of the previous partner
     */
    Partner findPartner( const Eigen::MatrixXd& x, const VectorDim& xc, int startElement ) const
    {
      int k             = std::clamp( startElement, 0, nBarElements - 1 );
      int lastDirection = 0;
      for ( int visited = 0; visited <= nBarElements; visited++ ) {
        const BarCoordinates xe = barElementCoordinates( x, k );
        // the projection onto the chord of the element decides the direction; the closest point iteration runs only
        // near the element, where its (extrapolated) polynomial is well conditioned
        const VectorDim chord    = xe.col( 1 ) - xe.col( 0 );
        const double    etaChord = 2.0 * ( xc - xe.col( 0 ) ).dot( chord ) / chord.squaredNorm() - 1.0;
        const bool      far      = std::abs( etaChord ) > 1.5 &&
                         ( ( etaChord < 0 && k > 0 ) || ( etaChord > 0 && k < nBarElements - 1 ) );
        const double eta = far ? etaChord : closestPoint( k, xe, xc, std::clamp( etaChord, -1.0, 1.0 ) );
        // a partner at a node shared by two collinear elements projects to eta = +-1 up to round-off on both: only a
        // clear excess moves to the neighbor, so that round-off is not taken for a kink
        const int direction = eta < -1.0 - etaTolerance && k > 0
                                ? -1
                                : ( eta > 1.0 + etaTolerance && k < nBarElements - 1 ? 1 : 0 );
        if ( direction == 0 )
          return { k, eta, false };
        // turning back: the point lies on the convex side of a kink, its partner is the kink node
        if ( direction == -lastDirection )
          return { k, direction < 0 ? -1.0 : 1.0, true };
        lastDirection = direction;
        k += direction;
      }
      throw std::logic_error( "EmbeddedLargeSlipBondElement: the partner search did not terminate" );
    }

    /// The slip of a channel point and its derivatives w.r.t. the local dofs (host, then the partner's bar element)
    struct LocalKinematics {
      Eigen::Vector3d slip = Eigen::Vector3d::Zero(); ///< [s_t, s_n1, (s_n2)]
      SlipJacobian    dSlip_dU;                       ///< B
      double          eta = 0;                        ///< partner parameter
    };

    /**
     * @brief The slip of a channel point for the local displacements uLocal (host nodes, then the nodes of bar element
     * k) and its first derivatives B, analytically (see the class documentation).
     * @param cp       the channel point
     * @param k        bar element of the partner
     * @param uLocal   local displacements
     * @param etaStart start of the closest point iteration, or the fixed partner at a kink node
     * @param etaFixed the partner is a kink node
     */
    LocalKinematics localKinematics( const ChannelPoint& cp,
                                     int                 k,
                                     const LocalVector&  uLocal,
                                     double              etaStart,
                                     bool                etaFixed ) const
    {
      LocalKinematics kin;

      const Eigen::Map< const Eigen::Matrix< double, nDim, nHostNodes > > uHost( uLocal.data() );
      const Eigen::Map< const BarCoordinates >                            uBar( uLocal.data() + nHostNodes * nDim );
      // positions relative to the first node of the bar element (round-off far from the origin, see closestPoint)
      const BarCoordinates xeGlobal = barElementCoordinates( windowCoordinates, k ) + uBar;
      const VectorDim      origin   = xeGlobal.col( 0 );
      const VectorDim      xc       = cp.X - origin + uHost * cp.hostN.transpose();
      const BarCoordinates xe       = xeGlobal.colwise() - origin;

      kin.eta                   = etaFixed ? etaStart : closestPoint( k, xe, xc, etaStart );
      const BarPoint  p         = barPoint( k, kin.eta );
      const VectorDim gap       = xe * p.N.transpose() - xc;
      const VectorDim dx_dEta   = xe * p.dN_dEta.transpose();
      const VectorDim d2x_dEta2 = xe * p.d2N_dEta2.transpose();
      const double    length    = dx_dEta.norm();
      const VectorDim t         = dx_dEta / length;
      const auto      frame     = EmbeddedBondKinematics::transportFrame< nDim >( cp.R, t );

      kin.slip( 0 ) = cp.S - p.S;
      for ( int i = 1; i < nDim; i++ )
        kin.slip( i ) = frame.R.col( i ).dot( gap );

      // derivatives at fixed eta: of the channel point, of the partner and of the bar tangent dx/deta
      using DimByLocal  = Eigen::Matrix< double, nDim, nLocalDofs >;
      DimByLocal dXc_dU = DimByLocal::Zero(), dXb_dU = DimByLocal::Zero(), dDx_dU = DimByLocal::Zero();
      for ( int a = 0; a < nHostNodes; a++ )
        dXc_dU.template block< nDim, nDim >( 0, a * nDim ) = cp.hostN( a ) * MatrixDim::Identity();
      for ( int b = 0; b < nBarNodes; b++ ) {
        dXb_dU.template block< nDim, nDim >( 0, ( nHostNodes + b ) * nDim ) = p.N( b ) * MatrixDim::Identity();
        dDx_dU.template block< nDim, nDim >( 0, ( nHostNodes + b ) * nDim ) = p.dN_dEta( b ) * MatrixDim::Identity();
      }

      // linearization of the closest point
      Eigen::Matrix< double, 1, nLocalDofs > dEta_dU = Eigen::Matrix< double, 1, nLocalDofs >::Zero();
      if ( !etaFixed )
        dEta_dU = -( dx_dEta.transpose() * ( dXb_dU - dXc_dU ) + gap.transpose() * dDx_dU ) /
                  ( dx_dEta.dot( dx_dEta ) + gap.dot( d2x_dEta2 ) );

      const DimByLocal dGap_dU     = dXb_dU - dXc_dU + dx_dEta * dEta_dU;
      const DimByLocal dTangent_dU = ( MatrixDim::Identity() - t * t.transpose() ) / length *
                                     ( dDx_dU + d2x_dEta2 * dEta_dU );

      kin.dSlip_dU.row( 0 ) = -p.dS_dEta * dEta_dU;
      for ( int i = 1; i < nDim; i++ )
        kin.dSlip_dU.row( i ) = frame.R.col( i ).transpose() * dGap_dU +
                                gap.transpose() * frame.dNormal_dTangent[i - 1] * dTangent_dU;
      return kin;
    }

    int getNumberOfRequiredStateVars() override
    {
      if ( !channelPoints[0].bondSlipLaw )
        throw std::logic_error( Marmot::Utilities::MakeString()
                                << "EmbeddedLargeSlipBondElement " << elLabel
                                << ": the bond-slip law must be assigned before the state variables" );
      return static_cast< int >( channelPoints.size() ) *
             ( QPStateVarManager::getNumberOfRequiredStateVarsQuadraturePointOnly() +
               channelPoints[0].bondSlipLaw->getNumberOfRequiredStateVars() );
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
                                     << "EmbeddedLargeSlipBondElement " << elLabel
                                     << ": requires the properties [perimeter, c_s, c_e, start is a bar end, end is a "
                                        "bar end, (number of channel points)]" );
      new ( &elementProperties )
        Eigen::Map< const Eigen::VectorXd >( property.elementProperties, property.nElementProperties );

      auto flag = [&]( int index ) {
        const double value = elementProperties[index];
        if ( value != 0.0 && value != 1.0 )
          throw std::invalid_argument( Marmot::Utilities::MakeString()
                                       << "EmbeddedLargeSlipBondElement " << elLabel << ": property " << index
                                       << " is a flag, 0 or 1, got " << value );
        return value == 1.0;
      };
      perimeter     = elementProperties[0];
      cStart        = elementProperties[1];
      cEnd          = elementProperties[2];
      startIsBarEnd = flag( 3 );
      endIsBarEnd   = flag( 4 );
      if ( !( perimeter > 0 ) )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << "EmbeddedLargeSlipBondElement " << elLabel << ": requires perimeter > 0" );
      if ( !( 0 <= cStart && cStart < cEnd && cEnd <= nBarElements ) )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << "EmbeddedLargeSlipBondElement " << elLabel
                                     << ": requires 0 <= c_s < c_e <= number of bar elements, got [" << cStart << ", "
                                     << cEnd << "]" );

      const int nPoints = property.nElementProperties > 5 ? static_cast< int >( std::lround( elementProperties[5] ) )
                                                          : defaultNumberOfChannelPoints();
      EmbeddedBondKinematics::gaussLobatto( nPoints ); // validates
      if ( channelPoints[0].bondSlipLaw && nPoints != static_cast< int >( channelPoints.size() ) )
        throw std::logic_error( Marmot::Utilities::MakeString()
                                << "EmbeddedLargeSlipBondElement " << elLabel
                                << ": the element properties (number of channel points) must be assigned before the "
                                   "bond-slip law" );
      channelPoints.resize( nPoints );
    }

    void assignProperty( const Marmot::Discretizations::MaterialSection& section ) override
    {
      for ( auto& cp : channelPoints )
        cp.bondSlipLaw = std::unique_ptr< Marmot::Materials::BondSlipLaw >(
          Marmot::Registry::BondSlipLawFactory::createBondSlipLaw( section.materialName,
                                                                   section.materialProperties,
                                                                   section.nMaterialProperties,
                                                                   elLabel ) );
    }

    void assignNodeCoordinates( const double* coordinates ) override
    {
      hostGeometry.assignNodeCoordinates( coordinates );
      new ( &windowCoordinates )
        Eigen::Map< const Eigen::MatrixXd >( coordinates + nHostNodes * nDim, nDim, nWindowNodes );
    }

    void initializeYourself() override
    {
      // reference arc lengths at the starts of the bar elements
      cumulativeLengths_.assign( nBarElements + 1, 0.0 );
      for ( int k = 0; k < nBarElements; k++ )
        cumulativeLengths_[k + 1] = arcLength( k, 1.0 ); // reads cumulativeLengths_[k], set before
      characteristicLength_ = std::max( 1e-12, cumulativeLengths_.back() / nBarElements );

      // the channel lies in one bar element
      const int k = std::min( static_cast< int >( std::floor( cStart ) ), nBarElements - 1 );
      if ( cEnd > k + 1 + 1e-12 )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << "EmbeddedLargeSlipBondElement " << elLabel
                                     << ": the channel [c_s, c_e] must lie in one bar element" );

      const auto&          rule   = EmbeddedBondKinematics::gaussLobatto( static_cast< int >( channelPoints.size() ) );
      const BarCoordinates Xe     = barElementCoordinates( windowCoordinates, k );
      const double         SStart = arcLength( k, 2.0 * ( cStart - k ) - 1.0 );
      const double         SEnd   = arcLength( k, 2.0 * ( cEnd - k ) - 1.0 );
      double               cumulativeWeight = 0.0;

      for ( size_t q = 0; q < channelPoints.size(); q++ ) {
        auto& cp = channelPoints[q];
        cp.c     = cStart + 0.5 * ( cEnd - cStart ) * ( rule[q].xi( 0 ) + 1 );
        // the weights (sum 2) tile the channel into the tributary parts of the points
        cp.tributaryStart = SStart + ( SEnd - SStart ) * cumulativeWeight / 2.0;
        cumulativeWeight += rule[q].weight;
        cp.tributaryEnd = SStart + ( SEnd - SStart ) * cumulativeWeight / 2.0;

        const double    eta     = 2.0 * ( cp.c - k ) - 1.0;
        const VectorDim dX_dEta = Xe * barDNdEta( eta ).transpose();
        cp.X                    = Xe * barN( eta ).transpose();
        cp.S                    = arcLength( k, eta );
        // dS = |dX/deta| deta, deta = (c_e - c_s) dxi on the [-1, 1] of the rule
        cp.perimeterTimesLength = perimeter * dX_dEta.norm() * ( cEnd - cStart ) * rule[q].weight;
        cp.R                    = EmbeddedBondKinematics::localFrame< nDim >( dX_dEta / dX_dEta.norm() );
        cp.hostXi               = hostGeometry.parentCoordinates( cp.X );
        if ( !HostGeometry::isInsideParentDomain( cp.hostXi, insideTolerance ) )
          throw std::invalid_argument( Marmot::Utilities::MakeString()
                                       << "EmbeddedLargeSlipBondElement " << elLabel
                                       << ": a channel point lies outside of the host element, parent coordinates "
                                       << cp.hostXi.transpose() );
        cp.hostN = hostGeometry.N( cp.hostXi );
      }
    }

    void setInitialConditions( StateTypes state, const double* ) override
    {
      if ( state != Marmot::Discretizations::Elements::Element::MarmotMaterialInitialization )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__ << ": invalid initial condition" );
      const int nQp = nStateVars_ / static_cast< int >( channelPoints.size() );
      for ( size_t i = 0; i < channelPoints.size(); i++ ) {
        QPStateVarManager sv( stateVars_ + i * nQp, nQp );
        sv.partner         = channelPoints[i].c;
        sv.active          = 1.0;
        sv.coveredFraction = 1.0;
        channelPoints[i].bondSlipLaw->initializeYourself( sv.materialStateVars.data(),
                                                          static_cast< int >( sv.materialStateVars.size() ) );
      }
    }

    /// the fade f(x) = 3 x^2 - 2 x^3 of the covered fraction x
    static double fade( double coveredFraction )
    {
      return coveredFraction * coveredFraction * ( 3.0 - 2.0 * coveredFraction );
    }

    /**
     * @brief The internal forces and (optionally) the tangent for the total displacements Q; the states hold the
     * values at the beginning of the increment and are updated in place.
     * @param Q      total element displacements
     * @param P      internal forces (added to)
     * @param K      tangent (added to), or nullptr
     * @param states the state variables of all channel points
     * @param time   time at the beginning of the increment
     * @param dT     time increment
     */
    void evaluate( const Eigen::Ref< const Eigen::VectorXd >& Q,
                   Eigen::Ref< Eigen::VectorXd >              P,
                   Eigen::Ref< Eigen::MatrixXd >*             K,
                   double*                                    states,
                   double                                     time,
                   double                                     dT ) const
    {
      const int nQp = nStateVars_ / static_cast< int >( channelPoints.size() );

      const Eigen::Map< const Eigen::Matrix< double, nDim, nHostNodes > > uHost( Q.data() );
      const Eigen::Map< const Eigen::MatrixXd > uWindow( Q.data() + nHostNodes * nDim, nDim, nWindowNodes );
      const Eigen::MatrixXd                     x = windowCoordinates + uWindow;

      for ( size_t i = 0; i < channelPoints.size(); i++ ) {
        const auto&       cp = channelPoints[i];
        QPStateVarManager sv( states + i * nQp, nQp );

        const VectorDim xc = cp.X + uHost * cp.hostN.transpose();
        // the search starts from the previous partner; 0: not searched yet (the states need no initialization), from
        // the element of the channel point
        const double  previousPartner = sv.partner > 0 ? sv.partner : cp.c;
        const Partner partner         = findPartner( x, xc, static_cast< int >( std::floor( previousPartner ) ) );
        if ( ( partner.element == 0 && partner.eta < -1.0 && !startIsBarEnd ) ||
             ( partner.element == nBarElements - 1 && partner.eta > 1.0 && !endIsBarEnd ) )
          throw std::runtime_error( Marmot::Utilities::MakeString()
                                    << "EmbeddedLargeSlipBondElement " << elLabel
                                    << ": the slip exceeds the window of bar elements (the window is too small)" );

        // the local dofs: the host, then the partner's bar element
        const auto                    idx = barElementNodes( partner.element );
        std::array< int, nLocalDofs > dofs;
        for ( int j = 0; j < nHostNodes * nDim; j++ )
          dofs[j] = j;
        for ( int b = 0; b < nBarNodes; b++ )
          for ( int d = 0; d < nDim; d++ )
            dofs[( nHostNodes + b ) * nDim + d] = ( nHostNodes + idx[b] ) * nDim + d;
        LocalVector uLocal;
        for ( int j = 0; j < nLocalDofs; j++ )
          uLocal( j ) = Q( dofs[j] );

        const LocalKinematics kin = localKinematics( cp, partner.element, uLocal, partner.eta, partner.atKink );
        const double          st  = kin.slip( 0 );

        // covered fraction of the tributary part: the bar ends moved by the tangential slip
        double coveredStart = cp.tributaryStart, coveredEnd = cp.tributaryEnd;
        double dStart_dSt = 0, dEnd_dSt = 0;
        if ( startIsBarEnd && cumulativeLengths_.front() + st > coveredStart ) {
          coveredStart = cumulativeLengths_.front() + st;
          dStart_dSt   = 1;
        }
        if ( endIsBarEnd && cumulativeLengths_.back() + st < coveredEnd ) {
          coveredEnd = cumulativeLengths_.back() + st;
          dEnd_dSt   = 1;
        }
        const double tributaryLength = cp.tributaryEnd - cp.tributaryStart;
        const double coveredRaw      = ( coveredEnd - coveredStart ) / tributaryLength;
        const double covered         = std::clamp( coveredRaw, 0.0, 1.0 );
        const double dCovered_dSt    = coveredRaw > 0 && coveredRaw < 1 ? ( dEnd_dSt - dStart_dSt ) / tributaryLength
                                                                        : 0.0;
        const double bondFade        = fade( covered );
        const double dFade_dSt       = 6.0 * covered * ( 1.0 - covered ) * dCovered_dSt;

        sv.partner         = partner.element + 0.5 * ( std::clamp( partner.eta, -1.0, 1.0 ) + 1.0 );
        sv.coveredFraction = covered;
        if ( bondFade <= 0.0 ) {
          // the bar has left the channel point: no bond; slip, bond stress and law history are kept for its return
          sv.active = 0.0;
          continue;
        }

        const double                          dissipationDensityOld = sv.dissipationDensity;
        Marmot::Materials::BondSlipLaw::State state{ sv.bondStress,
                                                     sv.elasticEnergyDensity,
                                                     sv.dissipationDensity,
                                                     sv.materialStateVars.data() };
        Eigen::Matrix3d                       dTau_dSlip;
        cp.bondSlipLaw->computeBondStress( state, dTau_dSlip, kin.slip, kin.slip - sv.slip, { time, dT } );

        sv.slip                 = kin.slip;
        sv.bondStress           = state.bondStress;
        sv.elasticEnergyDensity = state.elasticEnergyDensity;
        sv.dissipationDensity   = state.dissipation;
        // the energies of the point are those of its bonded part f: the forces, the reported elastic energy and the
        // dissipation increments are faded alike; the densities of the law are kept unfaded
        sv.dissipation += bondFade * cp.perimeterTimesLength * ( state.dissipation - dissipationDensityOld );
        sv.active = 1.0;

        const VectorDim   tau    = state.bondStress.template head< nDim >();
        const LocalVector BtTau  = kin.dSlip_dU.transpose() * tau;
        const LocalVector PLocal = bondFade * cp.perimeterTimesLength * BtTau;
        for ( int j = 0; j < nLocalDofs; j++ )
          P( dofs[j] ) += PLocal( j );

        if ( K ) {
          const LocalMatrix KLocal = cp.perimeterTimesLength *
                                     ( bondFade * kin.dSlip_dU.transpose() *
                                         dTau_dSlip.template topLeftCorner< nDim, nDim >() * kin.dSlip_dU +
                                       dFade_dSt * BtTau * kin.dSlip_dU.row( 0 ) +
                                       bondFade * geometricStiffness( cp, partner, uLocal, kin.eta, tau ) );
          for ( int a = 0; a < nLocalDofs; a++ )
            for ( int b = 0; b < nLocalDofs; b++ )
              ( *K )( dofs[a], dofs[b] ) += KLocal( a, b );
        }
      }
    }

    /**
     * @brief The geometric part of the tangent, \f$\partial (\boldsymbol{B}^T \boldsymbol{\tau}) / \partial
     * \boldsymbol{u}\f$ at fixed \f$\boldsymbol{\tau}\f$: the second derivatives of the closest point projection and of
     * the frame, by central differences of the analytic B over the local dofs (geometry only, no law evaluation).
     * @param cp      the channel point
     * @param partner its partner
     * @param uLocal  local displacements
     * @param eta     the partner parameter for uLocal (start of the closest point iterations)
     * @param tau     the bond stress, fixed
     */
    LocalMatrix geometricStiffness( const ChannelPoint& cp,
                                    const Partner&      partner,
                                    const LocalVector&  uLocal,
                                    double              eta,
                                    const VectorDim&    tau ) const
    {
      LocalMatrix  G;
      const double h = 1e-7 * characteristicLength_;
      LocalVector  u = uLocal;
      for ( int j = 0; j < nLocalDofs; j++ ) {
        u( j )                      = uLocal( j ) + h;
        const LocalVector BtTauPlus = localKinematics( cp, partner.element, u, eta, partner.atKink )
                                        .dSlip_dU.transpose() *
                                      tau;
        u( j )                       = uLocal( j ) - h;
        const LocalVector BtTauMinus = localKinematics( cp, partner.element, u, eta, partner.atKink )
                                         .dSlip_dU.transpose() *
                                       tau;
        u( j )     = uLocal( j );
        G.col( j ) = ( BtTauPlus - BtTauMinus ) / ( 2 * h );
      }
      return G;
    }

    void computeKernels( const double* QTotal_, const double*, double* Pe_, double* Ke_, double time, double dT )
      override
    {
      const Eigen::Map< const Eigen::VectorXd > QTotal( QTotal_, nDofs );
      Eigen::Map< Eigen::VectorXd >             Pe( Pe_, nDofs );
      Eigen::Map< Eigen::MatrixXd >             KeMap( Ke_, nDofs, nDofs );
      Eigen::Ref< Eigen::MatrixXd >             Ke( KeMap );
      evaluate( QTotal, Pe, &Ke, stateVars_, time, dT );
    }

    void computeKernelsExplicit( const double* QTotal_, const double*, double* Pe_, double time, double dT ) override
    {
      const Eigen::Map< const Eigen::VectorXd > QTotal( QTotal_, nDofs );
      Eigen::Map< Eigen::VectorXd >             Pe( Pe_, nDofs );
      evaluate( QTotal, Pe, nullptr, stateVars_, time, dT );
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

    /// the elastic energy of the bonded parts, f x density x perimeter x tributary length
    void computeInternalEnergy( double& internalEnergy ) override
    {
      const int nQp = nStateVars_ / static_cast< int >( channelPoints.size() );
      for ( size_t i = 0; i < channelPoints.size(); i++ ) {
        const QPStateVarManager sv( stateVars_ + i * nQp, nQp );
        if ( sv.active > 0.5 )
          internalEnergy += fade( sv.coveredFraction ) * sv.elasticEnergyDensity *
                            channelPoints[i].perimeterTimesLength;
      }
    }

    Marmot::Utilities::StateView getStateView( const std::string& stateName, int qpNumber ) override
    {
      const int         nQp  = nStateVars_ / static_cast< int >( channelPoints.size() );
      double*           base = stateVars_ + qpNumber * nQp;
      QPStateVarManager sv( base, nQp );
      if ( sv.contains( stateName ) )
        return sv.getStateView( stateName );
      return channelPoints[qpNumber].bondSlipLaw->getStateView( stateName, sv.materialStateVars.data() );
    }

    std::vector< double > getCoordinatesAtCenter() override
    {
      VectorDim X = VectorDim::Zero();
      for ( const auto& cp : channelPoints )
        X += cp.X / channelPoints.size();
      return std::vector< double >( X.data(), X.data() + nDim );
    }

    std::vector< std::vector< double > > getCoordinatesAtQuadraturePoints() override
    {
      std::vector< std::vector< double > > result;
      for ( const auto& cp : channelPoints )
        result.emplace_back( cp.X.data(), cp.X.data() + nDim );
      return result;
    }

    int getNumberOfQuadraturePoints() override { return static_cast< int >( channelPoints.size() ); }
  };

} // namespace Marmot::Discretizations::Elements
