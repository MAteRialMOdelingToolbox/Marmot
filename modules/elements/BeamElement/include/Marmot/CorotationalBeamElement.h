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
#include "Marmot/BeamElement.h"
#include <Eigen/Dense>
#include <cmath>
#include <numbers>

namespace Marmot::Discretizations::Elements {

  /**
   * @class Marmot::Discretizations::Elements::CorotationalBeamElement2D
   * @brief Co-rotational (finite rotation, small local strain) 2-node Euler-Bernoulli beam in the plane, with the fiber
   * section (BeamFiberSection) of BeamElement for any hypoelastic material of Marmot.
   *
   * Element name: @b B2D2CR. Node fields, dofs (\f$u_x, u_y, \theta_z\f$ per node), element properties
   * \f$[ n, y_1, A_1, \dots, y_n, A_n ]\f$, section integration (3 Gauss points along the axis, each with its own
   * BeamFiberSection) and quadrature point states are those of B2D2 (BeamElement<2, 2>); the states @b section
   * strains, @b axial strain and @b curvature are the local (co-rotated) ones.
   *
   * @par Formulation (Crisfield, Non-linear Finite Element Analysis of Solids and Structures, Vol. 1, ch. 7.4)
   * The rigid body motion of the chord is removed exactly: with the reference chord \f$\boldsymbol{X}_{21}\f$ (length
   * \f$L_0\f$, angle \f$\beta_0\f$) and the current chord \f$\boldsymbol{x}_{21} = \boldsymbol{X}_{21} +
   * \boldsymbol{u}_2 - \boldsymbol{u}_1\f$ (length \f$l_n\f$, angle \f$\beta\f$), the local deformational quantities
   * are
   * - the elongation \f$\bar u = l_n - L_0 = ( 2 \boldsymbol{X}_{21} \cdot \boldsymbol{d}_{21} + \boldsymbol{d}_{21}
   *   \cdot \boldsymbol{d}_{21} ) / ( l_n + L_0 )\f$ (exactly 0 for a rigid motion, no cancellation),
   * - the end rotations relative to the chord \f$\bar\theta_i = \theta_i - \alpha\f$, \f$\alpha = \beta - \beta_0\f$
   *   (evaluated with atan2 of \f$\sin\alpha, \cos\alpha\f$ and wrapped into \f$(-\pi, \pi]\f$, so arbitrarily large
   *   rotations, also several turns, are allowed; the local rotations must stay below \f$\pi\f$).
   *
   * In the local frame the geometrically linear B2D2 kinematics apply (\f$\varepsilon_0 = \bar u / L_0\f$, \f$\kappa
   * = v''\f$ of the cubic Hermite deflection with \f$v = 0\f$ at both ends and slopes \f$\bar\theta_i\f$); the section
   * gives the local forces \f$\boldsymbol{f}_l = [ N, M_1, M_2 ]\f$ (conjugate to \f$[ \bar u, \bar\theta_1,
   * \bar\theta_2 ]\f$) and the local tangent \f$\boldsymbol{K}_l\f$. With \f$\boldsymbol{r} = [ -c, -s, 0, c, s, 0
   * ]\f$, \f$\boldsymbol{z} = [ s, -c, 0, -s, c, 0 ]\f$ (\f$c = \cos\beta, s = \sin\beta\f$) and \f$\boldsymbol{B} =
   * [ \boldsymbol{r}; \boldsymbol{e}_3 - \boldsymbol{z} / l_n; \boldsymbol{e}_6 - \boldsymbol{z} / l_n ]\f$:
   * \f[ \boldsymbol{P} = \boldsymbol{B}^T \boldsymbol{f}_l, \quad \boldsymbol{K} = \boldsymbol{B}^T \boldsymbol{K}_l
   * \boldsymbol{B} + \frac{N}{l_n} \boldsymbol{z} \boldsymbol{z}^T + \frac{M_1 + M_2}{l_n^2} ( \boldsymbol{r}
   * \boldsymbol{z}^T + \boldsymbol{z} \boldsymbol{r}^T ), \f]
   * the consistent tangent (material + geometric). The section is incremental (hypoelastic materials): it is fed with
   * the increment of the local generalized strains between the last converged and the current state, i.e., the strains
   * are co-rotated with the chord.
   *
   * Validity: arbitrary rigid rotations and translations, large displacements; the local rotations relative to the
   * chord and the local strains must be small (refine the mesh where the beam bends strongly). Axial strain measure:
   * engineering strain of the chord.
   *
   * Body forces (dead loads b A): lumped to the nodes (\f$b A L_0 / 2\f$ each), no nodal moments, no load stiffness.
   * Not implemented: distributed (surface) loads, inertia.
   */
  class CorotationalBeamElement2D : public BeamElement< 2, 2 > {

  public:
    using Base     = BeamElement< 2, 2 >;
    using Vector6d = Eigen::Matrix< double, 6, 1 >;
    using Matrix6d = Eigen::Matrix< double, 6, 6 >;

    /// local generalized strains [eps0, kappa] of the local dofs [u, theta1, theta2], per quadrature point
    std::vector< Eigen::Matrix< double, 2, 3 > > localB;

    /// The reference chord
    Eigen::Vector2d X21 = Eigen::Vector2d::Zero();

    /// The current chord and the local deformational quantities for the global dofs
    struct Chord {
      double          ln;     ///< current length
      double          c, s;   ///< direction cosines of the current chord
      Eigen::Vector3d qLocal; ///< [elongation, theta1 - alpha, theta2 - alpha]
    };

    explicit CorotationalBeamElement2D( int elementID ) : Base( elementID ) {}

    Chord chord( const double* Q ) const
    {
      const Eigen::Vector2d d21( Q[3] - Q[0], Q[4] - Q[1] );
      const Eigen::Vector2d x21 = X21 + d21;
      Chord                 ch;
      ch.ln = x21.norm();
      if ( !( ch.ln > 1e-12 * length ) )
        throw Marmot::Utilities::StressUpdateFailed( Marmot::Utilities::MakeString()
                                                     << __PRETTY_FUNCTION__ << ": element " << elLabel
                                                     << " collapsed to zero length" );
      ch.c = x21( 0 ) / ch.ln;
      ch.s = x21( 1 ) / ch.ln;
      // rigid rotation of the chord, alpha = beta - beta0 in (-pi, pi]
      const double     c0 = X21( 0 ) / length, s0 = X21( 1 ) / length;
      const double     alpha = std::atan2( c0 * ch.s - s0 * ch.c, c0 * ch.c + s0 * ch.s );
      constexpr double twoPi = 2 * std::numbers::pi;
      ch.qLocal << ( 2 * X21.dot( d21 ) + d21.squaredNorm() ) / ( ch.ln + length ),
        std::remainder( Q[2] - alpha, twoPi ), std::remainder( Q[5] - alpha, twoPi );
      return ch;
    }

    void initializeYourself() override
    {
      Base::initializeYourself();
      X21 = coordinates.col( 1 ) - coordinates.col( 0 );
      localB.clear();
      for ( const auto& qp : qps ) {
        const double s                                 = 0.5 * ( qp.xi + 1 );
        const auto [N, dNds]                           = lagrange( s );
        const Eigen::Matrix< double, 1, nHermite > d2H = hermite( s, 2 );
        Eigen::Matrix< double, 2, 3 >              Bl  = Eigen::Matrix< double, 2, 3 >::Zero();
        Bl( 0, 0 )                                     = dNds( 1 ) / length; // = 1 / L0
        Bl( 1, 1 )                                     = d2H( 1 ) / length;  // slope part of node 1
        Bl( 1, 2 )                                     = d2H( 3 ) / length;  // slope part of node 2
        localB.push_back( Bl );
      }
    }

    void computeKernels( const double* QTotal_, const double* dQ_, double* Pe_, double* Ke_, double time, double dT )
      override
    {
      Eigen::Map< const Vector6d > QTotal( QTotal_ );
      Eigen::Map< const Vector6d > dQ( dQ_ );
      Eigen::Map< Vector6d >       Pe( Pe_ );
      Eigen::Map< Matrix6d >       Ke( Ke_ );

      const Vector6d QOld = QTotal - dQ;
      const Chord    ch   = chord( QTotal.data() );
      const Chord    chO  = chord( QOld.data() );
      // local increments; the rotations are wrapped individually, the difference stays small
      Eigen::Vector3d dqLocal = ch.qLocal - chO.qLocal;
      for ( int i = 1; i < 3; i++ )
        dqLocal( i ) = std::remainder( dqLocal( i ), 2 * std::numbers::pi );

      Eigen::Vector3d fl = Eigen::Vector3d::Zero();
      Eigen::Matrix3d Kl = Eigen::Matrix3d::Zero();
      for ( size_t i = 0; i < qps.size(); i++ ) {
        auto&                                qp = qps[i];
        const Eigen::Matrix< double, 2, 3 >& Bl = localB[i];
        VectorGen                            s;
        MatrixGen                            D;
        double                               energy, dissipation;
        qp.section.computeSection( qp.stateVars + nQpStateVarsOwn, s, D, energy, dissipation, Bl * dqLocal, time, dT );

        Eigen::Map< VectorGen >( qp.stateVars + idxSectionForces )  = s;
        Eigen::Map< VectorGen >( qp.stateVars + idxSectionStrains ) = Bl * ch.qLocal;
        qp.stateVars[idxElasticEnergy]                              = energy * qp.dV;
        qp.stateVars[idxDissipation]                                = dissipation * qp.dV;

        fl += Bl.transpose() * s * qp.dV;
        Kl += Bl.transpose() * D * Bl * qp.dV;
      }

      Vector6d r, z;
      r << -ch.c, -ch.s, 0, ch.c, ch.s, 0;
      z << ch.s, -ch.c, 0, -ch.s, ch.c, 0;
      Eigen::Matrix< double, 3, 6 > B;
      B.row( 0 ) = r.transpose();
      B.row( 1 ) = -z.transpose() / ch.ln;
      B.row( 2 ) = -z.transpose() / ch.ln;
      B( 1, 2 ) += 1;
      B( 2, 5 ) += 1;

      Pe += B.transpose() * fl;
      Ke += B.transpose() * Kl * B + fl( 0 ) / ch.ln * z * z.transpose() +
            ( fl( 1 ) + fl( 2 ) ) / ( ch.ln * ch.ln ) * ( r * z.transpose() + z * r.transpose() );
    }

    /// dead load b A per length, lumped: b A L0 / 2 at each node, no moments
    void computeBodyForce( double* P_, double*, const double* load, const double*, double, double ) override
    {
      Eigen::Map< Vector6d >              P( P_ );
      Eigen::Map< const Eigen::Vector2d > b( load );
      const Eigen::Vector2d               F = 0.5 * b * area() * length;
      P.segment< 2 >( 0 ) += F;
      P.segment< 2 >( 3 ) += F;
    }
  };

} // namespace Marmot::Discretizations::Elements
