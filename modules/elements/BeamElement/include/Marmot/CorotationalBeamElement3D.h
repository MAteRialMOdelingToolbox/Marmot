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
#include "autodiff/forward/dual/eigen.hpp"
#include <Eigen/Dense>
#include <autodiff/forward/dual/dual.hpp>
#include <cmath>
#include <type_traits>

namespace Marmot::Discretizations::Elements {

  /**
   * @class Marmot::Discretizations::Elements::CorotationalBeamElement3D
   * @brief Co-rotational (finite rotation, small local strain) 2-node Euler-Bernoulli beam in space, with the fiber
   * section (BeamFiberSection) of BeamElement for any hypoelastic material of Marmot.
   *
   * Element name: @b B3D2CR. Node fields, dofs (\f$u_x, u_y, u_z, \Psi_x, \Psi_y, \Psi_z\f$ per node), element
   * properties \f$[ v_x, v_y, v_z, J, n, y_1, z_1, A_1, \dots ]\f$, section integration (3 Gauss points along the axis,
   * each with its own BeamFiberSection) and quadrature point states are those of B3D2 (BeamElement<3, 2>); the states
   * @b section strains, @b axial strain, @b curvatures and @b twist rate are the local (co-rotated) ones.
   *
   * @par Rotation parametrization
   * The nodal rotation dofs are the components of the total rotation vector \f$\boldsymbol{\Psi}_i\f$ of node \f$i\f$
   * (the framework accumulates the rotation increments additively, which gives exactly the total vector); the nodal
   * triad is \f$\boldsymbol{R}_{g,i} = \exp( \hat{\boldsymbol{\Psi}}_i ) \boldsymbol{R}_0\f$ with the reference
   * element frame \f$\boldsymbol{R}_0 = [ \boldsymbol{e}_1, \boldsymbol{e}_2, \boldsymbol{e}_3 ]\f$ of B3D2 (axis and
   * orientation vector). Residual and tangent are the exact derivatives with respect to these dofs (additive rotation
   * vector formulation, valid for \f$| \boldsymbol{\Psi}_i | < 2 \pi\f$). External nodal moments are therefore
   * conjugate to \f$\delta \boldsymbol{\Psi}\f$, not to the spins; this is irrelevant for embedded fibers, which have
   * no moment loads. No element state beyond the section states is needed.
   *
   * @par Co-rotated frame and local deformations (Battini and Pacoste, Comput. Methods Appl. Mech. Engrg. 191 (2002)
   * 1755-1789; Crisfield, Non-linear Finite Element Analysis of Solids and Structures, Vol. 2, ch. 17)
   * With the current chord \f$\boldsymbol{x}_{21} = \boldsymbol{X}_{21} + \boldsymbol{u}_2 - \boldsymbol{u}_1\f$
   * (length \f$l_n\f$): \f$\boldsymbol{r}_1 = \boldsymbol{x}_{21} / l_n\f$, \f$\boldsymbol{p} = ( \boldsymbol{R}_{g,1}
   * \boldsymbol{e}_2 + \boldsymbol{R}_{g,2} \boldsymbol{e}_2 ) / 2\f$ (the mean of the second axes of the nodal triads,
   * in the frame of the reference triad), \f$\boldsymbol{r}_3 = \boldsymbol{r}_1 \times \boldsymbol{p} / | \cdot |\f$,
   * \f$\boldsymbol{r}_2 = \boldsymbol{r}_3 \times \boldsymbol{r}_1\f$ define the rigid rotation \f$\boldsymbol{R}_r =
   * [ \boldsymbol{r}_1, \boldsymbol{r}_2, \boldsymbol{r}_3 ]\f$. The local deformational quantities are the elongation
   * \f$\bar u = l_n - L_0 = ( 2 \boldsymbol{X}_{21} \cdot \boldsymbol{d}_{21} + \boldsymbol{d}_{21} \cdot
   * \boldsymbol{d}_{21} ) / ( l_n + L_0 )\f$ and the local nodal rotation vectors \f$\bar{\boldsymbol{\theta}}_i =
   * \log( \boldsymbol{R}_r^T \boldsymbol{R}_{g,i} )\f$ (components twist, rotation about the local y and z axes); they
   * vanish exactly for any rigid motion. In the local frame the geometrically linear B3D2 kinematics apply:
   * \f$\varepsilon_0 = \bar u / L_0\f$, \f$\kappa_y\f$ and \f$\kappa_z\f$ from the cubic Hermite deflections with zero
   * end values and end slopes \f$\bar\theta_{i,y}\f$, \f$\bar\theta_{i,z}\f$, \f$\chi = ( \bar\theta_{2,x} -
   * \bar\theta_{1,x} ) / L_0\f$. The section gives the local forces \f$\boldsymbol{f}_l\f$ (conjugate to
   * \f$\boldsymbol{q}_l = [ \bar u, \bar{\boldsymbol{\theta}}_1, \bar{\boldsymbol{\theta}}_2 ]\f$) and the local
   * tangent \f$\boldsymbol{K}_l\f$. With \f$\boldsymbol{G} = \partial \boldsymbol{q}_l / \partial \boldsymbol{Q}\f$:
   * \f[ \boldsymbol{P} = \boldsymbol{G}^T \boldsymbol{f}_l, \quad \boldsymbol{K} = \boldsymbol{G}^T \boldsymbol{K}_l
   * \boldsymbol{G} + \frac{\partial^2 ( \boldsymbol{f}_l \cdot \boldsymbol{q}_l )}{\partial \boldsymbol{Q}^2}
   * \Big|_{\boldsymbol{f}_l}, \f]
   * the consistent tangent (material + geometric). \f$\boldsymbol{G}\f$ and the geometric stiffness are evaluated
   * exactly by forward automatic differentiation (dual numbers) of the closed-form map \f$\boldsymbol{Q} \mapsto
   * \boldsymbol{q}_l\f$; the exponential and the logarithm use Taylor series near zero rotation, so the derivatives are
   * exact also in the reference configuration. As in B2D2CR the section is incremental: it is fed with the increment
   * of the local generalized strains between the last converged and the current state.
   *
   * Validity: arbitrary rigid rotations and translations (\f$| \boldsymbol{\Psi}_i | < 2 \pi\f$), large displacements;
   * the local rotations relative to the co-rotated frame (\f$< \pi\f$) and the local strains must be small. For
   * \f$\boldsymbol{R}_{g,i} = \boldsymbol{R}_0\f$ (no rotation) the tangent equals that of B3D2.
   *
   * A single element (or a chain) whose nodes are only coupled through their displacements (embedded fibers) has a
   * zero-energy twist about its own axis; the model must restrain it (e.g. a soft rotational spring).
   *
   * Body forces (dead loads b A): lumped to the nodes (\f$b A L_0 / 2\f$ each), no nodal moments, no load stiffness.
   * Not implemented: distributed (surface) loads, inertia.
   */
  class CorotationalBeamElement3D : public BeamElement< 3, 2 > {

  public:
    using Base      = BeamElement< 3, 2 >;
    using Vector12d = Eigen::Matrix< double, 12, 1 >;
    using Matrix12d = Eigen::Matrix< double, 12, 12 >;
    using Vector7d  = Eigen::Matrix< double, 7, 1 >;
    using Matrix7d  = Eigen::Matrix< double, 7, 7 >;
    using MatrixBl  = Eigen::Matrix< double, 4, 7 >;

    /// local generalized strains [eps0, kappa_y, kappa_z, chi] of the local dofs [u, theta1 (3), theta2 (3)], per
    /// quadrature point
    std::vector< MatrixBl > localB;

    /// The reference chord
    Eigen::Vector3d X21 = Eigen::Vector3d::Zero();
    /// The reference element frame, columns e1, e2, e3
    Eigen::Matrix3d R0 = Eigen::Matrix3d::Identity();

    explicit CorotationalBeamElement3D( int elementID ) : Base( elementID ) {}

    template < typename T >
    static double value( const T& x )
    {
      if constexpr ( std::is_arithmetic_v< T > )
        return x;
      else
        return static_cast< double >( autodiff::val( x ) );
    }

    /// The rotation matrix exp(skew(w)) (Rodrigues), with Taylor series near zero
    template < typename T >
    static Eigen::Matrix< T, 3, 3 > expSO3( const Eigen::Matrix< T, 3, 1 >& w )
    {
      using std::cos;
      using std::sin;
      using std::sqrt;
      const T t2 = w.dot( w );
      T       a, b;
      if ( value( t2 ) < 1e-8 ) {
        a = 1. - t2 / 6. + t2 * t2 / 120.;
        b = 0.5 - t2 / 24. + t2 * t2 / 720.;
      }
      else {
        const T t = sqrt( t2 );
        a         = sin( t ) / t;
        b         = ( 1. - cos( t ) ) / t2;
      }
      Eigen::Matrix< T, 3, 3 > K;
      K << T( 0. ), -w( 2 ), w( 1 ), w( 2 ), T( 0. ), -w( 0 ), -w( 1 ), w( 0 ), T( 0. );
      Eigen::Matrix< T, 3, 3 > R = Eigen::Matrix< T, 3, 3 >::Identity();
      R += a * K + b * ( K * K );
      return R;
    }

    /// The rotation vector of a rotation matrix (rotation angle < pi) via its quaternion, with a Taylor series near
    /// zero
    template < typename T >
    static Eigen::Matrix< T, 3, 1 > logSO3( const Eigen::Matrix< T, 3, 3 >& R )
    {
      using std::atan;
      using std::sqrt;
      const T trace = R( 0, 0 ) + R( 1, 1 ) + R( 2, 2 );
      if ( !( value( trace ) > -1 + 1e-12 ) )
        throw Marmot::Utilities::StressUpdateFailed( Marmot::Utilities::MakeString()
                                                     << __PRETTY_FUNCTION__ << ": local rotation reached pi" );
      const T                  qw = 0.5 * sqrt( 1. + trace );
      Eigen::Matrix< T, 3, 1 > v;
      v << R( 2, 1 ) - R( 1, 2 ), R( 0, 2 ) - R( 2, 0 ), R( 1, 0 ) - R( 0, 1 );
      v /= 4. * qw;
      const T s2 = v.dot( v );
      const T x  = s2 / ( qw * qw );
      T       g;
      if ( value( x ) < 1e-8 )
        g = 2. / qw * ( 1. - x / 3. + x * x / 5. );
      else {
        const T s = sqrt( s2 );
        g         = 2. * atan( s / qw ) / s;
      }
      return g * v;
    }

    /// The local deformations [elongation, theta1 (3), theta2 (3)] of the global dofs Q (12)
    template < typename T >
    Eigen::Matrix< T, Eigen::Dynamic, 1 > localDeformations( const Eigen::Matrix< T, Eigen::Dynamic, 1 >& Q ) const
    {
      using std::sqrt;
      using V3      = Eigen::Matrix< T, 3, 1 >;
      using M3      = Eigen::Matrix< T, 3, 3 >;
      const V3 u1   = Q.template segment< 3 >( 0 );
      const V3 w1   = Q.template segment< 3 >( 3 );
      const V3 u2   = Q.template segment< 3 >( 6 );
      const V3 w2   = Q.template segment< 3 >( 9 );
      const V3 X    = X21.template cast< T >();
      const V3 d21  = u2 - u1;
      const V3 x21  = X + d21;
      const T  ln   = sqrt( x21.dot( x21 ) );
      const T  ubar = ( 2. * X.dot( d21 ) + d21.dot( d21 ) ) / ( ln + length );
      const M3 R0T  = R0.template cast< T >();
      const M3 Rg1  = expSO3< T >( w1 ) * R0T;
      const M3 Rg2  = expSO3< T >( w2 ) * R0T;
      const V3 r1   = x21 / ln;
      const V3 p    = 0.5 * ( Rg1.col( 1 ) + Rg2.col( 1 ) );
      V3       r3   = r1.cross( p );
      r3            = r3 / sqrt( r3.dot( r3 ) );
      const V3 r2   = r3.cross( r1 );
      M3       Rr;
      Rr.col( 0 ) = r1, Rr.col( 1 ) = r2, Rr.col( 2 ) = r3;
      const M3                              R1 = Rr.transpose() * Rg1;
      const M3                              R2 = Rr.transpose() * Rg2;
      Eigen::Matrix< T, Eigen::Dynamic, 1 > q( 7 );
      q( 0 )                       = ubar;
      q.template segment< 3 >( 1 ) = logSO3< T >( R1 );
      q.template segment< 3 >( 4 ) = logSO3< T >( R2 );
      return q;
    }

    Vector7d localDeformations( const Vector12d& Q ) const
    {
      return localDeformations< double >( Eigen::VectorXd( Q ) );
    }

    void initializeYourself() override
    {
      Base::initializeYourself();
      X21 = coordinates.col( 1 ) - coordinates.col( 0 );
      R0  = R.transpose();
      localB.clear();
      for ( const auto& qp : qps ) {
        const double                               s   = 0.5 * ( qp.xi + 1 );
        const Eigen::Matrix< double, 1, nHermite > d2H = hermite( s, 2 );
        MatrixBl                                   Bl  = MatrixBl::Zero();
        Bl( 0, 0 )                                     = 1 / length;        // eps0 = u / L0
        Bl( 1, 2 )                                     = d2H( 1 ) / length; // kappa_y: slopes theta_y of node 1
        Bl( 1, 5 )                                     = d2H( 3 ) / length; //          and node 2
        Bl( 2, 3 )                                     = d2H( 1 ) / length; // kappa_z: slopes theta_z
        Bl( 2, 6 )                                     = d2H( 3 ) / length;
        Bl( 3, 1 )                                     = -1 / length;       // chi = (theta_2x - theta_1x) / L0
        Bl( 3, 4 )                                     = 1 / length;
        localB.push_back( Bl );
      }
    }

    void computeKernels( const double* QTotal_, const double* dQ_, double* Pe_, double* Ke_, double time, double dT )
      override
    {
      Eigen::Map< const Vector12d > QTotal( QTotal_ );
      Eigen::Map< const Vector12d > dQ( dQ_ );
      Eigen::Map< Vector12d >       Pe( Pe_ );
      Eigen::Map< Matrix12d >       Ke( Ke_ );

      const Vector12d QOld = QTotal - dQ;
      const Vector7d  q    = localDeformations( Vector12d( QTotal ) );
      const Vector7d  dq   = q - localDeformations( QOld );

      Vector7d fl = Vector7d::Zero();
      Matrix7d Kl = Matrix7d::Zero();
      for ( size_t i = 0; i < qps.size(); i++ ) {
        auto&           qp = qps[i];
        const MatrixBl& Bl = localB[i];
        VectorGen       s;
        MatrixGen       D;
        double          energy, dissipation;
        qp.section.computeSection( qp.stateVars + nQpStateVarsOwn, s, D, energy, dissipation, Bl * dq, time, dT );

        Eigen::Map< VectorGen >( qp.stateVars + idxSectionForces )  = s;
        Eigen::Map< VectorGen >( qp.stateVars + idxSectionStrains ) = Bl * q;
        qp.stateVars[idxElasticEnergy]                              = energy * qp.dV;
        qp.stateVars[idxDissipation]                                = dissipation * qp.dV;

        fl += Bl.transpose() * s * qp.dV;
        Kl += Bl.transpose() * D * Bl * qp.dV;
      }

      using autodiff::at;
      using autodiff::wrt;
      // G = d q / d Q
      autodiff::VectorXdual Qd = QTotal.cast< autodiff::dual >();
      autodiff::VectorXdual F;
      const Eigen::MatrixXd
        G = autodiff::jacobian( [this]( const autodiff::VectorXdual& x )
                                  -> autodiff::VectorXdual { return localDeformations< autodiff::dual >( x ); },
                                wrt( Qd ),
                                at( Qd ),
                                F );

      // geometric stiffness: the Hessian of f_l . q_l(Q) at fixed f_l
      autodiff::VectorXdual2nd Qh = QTotal.cast< autodiff::dual2nd >();
      autodiff::dual2nd        phi;
      Eigen::VectorXd          g;
      const Eigen::MatrixXd    H = autodiff::hessian(
        [this, &fl]( const autodiff::VectorXdual2nd& x ) -> autodiff::dual2nd {
          const autodiff::VectorXdual2nd ql = localDeformations< autodiff::dual2nd >( x );
          autodiff::dual2nd              r  = 0.;
          for ( int k = 0; k < 7; k++ )
            r += fl( k ) * ql( k );
          return r;
        },
        wrt( Qh ),
        at( Qh ),
        phi,
        g );

      Pe += G.transpose() * fl;
      Ke += G.transpose() * Kl * G + H;
    }

    void computeKernelsExplicit( const double* QTotal, const double* dQ, double* Pe, double time, double dT ) override
    {
      Matrix12d K = Matrix12d::Zero();
      computeKernels( QTotal, dQ, Pe, K.data(), time, dT );
    }

    /// dead load b A per length, lumped: b A L0 / 2 at each node, no moments
    void computeBodyForce( double* P_, double*, const double* load, const double*, double, double ) override
    {
      Eigen::Map< Vector12d >             P( P_ );
      Eigen::Map< const Eigen::Vector3d > b( load );
      const Eigen::Vector3d               Fb = 0.5 * b * area() * length;
      P.segment< 3 >( 0 ) += Fb;
      P.segment< 3 >( 6 ) += Fb;
    }
  };

} // namespace Marmot::Discretizations::Elements
