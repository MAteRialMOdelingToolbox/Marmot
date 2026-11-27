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
#include "Marmot/MarmotAutomaticDifferentiationTypes.h"
#include <Eigen/Dense>
#include <cmath>

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
   * conjugate to \f$\delta \boldsymbol{\Psi}\f$, not to the spins: a nodal moment load is a physical moment only for
   * small rotations. No element state beyond the section states is needed.
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
   * tangent \f$\boldsymbol{K}_l\f$. With \f$\partial \boldsymbol{q}_l / \partial \boldsymbol{Q}\f$ (code:
   * dLocalDeformations_dQ):
   * \f[ \boldsymbol{P} = \left( \frac{\partial \boldsymbol{q}_l}{\partial \boldsymbol{Q}} \right)^T \boldsymbol{f}_l,
   * \quad \boldsymbol{K} = \left( \frac{\partial \boldsymbol{q}_l}{\partial \boldsymbol{Q}} \right)^T \boldsymbol{K}_l
   * \frac{\partial \boldsymbol{q}_l}{\partial \boldsymbol{Q}} + \frac{\partial^2 ( \boldsymbol{f}_l \cdot
   * \boldsymbol{q}_l )}{\partial \boldsymbol{Q}^2} \Big|_{\boldsymbol{f}_l}, \f]
   * the consistent tangent: material stiffness plus geometric stiffness (code: geometricStiffness, the Hessian of the
   * work of the fixed local forces on the local deformations). Both derivatives are evaluated exactly by forward
   * automatic differentiation (dual numbers) of the closed-form map \f$\boldsymbol{Q} \mapsto \boldsymbol{q}_l\f$;
   * the exponential and the logarithm use Taylor series near zero rotation, so the derivatives are exact also in the
   * reference configuration. As in B2D2CR the section is incremental: it is fed with the increment of the local
   * generalized strains between the last converged and the current state. The steps in computeKernels(), named as in
   * CorotationalBeamElement2D: (1) rigid frame and (2) local deformations (localDeformations()), (3) section
   * (BeamFiberSection, per quadrature point), (4) transformation to the global dofs and geometric stiffness.
   *
   * @par Relation to Battini and Pacoste (2002)
   * The local deformations \f$\boldsymbol{q}_l = [ \bar u, \bar{\boldsymbol{\theta}}_1, \bar{\boldsymbol{\theta}}_2
   * ]\f$, the co-rotated frame and the local element are theirs. Their matrices (the projector \f$\boldsymbol{P}\f$,
   * \f$\boldsymbol{G}\f$ = the derivative of the rigid rotation of the frame, \f$\boldsymbol{E}\f$, the matrix
   * \f$\boldsymbol{B}\f$ of the local deformations and the terms of the geometric stiffness) are formed for variations
   * of the global rotations by spins \f$\delta\boldsymbol{w}\f$, while the dofs here are additive rotation vectors,
   * \f$\delta\boldsymbol{w}_i = \boldsymbol{T}_s( \boldsymbol{\Psi}_i ) \delta\boldsymbol{\Psi}_i\f$ with the
   * tangential transformation \f$\boldsymbol{T}_s\f$ of the rotation vector. Hence \f$\partial \boldsymbol{q}_l /
   * \partial \boldsymbol{Q} = \boldsymbol{B}_{\mathrm{BP}} \, \mathrm{diag}( \boldsymbol{I}, \boldsymbol{T}_s(
   * \boldsymbol{\Psi}_1 ), \boldsymbol{I}, \boldsymbol{T}_s( \boldsymbol{\Psi}_2 ) )\f$, and the geometric stiffness
   * here contains, besides their terms, the derivative of \f$\boldsymbol{T}_s\f$. Their matrices are not formed:
   * automatic differentiation gives the exact derivatives for the additive parametrization directly, and \f$\partial^2
   * ( \boldsymbol{f}_l \cdot \boldsymbol{q}_l ) / \partial \boldsymbol{Q}^2\f$ is the compact statement of their
   * geometric stiffness in these dofs. At \f$\boldsymbol{\Psi}_i = \mathbf{0}\f$, \f$\boldsymbol{T}_s =
   * \boldsymbol{I}\f$ and the two coincide.
   *
   * Validity: arbitrary rigid rotations and translations (\f$| \boldsymbol{\Psi}_i | < 2 \pi\f$), large displacements;
   * the local rotations relative to the co-rotated frame (\f$< \pi\f$) and the local strains must be small. For
   * \f$\boldsymbol{R}_{g,i} = \boldsymbol{R}_0\f$ (no rotation) the tangent equals that of B3D2.
   *
   * A single element (or a chain) whose nodes are coupled to the rest of the model only through their displacements
   * has a zero-energy twist about its own axis; the model must restrain it (e.g. a soft rotational spring).
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

    /// The value of a scalar or of a dual number (for the branch decisions inside the differentiated functions)
    template < typename T >
    static double value( const T& x )
    {
      return static_cast< double >( x );
    }

    /// The rotation matrix exp(skew(w)) (Rodrigues), with Taylor series near zero: below |w|^2 = 1e-8 the truncation
    /// error of the series is O(|w|^12) < 1e-24, far below round-off, and the series (unlike sin(t)/t) has exact
    /// derivatives at w = 0, which the automatic differentiation of the tangent needs in the reference configuration
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
    /// zero (below x = 1e-8 the truncation error is O(x^3) < 1e-24, with exact derivatives at zero rotation, as in
    /// expSO3); throws (cutback) when the rotation reaches pi, where the logarithm is not unique
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

    /**
     * @brief (1) The rigid (co-rotated) frame and (2) the local deformations of the global dofs.
     * @tparam T double or a dual number type (automatic differentiation of the map)
     * @param[in] Q The global dofs (12): per node the displacement and the rotation vector.
     * @return The local deformations [elongation, theta1 (3), theta2 (3)].
     */
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
      // undefined frame: the second axes of the end nodes cancel (relative twist pi) or align with the chord; r1 and
      // p are of unit size, so |r3|^2 < 1e-24 means they are parallel to round-off
      if ( !( value( r3.dot( r3 ) ) > 1e-24 ) )
        throw Marmot::Utilities::StressUpdateFailed( Marmot::Utilities::MakeString()
                                                     << __PRETTY_FUNCTION__ << ": element " << elLabel
                                                     << ": the co-rotated frame is undefined (relative twist pi)" );
      r3          = r3 / sqrt( r3.dot( r3 ) );
      const V3 r2 = r3.cross( r1 );
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

    /// The local deformations (double)
    /// @param[in] Q The global dofs (12).
    /// @return The local deformations [elongation, theta1 (3), theta2 (3)].
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

    /**
     * @brief (2) local deformations and (3) the sections: the local forces and the local tangent; updates the states.
     * @param[in]  QTotal The global dofs at the end of the increment.
     * @param[in]  dQ     The increment of the global dofs.
     * @param[out] fl     The local forces, conjugate to the local deformations.
     * @param[out] Kl     The local tangent.
     * @param[in]  time   (Pseudo-)time.
     * @param[in]  dT     (Pseudo-)time increment.
     */
    void localForces( const Vector12d& QTotal, const Vector12d& dQ, Vector7d& fl, Matrix7d& Kl, double time, double dT )
    {
      const Vector7d q  = localDeformations( QTotal );
      const Vector7d dq = q - localDeformations( Vector12d( QTotal - dQ ) );

      fl.setZero();
      Kl.setZero();
      for ( size_t i = 0; i < qps.size(); i++ ) {
        auto&           qp = qps[i];
        const MatrixBl& Bl = localB[i];
        VectorGen       sectionForces;
        MatrixGen       dSectionForces_dSectionStrains;
        double          energy, dissipation;
        qp.section.computeSection( qp.state->sectionState,
                                   sectionForces,
                                   dSectionForces_dSectionStrains,
                                   energy,
                                   dissipation,
                                   Bl * dq,
                                   time,
                                   dT );

        qp.state->sectionForces  = sectionForces;
        qp.state->sectionStrains = Bl * q;
        qp.state->elasticEnergy  = energy * qp.dV;
        qp.state->dissipation    = dissipation * qp.dV;

        fl += Bl.transpose() * sectionForces * qp.dV;
        Kl += Bl.transpose() * dSectionForces_dSectionStrains * Bl * qp.dV;
      }
    }

    /// d q_l / d Q (7 x 12), exact by forward automatic differentiation of localDeformations(), all 12 directions
    /// in one evaluation
    /// @param[in] QTotal The global dofs.
    /// @return The Jacobian of the local deformations.
    Eigen::MatrixXd dLocalDeformations_dQ( const Vector12d& QTotal ) const
    {
      using namespace Marmot::Math::Differentiation::Automatic;
      using scalar                                  = JacobianDual< 12 >;
      Eigen::Matrix< scalar, Eigen::Dynamic, 1 > Qd = QTotal.cast< scalar >();
      for ( int j = 0; j < 12; j++ )
        seedDirection( Qd( j ), j );
      const Eigen::Matrix< scalar, Eigen::Dynamic, 1 > ql = localDeformations< scalar >( Qd );
      Eigen::MatrixXd                                  J( ql.size(), 12 );
      for ( int i = 0; i < ql.size(); i++ )
        for ( int j = 0; j < 12; j++ )
          J( i, j ) = partialDerivative( ql( i ), j );
      return J;
    }

    /// The geometric stiffness d^2 (f_l . q_l) / d Q^2 at fixed local forces, exact by automatic differentiation,
    /// the whole Hessian in one evaluation
    /// @param[in] QTotal The global dofs.
    /// @param[in] fl     The local forces.
    /// @return The geometric stiffness (12 x 12).
    Eigen::MatrixXd geometricStiffness( const Vector12d& QTotal, const Vector7d& fl ) const
    {
      using namespace Marmot::Math::Differentiation::Automatic;
      using scalar                                  = HessianDual< 12 >;
      Eigen::Matrix< scalar, Eigen::Dynamic, 1 > Qh = QTotal.cast< scalar >();
      for ( int j = 0; j < 12; j++ )
        seedDirection( Qh( j ), j );
      const Eigen::Matrix< scalar, Eigen::Dynamic, 1 > ql                = localDeformations< scalar >( Qh );
      scalar                                           workOfLocalForces = 0.;
      for ( int k = 0; k < 7; k++ )
        workOfLocalForces += fl( k ) * ql( k );
      Eigen::MatrixXd H( 12, 12 );
      for ( int i = 0; i < 12; i++ )
        for ( int j = 0; j < 12; j++ )
          H( i, j ) = secondPartialDerivative( workOfLocalForces, i, j );
      return H;
    }

    void computeKernels( const double* QTotal_, const double* dQ_, double* Pe_, double* Ke_, double time, double dT )
      override
    {
      Eigen::Map< const Vector12d > QTotal( QTotal_ );
      Eigen::Map< const Vector12d > dQ( dQ_ );
      Eigen::Map< Vector12d >       Pe( Pe_ );
      Eigen::Map< Matrix12d >       Ke( Ke_ );

      // (1) rigid frame, (2) local deformations, (3) section
      Vector7d fl;
      Matrix7d Kl;
      localForces( QTotal, dQ, fl, Kl, time, dT );

      // (4) transformation to the global dofs and geometric stiffness
      const Eigen::MatrixXd dq_dQ = dLocalDeformations_dQ( QTotal );
      Pe += dq_dQ.transpose() * fl;
      Ke += dq_dQ.transpose() * Kl * dq_dQ + geometricStiffness( QTotal, fl );
    }

    /// The internal forces only: the geometric stiffness (the Hessian) is not needed
    void computeKernelsExplicit( const double* QTotal_,
                                 const double* dQ_,
                                 double*       Pe_,
                                 double        time,
                                 double        dT ) override
    {
      Eigen::Map< const Vector12d > QTotal( QTotal_ );
      Eigen::Map< const Vector12d > dQ( dQ_ );
      Eigen::Map< Vector12d >       Pe( Pe_ );
      Vector7d                      fl;
      Matrix7d                      Kl;
      localForces( QTotal, dQ, fl, Kl, time, dT );
      Pe += dLocalDeformations_dQ( QTotal ).transpose() * fl;
    }

    /// dead load b A per length, lumped: b A L0 / 2 at each node, no moments; it acts on the axis (the moment of a
    /// section whose centroid is off the axis is neglected, unlike in the geometrically linear element)
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
