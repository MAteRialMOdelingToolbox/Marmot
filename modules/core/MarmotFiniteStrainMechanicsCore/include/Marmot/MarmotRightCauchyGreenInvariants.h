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
#include "Marmot/MarmotFastorTensorBasics.h"
#include <array>
#include <cstddef>
#include <tuple>
#include <type_traits>

/** @brief Invariants of the right Cauchy-Green tensor \f$\boldsymbol C\f$ and their derivatives w.r.t.
 * \f$\boldsymbol C\f$, and the chain rule that turns the scalar derivatives of an isotropic function of these
 * invariants into \f$\partial\Psi/\partial\boldsymbol C\f$ and
 * \f$\partial^2\Psi/\partial\boldsymbol C\partial\boldsymbol C\f$.
 *
 * Two sets of three invariants \f$\boldsymbol x\f$ are available, the principal invariants
 * \f$(I_1, I_2, J)\f$ (#computeInvariants) and the isochoric invariants \f$(\bar I_1, \bar I_2, J)\f$
 * (#computeIsochoricInvariants), with
 * \f[
 *   I_1 = \mathrm{tr}\,\boldsymbol C, \qquad
 *   I_2 = \tfrac{1}{2}\left(I_1^2 - \mathrm{tr}\,\boldsymbol C^2\right), \qquad
 *   J = \sqrt{\det\boldsymbol C}, \qquad
 *   \bar I_1 = J^{-2/3} I_1, \qquad \bar I_2 = J^{-4/3} I_2.
 * \f]
 * All derivatives treat \f$\boldsymbol C\f$ as a general (not necessarily symmetric) tensor, i.e. they are the
 * derivatives w.r.t. the nine independent components \f$C_{ij}\f$, as obtained from automatic or complex-step
 * differentiation.
 *
 * A model only provides a scalar function with its scalar derivatives (#ScalarDerivatives); the tensor algebra
 * is written once, in the two functions computing the invariants and in #composeWithInvariants:
 * @code
 * using namespace Marmot::ContinuumMechanics::Invariants::RightCauchyGreen;
 * const auto inv = computeIsochoricInvariants< 2 >( C );
 * ScalarDerivatives< double, 2 > f;      // f( Ibar1, Ibar2, J ), its gradient and Hessian
 * ...
 * const auto [psi, dPsi_dC, d2Psi_dCdC] = composeWithInvariants( f, inv );
 * @endcode
 */
namespace Marmot::ContinuumMechanics::Invariants::RightCauchyGreen {

  /// @brief Position of \f$I_1\f$ or \f$\bar I_1\f$ in the invariant arrays.
  inline constexpr std::size_t indexI1 = 0;
  /// @brief Position of \f$I_2\f$ or \f$\bar I_2\f$ in the invariant arrays.
  inline constexpr std::size_t indexI2 = 1;
  /// @brief Position of \f$J\f$ in the invariant arrays.
  inline constexpr std::size_t indexJ = 2;

  /// @brief Placeholder for a derivative that is not computed at the requested order.
  struct NotComputed {};

  /** @brief Value and scalar derivatives of a function \f$f(x_1, x_2, x_3)\f$ of three invariants.
   *
   * Only the entries up to @p Order are meaningful: @c f for every order, @c df for @p Order >= 1 and @c d2f for
   * @p Order >= 2. Entries that are not set are zero.
   *
   * @tparam T Scalar type, e.g. double, std::complex<double>, autodiff scalar.
   * @tparam Order Highest derivative that is filled (0, 1 or 2).
   */
  template < typename T, int Order >
  struct ScalarDerivatives {
    /// @brief Function value.
    T f{};
    /// @brief Gradient \f$\partial f/\partial x_a\f$.
    std::array< T, 3 > df{};
    /// @brief Hessian \f$\partial^2 f/\partial x_a\partial x_b\f$.
    std::array< std::array< T, 3 >, 3 > d2f{};
  };

  /** @brief Three invariants \f$\boldsymbol x\f$ of \f$\boldsymbol C\f$ and their derivatives w.r.t.
   * \f$\boldsymbol C\f$ up to @p Order.
   *
   * @tparam T Scalar type, e.g. double, std::complex<double>, autodiff scalar.
   * @tparam Order Highest derivative that is available (0, 1 or 2).
   */
  template < typename T, int Order >
  struct InvariantsAndDerivatives {
    /// @brief The invariants, \f$(I_1, I_2, J)\f$ or \f$(\bar I_1, \bar I_2, J)\f$.
    std::array< T, 3 > x;
    /// @brief First derivatives \f$\partial x_a/\partial\boldsymbol C\f$ (@p Order >= 1).
    std::conditional_t< ( Order >= 1 ), std::array< Math::Tensors::FastorTensors::Tensor33t< T >, 3 >, NotComputed >
      dx_dC;
    /// @brief Second derivatives \f$\partial^2 x_a/\partial\boldsymbol C\partial\boldsymbol C\f$ (@p Order >= 2).
    std::conditional_t< ( Order >= 2 ), std::array< Math::Tensors::FastorTensors::Tensor3333t< T >, 3 >, NotComputed >
      d2x_dCdC;
  };

  /** @brief Second principal invariant \f$I_2 = \tfrac{1}{2}\left(I_1^2 - \mathrm{tr}\,\boldsymbol C^2\right)\f$,
   * evaluated as the sum of the principal 2x2 minors of \f$\boldsymbol C\f$.
   *
   * @tparam T Scalar type, e.g. double, std::complex<double>, autodiff scalar.
   * @param C Right Cauchy-Green tensor.
   * @return \f$I_2\f$.
   */
  template < typename T >
  T secondInvariant( const Math::Tensors::FastorTensors::Tensor33t< T >& C )
  {
    return C( 0, 0 ) * C( 1, 1 ) + C( 1, 1 ) * C( 2, 2 ) + C( 0, 0 ) * C( 2, 2 ) - C( 0, 1 ) * C( 1, 0 ) -
           C( 1, 2 ) * C( 2, 1 ) - C( 0, 2 ) * C( 2, 0 );
  }

  /** @brief Compute the principal invariants \f$(I_1, I_2, J)\f$ of \f$\boldsymbol C\f$ and their derivatives
   * w.r.t. \f$\boldsymbol C\f$ up to @p Order,
   * \f[
   *   \frac{\partial I_1}{\partial\boldsymbol C} = \boldsymbol I, \qquad
   *   \frac{\partial I_2}{\partial\boldsymbol C} = I_1\boldsymbol I - \boldsymbol C^T, \qquad
   *   \frac{\partial J}{\partial\boldsymbol C} = \frac{J}{2}\boldsymbol C^{-T},
   * \f]
   * \f[
   *   \frac{\partial^2 I_2}{\partial C_{ij}\partial C_{kl}} = \delta_{ij}\delta_{kl} - \delta_{il}\delta_{jk},
   *   \qquad
   *   \frac{\partial^2 J}{\partial C_{ij}\partial C_{kl}} = \frac{J}{4}C^{-1}_{ji}C^{-1}_{lk} -
   *   \frac{J}{2}C^{-1}_{jk}C^{-1}_{li}.
   * \f]
   *
   * @tparam Order Highest derivative to compute (0, 1 or 2).
   * @tparam T Scalar type, e.g. double, std::complex<double>, autodiff scalar.
   * @param C Right Cauchy-Green tensor.
   * @param withI2 Pass @c false for a function that does not depend on \f$I_2\f$, to skip its evaluation;
   * \f$I_2\f$ and its derivatives are then not meaningful.
   * @return The invariants and their derivatives.
   */
  template < int Order, typename T >
  InvariantsAndDerivatives< T, Order > computeInvariants( const Math::Tensors::FastorTensors::Tensor33t< T >& C,
                                                          const bool withI2 = true )
  {
    static_assert( Order >= 0 && Order <= 2, "computeInvariants: only derivatives up to the second order" );
    using namespace Fastor;

    InvariantsAndDerivatives< T, Order > inv;

    const T I1     = trace( C );
    const T J      = sqrt( determinant( C ) );
    inv.x[indexI1] = I1;
    inv.x[indexI2] = withI2 ? secondInvariant( C ) : T( 0. );
    inv.x[indexJ]  = J;

    if constexpr ( Order >= 1 ) {
      const Math::Tensors::FastorTensors::Tensor33t< T > CInv  = inverse( C );
      const T                                            halfJ = 0.5 * J;
      const auto                                         delta = []( int a, int b ) { return a == b ? 1. : 0.; };

      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ ) {
          inv.dx_dC[indexI1]( i, j ) = T( delta( i, j ) );
          inv.dx_dC[indexI2]( i, j ) = I1 * delta( i, j ) - C( j, i );
          inv.dx_dC[indexJ]( i, j )  = halfJ * CInv( j, i );
        }

      if constexpr ( Order >= 2 ) {
        const T quarterJ = 0.25 * J;
        for ( int i = 0; i < 3; i++ )
          for ( int j = 0; j < 3; j++ )
            for ( int k = 0; k < 3; k++ )
              for ( int l = 0; l < 3; l++ ) {
                inv.d2x_dCdC[indexI1]( i, j, k, l ) = T( 0. );
                inv.d2x_dCdC[indexI2]( i, j, k, l ) = T( delta( i, j ) * delta( k, l ) -
                                                         delta( i, l ) * delta( j, k ) );
                inv.d2x_dCdC[indexJ]( i, j, k, l )  = quarterJ * CInv( j, i ) * CInv( l, k ) -
                                                     halfJ * CInv( j, k ) * CInv( l, i );
              }
      }
    }

    return inv;
  }

  /** @brief Compute the isochoric invariants \f$(\bar I_1, \bar I_2, J)\f$ of \f$\boldsymbol C\f$, i.e. the
   * invariants of \f$\bar{\boldsymbol C} = J^{-2/3}\boldsymbol C\f$ together with \f$J\f$, and their derivatives
   * w.r.t. \f$\boldsymbol C\f$ up to @p Order,
   * \f[
   *   \frac{\partial\bar I_1}{\partial\boldsymbol C} = J^{-2/3}\left(\boldsymbol I -
   *   \frac{I_1}{3}\boldsymbol C^{-T}\right), \qquad
   *   \frac{\partial\bar I_2}{\partial\boldsymbol C} = J^{-4/3}\left(I_1\boldsymbol I - \boldsymbol C^T -
   *   \frac{2 I_2}{3}\boldsymbol C^{-T}\right).
   * \f]
   * The derivatives of \f$\bar I_1\f$ and \f$\bar I_2\f$ are evaluated in this form, in which they vanish
   * exactly at \f$\boldsymbol C = \lambda\boldsymbol I\f$; assembling them from the derivatives of
   * \f$(I_1, I_2, J)\f$ would leave the round-off of two large cancelling terms.
   *
   * @tparam Order Highest derivative to compute (0, 1 or 2).
   * @tparam T Scalar type, e.g. double, std::complex<double>, autodiff scalar.
   * @param C Right Cauchy-Green tensor.
   * @param withI2 Pass @c false for a function that does not depend on \f$\bar I_2\f$, to skip its
   * evaluation; \f$\bar I_2\f$ and its derivatives are then not meaningful.
   * @return The isochoric invariants and their derivatives.
   */
  template < int Order, typename T >
  InvariantsAndDerivatives< T, Order > computeIsochoricInvariants(
    const Math::Tensors::FastorTensors::Tensor33t< T >& C,
    const bool                                          withI2 = true )
  {
    static_assert( Order >= 0 && Order <= 2, "computeIsochoricInvariants: only derivatives up to the second order" );
    using namespace Fastor;

    InvariantsAndDerivatives< T, Order > inv;

    const T I1     = trace( C );
    const T I2     = withI2 ? secondInvariant( C ) : T( 0. );
    const T J      = sqrt( determinant( C ) );
    const T Jm23   = pow( J, -2. / 3. );
    const T Jm43   = Jm23 * Jm23;
    const T Ibar1  = I1 * Jm23;
    const T Ibar2  = I2 * Jm43;
    inv.x[indexI1] = Ibar1;
    inv.x[indexI2] = Ibar2;
    inv.x[indexJ]  = J;

    if constexpr ( Order >= 1 ) {
      const Math::Tensors::FastorTensors::Tensor33t< T > CInv  = inverse( C );
      const T                                            halfJ = 0.5 * J;
      const auto                                         delta = []( int a, int b ) { return a == b ? 1. : 0.; };

      auto& dIbar1 = inv.dx_dC[indexI1];
      auto& dIbar2 = inv.dx_dC[indexI2];
      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ ) {
          dIbar1( i, j )            = Jm23 * ( delta( i, j ) - I1 / 3. * CInv( j, i ) );
          dIbar2( i, j )            = Jm43 * ( I1 * delta( i, j ) - C( j, i ) - 2. / 3. * I2 * CInv( j, i ) );
          inv.dx_dC[indexJ]( i, j ) = halfJ * CInv( j, i );
        }

      if constexpr ( Order >= 2 ) {
        const T quarterJ = 0.25 * J;
        for ( int i = 0; i < 3; i++ )
          for ( int j = 0; j < 3; j++ )
            for ( int k = 0; k < 3; k++ )
              for ( int l = 0; l < 3; l++ ) {
                // C^{-T}_ij = CInv( j, i ) and d C^{-T}_ij / d C_kl = - CInv( j, k ) * CInv( l, i )
                const T CInvT_ij    = CInv( j, i );
                const T CInvT_kl    = CInv( l, k );
                const T dCInvT_ijkl = -CInv( j, k ) * CInv( l, i );

                inv.d2x_dCdC[indexI1]( i, j, k, l ) = -1. / 3. *
                                                        ( dIbar1( i, j ) * CInvT_kl + CInvT_ij * dIbar1( k, l ) ) -
                                                      Ibar1 / 9. * CInvT_ij * CInvT_kl - Ibar1 / 3. * dCInvT_ijkl;

                inv.d2x_dCdC[indexI2]( i, j, k, l ) = -2. / 3. * dIbar2( i, j ) * CInvT_kl +
                                                      Jm43 * ( delta( i, j ) * delta( k, l ) -
                                                               delta( i, l ) * delta( j, k ) -
                                                               2. / 3. * ( I1 * delta( k, l ) - C( l, k ) ) * CInvT_ij -
                                                               2. / 3. * I2 * dCInvT_ijkl );

                inv.d2x_dCdC[indexJ]( i, j, k, l ) = quarterJ * CInvT_ij * CInvT_kl + halfJ * dCInvT_ijkl;
              }
      }
    }

    return inv;
  }

  /** @brief Chain rule from the scalar derivatives of \f$f(\boldsymbol x)\f$ to the derivatives of
   * \f$\Psi(\boldsymbol C) = f(\boldsymbol x(\boldsymbol C))\f$ w.r.t. \f$\boldsymbol C\f$,
   * \f[
   *   \frac{\partial\Psi}{\partial\boldsymbol C} = \sum_a \frac{\partial f}{\partial x_a}
   *   \frac{\partial x_a}{\partial\boldsymbol C}, \qquad
   *   \frac{\partial^2\Psi}{\partial\boldsymbol C\partial\boldsymbol C} = \sum_{a,b}
   *   \frac{\partial^2 f}{\partial x_a\partial x_b}
   *   \frac{\partial x_a}{\partial\boldsymbol C}\otimes\frac{\partial x_b}{\partial\boldsymbol C} + \sum_a
   *   \frac{\partial f}{\partial x_a}\frac{\partial^2 x_a}{\partial\boldsymbol C\partial\boldsymbol C}.
   * \f]
   *
   * @tparam T Scalar type, e.g. double, std::complex<double>, autodiff scalar.
   * @tparam Order Highest derivative to compute (0, 1 or 2).
   * @param f \f$f\f$ and its derivatives w.r.t. the invariants in @p inv.
   * @param inv Invariants of \f$\boldsymbol C\f$ and their derivatives, from #computeInvariants or
   * #computeIsochoricInvariants.
   * @return \f$\Psi\f$ for @p Order 0, the tuple \f$(\Psi, \partial\Psi/\partial\boldsymbol C)\f$ for @p Order 1
   * and the tuple \f$(\Psi, \partial\Psi/\partial\boldsymbol C, \partial^2\Psi/\partial\boldsymbol
   * C\partial\boldsymbol C)\f$ for @p Order 2.
   */
  template < typename T, int Order >
  auto composeWithInvariants( const ScalarDerivatives< T, Order >& f, const InvariantsAndDerivatives< T, Order >& inv )
  {
    static_assert( Order >= 0 && Order <= 2, "composeWithInvariants: only derivatives up to the second order" );

    if constexpr ( Order == 0 ) {
      return f.f;
    }
    else {
      Math::Tensors::FastorTensors::Tensor33t< T > dPsi_dC;
      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ ) {
          T sum = T( 0. );
          for ( std::size_t a = 0; a < 3; a++ )
            sum += f.df[a] * inv.dx_dC[a]( i, j );
          dPsi_dC( i, j ) = sum;
        }

      if constexpr ( Order == 1 ) {
        return std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T > >{ f.f, dPsi_dC };
      }
      else {
        // sum_b d2f_ab dx_b/dC, so that the dyadic part needs one product per entry and invariant
        std::array< Math::Tensors::FastorTensors::Tensor33t< T >, 3 > d2f_dx;
        for ( std::size_t a = 0; a < 3; a++ )
          for ( int k = 0; k < 3; k++ )
            for ( int l = 0; l < 3; l++ ) {
              T sum = T( 0. );
              for ( std::size_t b = 0; b < 3; b++ )
                sum += f.d2f[a][b] * inv.dx_dC[b]( k, l );
              d2f_dx[a]( k, l ) = sum;
            }

        Math::Tensors::FastorTensors::Tensor3333t< T > d2Psi_dCdC;
        for ( int i = 0; i < 3; i++ )
          for ( int j = 0; j < 3; j++ )
            for ( int k = 0; k < 3; k++ )
              for ( int l = 0; l < 3; l++ ) {
                T sum = T( 0. );
                for ( std::size_t a = 0; a < 3; a++ )
                  sum += inv.dx_dC[a]( i, j ) * d2f_dx[a]( k, l ) + f.df[a] * inv.d2x_dCdC[a]( i, j, k, l );
                d2Psi_dCdC( i, j, k, l ) = sum;
              }

        return std::tuple< T,
                           Math::Tensors::FastorTensors::Tensor33t< T >,
                           Math::Tensors::FastorTensors::Tensor3333t< T > >{ f.f, dPsi_dC, d2Psi_dCdC };
      }
    }
  }

} // namespace Marmot::ContinuumMechanics::Invariants::RightCauchyGreen
