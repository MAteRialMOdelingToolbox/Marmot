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
   * @c f is available for every order, @c df for @p Order >= 1 and @c d2f for @p Order >= 2. Entries that are
   * not set are zero.
   *
   * @tparam T Scalar type, e.g. double, std::complex<double> or an automatic differentiation number.
   * @tparam Order Highest derivative that is available (0, 1 or 2).
   */
  template < typename T, int Order >
  struct ScalarDerivatives {
    /// @brief Function value.
    T f{};
    /// @brief Gradient \f$\partial f/\partial x_a\f$ (@p Order >= 1).
    std::conditional_t< ( Order >= 1 ), std::array< T, 3 >, NotComputed > df{};
    /// @brief Hessian \f$\partial^2 f/\partial x_a\partial x_b\f$ (@p Order >= 2).
    std::conditional_t< ( Order >= 2 ), std::array< std::array< T, 3 >, 3 >, NotComputed > d2f{};
  };

  /** @brief Second derivatives of three invariants w.r.t. \f$\boldsymbol C\f$ in compact form.
   *
   * The second derivative of each invariant \f$x_a\f$ of #computeInvariants and #computeIsochoricInvariants is a
   * combination of dyads of the first derivatives and of two fixed fourth-order tensors,
   * \f[
   *   \frac{\partial^2 x_a}{\partial C_{ij}\partial C_{kl}} = \sum_{b,c} \Gamma_{a,bc}\,
   *   \frac{\partial x_b}{\partial C_{ij}}\frac{\partial x_c}{\partial C_{kl}}
   *   + \beta_a\, C^{-1}_{li} C^{-1}_{jk}
   *   + \gamma_a \left(\delta_{ij}\delta_{kl} - \delta_{il}\delta_{jk}\right),
   * \f]
   * so only the scalars \f$\Gamma_{a,bc}\f$, \f$\beta_a\f$ and \f$\gamma_a\f$ are stored. #composeWithInvariants
   * assembles \f$\partial^2\Psi/\partial\boldsymbol C\partial\boldsymbol C\f$ from them in one pass;
   * #secondDerivative returns the full tensor of one invariant.
   *
   * @tparam T Scalar type, e.g. double, std::complex<double> or an automatic differentiation number.
   */
  template < typename T >
  struct SecondDerivativeCoefficients {
    /// @brief \f$\Gamma_{a,bc}\f$, coefficients of the dyads of the first derivatives.
    std::array< std::array< std::array< T, 3 >, 3 >, 3 > dyads{};
    /// @brief \f$\beta_a\f$, coefficient of \f$C^{-1}_{li} C^{-1}_{jk}\f$.
    std::array< T, 3 > inverse{};
    /// @brief \f$\gamma_a\f$, coefficient of \f$\delta_{ij}\delta_{kl} - \delta_{il}\delta_{jk}\f$.
    std::array< T, 3 > identity{};
  };

  /** @brief Three invariants \f$\boldsymbol x\f$ of \f$\boldsymbol C\f$ and their derivatives w.r.t.
   * \f$\boldsymbol C\f$ up to @p Order.
   *
   * @tparam T Scalar type, e.g. double, std::complex<double> or an automatic differentiation number.
   * @tparam Order Highest derivative that is available (0, 1 or 2).
   */
  template < typename T, int Order >
  struct InvariantsAndDerivatives {
    /// @brief The invariants, \f$(I_1, I_2, J)\f$ or \f$(\bar I_1, \bar I_2, J)\f$.
    std::array< T, 3 > x;
    /// @brief First derivatives \f$\partial x_a/\partial\boldsymbol C\f$ (@p Order >= 1).
    std::conditional_t< ( Order >= 1 ), std::array< Math::Tensors::FastorTensors::Tensor33t< T >, 3 >, NotComputed >
      dx_dC;
    /// @brief \f$\boldsymbol C^{-T}\f$ (@p Order >= 2).
    std::conditional_t< ( Order >= 2 ), Math::Tensors::FastorTensors::Tensor33t< T >, NotComputed > CInvT;
    /// @brief Second derivatives \f$\partial^2 x_a/\partial\boldsymbol C\partial\boldsymbol C\f$ in compact form
    /// (@p Order >= 2), see #secondDerivative for the full tensor.
    std::conditional_t< ( Order >= 2 ), SecondDerivativeCoefficients< T >, NotComputed > d2x_dCdC;
  };

  /// @cond
  namespace Detail {
    /* Cofactors cof( i, j ) of C, so that det C = sum_j C( 0, j ) cof( 0, j ), C^-T = cof / det C and
     * I2 = cof( 0, 0 ) + cof( 1, 1 ) + cof( 2, 2 ). Every order of the invariants evaluates these same expressions,
     * so the energy does not depend on the requested order, not even in the last bit. */
    template < typename T >
    T cofactor( const Math::Tensors::FastorTensors::Tensor33t< T >& C, const int i, const int j )
    {
      const int i1 = ( i + 1 ) % 3, i2 = ( i + 2 ) % 3;
      const int j1 = ( j + 1 ) % 3, j2 = ( j + 2 ) % 3;
      return C( i1, j1 ) * C( i2, j2 ) - C( i1, j2 ) * C( i2, j1 );
    }

    template < typename T >
    T determinantFromFirstRowCofactors( const Math::Tensors::FastorTensors::Tensor33t< T >& C,
                                        const T&                                            cof00,
                                        const T&                                            cof01,
                                        const T&                                            cof02 )
    {
      return C( 0, 0 ) * cof00 + C( 0, 1 ) * cof01 + C( 0, 2 ) * cof02;
    }

    /* The quantities shared by both sets of invariants: I1, I2, det C, J and, for Order >= 1, C^-T. */
    template < typename T, int Order >
    struct Basics {
      T                                                                                               I1, I2, detC, J;
      std::conditional_t< ( Order >= 1 ), Math::Tensors::FastorTensors::Tensor33t< T >, NotComputed > CInvT;
    };

    template < int Order, typename T >
    Basics< T, Order > computeBasics( const Math::Tensors::FastorTensors::Tensor33t< T >& C, const bool withI2 )
    {
      Basics< T, Order > b;
      b.I1 = C( 0, 0 ) + C( 1, 1 ) + C( 2, 2 );

      if constexpr ( Order == 0 ) {
        const T cof00 = cofactor( C, 0, 0 );
        b.detC        = determinantFromFirstRowCofactors( C, cof00, cofactor( C, 0, 1 ), cofactor( C, 0, 2 ) );
        b.I2          = withI2 ? T( cof00 + cofactor( C, 1, 1 ) + cofactor( C, 2, 2 ) ) : T( 0. );
      }
      else if constexpr ( std::is_floating_point_v< T > ) {
        for ( int i = 0; i < 3; i++ )
          for ( int j = 0; j < 3; j++ )
            b.CInvT( i, j ) = cofactor( C, i, j );
        b.detC = determinantFromFirstRowCofactors( C, b.CInvT( 0, 0 ), b.CInvT( 0, 1 ), b.CInvT( 0, 2 ) );
        b.I2   = withI2 ? T( b.CInvT( 0, 0 ) + b.CInvT( 1, 1 ) + b.CInvT( 2, 2 ) ) : T( 0. );

        const T invDetC = 1. / b.detC;
        for ( int i = 0; i < 3; i++ )
          for ( int j = 0; j < 3; j++ )
            b.CInvT( i, j ) *= invDetC;
      }
      else {
        // complex and dual numbers: the inverse of Fastor is faster than the cofactors written out
        const T cof00 = cofactor( C, 0, 0 );
        b.detC        = determinantFromFirstRowCofactors( C, cof00, cofactor( C, 0, 1 ), cofactor( C, 0, 2 ) );
        b.I2          = withI2 ? T( cof00 + cofactor( C, 1, 1 ) + cofactor( C, 2, 2 ) ) : T( 0. );
        b.CInvT       = Fastor::transpose( Fastor::inverse( C ) );
      }

      b.J = sqrt( b.detC );
      return b;
    }
  } // namespace Detail
  /// @endcond

  /** @brief Second principal invariant \f$I_2 = \tfrac{1}{2}\left(I_1^2 - \mathrm{tr}\,\boldsymbol C^2\right)\f$,
   * evaluated as the sum of the principal 2x2 minors of \f$\boldsymbol C\f$.
   *
   * @tparam T Scalar type, e.g. double, std::complex<double> or an automatic differentiation number.
   * @param C Right Cauchy-Green tensor.
   * @return \f$I_2\f$.
   */
  template < typename T >
  T secondInvariant( const Math::Tensors::FastorTensors::Tensor33t< T >& C )
  {
    return Detail::cofactor( C, 0, 0 ) + Detail::cofactor( C, 1, 1 ) + Detail::cofactor( C, 2, 2 );
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
   * @tparam T Scalar type, e.g. double, std::complex<double> or an automatic differentiation number.
   * @param C Right Cauchy-Green tensor.
   * @param withI2 Pass @c false for a function that does not depend on \f$I_2\f$, to skip its evaluation;
   * \f$I_2\f$ and its first derivative are then zero and its second derivative is not meaningful.
   * @return The invariants and their derivatives.
   */
  template < int Order, typename T >
  InvariantsAndDerivatives< T, Order > computeInvariants( const Math::Tensors::FastorTensors::Tensor33t< T >& C,
                                                          const bool withI2 = true )
  {
    static_assert( Order >= 0 && Order <= 2, "computeInvariants: only derivatives up to the second order" );

    const auto b = Detail::computeBasics< Order >( C, withI2 );

    InvariantsAndDerivatives< T, Order > inv;
    inv.x[indexI1] = b.I1;
    inv.x[indexI2] = b.I2;
    inv.x[indexJ]  = b.J;

    if constexpr ( Order >= 1 ) {
      const T halfJ = 0.5 * b.J;
      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ ) {
          const double delta_ij      = i == j ? 1. : 0.;
          inv.dx_dC[indexI1]( i, j ) = T( delta_ij );
          inv.dx_dC[indexI2]( i, j ) = withI2 ? T( b.I1 * delta_ij - C( j, i ) ) : T( 0. );
          inv.dx_dC[indexJ]( i, j )  = halfJ * b.CInvT( i, j );
        }

      if constexpr ( Order >= 2 ) {
        // with C^-T = 2 / J dJ/dC: J / 4 C^-T (x) C^-T = 1 / J dJ/dC (x) dJ/dC
        inv.CInvT                                  = b.CInvT;
        inv.d2x_dCdC.identity[indexI2]             = T( 1. );
        inv.d2x_dCdC.dyads[indexJ][indexJ][indexJ] = 1. / b.J;
        inv.d2x_dCdC.inverse[indexJ]               = -halfJ;
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
   * @tparam T Scalar type, e.g. double, std::complex<double> or an automatic differentiation number.
   * @param C Right Cauchy-Green tensor.
   * @param withI2 Pass @c false for a function that does not depend on \f$\bar I_2\f$, to skip its
   * evaluation; \f$\bar I_2\f$ and its first derivative are then zero and its second derivative is not
   * meaningful.
   * @return The isochoric invariants and their derivatives.
   */
  template < int Order, typename T >
  InvariantsAndDerivatives< T, Order > computeIsochoricInvariants(
    const Math::Tensors::FastorTensors::Tensor33t< T >& C,
    const bool                                          withI2 = true )
  {
    static_assert( Order >= 0 && Order <= 2, "computeIsochoricInvariants: only derivatives up to the second order" );

    const auto b = Detail::computeBasics< Order >( C, withI2 );

    const T Jm23  = pow( b.J, -2. / 3. );
    const T Jm43  = Jm23 * Jm23;
    const T Ibar1 = b.I1 * Jm23;
    const T Ibar2 = b.I2 * Jm43;

    InvariantsAndDerivatives< T, Order > inv;
    inv.x[indexI1] = Ibar1;
    inv.x[indexI2] = Ibar2;
    inv.x[indexJ]  = b.J;

    if constexpr ( Order >= 1 ) {
      const T halfJ       = 0.5 * b.J;
      const T I1Third     = b.I1 / 3.;
      const T twoThirdsI2 = 2. / 3. * b.I2;
      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ ) {
          const double delta_ij      = i == j ? 1. : 0.;
          inv.dx_dC[indexI1]( i, j ) = Jm23 * ( delta_ij - I1Third * b.CInvT( i, j ) );
          inv.dx_dC[indexI2]( i, j ) = withI2
                                         ? T( Jm43 * ( b.I1 * delta_ij - C( j, i ) - twoThirdsI2 * b.CInvT( i, j ) ) )
                                         : T( 0. );
          inv.dx_dC[indexJ]( i, j )  = halfJ * b.CInvT( i, j );
        }

      if constexpr ( Order >= 2 ) {
        // dyads with C^-T are written with dJ/dC: C^-T = s dJ/dC, s = 2 / J
        const T s  = 2. / b.J;
        auto&   d2 = inv.d2x_dCdC;
        inv.CInvT  = b.CInvT;

        // d2Ibar1 = -1/3 ( dIbar1 (x) C^-T + C^-T (x) dIbar1 ) - Ibar1/9 C^-T (x) C^-T - Ibar1/3 dC^-T/dC
        d2.dyads[indexI1][indexI1][indexJ] = -s / 3.;
        d2.dyads[indexI1][indexJ][indexI1] = -s / 3.;
        d2.dyads[indexI1][indexJ][indexJ]  = -Ibar1 / 9. * s * s;
        d2.inverse[indexI1]                = Ibar1 / 3.;

        // d2Ibar2 = -2/3 ( dIbar2 (x) C^-T + C^-T (x) dIbar2 ) - 4/9 Ibar2 C^-T (x) C^-T - 2/3 Ibar2 dC^-T/dC
        //           + J^-4/3 d2I2
        d2.dyads[indexI2][indexI2][indexJ] = -2. / 3. * s;
        d2.dyads[indexI2][indexJ][indexI2] = -2. / 3. * s;
        d2.dyads[indexI2][indexJ][indexJ]  = -4. / 9. * Ibar2 * s * s;
        d2.inverse[indexI2]                = 2. / 3. * Ibar2;
        d2.identity[indexI2]               = Jm43;

        d2.dyads[indexJ][indexJ][indexJ] = 1. / b.J;
        d2.inverse[indexJ]               = -halfJ;
      }
    }

    return inv;
  }

  /// @cond
  namespace Detail {
    /* G_ij = sum_a c_a dx_a,ij */
    template < typename T >
    void assembleSecondOrder( Math::Tensors::FastorTensors::Tensor33t< T >&                        G,
                              const std::array< Math::Tensors::FastorTensors::Tensor33t< T >, 3 >& dx_dC,
                              const std::array< T, 3 >&                                            c )
    {
      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          G( i, j ) = c[0] * dx_dC[0]( i, j ) + c[1] * dx_dC[1]( i, j ) + c[2] * dx_dC[2]( i, j );
    }

    /* H_ijkl = cDyads_bc dx_b,ij dx_c,kl + cInverse C^-T_il C^-T_kj + cIdentity ( d_ij d_kl - d_il d_jk ) */
    template < typename T >
    void assembleFourthOrder( Math::Tensors::FastorTensors::Tensor3333t< T >&                      H,
                              const std::array< Math::Tensors::FastorTensors::Tensor33t< T >, 3 >& dx_dC,
                              const Math::Tensors::FastorTensors::Tensor33t< T >&                  CInvT,
                              const std::array< std::array< T, 3 >, 3 >&                           cDyads,
                              const T&                                                             cInverse,
                              const T&                                                             cIdentity )
    {
      // sum_c cDyads_bc dx_c, so that the dyadic part needs one product per entry and invariant
      std::array< Math::Tensors::FastorTensors::Tensor33t< T >, 3 > right;
      for ( std::size_t b = 0; b < 3; b++ )
        for ( int k = 0; k < 3; k++ )
          for ( int l = 0; l < 3; l++ )
            right[b]( k, l ) = cDyads[b][0] * dx_dC[0]( k, l ) + cDyads[b][1] * dx_dC[1]( k, l ) +
                               cDyads[b][2] * dx_dC[2]( k, l );

      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          for ( int k = 0; k < 3; k++ ) {
            const T cInverseCInvT_kj = cInverse * CInvT( k, j );
            for ( int l = 0; l < 3; l++ )
              H( i, j, k, l ) = dx_dC[0]( i, j ) * right[0]( k, l ) + dx_dC[1]( i, j ) * right[1]( k, l ) +
                                dx_dC[2]( i, j ) * right[2]( k, l ) + CInvT( i, l ) * cInverseCInvT_kj;
          }

      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ ) {
          H( i, i, j, j ) += cIdentity;
          H( i, j, j, i ) -= cIdentity;
        }
    }
  } // namespace Detail
  /// @endcond

  /** @brief Full second derivative \f$\partial^2 x_a/\partial\boldsymbol C\partial\boldsymbol C\f$ of one
   * invariant, from the compact form in @p inv.
   *
   * @tparam T Scalar type, e.g. double, std::complex<double> or an automatic differentiation number.
   * @param inv Invariants of \f$\boldsymbol C\f$ and their derivatives, from #computeInvariants or
   * #computeIsochoricInvariants.
   * @param a Index of the invariant (#indexI1, #indexI2 or #indexJ).
   * @return \f$\partial^2 x_a/\partial\boldsymbol C\partial\boldsymbol C\f$.
   */
  template < typename T >
  Math::Tensors::FastorTensors::Tensor3333t< T > secondDerivative( const InvariantsAndDerivatives< T, 2 >& inv,
                                                                   const std::size_t                       a )
  {
    Math::Tensors::FastorTensors::Tensor3333t< T > d2x_dCdC;
    Detail::assembleFourthOrder( d2x_dCdC,
                                 inv.dx_dC,
                                 inv.CInvT,
                                 inv.d2x_dCdC.dyads[a],
                                 inv.d2x_dCdC.inverse[a],
                                 inv.d2x_dCdC.identity[a] );
    return d2x_dCdC;
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
   * The second derivatives of the invariants enter in their compact form (#SecondDerivativeCoefficients), so
   * the fourth-order tensor is assembled in one pass.
   *
   * @tparam T Scalar type, e.g. double, std::complex<double> or an automatic differentiation number.
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
    else if constexpr ( Order == 1 ) {
      std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T > > result;
      std::get< 0 >( result ) = f.f;
      Detail::assembleSecondOrder( std::get< 1 >( result ), inv.dx_dC, f.df );
      return result;
    }
    else {
      // filled in place: the result is not copied on the way out
      std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T >, Math::Tensors::FastorTensors::Tensor3333t< T > >
        result;
      std::get< 0 >( result ) = f.f;
      Detail::assembleSecondOrder( std::get< 1 >( result ), inv.dx_dC, f.df );

      // d2f_bc + sum_a df_a Gamma_a,bc, and the coefficients of the two fixed fourth-order tensors
      std::array< std::array< T, 3 >, 3 > cDyads    = f.d2f;
      T                                   cInverse  = T( 0. );
      T                                   cIdentity = T( 0. );
      for ( std::size_t a = 0; a < 3; a++ ) {
        for ( std::size_t b = 0; b < 3; b++ )
          for ( std::size_t c = 0; c < 3; c++ )
            cDyads[b][c] += f.df[a] * inv.d2x_dCdC.dyads[a][b][c];
        cInverse += f.df[a] * inv.d2x_dCdC.inverse[a];
        cIdentity += f.df[a] * inv.d2x_dCdC.identity[a];
      }
      Detail::assembleFourthOrder( std::get< 2 >( result ), inv.dx_dC, inv.CInvT, cDyads, cInverse, cIdentity );
      return result;
    }
  }

} // namespace Marmot::ContinuumMechanics::Invariants::RightCauchyGreen
