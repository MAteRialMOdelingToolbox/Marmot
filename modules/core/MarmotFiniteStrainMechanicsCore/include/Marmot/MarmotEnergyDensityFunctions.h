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
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotRightCauchyGreenInvariants.h"
#include <array>
#include <concepts>
#include <span>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>

namespace Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions {

  /** @brief Hyperelastic Energy Density Function Wa acc. Pence & Gou (2015), Eq. (2.11)
   *
   *  The energy density function \f$W_a\f$ is given as
   *  \f[
   *    W_a = \frac{G}{2} (I_1 - 3) + \left(\frac{K}{2} - \frac{G}{3}\right) (J - 1)^2 - G \ln(J)
   *  \f]
   *  where \f$ I_1 = \text{tr}(\boldsymbol{C}) \f$ is the first invariant of the right Cauchy-Green tensor
   *  \f$ \boldsymbol{C} = \boldsymbol{F}^T \boldsymbol{F} \f$, \f$ J = \sqrt{\det(\boldsymbol{C})} =
   * \det(\boldsymbol{F}) \f$ is the determinant of the deformation gradient, and \f$ K, G \f$ are the bulk and shear
   * modulus, respectively.
   *
   * @tparam T Scalar type, e.g. double, float, etc.
   * @param C Right Cauchy-Green tensor
   * @param K Bulk modulus
   * @param G Shear modulus
   * @return Energy density
   */
  template < typename T >
  T PenceGouPotentialA( const Math::Tensors::FastorTensors::Tensor33t< T >& C, const double K, const double G )
  {
    using namespace Fastor;
    using namespace Marmot::Math::Tensors::FastorIndices;

    const T J  = sqrt( determinant( C ) );
    const T I1 = trace( C );

    T res = G / 2. * ( I1 - 3. ) + ( K / 2. - G / 3. ) * pow( J - 1, 2 ) - G * log( J );

    return res;
  }

  /** @brief Hyperelastic Energy Density Function Wb acc. Pence & Gou (2015), Eq. (2.12)
   *
   *  The energy density function \f$W_b\f$ is given as
   *  \f[
   *    W_b = \frac{K}{8} \left(J - \frac{1}{J}\right)^2 + \frac{G}{2} \left(I_1 J^{-\frac{2}{3}} - 3\right)
   *  \f]
   *  where \f$ I_1 = \text{tr}(\boldsymbol{C}) \f$ is the first invariant of the right Cauchy-Green tensor
   *  \f$ \boldsymbol{C} = \boldsymbol{F}^T \boldsymbol{F} \f$, \f$ J = \sqrt{\det(\boldsymbol{C})} =
   * \det(\boldsymbol{F}) \f$ is the determinant of the deformation gradient, and \f$ K, G \f$ are the bulk and shear
   * modulus, respectively.
   *
   * @tparam T Scalar type, e.g. double, float, etc.
   * @param C Right Cauchy-Green tensor
   * @param K Bulk modulus
   * @param G Shear modulus
   * @return Energy density
   */
  template < typename T >
  T PenceGouPotentialB( const Math::Tensors::FastorTensors::Tensor33t< T >& C, const double K, const double G )
  {
    using namespace Fastor;
    using namespace Marmot::Math::Tensors::FastorIndices;

    const T detC = determinant( C );
    const T I1   = trace( C );

    T res = K / 8. * ( detC + 1. / detC - 2. ) + G / 2. * ( I1 * pow( detC, -1. / 3 ) - 3. );

    return res;
  }

  /** @brief Hyperelastic Energy Density Function Wc acc. Pence & Gou (2015), Eq. (2.13)
   *
   *  The energy density function \f$W_c\f$ is given as
   *  \f[
   *    W_c = \frac{G}{2} (I_1 - 3) + \frac{3 G^2}{3 K - 2 G} \left(J^{\frac{2}{3} - \frac{K}{G}} - 1\right)
   *  \f]
   *  where \f$ I_1 = \text{tr}(\boldsymbol{C}) \f$ is the first invariant of the right Cauchy-Green tensor
   *  \f$ \boldsymbol{C} = \boldsymbol{F}^T \boldsymbol{F} \f$, \f$ J = \sqrt{\det(\boldsymbol{C})} =
   * \det(\boldsymbol{F}) \f$ is the determinant of the deformation gradient, and \f$ K, G \f$ are the bulk and shear
   * modulus, respectively.
   *
   * @tparam T Scalar type, e.g. double, float, etc.
   * @param C Right Cauchy-Green tensor
   * @param K Bulk modulus
   * @param G Shear modulus
   * @return Energy density
   */
  template < typename T >
  T PenceGouPotentialC( const Math::Tensors::FastorTensors::Tensor33t< T >& C, const double K, const double G )
  {
    using namespace Fastor;
    using namespace Marmot::Math::Tensors::FastorIndices;

    const T J  = sqrt( determinant( C ) );
    const T I1 = trace( C );

    T res = G / 2. * ( I1 - 3. ) + 3. * G * G / ( 3. * K - 2. * G ) * ( pow( J, 2. / 3 - K / G ) - 1 );

    return res;
  }

  /** @cond */
  namespace Detail {
    /* Coupled::StandardNeoHooke as a function of (I1, J), with its scalar derivatives w.r.t. (I1, I2, J):
     * G/2 ( I1 - 3 ) - ( G + lambda/2 ) ln J + lambda/4 ( J^2 - 1 ), lambda = K - 2/3 G. */
    template < int Order, typename T >
    Invariants::RightCauchyGreen::ScalarDerivatives< T, Order > standardNeoHookeOfInvariants(
      const std::array< T, 3 >& x,
      const double              K,
      const double              G )
    {
      [[maybe_unused]] constexpr std::size_t iI1 = Invariants::RightCauchyGreen::indexI1;
      [[maybe_unused]] constexpr std::size_t iJ  = Invariants::RightCauchyGreen::indexJ;

      const double lambda = K - 2.0 / 3.0 * G;
      const double cLnJ   = G + lambda / 2.;
      const T&     I1     = x[Invariants::RightCauchyGreen::indexI1];
      const T&     J      = x[Invariants::RightCauchyGreen::indexJ];

      Invariants::RightCauchyGreen::ScalarDerivatives< T, Order > f;
      f.f = G / 2. * ( I1 - 3. ) - cLnJ * log( J ) + lambda / 4. * ( J * J - 1. );
      if constexpr ( Order >= 1 ) {
        f.df[iI1] = T( G / 2. );
        f.df[iJ]  = lambda / 2. * J - cLnJ / J;
      }
      if constexpr ( Order >= 2 )
        f.d2f[iJ][iJ] = lambda / 2. + cLnJ / ( J * J );
      return f;
    }
  } // namespace Detail
  /** @endcond */

  /** @cond */
  namespace Detail {
    /* Every isochoric potential is one scalar function of the isochoric invariants (Ibar1, Ibar2) and every
     * volumetric potential one scalar function of J, each with its scalar derivatives up to the requested order.
     * The plain, FirstOrder and SecondOrder overloads w.r.t. C all funnel through isochoricPotential() and
     * volumetricPotential(), which leave the tensor algebra to Invariants::RightCauchyGreen, so the overloads
     * cannot drift apart. */

    template < typename T, int Order >
    using ScalarDerivatives = Invariants::RightCauchyGreen::ScalarDerivatives< T, Order >;

    /* Psi( C ) and its derivatives w.r.t. C up to Order for a scalar function
     * fBar: (Ibar1, Ibar2, J) -> ScalarDerivatives w.r.t. (Ibar1, Ibar2, J). UsesIbar2 = false skips the evaluation
     * of Ibar2 for a function of Ibar1 alone. */
    template < int Order, bool UsesIbar2, typename T, typename ScalarFunction >
    auto isochoricPotential( const Math::Tensors::FastorTensors::Tensor33t< T >& C, const ScalarFunction& fBar )
    {
      const auto inv = Invariants::RightCauchyGreen::computeIsochoricInvariants< Order >( C, UsesIbar2 );
      const ScalarDerivatives< T, Order > f = fBar( inv.x );
      return Invariants::RightCauchyGreen::composeWithInvariants( f, inv );
    }

    /* U( C ) and its derivatives w.r.t. C up to Order for a scalar function
     * U: J -> ScalarDerivatives w.r.t. (I1, I2, J), of which only the entries of J are set. */
    template < int Order, typename T, typename ScalarFunction >
    auto volumetricPotential( const Math::Tensors::FastorTensors::Tensor33t< T >& C, const ScalarFunction& U )
    {
      const auto                          inv = Invariants::RightCauchyGreen::computeInvariants< Order >( C, false );
      const ScalarDerivatives< T, Order > f   = U( inv.x[Invariants::RightCauchyGreen::indexJ] );
      return Invariants::RightCauchyGreen::composeWithInvariants( f, inv );
    }

    template < int Order, typename T >
    ScalarDerivatives< T, Order > neoHookeIsochoric( const std::array< T, 3 >& xBar, const double mu )
    {
      ScalarDerivatives< T, Order > f;
      f.f = mu / 2. * ( xBar[0] - 3. );
      if constexpr ( Order >= 1 )
        f.df[0] = T( mu / 2. );
      return f;
    }

    template < int Order, typename T >
    ScalarDerivatives< T, Order > yeohIsochoric( const std::array< T, 3 >& xBar,
                                                 const double              C10,
                                                 const double              C20,
                                                 const double              C30 )
    {
      const T                       Ibar1m3 = xBar[0] - 3.;
      ScalarDerivatives< T, Order > f;
      f.f = C10 * Ibar1m3 + C20 * Ibar1m3 * Ibar1m3 + C30 * Ibar1m3 * Ibar1m3 * Ibar1m3;
      if constexpr ( Order >= 1 )
        f.df[0] = C10 + 2. * C20 * Ibar1m3 + 3. * C30 * Ibar1m3 * Ibar1m3;
      if constexpr ( Order >= 2 )
        f.d2f[0][0] = 2. * C20 + 6. * C30 * Ibar1m3;
      return f;
    }

    template < int Order, typename T >
    ScalarDerivatives< T, Order > mooneyRivlinIsochoric( const std::array< T, 3 >& xBar,
                                                         const double              C10,
                                                         const double              C01 )
    {
      ScalarDerivatives< T, Order > f;
      f.f = C10 * ( xBar[0] - 3. ) + C01 * ( xBar[1] - 3. );
      if constexpr ( Order >= 1 ) {
        f.df[0] = T( C10 );
        f.df[1] = T( C01 );
      }
      return f;
    }

    /** @brief Arruda-Boyce 8-chain potential as a function of the isochoric first invariant
     * \f$\bar I_1\f$ alone.
     *
     * Closed-form Cohen (1991) Pade approximation to the inverse Langevin function,
     * integrated in \f$\bar I_1\f$. Regularized near the locking limit
     * (\f$\bar I_1 \to 3\lambda_L^2\f$, i.e. the average chain stretch approaching
     * the locking stretch): below a small positive floor \f$w_f\f$ of the argument
     * \f$w\f$, \f$\ln w\f$ is continued linearly (\f$C^1\f$) about \f$w_f\f$, so the
     * energy and its derivative stay consistent and any complex-step or dual-number
     * perturbation carried in \f$\bar I_1\f$ survives the regularization.
     * The floor is scaled relative to \f$w_0 = 1 - 1/\lambda_L^2\f$ so that the
     * reference state (\f$\bar I_1 = 3, w = w_0\f$) remains strictly unregularized
     * and exactly zero-energy even for locking stretches arbitrarily close to 1.
     *
     * @tparam Order Highest scalar derivative to compute.
     * @tparam T Scalar type, e.g. double, std::complex<double> or dual3rd.
     * @param xBar Isochoric invariants \f$(\bar I_1, \bar I_2, J)\f$.
     * @param mu Shear-modulus-like parameter (small-strain limit -> standard
     * neo-Hookean shear modulus).
     * @param lambdaL Locking stretch (must be > 1).
     * @return \f$\Psi_{\rm iso}\f$ and its derivatives w.r.t. \f$(\bar I_1, \bar I_2, J)\f$.
     */
    template < int Order, typename T >
    ScalarDerivatives< T, Order > arrudaBoyce8ChainIsochoric( const std::array< T, 3 >& xBar,
                                                              const double              mu,
                                                              const double              lambdaL )
    {
      const T&     Ibar1    = xBar[0];
      const double lambdaL2 = lambdaL * lambdaL;
      const double w0       = 1. - 1. / lambdaL2;
      const T      w        = w0 - ( Ibar1 - 3. ) / ( 3. * lambdaL2 );

      // w0 < 1 for lambdaL > 1, so the floor always lies below the reference state w = w0
      const double wFloor = 1e-2 * w0;

      ScalarDerivatives< T, Order > f;
      if ( Math::Tensors::makeReal( w ) > wFloor ) {
        f.f = mu / 6. * ( Ibar1 - 3. ) - mu * lambdaL2 * log( w / w0 );
        if constexpr ( Order >= 1 )
          f.df[0] = mu / 6. + mu / 3. / w;
        if constexpr ( Order >= 2 )
          f.d2f[0][0] = mu / ( 9. * lambdaL2 ) / ( w * w );
        return f;
      }

      // C1 linear continuation of ln(w) about wFloor: d/dIbar1 of the energy is the constant below
      f.f = mu / 6. * ( Ibar1 - 3. ) - mu * lambdaL2 * ( std::log( wFloor / w0 ) + ( w - wFloor ) / wFloor );
      if constexpr ( Order >= 1 )
        f.df[0] = T( mu / 6. + mu / 3. / wFloor );
      return f;
    }

    /* Volumetric strain energy functions U(J) acc. Doll & Schweizerhof (2000), J. Appl. Mech. 67, 17-21,
     * https://doi.org/10.1115/1.321146, as scalar functions of J with dU/dJ and d2U/dJ2. */

    template < int Order, typename T >
    ScalarDerivatives< T, Order > sussmanBatheVolumetric( const T& J, const double kappa )
    {
      ScalarDerivatives< T, Order >          U;
      [[maybe_unused]] constexpr std::size_t iJ = Invariants::RightCauchyGreen::indexJ;
      U.f                                       = kappa / 2. * ( J - 1. ) * ( J - 1. );
      if constexpr ( Order >= 1 )
        U.df[iJ] = kappa * ( J - 1. );
      if constexpr ( Order >= 2 )
        U.d2f[iJ][iJ] = T( kappa );
      return U;
    }

    template < int Order, typename T >
    ScalarDerivatives< T, Order > penaltyVolumetric( const T& J, const double kappa )
    {
      const T                                lnJ = log( J );
      ScalarDerivatives< T, Order >          U;
      [[maybe_unused]] constexpr std::size_t iJ = Invariants::RightCauchyGreen::indexJ;
      U.f                                       = kappa / 2. * lnJ * lnJ;
      if constexpr ( Order >= 1 )
        U.df[iJ] = kappa * lnJ / J;
      if constexpr ( Order >= 2 )
        U.d2f[iJ][iJ] = kappa * ( 1. - lnJ ) / ( J * J );
      return U;
    }

    template < int Order, typename T >
    ScalarDerivatives< T, Order > simoTaylorVolumetric( const T& J, const double kappa )
    {
      const T                                lnJ = log( J );
      ScalarDerivatives< T, Order >          U;
      [[maybe_unused]] constexpr std::size_t iJ = Invariants::RightCauchyGreen::indexJ;
      U.f                                       = kappa / 4. * ( ( J - 1. ) * ( J - 1. ) + lnJ * lnJ );
      if constexpr ( Order >= 1 )
        U.df[iJ] = kappa / 2. * ( ( J - 1. ) + lnJ / J );
      if constexpr ( Order >= 2 )
        U.d2f[iJ][iJ] = kappa / 2. * ( 1. + ( 1. - lnJ ) / ( J * J ) );
      return U;
    }

    template < int Order, typename T >
    ScalarDerivatives< T, Order > ogdenVolumetric( const T& J, const double kappa, const double theta )
    {
      const T                                JmTheta = pow( J, -theta );
      ScalarDerivatives< T, Order >          U;
      [[maybe_unused]] constexpr std::size_t iJ = Invariants::RightCauchyGreen::indexJ;
      U.f                                       = kappa / ( theta * theta ) * ( theta * log( J ) + JmTheta - 1. );
      if constexpr ( Order >= 1 )
        U.df[iJ] = kappa / theta * ( 1. - JmTheta ) / J;
      if constexpr ( Order >= 2 )
        U.d2f[iJ][iJ] = kappa / theta * ( ( theta + 1. ) * JmTheta - 1. ) / ( J * J );
      return U;
    }

    template < int Order, typename T >
    ScalarDerivatives< T, Order > liuHofstetterMangVolumetric( const T& J, const double kappa )
    {
      const T                                lnJ = log( J );
      ScalarDerivatives< T, Order >          U;
      [[maybe_unused]] constexpr std::size_t iJ = Invariants::RightCauchyGreen::indexJ;
      U.f                                       = kappa * ( J * lnJ - J + 1. );
      if constexpr ( Order >= 1 )
        U.df[iJ] = kappa * lnJ;
      if constexpr ( Order >= 2 )
        U.d2f[iJ][iJ] = kappa / J;
      return U;
    }

    template < int Order, typename T >
    ScalarDerivatives< T, Order > generalizedDollSchweizerhofVolumetric( const T&                        J,
                                                                         const double                    kappa,
                                                                         const std::span< const double > alphas,
                                                                         const std::span< const double > betas )
    {
      if ( alphas.size() != betas.size() || alphas.empty() )
        throw std::invalid_argument( "Generalized Doll-Schweizerhof volumetric potential: alphas and betas must be "
                                     "non-empty and of equal size." );

      double sumAlphaPlusBeta = 0.;
      for ( std::size_t i = 0; i < alphas.size(); i++ )
        sumAlphaPlusBeta += alphas[i] + betas[i];
      const double scale = kappa / sumAlphaPlusBeta;

      ScalarDerivatives< T, Order >          U;
      [[maybe_unused]] constexpr std::size_t iJ = Invariants::RightCauchyGreen::indexJ;
      for ( std::size_t i = 0; i < alphas.size(); i++ ) {
        const double a = alphas[i];
        const double b = betas[i];
        // all powers of J from J^a and J^-b
        const T powJa      = pow( J, a );
        const T powJMinusB = pow( J, -b );
        U.f += scale * ( powJa * J / ( a + 1. ) + powJMinusB * J / ( b - 1. ) - 1. / ( a + 1. ) - 1. / ( b - 1. ) );
        if constexpr ( Order >= 1 )
          U.df[iJ] += scale * ( powJa - powJMinusB );
        if constexpr ( Order >= 2 )
          U.d2f[iJ][iJ] += scale * ( a * powJa + b * powJMinusB ) / J;
      }
      return U;
    }

    template < int Order, typename T >
    ScalarDerivatives< T, Order > dollSchweizerhofExponentialVolumetric( const T& J, const double kappa )
    {
      const T                                expJm1 = exp( J - 1. );
      ScalarDerivatives< T, Order >          U;
      [[maybe_unused]] constexpr std::size_t iJ = Invariants::RightCauchyGreen::indexJ;
      U.f                                       = kappa / 2. * ( expJm1 - log( J ) - 1. );
      if constexpr ( Order >= 1 )
        U.df[iJ] = kappa / 2. * ( expJm1 - 1. / J );
      if constexpr ( Order >= 2 )
        U.d2f[iJ][iJ] = kappa / 2. * ( expJm1 + 1. / ( J * J ) );
      return U;
    }

    template < int Order, typename T >
    ScalarDerivatives< T, Order > dollSchweizerhofLogarithmicVolumetric( const T& J, const double kappa )
    {
      const T                                lnJ = log( J );
      ScalarDerivatives< T, Order >          U;
      [[maybe_unused]] constexpr std::size_t iJ = Invariants::RightCauchyGreen::indexJ;
      U.f                                       = kappa / 2. * ( J - 1. ) * lnJ;
      if constexpr ( Order >= 1 )
        U.df[iJ] = kappa / 2. * ( lnJ + ( J - 1. ) / J );
      if constexpr ( Order >= 2 )
        U.d2f[iJ][iJ] = kappa / 2. * ( 1. / J + 1. / ( J * J ) );
      return U;
    }
  } // namespace Detail
  /** @endcond */

  /** @brief Scalar function of the isochoric invariants \f$(\bar I_1, \bar I_2, J)\f$ with its scalar derivatives,
   * e.g. Isochoric::Yeoh: the isochoric part of the single-pass @c DecoupledPotential. A user-defined type qualifies
   * with a static member @c usesIbar2 and a member template
   * @c scalarDerivatives<Order>( xBar ) returning Invariants::RightCauchyGreen::ScalarDerivatives. */
  template < typename Part >
  concept IsochoricScalarFunction = requires( const Part& part, const std::array< double, 3 >& xBar ) {
    {
      Part::usesIbar2
    } -> std::convertible_to< bool >;
    {
      part.template scalarDerivatives< 0 >( xBar )
    } -> std::same_as< Invariants::RightCauchyGreen::ScalarDerivatives< double, 0 > >;
  };

  /** @brief Scalar function \f$U(J)\f$ with its scalar derivatives, e.g. Volumetric::SimoTaylor: the volumetric part
   * of the single-pass @c DecoupledPotential. A user-defined type qualifies with a member template
   * @c scalarDerivatives<Order>( J ) returning Invariants::RightCauchyGreen::ScalarDerivatives with the entries of
   * \f$J\f$ set. */
  template < typename Part >
  concept VolumetricScalarFunction = requires( const Part& part, const double& J ) {
    {
      part.template scalarDerivatives< 0 >( J )
    } -> std::same_as< Invariants::RightCauchyGreen::ScalarDerivatives< double, 0 > >;
  };

  namespace Isochoric {
    /** @brief Scalar function of the isochoric neo-Hooke potential, see #NeoHookePotential; for
     * Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::DecoupledPotential. */
    struct NeoHooke {
      double mu; ///< Shear modulus.

      /// Whether the function depends on the second isochoric invariant.
      static constexpr bool usesIbar2 = false;

      /** @brief Value and derivatives up to @c Order w.r.t. \f$(\bar I_1, \bar I_2, J)\f$.
       * @param xBar Isochoric invariants \f$(\bar I_1, \bar I_2, J)\f$.
       * @return Value and scalar derivatives. */
      template < int Order, typename T >
      Invariants::RightCauchyGreen::ScalarDerivatives< T, Order > scalarDerivatives(
        const std::array< T, 3 >& xBar ) const
      {
        return Detail::neoHookeIsochoric< Order >( xBar, mu );
      }
    };

    /** @brief Scalar function of the isochoric Yeoh potential, see #YeohPotential; for
     * Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::DecoupledPotential. */
    struct Yeoh {
      double C10; ///< Parameter \f$C_{10}\f$.
      double C20; ///< Parameter \f$C_{20}\f$.
      double C30; ///< Parameter \f$C_{30}\f$.

      /// Whether the function depends on the second isochoric invariant.
      static constexpr bool usesIbar2 = false;

      /** @brief Value and derivatives up to @c Order w.r.t. \f$(\bar I_1, \bar I_2, J)\f$.
       * @param xBar Isochoric invariants \f$(\bar I_1, \bar I_2, J)\f$.
       * @return Value and scalar derivatives. */
      template < int Order, typename T >
      Invariants::RightCauchyGreen::ScalarDerivatives< T, Order > scalarDerivatives(
        const std::array< T, 3 >& xBar ) const
      {
        return Detail::yeohIsochoric< Order >( xBar, C10, C20, C30 );
      }
    };

    /** @brief Scalar function of the isochoric Mooney-Rivlin potential, see #MooneyRivlinPotential; for
     * Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::DecoupledPotential. */
    struct MooneyRivlin {
      double C10; ///< Parameter \f$C_{10}\f$.
      double C01; ///< Parameter \f$C_{01}\f$.

      /// Whether the function depends on the second isochoric invariant.
      static constexpr bool usesIbar2 = true;

      /** @brief Value and derivatives up to @c Order w.r.t. \f$(\bar I_1, \bar I_2, J)\f$.
       * @param xBar Isochoric invariants \f$(\bar I_1, \bar I_2, J)\f$.
       * @return Value and scalar derivatives. */
      template < int Order, typename T >
      Invariants::RightCauchyGreen::ScalarDerivatives< T, Order > scalarDerivatives(
        const std::array< T, 3 >& xBar ) const
      {
        return Detail::mooneyRivlinIsochoric< Order >( xBar, C10, C01 );
      }
    };

    /** @brief Scalar function of the isochoric Arruda-Boyce 8-chain potential, see #ArrudaBoyce8ChainPotential; for
     * Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::DecoupledPotential. */
    struct ArrudaBoyce8Chain {
      double mu;      ///< Shear-modulus-like parameter.
      double lambdaL; ///< Locking stretch.

      /// Whether the function depends on the second isochoric invariant.
      static constexpr bool usesIbar2 = false;

      /** @brief Value and derivatives up to @c Order w.r.t. \f$(\bar I_1, \bar I_2, J)\f$.
       * @param xBar Isochoric invariants \f$(\bar I_1, \bar I_2, J)\f$.
       * @return Value and scalar derivatives. */
      template < int Order, typename T >
      Invariants::RightCauchyGreen::ScalarDerivatives< T, Order > scalarDerivatives(
        const std::array< T, 3 >& xBar ) const
      {
        return Detail::arrudaBoyce8ChainIsochoric< Order >( xBar, mu, lambdaL );
      }
    };
  } // namespace Isochoric

  namespace Volumetric {
    /** @brief Scalar function \f$U(J)\f$ of the penalty volumetric potential, see #PenaltyPotential; for
     * Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::DecoupledPotential. */
    struct Penalty {
      double kappa; ///< Bulk modulus.

      /** @brief Value and derivatives up to @c Order w.r.t. \f$J\f$, stored in the entries of \f$J\f$.
       * @param J Determinant of the deformation gradient.
       * @return Value and scalar derivatives. */
      template < int Order, typename T >
      Invariants::RightCauchyGreen::ScalarDerivatives< T, Order > scalarDerivatives( const T& J ) const
      {
        return Detail::penaltyVolumetric< Order >( J, kappa );
      }
    };

    /** @brief Scalar function \f$U(J)\f$ of the Sussman-Bathe volumetric potential, see #SussmanBathePotential; for
     * Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::DecoupledPotential. */
    struct SussmanBathe {
      double kappa; ///< Bulk modulus.

      /** @brief Value and derivatives up to @c Order w.r.t. \f$J\f$, stored in the entries of \f$J\f$.
       * @param J Determinant of the deformation gradient.
       * @return Value and scalar derivatives. */
      template < int Order, typename T >
      Invariants::RightCauchyGreen::ScalarDerivatives< T, Order > scalarDerivatives( const T& J ) const
      {
        return Detail::sussmanBatheVolumetric< Order >( J, kappa );
      }
    };

    /** @brief Scalar function \f$U(J)\f$ of the Simo-Taylor volumetric potential, see #SimoTaylorPotential; for
     * Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::DecoupledPotential. */
    struct SimoTaylor {
      double kappa; ///< Bulk modulus.

      /** @brief Value and derivatives up to @c Order w.r.t. \f$J\f$, stored in the entries of \f$J\f$.
       * @param J Determinant of the deformation gradient.
       * @return Value and scalar derivatives. */
      template < int Order, typename T >
      Invariants::RightCauchyGreen::ScalarDerivatives< T, Order > scalarDerivatives( const T& J ) const
      {
        return Detail::simoTaylorVolumetric< Order >( J, kappa );
      }
    };

    /** @brief Scalar function \f$U(J)\f$ of the Ogden volumetric potential, see #OgdenPotential; for
     * Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::DecoupledPotential. */
    struct Ogden {
      double kappa; ///< Bulk modulus.
      double theta; ///< Exponent \f$\theta\f$.

      /** @brief Value and derivatives up to @c Order w.r.t. \f$J\f$, stored in the entries of \f$J\f$.
       * @param J Determinant of the deformation gradient.
       * @return Value and scalar derivatives. */
      template < int Order, typename T >
      Invariants::RightCauchyGreen::ScalarDerivatives< T, Order > scalarDerivatives( const T& J ) const
      {
        return Detail::ogdenVolumetric< Order >( J, kappa, theta );
      }
    };

    /** @brief Scalar function \f$U(J)\f$ of the Liu-Hofstetter-Mang volumetric potential, see
     * #LiuHofstetterMangPotential; for
     * Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::DecoupledPotential. */
    struct LiuHofstetterMang {
      double kappa; ///< Bulk modulus.

      /** @brief Value and derivatives up to @c Order w.r.t. \f$J\f$, stored in the entries of \f$J\f$.
       * @param J Determinant of the deformation gradient.
       * @return Value and scalar derivatives. */
      template < int Order, typename T >
      Invariants::RightCauchyGreen::ScalarDerivatives< T, Order > scalarDerivatives( const T& J ) const
      {
        return Detail::liuHofstetterMangVolumetric< Order >( J, kappa );
      }
    };

    /** @brief Scalar function \f$U(J)\f$ of the two-parameter Doll-Schweizerhof volumetric potential, see
     * #DollSchweizerhofPotential; for
     * Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::DecoupledPotential. */
    struct DollSchweizerhof {
      double kappa; ///< Bulk modulus.
      double alpha; ///< Shape parameter \f$\alpha > 0\f$.
      double beta;  ///< Shape parameter \f$\beta > 1\f$.

      /** @brief Value and derivatives up to @c Order w.r.t. \f$J\f$, stored in the entries of \f$J\f$.
       * @param J Determinant of the deformation gradient.
       * @return Value and scalar derivatives. */
      template < int Order, typename T >
      Invariants::RightCauchyGreen::ScalarDerivatives< T, Order > scalarDerivatives( const T& J ) const
      {
        const double alphas[] = { alpha };
        const double betas[]  = { beta };
        return Detail::generalizedDollSchweizerhofVolumetric< Order >( J, kappa, alphas, betas );
      }
    };

    /** @brief Scalar function \f$U(J)\f$ of the generalized Doll-Schweizerhof volumetric potential, see
     * #GeneralizedDollSchweizerhofPotential; for
     * Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::DecoupledPotential. */
    struct GeneralizedDollSchweizerhof {
      double                    kappa;  ///< Bulk modulus.
      std::span< const double > alphas; ///< Shape parameters \f$\alpha_i > 0\f$ (not owned).
      std::span< const double > betas;  ///< Shape parameters \f$\beta_i > 1\f$ (not owned).

      /** @brief Value and derivatives up to @c Order w.r.t. \f$J\f$, stored in the entries of \f$J\f$.
       * @param J Determinant of the deformation gradient.
       * @return Value and scalar derivatives. */
      template < int Order, typename T >
      Invariants::RightCauchyGreen::ScalarDerivatives< T, Order > scalarDerivatives( const T& J ) const
      {
        return Detail::generalizedDollSchweizerhofVolumetric< Order >( J, kappa, alphas, betas );
      }
    };

    /** @brief Scalar function \f$U(J)\f$ of the exponential Doll-Schweizerhof volumetric potential, see
     * #DollSchweizerhofExponentialPotential; for
     * Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::DecoupledPotential. */
    struct DollSchweizerhofExponential {
      double kappa; ///< Bulk modulus.

      /** @brief Value and derivatives up to @c Order w.r.t. \f$J\f$, stored in the entries of \f$J\f$.
       * @param J Determinant of the deformation gradient.
       * @return Value and scalar derivatives. */
      template < int Order, typename T >
      Invariants::RightCauchyGreen::ScalarDerivatives< T, Order > scalarDerivatives( const T& J ) const
      {
        return Detail::dollSchweizerhofExponentialVolumetric< Order >( J, kappa );
      }
    };

    /** @brief Scalar function \f$U(J)\f$ of the logarithmic Doll-Schweizerhof volumetric potential, see
     * #DollSchweizerhofLogarithmicPotential; for
     * Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::DecoupledPotential. */
    struct DollSchweizerhofLogarithmic {
      double kappa; ///< Bulk modulus.

      /** @brief Value and derivatives up to @c Order w.r.t. \f$J\f$, stored in the entries of \f$J\f$.
       * @param J Determinant of the deformation gradient.
       * @return Value and scalar derivatives. */
      template < int Order, typename T >
      Invariants::RightCauchyGreen::ScalarDerivatives< T, Order > scalarDerivatives( const T& J ) const
      {
        return Detail::dollSchweizerhofLogarithmicVolumetric< Order >( J, kappa );
      }
    };
  } // namespace Volumetric

  /** @brief Scalar function of the principal invariants \f$(I_1, I_2, J)\f$ with its scalar derivatives, e.g.
   * Coupled::StandardNeoHooke: a potential that is not the sum of an isochoric and a volumetric part, for
   * @c CoupledPotential. A user-defined type qualifies with a static member @c usesI2 and a member template
   * @c scalarDerivatives<Order>( x ) returning Invariants::RightCauchyGreen::ScalarDerivatives. */
  template < typename Potential >
  concept PrincipalInvariantScalarFunction = requires( const Potential& potential, const std::array< double, 3 >& x ) {
    {
      Potential::usesI2
    } -> std::convertible_to< bool >;
    {
      potential.template scalarDerivatives< 0 >( x )
    } -> std::same_as< Invariants::RightCauchyGreen::ScalarDerivatives< double, 0 > >;
  };

  /** @brief Hyperelastic potentials as scalar functions of the principal invariants \f$(I_1, I_2, J)\f$ of
   * \f$\boldsymbol C\f$, which do not split into an isochoric and a volumetric part; evaluate them with
   * Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::CoupledPotential. */
  namespace Coupled {
    /** @brief Standard compressible Neo-Hooke potential
     *  \f[
     *    \Psi = \frac{G}{2} \left(I_1 - 3 - 2 \ln J\right) + \frac{\lambda}{2} \left(\frac{J^2 - 1}{2} - \ln J\right),
     *    \qquad \lambda = K - \frac{2}{3} G,
     *  \f]
     *  with \f$ I_1 = \text{tr}(\boldsymbol{C}) \f$ and \f$ J = \sqrt{\det(\boldsymbol{C})} \f$. It depends on
     *  \f$I_1\f$ and not on the isochoric invariant \f$\bar I_1\f$. */
    struct StandardNeoHooke {
      double K; ///< Bulk modulus.
      double G; ///< Shear modulus.

      /// Whether the function depends on the second invariant.
      static constexpr bool usesI2 = false;

      /** @brief Value and derivatives up to @c Order w.r.t. \f$(I_1, I_2, J)\f$.
       * @param x Principal invariants \f$(I_1, I_2, J)\f$.
       * @return Value and scalar derivatives. */
      template < int Order, typename T >
      Invariants::RightCauchyGreen::ScalarDerivatives< T, Order > scalarDerivatives( const std::array< T, 3 >& x ) const
      {
        return Detail::standardNeoHookeOfInvariants< Order >( x, K, G );
      }
    };
  } // namespace Coupled

  /** @cond */
  namespace Detail {
    /* Psi( I1, I2, J ) and its derivatives w.r.t. C up to Order. */
    template < int Order, typename T, PrincipalInvariantScalarFunction Potential >
    auto coupledPotential( const Math::Tensors::FastorTensors::Tensor33t< T >& C, const Potential& potential )
    {
      const auto inv = Invariants::RightCauchyGreen::computeInvariants< Order >( C, Potential::usesI2 );
      const ScalarDerivatives< T, Order > f = potential.template scalarDerivatives< Order >( inv.x );
      return Invariants::RightCauchyGreen::composeWithInvariants( f, inv );
    }
  } // namespace Detail
  /** @endcond */

  /** @cond */
  namespace Detail {
    /* Psi_iso + U and its derivatives w.r.t. C up to Order: the scalar derivatives of both parts are added
     * before the chain rule, so the invariants of C are computed once and the tensors assembled once. */
    template < int Order, typename T, IsochoricScalarFunction IsochoricPart, VolumetricScalarFunction VolumetricPart >
    auto decoupledPotential( const Math::Tensors::FastorTensors::Tensor33t< T >& C,
                             const IsochoricPart&                                isochoricPart,
                             const VolumetricPart&                               volumetricPart )
    {
      [[maybe_unused]] constexpr std::size_t iJ = Invariants::RightCauchyGreen::indexJ;

      const auto inv = Invariants::RightCauchyGreen::computeIsochoricInvariants< Order >( C, IsochoricPart::usesIbar2 );
      ScalarDerivatives< T, Order >       f = isochoricPart.template scalarDerivatives< Order >( inv.x );
      const ScalarDerivatives< T, Order > U = volumetricPart.template scalarDerivatives< Order >(
        inv.x[Invariants::RightCauchyGreen::indexJ] );
      f.f += U.f;
      if constexpr ( Order >= 1 )
        f.df[iJ] += U.df[iJ];
      if constexpr ( Order >= 2 )
        f.d2f[iJ][iJ] += U.d2f[iJ][iJ];
      return Invariants::RightCauchyGreen::composeWithInvariants( f, inv );
    }
  } // namespace Detail
  /** @endcond */

  /** @brief Purely isochoric (volume-preserving) hyperelastic energy density functions
   * \f$\Psi_{\rm iso}(\bar{\boldsymbol C})\f$, \f$\bar{\boldsymbol C} = J^{-2/3}\boldsymbol C\f$.
   *
   * They carry no volumetric stiffness; combine them with a potential from
   * Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::Volumetric via #DecoupledPotential or
   * #makeDecoupledPotential.
   */
  namespace Isochoric {

    /** @brief Isochoric Mooney-Rivlin hyperelastic energy density function
     * (no volumetric term; combine with a Volumetric potential, see #DecoupledPotential).
     *
     *  The isochoric energy density is
     *  \f[
     *    W_{\rm iso} = C_{10} (\bar{I}_1 - 3) + C_{01} (\bar{I}_2 - 3)
     *  \f]
     *  where \f$ \bar{I}_1 = I_1 J^{-2/3} \f$ and \f$ \bar{I}_2 = I_2 J^{-4/3} \f$ are the first
     *  and second invariant of the isochoric right Cauchy-Green tensor, respectively, \f$ I_1 =
     *  \text{tr}(\boldsymbol{C}) \f$ and \f$ I_2 = 0.5 (I_1^2 - \text{tr}(\boldsymbol{C}^2)) \f$ are the first and
     *  second invariant of the right Cauchy-Green tensor \f$ \boldsymbol{C} = \boldsymbol{F}^T \boldsymbol{F} \f$, \f$
     * J = \sqrt{\det(\boldsymbol{C})} = \det(\boldsymbol{F}) \f$ is the determinant of the deformation gradient.
     *
     * @tparam T Scalar type, e.g. double, float, etc.
     * @param C Right Cauchy-Green tensor
     * @param C10 Mooney-Rivlin material parameter C10
     * @param C01 Mooney-Rivlin material parameter C01
     * @return Energy density (isochoric part only).
     */
    template < typename T >
    T MooneyRivlinPotential( const Math::Tensors::FastorTensors::Tensor33t< T >& C, const double C10, const double C01 )
    {
      return Detail::isochoricPotential< 0, true >( C, [C10, C01]( const std::array< T, 3 >& xBar ) {
        return Detail::mooneyRivlinIsochoric< 0 >( xBar, C10, C01 );
      } );
    }

    /** @brief Isochoric Yeoh hyperelastic energy density function (no volumetric
     * term; combine with a Volumetric potential, see #DecoupledPotential).
     *
     *  The isochoric energy density is
     *  \f[
     *    W_{\rm iso} = C_{10} (\bar{I}_1 - 3) + C_{20} (\bar{I}_1 - 3)^2 + C_{30} (\bar{I}_1 - 3)^3
     *  \f]
     *  where \f$\bar{I}_1 = I_1 J^{-2/3}\f$, \f$I_1 = \mathrm{tr}(\boldsymbol{C})\f$, and
     *  \f$J = \sqrt{\det(\boldsymbol{C})}\f$.
     *
     * @tparam T Scalar type, e.g. double, float or an automatic differentiation number.
     * @param C Right Cauchy-Green tensor.
     * @param C10 Yeoh material parameter.
     * @param C20 Yeoh material parameter.
     * @param C30 Yeoh material parameter.
     * @return Energy density (isochoric part only).
     */
    template < typename T >
    T YeohPotential( const Math::Tensors::FastorTensors::Tensor33t< T >& C,
                     const double                                        C10,
                     const double                                        C20,
                     const double                                        C30 )
    {
      return Detail::isochoricPotential< 0, false >( C, [C10, C20, C30]( const std::array< T, 3 >& xBar ) {
        return Detail::yeohIsochoric< 0 >( xBar, C10, C20, C30 );
      } );
    }

    /** @brief Isochoric neo-Hookean hyperelastic energy density function
     * (no volumetric term; combine with Volumetric::PenaltyPotential
     * for the shared volumetric term),
     * acc. Pence & Gou (2015).
     *
     * The isochoric energy density is
     * \f[
     *   W_{\rm iso} = \frac{\mu}{2}\left(\bar I_1 - 3\right)
     * \f]
     * where \f$\bar I_1 = I_1 J^{-2/3}\f$, \f$I_1 = \mathrm{tr}(\boldsymbol C)\f$, and
     * \f$J = \sqrt{\det\boldsymbol C}\f$.
     *
     * @tparam T Scalar type, e.g. double, float or an automatic differentiation number.
     * @param C Right Cauchy-Green tensor.
     * @param mu Shear modulus.
     * @return Energy density (isochoric part only).
     */
    template < typename T >
    T NeoHookePotential( const Math::Tensors::FastorTensors::Tensor33t< T >& C, const double mu )
    {
      return Detail::isochoricPotential< 0, false >( C, [mu]( const std::array< T, 3 >& xBar ) {
        return Detail::neoHookeIsochoric< 0 >( xBar, mu );
      } );
    }

    /** @brief Isochoric Arruda-Boyce 8-chain hyperelastic energy density function (isochoric
     * part only -- no volumetric term, since the classical 8-chain derivation is purely
     * isochoric; combine with e.g. Volumetric::PenaltyPotential),
     * via the closed-form Cohen (1991) Pade approximation to the inverse Langevin
     * function.
     *
     * The isochoric energy density is
     * \f[
     *   W_{AB} = \frac{\mu}{6}\left(\bar I_1 - 3\right) - \mu \lambda_L^2 \ln\left(
     *   \frac{1 - \bar I_1/(3\lambda_L^2)}{1 - 1/\lambda_L^2} \right)
     * \f]
     * where \f$\bar I_1 = I_1 J^{-2/3}\f$, \f$I_1 = \mathrm{tr}(\boldsymbol C)\f$,
     * \f$J=\sqrt{\det\boldsymbol C}\f$, \f$\mu\f$ is the shear-modulus-like parameter,
     * and \f$\lambda_L\f$ is the locking stretch. Reduces exactly to the isochoric
     * part of Coupled::StandardNeoHooke as \f$\lambda_L\to\infty\f$, and is exactly
     * stress-free (zero energy and zero gradient) at \f$\boldsymbol C=\boldsymbol I\f$
     * for any finite \f$\lambda_L\f$, since \f$\bar I_1\f$ is invariant under
     * \f$\boldsymbol C \to \lambda\boldsymbol C\f$ for any scalar \f$\lambda\f$.
     *
     * @tparam T Scalar type, e.g. double, float or an automatic differentiation number.
     * @param C Right Cauchy-Green tensor.
     * @param mu Shear-modulus-like parameter.
     * @param lambdaL Locking stretch (must be > 1; the model has a genuine
     * physical/numerical singularity as the average chain stretch approaches
     * \f$\lambda_L\f$).
     * @return Energy density (isochoric part only).
     */
    template < typename T >
    T ArrudaBoyce8ChainPotential( const Math::Tensors::FastorTensors::Tensor33t< T >& C,
                                  const double                                        mu,
                                  const double                                        lambdaL )
    {
      return Detail::isochoricPotential< 0, false >( C, [mu, lambdaL]( const std::array< T, 3 >& xBar ) {
        return Detail::arrudaBoyce8ChainIsochoric< 0 >( xBar, mu, lambdaL );
      } );
    }
    namespace FirstOrder {

      /** @brief Isochoric neo-Hookean potential (see the plain
       * Isochoric::NeoHookePotential, isochoric part only) and its first derivative w.r.t. C,
       * \f[
       *   \frac{\partial W_{\rm iso}}{\partial \boldsymbol C} = \frac{\mu}{2}\frac{\partial \bar I_1}{\partial
       *   \boldsymbol C}, \qquad \frac{\partial \bar I_1}{\partial \boldsymbol C} = J^{-2/3}\boldsymbol I -
       *   \frac{\bar I_1}{3}\boldsymbol C^{-1}.
       * \f]
       *
       * @tparam T Scalar type, e.g. double, std::complex<double>.
       * @param C Right Cauchy-Green tensor.
       * @param mu Shear modulus.
       * @return A tuple containing the (isochoric) energy density and its first derivative w.r.t. C.
       */
      template < typename T >
      std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T > > NeoHookePotential(
        const Math::Tensors::FastorTensors::Tensor33t< T >& C,
        const double                                        mu )
      {
        return Detail::isochoricPotential< 1, false >( C, [mu]( const std::array< T, 3 >& xBar ) {
          return Detail::neoHookeIsochoric< 1 >( xBar, mu );
        } );
      }

      /** @brief Isochoric Yeoh potential (see the plain Isochoric::YeohPotential, isochoric part
       * only) and its first derivative w.r.t. C.
       *
       * @tparam T Scalar type, e.g. double, std::complex<double>.
       * @param C Right Cauchy-Green tensor.
       * @param C10 Yeoh material parameter.
       * @param C20 Yeoh material parameter.
       * @param C30 Yeoh material parameter.
       * @return A tuple containing the (isochoric) energy density and its first derivative w.r.t. C.
       */
      template < typename T >
      std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T > > YeohPotential(
        const Math::Tensors::FastorTensors::Tensor33t< T >& C,
        const double                                        C10,
        const double                                        C20,
        const double                                        C30 )
      {
        return Detail::isochoricPotential< 1, false >( C, [C10, C20, C30]( const std::array< T, 3 >& xBar ) {
          return Detail::yeohIsochoric< 1 >( xBar, C10, C20, C30 );
        } );
      }

      /** @brief Isochoric Mooney-Rivlin potential (see the plain Isochoric::MooneyRivlinPotential,
       * isochoric part only) and its first derivative w.r.t. C.
       *
       * @tparam T Scalar type, e.g. double, std::complex<double>.
       * @param C Right Cauchy-Green tensor.
       * @param C10 Mooney-Rivlin material parameter.
       * @param C01 Mooney-Rivlin material parameter.
       * @return A tuple containing the (isochoric) energy density and its first derivative w.r.t. C.
       */
      template < typename T >
      std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T > > MooneyRivlinPotential(
        const Math::Tensors::FastorTensors::Tensor33t< T >& C,
        const double                                        C10,
        const double                                        C01 )
      {
        return Detail::isochoricPotential< 1, true >( C, [C10, C01]( const std::array< T, 3 >& xBar ) {
          return Detail::mooneyRivlinIsochoric< 1 >( xBar, C10, C01 );
        } );
      }

      /** @brief Isochoric Arruda-Boyce 8-chain hyperelastic energy density function (isochoric
       * part only) and its first derivative w.r.t. C -- see the plain
       * Isochoric::ArrudaBoyce8ChainPotential for the formula and its reduction/stress-free
       * properties. This overload additionally returns
       * \f[
       *   \frac{\partial W_{AB}}{\partial \boldsymbol C} = \frac{\partial
       *   W_{AB}}{\partial \bar I_1} \frac{\partial \bar I_1}{\partial \boldsymbol C},
       *   \qquad \frac{\partial \bar I_1}{\partial \boldsymbol C} =
       *   J^{-2/3}\boldsymbol I - \frac{\bar I_1}{3}\boldsymbol C^{-1}
       * \f]
       * (the standard isochoric-invariant derivative identity -- note \f$J^{-2/3}\f$
       * multiplies ONLY the \f$\boldsymbol I\f$ term, not the \f$\boldsymbol C^{-1}\f$
       * term, since the latter's scaling is already absorbed into \f$\bar I_1\f$
       * itself), sharing the exact same closed-form \f$\Psi_{iso}(\bar I_1)\f$/
       * \f$\partial\Psi_{iso}/\partial \bar I_1\f$ evaluation as the plain overload
       * (both funnel through one shared internal helper), so the two cannot drift
       * apart.
       *
       * @tparam T Scalar type, e.g. double, std::complex<double>.
       * @param C Right Cauchy-Green tensor.
       * @param mu Shear-modulus-like parameter.
       * @param lambdaL Locking stretch.
       * @return A tuple containing the (isochoric) energy density and its first
       * derivative w.r.t. C.
       */
      template < typename T >
      std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T > > ArrudaBoyce8ChainPotential(
        const Math::Tensors::FastorTensors::Tensor33t< T >& C,
        const double                                        mu,
        const double                                        lambdaL )
      {
        return Detail::isochoricPotential< 1, false >( C, [mu, lambdaL]( const std::array< T, 3 >& xBar ) {
          return Detail::arrudaBoyce8ChainIsochoric< 1 >( xBar, mu, lambdaL );
        } );
      }
    } // namespace FirstOrder

    namespace SecondOrder {

      /** @brief Isochoric neo-Hookean potential (see the plain #NeoHookePotential) and its first and second derivative
       * w.r.t. C.
       *
       * @tparam T Scalar type, e.g. double, std::complex<double>.
       * @param C Right Cauchy-Green tensor.
       * @param mu Shear modulus.
       * @return A tuple containing the (isochoric) energy density and its first and second derivative w.r.t. C.
       */
      template < typename T >
      std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T >, Math::Tensors::FastorTensors::Tensor3333t< T > > NeoHookePotential(
        const Math::Tensors::FastorTensors::Tensor33t< T >& C,
        const double                                        mu )
      {
        return Detail::isochoricPotential< 2, false >( C, [mu]( const std::array< T, 3 >& xBar ) {
          return Detail::neoHookeIsochoric< 2 >( xBar, mu );
        } );
      }

      /** @brief Isochoric Yeoh potential (see the plain #YeohPotential) and its first and second derivative w.r.t. C.
       *
       * @tparam T Scalar type, e.g. double, std::complex<double>.
       * @param C Right Cauchy-Green tensor.
       * @param C10 Yeoh material parameter.
       * @param C20 Yeoh material parameter.
       * @param C30 Yeoh material parameter.
       * @return A tuple containing the (isochoric) energy density and its first and second derivative w.r.t. C.
       */
      template < typename T >
      std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T >, Math::Tensors::FastorTensors::Tensor3333t< T > > YeohPotential(
        const Math::Tensors::FastorTensors::Tensor33t< T >& C,
        const double                                        C10,
        const double                                        C20,
        const double                                        C30 )
      {
        return Detail::isochoricPotential< 2, false >( C, [C10, C20, C30]( const std::array< T, 3 >& xBar ) {
          return Detail::yeohIsochoric< 2 >( xBar, C10, C20, C30 );
        } );
      }

      /** @brief Isochoric Mooney-Rivlin potential (see the plain #MooneyRivlinPotential) and its first and second
       * derivative w.r.t. C.
       *
       * @tparam T Scalar type, e.g. double, std::complex<double>.
       * @param C Right Cauchy-Green tensor.
       * @param C10 Mooney-Rivlin material parameter.
       * @param C01 Mooney-Rivlin material parameter.
       * @return A tuple containing the (isochoric) energy density and its first and second derivative w.r.t. C.
       */
      template < typename T >
      std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T >, Math::Tensors::FastorTensors::Tensor3333t< T > > MooneyRivlinPotential(
        const Math::Tensors::FastorTensors::Tensor33t< T >& C,
        const double                                        C10,
        const double                                        C01 )
      {
        return Detail::isochoricPotential< 2, true >( C, [C10, C01]( const std::array< T, 3 >& xBar ) {
          return Detail::mooneyRivlinIsochoric< 2 >( xBar, C10, C01 );
        } );
      }

      /** @brief Isochoric Arruda-Boyce 8-chain potential (see the plain #ArrudaBoyce8ChainPotential) and its first and
       * second derivative w.r.t. C.
       *
       * @tparam T Scalar type, e.g. double, std::complex<double>.
       * @param C Right Cauchy-Green tensor.
       * @param mu Shear-modulus-like parameter.
       * @param lambdaL Locking stretch (must be > 1).
       * @return A tuple containing the (isochoric) energy density and its first and second derivative w.r.t. C.
       */
      template < typename T >
      std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T >, Math::Tensors::FastorTensors::Tensor3333t< T > > ArrudaBoyce8ChainPotential(
        const Math::Tensors::FastorTensors::Tensor33t< T >& C,
        const double                                        mu,
        const double                                        lambdaL )
      {
        return Detail::isochoricPotential< 2, false >( C, [mu, lambdaL]( const std::array< T, 3 >& xBar ) {
          return Detail::arrudaBoyce8ChainIsochoric< 2 >( xBar, mu, lambdaL );
        } );
      }
    } // namespace SecondOrder
  }   // namespace Isochoric

  /** @brief Purely volumetric strain energy functions \f$U(J)\f$, \f$J = \sqrt{\det\boldsymbol C}\f$.
   *
   * Combine them with a potential from
   * Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::Isochoric via #DecoupledPotential or
   * #makeDecoupledPotential.
   */
  namespace Volumetric {

    /** @brief Volumetric penalty potential (logarithmic, quadratic in \f$\ln J\f$),
     *
     * \f[
     *   W_{\rm vol} = \frac{\kappa}{8}\left(\ln I_3\right)^2, \qquad I_3 = J^2 = \det\boldsymbol C
     * \f]
     *
     * This is \f$U_2 = \frac{\kappa}{2}(\ln J)^2\f$ in Table 1 of Doll & Schweizerhof (2000),
     * https://doi.org/10.1115/1.321146. Its convexity is lost for \f$J > e\f$, and the volumetric stress vanishes as
     * \f$J\to\infty\f$; for large volumetric expansions consider e.g. #OgdenPotential or
     * #DollSchweizerhofPotential instead.
     *
     * @tparam T Scalar type, e.g. double, float or an automatic differentiation number.
     * @param C Right Cauchy-Green tensor.
     * @param kappa Bulk modulus.
     * @return Energy density (volumetric part only).
     */
    template < typename T >
    T PenaltyPotential( const Math::Tensors::FastorTensors::Tensor33t< T >& C, const double kappa )
    {
      return Detail::volumetricPotential< 0 >( C, [kappa]( const T& J ) {
        return Detail::penaltyVolumetric< 0 >( J, kappa );
      } );
    }

    /** @brief Quadratic volumetric strain energy function, \f$U_1\f$ in Table 1 of Doll & Schweizerhof (2000),
     * https://doi.org/10.1115/1.321146 (Sussman & Bathe 1987, Simo 1988, among others),
     * \f[
     *   U_1(J) = \frac{\kappa}{2}(J-1)^2, \qquad J = \sqrt{\det\boldsymbol C}.
     * \f]
     * Note that \f$U_1\f$ stays bounded as \f$J\to 0\f$ (violating condition (9) of the paper), so it should not be
     * used for large compression.
     *
     * @tparam T Scalar type, e.g. double, std::complex<double> or an automatic differentiation number.
     * @param C Right Cauchy-Green tensor.
     * @param kappa Bulk modulus.
     * @return Energy density (volumetric part only).
     */
    template < typename T >
    T SussmanBathePotential( const Math::Tensors::FastorTensors::Tensor33t< T >& C, const double kappa )
    {
      return Detail::volumetricPotential< 0 >( C, [kappa]( const T& J ) {
        return Detail::sussmanBatheVolumetric< 0 >( J, kappa );
      } );
    }

    /** @brief Volumetric strain energy function \f$U_3\f$ in Table 1 of Doll & Schweizerhof (2000),
     * https://doi.org/10.1115/1.321146 (Simo & Taylor 1982), the mean of \f$U_1\f$
     * (#SussmanBathePotential) and \f$U_2\f$ (#PenaltyPotential),
     * \f[
     *   U_3(J) = \frac{\kappa}{4}\left[(J-1)^2 + (\ln J)^2\right], \qquad J = \sqrt{\det\boldsymbol C}.
     * \f]
     *
     * @tparam T Scalar type, e.g. double, std::complex<double> or an automatic differentiation number.
     * @param C Right Cauchy-Green tensor.
     * @param kappa Bulk modulus.
     * @return Energy density (volumetric part only).
     */
    template < typename T >
    T SimoTaylorPotential( const Math::Tensors::FastorTensors::Tensor33t< T >& C, const double kappa )
    {
      return Detail::volumetricPotential< 0 >( C, [kappa]( const T& J ) {
        return Detail::simoTaylorVolumetric< 0 >( J, kappa );
      } );
    }

    /** @brief Ogden (1972) volumetric strain energy function, \f$U_4\f$ in Table 1 of Doll & Schweizerhof (2000),
     * https://doi.org/10.1115/1.321146,
     * \f[
     *   U_4(J) = \frac{\kappa}{\theta^2}\left(\theta\ln J + J^{-\theta} - 1\right), \qquad J = \sqrt{\det\boldsymbol
     * C}. \f] All physical conditions of the paper hold only for \f$\theta < -1\f$; the popular choice \f$\theta=-2\f$
     * gives \f$\frac{\kappa}{4}(J^2 - 1 - 2\ln J)\f$ (Simo & Taylor 1991), and Ogden's rubber fit \f$\theta=9\f$ is
     * admissible for nearly incompressible applications. \f$\theta\neq0\f$ is required.
     *
     * @tparam T Scalar type, e.g. double, std::complex<double> or an automatic differentiation number.
     * @param C Right Cauchy-Green tensor.
     * @param kappa Bulk modulus.
     * @param theta Shape parameter \f$\theta\f$ (\f$\neq 0\f$).
     * @return Energy density (volumetric part only).
     */
    template < typename T >
    T OgdenPotential( const Math::Tensors::FastorTensors::Tensor33t< T >& C, const double kappa, const double theta )
    {
      return Detail::volumetricPotential< 0 >( C, [kappa, theta]( const T& J ) {
        return Detail::ogdenVolumetric< 0 >( J, kappa, theta );
      } );
    }

    /** @brief Volumetric strain energy function \f$U_5\f$ in Table 1 of Doll & Schweizerhof (2000),
     * https://doi.org/10.1115/1.321146 (Liu, Hofstetter & Mang 1992, 1994),
     * \f[
     *   U_5(J) = \kappa\left(J\ln J - J + 1\right), \qquad J = \sqrt{\det\boldsymbol C}.
     * \f]
     * Note that \f$U_5\f$ stays bounded as \f$J\to 0\f$ (violating condition (9) of the paper), so it should not be
     * used for large compression.
     *
     * @tparam T Scalar type, e.g. double, std::complex<double> or an automatic differentiation number.
     * @param C Right Cauchy-Green tensor.
     * @param kappa Bulk modulus.
     * @return Energy density (volumetric part only).
     */
    template < typename T >
    T LiuHofstetterMangPotential( const Math::Tensors::FastorTensors::Tensor33t< T >& C, const double kappa )
    {
      return Detail::volumetricPotential< 0 >( C, [kappa]( const T& J ) {
        return Detail::liuHofstetterMangVolumetric< 0 >( J, kappa );
      } );
    }

    /** @brief Two-parameter volumetric strain energy function \f$U_6\f$ acc. Doll & Schweizerhof (2000),
     * https://doi.org/10.1115/1.321146, Eqs. (12)-(13),
     * \f[
     *   U_6(J) = \frac{\kappa}{\alpha+\beta}\left[\frac{J^{\alpha+1}}{\alpha+1} + \frac{J^{-(\beta-1)}}{\beta-1}\right]
     *   - \frac{\kappa}{(\alpha+1)(\beta-1)}, \qquad
     *   \frac{\partial U_6}{\partial J} = \kappa\,\frac{J^\alpha - J^{-\beta}}{\alpha+\beta},
     * \f]
     * with \f$\alpha > 0\f$ and \f$\beta > 1\f$, which fulfills all physical conditions of the paper.
     * Special cases: \f$\beta=\alpha\f$ gives \f$p(J) = -p(J^{-1})\f$, \f$\beta=\alpha+2\f$ gives
     * \f$U_6(J)=U_6(J^{-1})\f$; \f$\alpha=1,\ \beta\to1\f$ approaches \f$U_4\f$ with \f$\theta=-2\f$, and
     * \f$\alpha\to0\f$ recovers Murnaghan's pressure formula.
     *
     * @tparam T Scalar type, e.g. double, std::complex<double> or an automatic differentiation number.
     * @param C Right Cauchy-Green tensor.
     * @param kappa Bulk modulus.
     * @param alpha Shape parameter \f$\alpha > 0\f$.
     * @param beta Shape parameter \f$\beta > 1\f$.
     * @return Energy density (volumetric part only).
     */
    template < typename T >
    T DollSchweizerhofPotential( const Math::Tensors::FastorTensors::Tensor33t< T >& C,
                                 const double                                        kappa,
                                 const double                                        alpha,
                                 const double                                        beta )
    {
      const double alphas[] = { alpha };
      const double betas[]  = { beta };
      return Detail::volumetricPotential< 0 >( C, [kappa, &alphas, &betas]( const T& J ) {
        return Detail::generalizedDollSchweizerhofVolumetric< 0 >( J, kappa, alphas, betas );
      } );
    }

    /** @brief Generalized volumetric strain energy function \f$U_{6\,\rm gen}\f$ acc. Doll & Schweizerhof (2000),
     * https://doi.org/10.1115/1.321146, Eq. (15): an additive composition of \f$n\f$ terms of type
     * #DollSchweizerhofPotential,
     * \f[
     *   \frac{\partial U_{6\,\rm gen}}{\partial J} = \kappa\,
     *   \frac{\sum_{i=1}^n \left(J^{\alpha_i} - J^{-\beta_i}\right)}{\sum_{i=1}^n (\alpha_i+\beta_i)}, \qquad
     *   U_{6\,\rm gen}(J) = \frac{\kappa}{\sum_i(\alpha_i+\beta_i)} \sum_{i=1}^n \left[
     *   \frac{J^{\alpha_i+1}-1}{\alpha_i+1} + \frac{J^{-(\beta_i-1)}-1}{\beta_i-1}\right],
     * \f]
     * with \f$\alpha_i > 0\f$ and \f$\beta_i > 1\f$.
     *
     * @tparam T Scalar type, e.g. double, std::complex<double> or an automatic differentiation number.
     * @param C Right Cauchy-Green tensor.
     * @param kappa Bulk modulus.
     * @param alphas Shape parameters \f$\alpha_i > 0\f$.
     * @param betas Shape parameters \f$\beta_i > 1\f$ (same, non-zero size as alphas).
     * @return Energy density (volumetric part only).
     */
    template < typename T >
    T GeneralizedDollSchweizerhofPotential( const Math::Tensors::FastorTensors::Tensor33t< T >& C,
                                            const double                                        kappa,
                                            const std::span< const double >                     alphas,
                                            const std::span< const double >                     betas )
    {
      return Detail::volumetricPotential< 0 >( C, [kappa, alphas, betas]( const T& J ) {
        return Detail::generalizedDollSchweizerhofVolumetric< 0 >( J, kappa, alphas, betas );
      } );
    }

    /** @brief Exponential volumetric strain energy function \f$U_7\f$ acc. Doll & Schweizerhof (2000),
     * https://doi.org/10.1115/1.321146, Eq. (16)\f$_1\f$,
     * \f[
     *   U_7(J) = \frac{\kappa}{2}\left(e^{J-1} - \ln J - 1\right), \qquad J = \sqrt{\det\boldsymbol C}.
     * \f]
     *
     * @tparam T Scalar type, e.g. double, std::complex<double> or an automatic differentiation number.
     * @param C Right Cauchy-Green tensor.
     * @param kappa Bulk modulus.
     * @return Energy density (volumetric part only).
     */
    template < typename T >
    T DollSchweizerhofExponentialPotential( const Math::Tensors::FastorTensors::Tensor33t< T >& C, const double kappa )
    {
      return Detail::volumetricPotential< 0 >( C, [kappa]( const T& J ) {
        return Detail::dollSchweizerhofExponentialVolumetric< 0 >( J, kappa );
      } );
    }

    /** @brief Logarithmic volumetric strain energy function \f$U_8\f$ acc. Doll & Schweizerhof (2000),
     * https://doi.org/10.1115/1.321146, Eq. (16)\f$_2\f$,
     * \f[
     *   U_8(J) = \frac{\kappa}{2}(J-1)\ln J, \qquad J = \sqrt{\det\boldsymbol C}.
     * \f]
     *
     * @tparam T Scalar type, e.g. double, std::complex<double> or an automatic differentiation number.
     * @param C Right Cauchy-Green tensor.
     * @param kappa Bulk modulus.
     * @return Energy density (volumetric part only).
     */
    template < typename T >
    T DollSchweizerhofLogarithmicPotential( const Math::Tensors::FastorTensors::Tensor33t< T >& C, const double kappa )
    {
      return Detail::volumetricPotential< 0 >( C, [kappa]( const T& J ) {
        return Detail::dollSchweizerhofLogarithmicVolumetric< 0 >( J, kappa );
      } );
    }
    namespace FirstOrder {

      /** @brief Shared volumetric penalty potential (see the plain
       * #PenaltyPotential) and its first derivative w.r.t. C,
       * \f[
       *   \frac{\partial W_{\rm vol}}{\partial \boldsymbol C} = \frac{\kappa}{4}\ln(I_3)\,\boldsymbol C^{-1}.
       * \f]
       *
       * @tparam T Scalar type, e.g. double, std::complex<double>.
       * @param C Right Cauchy-Green tensor.
       * @param kappa Bulk modulus.
       * @return A tuple containing the (volumetric) energy density and its first derivative w.r.t. C.
       */
      template < typename T >
      std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T > > PenaltyPotential(
        const Math::Tensors::FastorTensors::Tensor33t< T >& C,
        const double                                        kappa )
      {
        return Detail::volumetricPotential< 1 >( C, [kappa]( const T& J ) {
          return Detail::penaltyVolumetric< 1 >( J, kappa );
        } );
      }

      /** @brief Volumetric strain energy function \f$U_1\f$ (see the plain
       * #SussmanBathePotential) and its first derivative w.r.t. C,
       * \f[
       *   \frac{\partial U_1}{\partial \boldsymbol C} = \frac{\partial U_1}{\partial J}\,\frac{J}{2}\boldsymbol C^{-T},
       *   \qquad \frac{\partial U_1}{\partial J} = \kappa(J-1).
       * \f]
       *
       * @tparam T Scalar type, e.g. double, std::complex<double>.
       * @param C Right Cauchy-Green tensor.
       * @param kappa Bulk modulus.
       * @return A tuple containing the (volumetric) energy density and its first derivative w.r.t. C.
       */
      template < typename T >
      std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T > > SussmanBathePotential(
        const Math::Tensors::FastorTensors::Tensor33t< T >& C,
        const double                                        kappa )
      {
        return Detail::volumetricPotential< 1 >( C, [kappa]( const T& J ) {
          return Detail::sussmanBatheVolumetric< 1 >( J, kappa );
        } );
      }

      /** @brief Volumetric strain energy function \f$U_3\f$ (see the plain
       * #SimoTaylorPotential) and its first derivative w.r.t. C,
       * \f[
       *   \frac{\partial U_3}{\partial \boldsymbol C} = \frac{\partial U_3}{\partial J}\,\frac{J}{2}\boldsymbol C^{-T},
       *   \qquad \frac{\partial U_3}{\partial J} = \frac{\kappa}{2}\left(J-1+\frac{\ln J}{J}\right).
       * \f]
       *
       * @tparam T Scalar type, e.g. double, std::complex<double>.
       * @param C Right Cauchy-Green tensor.
       * @param kappa Bulk modulus.
       * @return A tuple containing the (volumetric) energy density and its first derivative w.r.t. C.
       */
      template < typename T >
      std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T > > SimoTaylorPotential(
        const Math::Tensors::FastorTensors::Tensor33t< T >& C,
        const double                                        kappa )
      {
        return Detail::volumetricPotential< 1 >( C, [kappa]( const T& J ) {
          return Detail::simoTaylorVolumetric< 1 >( J, kappa );
        } );
      }

      /** @brief Volumetric strain energy function \f$U_4\f$ (see the plain
       * #OgdenPotential) and its first derivative w.r.t. C,
       * \f[
       *   \frac{\partial U_4}{\partial \boldsymbol C} = \frac{\partial U_4}{\partial J}\,\frac{J}{2}\boldsymbol C^{-T},
       *   \qquad \frac{\partial U_4}{\partial J} = \frac{\kappa}{\theta}\,\frac{1-J^{-\theta}}{J}.
       * \f]
       *
       * @tparam T Scalar type, e.g. double, std::complex<double>.
       * @param C Right Cauchy-Green tensor.
       * @param kappa Bulk modulus.
       * @param theta Shape parameter \f$\theta\f$ (\f$\neq 0\f$).
       * @return A tuple containing the (volumetric) energy density and its first derivative w.r.t. C.
       */
      template < typename T >
      std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T > > OgdenPotential(
        const Math::Tensors::FastorTensors::Tensor33t< T >& C,
        const double                                        kappa,
        const double                                        theta )
      {
        return Detail::volumetricPotential< 1 >( C, [kappa, theta]( const T& J ) {
          return Detail::ogdenVolumetric< 1 >( J, kappa, theta );
        } );
      }

      /** @brief Volumetric strain energy function \f$U_5\f$ (see the plain
       * #LiuHofstetterMangPotential) and its first derivative w.r.t. C,
       * \f[
       *   \frac{\partial U_5}{\partial \boldsymbol C} = \frac{\partial U_5}{\partial J}\,\frac{J}{2}\boldsymbol C^{-T},
       *   \qquad \frac{\partial U_5}{\partial J} = \kappa\ln J.
       * \f]
       *
       * @tparam T Scalar type, e.g. double, std::complex<double>.
       * @param C Right Cauchy-Green tensor.
       * @param kappa Bulk modulus.
       * @return A tuple containing the (volumetric) energy density and its first derivative w.r.t. C.
       */
      template < typename T >
      std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T > > LiuHofstetterMangPotential(
        const Math::Tensors::FastorTensors::Tensor33t< T >& C,
        const double                                        kappa )
      {
        return Detail::volumetricPotential< 1 >( C, [kappa]( const T& J ) {
          return Detail::liuHofstetterMangVolumetric< 1 >( J, kappa );
        } );
      }

      /** @brief Volumetric strain energy function \f$U_6\f$ (see the plain
       * #DollSchweizerhofPotential) and its first derivative w.r.t. C,
       * \f[
       *   \frac{\partial U_6}{\partial \boldsymbol C} = \frac{\partial U_6}{\partial J}\,\frac{J}{2}\boldsymbol C^{-T},
       *   \qquad \frac{\partial U_6}{\partial J} = \kappa\,\frac{J^\alpha - J^{-\beta}}{\alpha+\beta}.
       * \f]
       *
       * @tparam T Scalar type, e.g. double, std::complex<double>.
       * @param C Right Cauchy-Green tensor.
       * @param kappa Bulk modulus.
       * @param alpha Shape parameter \f$\alpha > 0\f$.
       * @param beta Shape parameter \f$\beta > 1\f$.
       * @return A tuple containing the (volumetric) energy density and its first derivative w.r.t. C.
       */
      template < typename T >
      std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T > > DollSchweizerhofPotential(
        const Math::Tensors::FastorTensors::Tensor33t< T >& C,
        const double                                        kappa,
        const double                                        alpha,
        const double                                        beta )
      {
        const double alphas[] = { alpha };
        const double betas[]  = { beta };
        return Detail::volumetricPotential< 1 >( C, [kappa, &alphas, &betas]( const T& J ) {
          return Detail::generalizedDollSchweizerhofVolumetric< 1 >( J, kappa, alphas, betas );
        } );
      }

      /** @brief Volumetric strain energy function \f$U_{6\,\rm gen}\f$ (see the plain
       * #GeneralizedDollSchweizerhofPotential) and its first derivative w.r.t. C,
       * \f[
       *   \frac{\partial U_{6\,\rm gen}}{\partial \boldsymbol C} = \frac{\partial U_{6\,\rm gen}}{\partial
       * J}\,\frac{J}{2}\boldsymbol C^{-T}, \qquad \frac{\partial U_{6\,\rm gen}}{\partial J} = \kappa\,\frac{\sum_i
       * (J^{\alpha_i} - J^{-\beta_i})}{\sum_i(\alpha_i+\beta_i)}. \f]
       *
       * @tparam T Scalar type, e.g. double, std::complex<double>.
       * @param C Right Cauchy-Green tensor.
       * @param kappa Bulk modulus.
       * @param alphas Shape parameters \f$\alpha_i > 0\f$.
       * @param betas Shape parameters \f$\beta_i > 1\f$ (same, non-zero size as alphas).
       * @return A tuple containing the (volumetric) energy density and its first derivative w.r.t. C.
       */
      template < typename T >
      std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T > > GeneralizedDollSchweizerhofPotential(
        const Math::Tensors::FastorTensors::Tensor33t< T >& C,
        const double                                        kappa,
        const std::span< const double >                     alphas,
        const std::span< const double >                     betas )
      {
        return Detail::volumetricPotential< 1 >( C, [kappa, alphas, betas]( const T& J ) {
          return Detail::generalizedDollSchweizerhofVolumetric< 1 >( J, kappa, alphas, betas );
        } );
      }

      /** @brief Volumetric strain energy function \f$U_7\f$ (see the plain
       * #DollSchweizerhofExponentialPotential) and its first derivative w.r.t. C,
       * \f[
       *   \frac{\partial U_7}{\partial \boldsymbol C} = \frac{\partial U_7}{\partial J}\,\frac{J}{2}\boldsymbol C^{-T},
       *   \qquad \frac{\partial U_7}{\partial J} = \frac{\kappa}{2}\left(e^{J-1} - \frac{1}{J}\right).
       * \f]
       *
       * @tparam T Scalar type, e.g. double, std::complex<double>.
       * @param C Right Cauchy-Green tensor.
       * @param kappa Bulk modulus.
       * @return A tuple containing the (volumetric) energy density and its first derivative w.r.t. C.
       */
      template < typename T >
      std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T > > DollSchweizerhofExponentialPotential(
        const Math::Tensors::FastorTensors::Tensor33t< T >& C,
        const double                                        kappa )
      {
        return Detail::volumetricPotential< 1 >( C, [kappa]( const T& J ) {
          return Detail::dollSchweizerhofExponentialVolumetric< 1 >( J, kappa );
        } );
      }

      /** @brief Volumetric strain energy function \f$U_8\f$ (see the plain
       * #DollSchweizerhofLogarithmicPotential) and its first derivative w.r.t. C,
       * \f[
       *   \frac{\partial U_8}{\partial \boldsymbol C} = \frac{\partial U_8}{\partial J}\,\frac{J}{2}\boldsymbol C^{-T},
       *   \qquad \frac{\partial U_8}{\partial J} = \frac{\kappa}{2}\left(\ln J + \frac{J-1}{J}\right).
       * \f]
       *
       * @tparam T Scalar type, e.g. double, std::complex<double>.
       * @param C Right Cauchy-Green tensor.
       * @param kappa Bulk modulus.
       * @return A tuple containing the (volumetric) energy density and its first derivative w.r.t. C.
       */
      template < typename T >
      std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T > > DollSchweizerhofLogarithmicPotential(
        const Math::Tensors::FastorTensors::Tensor33t< T >& C,
        const double                                        kappa )
      {
        return Detail::volumetricPotential< 1 >( C, [kappa]( const T& J ) {
          return Detail::dollSchweizerhofLogarithmicVolumetric< 1 >( J, kappa );
        } );
      }
    } // namespace FirstOrder

    namespace SecondOrder {

      /** @brief Volumetric penalty potential (see the plain #PenaltyPotential) and its first and second derivative
       * w.r.t. C.
       *
       * @tparam T Scalar type, e.g. double, std::complex<double>.
       * @param C Right Cauchy-Green tensor.
       * @param kappa Bulk modulus.
       * @return A tuple containing the (volumetric) energy density and its first and second derivative w.r.t. C.
       */
      template < typename T >
      std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T >, Math::Tensors::FastorTensors::Tensor3333t< T > > PenaltyPotential(
        const Math::Tensors::FastorTensors::Tensor33t< T >& C,
        const double                                        kappa )
      {
        return Detail::volumetricPotential< 2 >( C, [kappa]( const T& J ) {
          return Detail::penaltyVolumetric< 2 >( J, kappa );
        } );
      }

      /** @brief Volumetric strain energy function \f$U_1\f$ (see the plain #SussmanBathePotential) and its first and
       * second derivative w.r.t. C.
       *
       * @tparam T Scalar type, e.g. double, std::complex<double>.
       * @param C Right Cauchy-Green tensor.
       * @param kappa Bulk modulus.
       * @return A tuple containing the (volumetric) energy density and its first and second derivative w.r.t. C.
       */
      template < typename T >
      std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T >, Math::Tensors::FastorTensors::Tensor3333t< T > > SussmanBathePotential(
        const Math::Tensors::FastorTensors::Tensor33t< T >& C,
        const double                                        kappa )
      {
        return Detail::volumetricPotential< 2 >( C, [kappa]( const T& J ) {
          return Detail::sussmanBatheVolumetric< 2 >( J, kappa );
        } );
      }

      /** @brief Volumetric strain energy function \f$U_3\f$ (see the plain #SimoTaylorPotential) and its first and
       * second derivative w.r.t. C.
       *
       * @tparam T Scalar type, e.g. double, std::complex<double>.
       * @param C Right Cauchy-Green tensor.
       * @param kappa Bulk modulus.
       * @return A tuple containing the (volumetric) energy density and its first and second derivative w.r.t. C.
       */
      template < typename T >
      std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T >, Math::Tensors::FastorTensors::Tensor3333t< T > > SimoTaylorPotential(
        const Math::Tensors::FastorTensors::Tensor33t< T >& C,
        const double                                        kappa )
      {
        return Detail::volumetricPotential< 2 >( C, [kappa]( const T& J ) {
          return Detail::simoTaylorVolumetric< 2 >( J, kappa );
        } );
      }

      /** @brief Volumetric strain energy function \f$U_4\f$ (see the plain #OgdenPotential) and its first and second
       * derivative w.r.t. C.
       *
       * @tparam T Scalar type, e.g. double, std::complex<double>.
       * @param C Right Cauchy-Green tensor.
       * @param kappa Bulk modulus.
       * @param theta Shape parameter \f$\theta\f$ (\f$\neq 0\f$).
       * @return A tuple containing the (volumetric) energy density and its first and second derivative w.r.t. C.
       */
      template < typename T >
      std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T >, Math::Tensors::FastorTensors::Tensor3333t< T > > OgdenPotential(
        const Math::Tensors::FastorTensors::Tensor33t< T >& C,
        const double                                        kappa,
        const double                                        theta )
      {
        return Detail::volumetricPotential< 2 >( C, [kappa, theta]( const T& J ) {
          return Detail::ogdenVolumetric< 2 >( J, kappa, theta );
        } );
      }

      /** @brief Volumetric strain energy function \f$U_5\f$ (see the plain #LiuHofstetterMangPotential) and its first
       * and second derivative w.r.t. C.
       *
       * @tparam T Scalar type, e.g. double, std::complex<double>.
       * @param C Right Cauchy-Green tensor.
       * @param kappa Bulk modulus.
       * @return A tuple containing the (volumetric) energy density and its first and second derivative w.r.t. C.
       */
      template < typename T >
      std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T >, Math::Tensors::FastorTensors::Tensor3333t< T > > LiuHofstetterMangPotential(
        const Math::Tensors::FastorTensors::Tensor33t< T >& C,
        const double                                        kappa )
      {
        return Detail::volumetricPotential< 2 >( C, [kappa]( const T& J ) {
          return Detail::liuHofstetterMangVolumetric< 2 >( J, kappa );
        } );
      }

      /** @brief Volumetric strain energy function \f$U_6\f$ (see the plain #DollSchweizerhofPotential) and its first
       * and second derivative w.r.t. C.
       *
       * @tparam T Scalar type, e.g. double, std::complex<double>.
       * @param C Right Cauchy-Green tensor.
       * @param kappa Bulk modulus.
       * @param alpha Shape parameter \f$\alpha > 0\f$.
       * @param beta Shape parameter \f$\beta > 1\f$.
       * @return A tuple containing the (volumetric) energy density and its first and second derivative w.r.t. C.
       */
      template < typename T >
      std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T >, Math::Tensors::FastorTensors::Tensor3333t< T > > DollSchweizerhofPotential(
        const Math::Tensors::FastorTensors::Tensor33t< T >& C,
        const double                                        kappa,
        const double                                        alpha,
        const double                                        beta )
      {
        const double alphas[] = { alpha };
        const double betas[]  = { beta };
        return Detail::volumetricPotential< 2 >( C, [kappa, &alphas, &betas]( const T& J ) {
          return Detail::generalizedDollSchweizerhofVolumetric< 2 >( J, kappa, alphas, betas );
        } );
      }

      /** @brief Volumetric strain energy function \f$U_{6\,\rm gen}\f$ (see the plain
       * #GeneralizedDollSchweizerhofPotential) and its first and second derivative w.r.t. C.
       *
       * @tparam T Scalar type, e.g. double, std::complex<double>.
       * @param C Right Cauchy-Green tensor.
       * @param kappa Bulk modulus.
       * @param alphas Shape parameters \f$\alpha_i > 0\f$.
       * @param betas Shape parameters \f$\beta_i > 1\f$ (same, non-zero size as alphas).
       * @return A tuple containing the (volumetric) energy density and its first and second derivative w.r.t. C.
       */
      template < typename T >
      std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T >, Math::Tensors::FastorTensors::Tensor3333t< T > > GeneralizedDollSchweizerhofPotential(
        const Math::Tensors::FastorTensors::Tensor33t< T >& C,
        const double                                        kappa,
        const std::span< const double >                     alphas,
        const std::span< const double >                     betas )
      {
        return Detail::volumetricPotential< 2 >( C, [kappa, alphas, betas]( const T& J ) {
          return Detail::generalizedDollSchweizerhofVolumetric< 2 >( J, kappa, alphas, betas );
        } );
      }

      /** @brief Volumetric strain energy function \f$U_7\f$ (see the plain #DollSchweizerhofExponentialPotential) and
       * its first and second derivative w.r.t. C.
       *
       * @tparam T Scalar type, e.g. double, std::complex<double>.
       * @param C Right Cauchy-Green tensor.
       * @param kappa Bulk modulus.
       * @return A tuple containing the (volumetric) energy density and its first and second derivative w.r.t. C.
       */
      template < typename T >
      std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T >, Math::Tensors::FastorTensors::Tensor3333t< T > > DollSchweizerhofExponentialPotential(
        const Math::Tensors::FastorTensors::Tensor33t< T >& C,
        const double                                        kappa )
      {
        return Detail::volumetricPotential< 2 >( C, [kappa]( const T& J ) {
          return Detail::dollSchweizerhofExponentialVolumetric< 2 >( J, kappa );
        } );
      }

      /** @brief Volumetric strain energy function \f$U_8\f$ (see the plain #DollSchweizerhofLogarithmicPotential) and
       * its first and second derivative w.r.t. C.
       *
       * @tparam T Scalar type, e.g. double, std::complex<double>.
       * @param C Right Cauchy-Green tensor.
       * @param kappa Bulk modulus.
       * @return A tuple containing the (volumetric) energy density and its first and second derivative w.r.t. C.
       */
      template < typename T >
      std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T >, Math::Tensors::FastorTensors::Tensor3333t< T > > DollSchweizerhofLogarithmicPotential(
        const Math::Tensors::FastorTensors::Tensor33t< T >& C,
        const double                                        kappa )
      {
        return Detail::volumetricPotential< 2 >( C, [kappa]( const T& J ) {
          return Detail::dollSchweizerhofLogarithmicVolumetric< 2 >( J, kappa );
        } );
      }
    } // namespace SecondOrder
  }   // namespace Volumetric

  /** @cond */
  namespace Detail {
    template < typename R >
    struct IsTuple : std::false_type {};

    template < typename... Ts >
    struct IsTuple< std::tuple< Ts... > > : std::true_type {};

    /* Sum of two potential results of the same type: a scalar energy (plain potentials), or a tuple of energy and
     * derivatives (FirstOrder/SecondOrder potentials), summed entry by entry. */
    template < typename R >
    R sumPotentialResults( const R& a, const R& b )
    {
      if constexpr ( IsTuple< R >::value ) {
        return [&]< std::size_t... i >( std::index_sequence< i... > ) {
          return R( std::tuple_element_t< i, R >( std::get< i >( a ) + std::get< i >( b ) )... );
        }( std::make_index_sequence< std::tuple_size_v< R > >{} );
      }
      else
        return R( a + b );
    }
  } // namespace Detail
  /** @endcond */

  /** @brief Decoupled hyperelastic potential
   * \f$\Psi(\boldsymbol C) = \Psi_{\rm iso}(\boldsymbol C) + U(\boldsymbol C)\f$, assembled from an isochoric part
   * (e.g. from Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::Isochoric) and a volumetric part (e.g.
   * from Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::Volumetric).
   *
   * Both parts are callables of \f$\boldsymbol C\f$ that return the same type: either the energy alone (plain
   * potentials) or a tuple of the energy and its derivatives w.r.t. \f$\boldsymbol C\f$ (e.g. the @c FirstOrder
   * potentials), which are then summed entry by entry. Example:
   * @code
   * using namespace Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions;
   * const auto [psi, dPsi_dC] = DecoupledPotential(
   *   C,
   *   [&]( const auto& C_ ) { return Isochoric::FirstOrder::NeoHookePotential( C_, mu ); },
   *   [&]( const auto& C_ ) { return Volumetric::FirstOrder::PenaltyPotential( C_, kappa ); } );
   * @endcode
   *
   * @par Single pass
   * The parts can instead be given as scalar functions of the invariants (IsochoricScalarFunction and
   * VolumetricScalarFunction, e.g. Isochoric::Yeoh and Volumetric::SimoTaylor). The scalar derivatives of both
   * parts are then added before the chain rule, so the invariants of \f$\boldsymbol C\f$ are computed once and the
   * derivatives assembled once; the result is the energy, or with FirstOrder::DecoupledPotential and
   * SecondOrder::DecoupledPotential the tuple with its derivatives:
   * @code
   * const double psi = DecoupledPotential( C, Isochoric::Yeoh{ C10, C20, C30 }, Volumetric::SimoTaylor{ kappa } );
   * const auto [psi_, dPsi_dC, d2Psi_dCdC] = SecondOrder::DecoupledPotential( C,
   *                                                                           Isochoric::Yeoh{ C10, C20, C30 },
   *                                                                           Volumetric::SimoTaylor{ kappa } );
   * @endcode
   *
   * @tparam T Scalar type, e.g. double, std::complex<double> or an automatic differentiation number.
   * @tparam IsochoricPart Callable \f$\boldsymbol C \mapsto \Psi_{\rm iso}\f$ (and derivatives).
   * @tparam VolumetricPart Callable \f$\boldsymbol C \mapsto U\f$ (and derivatives).
   * @param C Right Cauchy-Green tensor.
   * @param isochoricPart Isochoric part of the potential.
   * @param volumetricPart Volumetric part of the potential.
   * @return Sum of both parts, of the type returned by the parts.
   */
  template < typename T, typename IsochoricPart, typename VolumetricPart >
  auto DecoupledPotential( const Math::Tensors::FastorTensors::Tensor33t< T >& C,
                           const IsochoricPart&                                isochoricPart,
                           const VolumetricPart&                               volumetricPart )
  {
    return Detail::sumPotentialResults( isochoricPart( C ), volumetricPart( C ) );
  }

  /** @cond */
  /* Single-pass overload of DecoupledPotential for two scalar functions, documented with the generic overload
   * (Doxygen merges the two, as they differ in their constraints only). */
  template < typename T, IsochoricScalarFunction IsochoricPart, VolumetricScalarFunction VolumetricPart >
  auto DecoupledPotential( const Math::Tensors::FastorTensors::Tensor33t< T >& C,
                           const IsochoricPart&                                isochoricPart,
                           const VolumetricPart&                               volumetricPart )
  {
    return Detail::decoupledPotential< 0 >( C, isochoricPart, volumetricPart );
  }
  /** @endcond */

  /** @brief Bundles an isochoric and a volumetric part into one reusable decoupled potential, i.e. a callable
   * \f$\boldsymbol C \mapsto \Psi_{\rm iso}(\boldsymbol C) + U(\boldsymbol C)\f$ that evaluates #DecoupledPotential.
   *
   * Since the returned callable is generic in the scalar type, it can be handed directly to automatic or
   * complex-step differentiation, or stored in a @c std::function:
   * @code
   * using namespace Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions;
   * const auto psi = makeDecoupledPotential( [=]( const auto& C ) { return Isochoric::NeoHookePotential( C, mu ); },
   *                                          [=]( const auto& C ) { return Volumetric::PenaltyPotential( C, kappa ); }
   * );
   * @endcode
   *
   * @tparam IsochoricPart Callable \f$\boldsymbol C \mapsto \Psi_{\rm iso}\f$ (and derivatives).
   * @tparam VolumetricPart Callable \f$\boldsymbol C \mapsto U\f$ (and derivatives).
   * @param isochoricPart Isochoric part of the potential.
   * @param volumetricPart Volumetric part of the potential.
   * @return Callable evaluating the sum of both parts.
   */
  template < typename IsochoricPart, typename VolumetricPart >
  auto makeDecoupledPotential( IsochoricPart isochoricPart, VolumetricPart volumetricPart )
  {
    return [isochoricPart, volumetricPart]< typename T >( const Math::Tensors::FastorTensors::Tensor33t< T >& C ) {
      return DecoupledPotential( C, isochoricPart, volumetricPart );
    };
  }

  /**
   * @brief Mooney-Rivlin Hyperelastic Energy Density Function
   * The energy density function \f$W_{MR}\f$ is given as
   * \f[
   *   W_{MR} = C_1 (\bar{I}_1 - 3) + C_2 (\bar{I}_2 - 3) + \frac{1}{D_1}\left( \frac{J^2 - 1}{2} - \ln J \right)
   * \f]
   * where \f$ \bar{I}_1 = I_1 J^{-\frac{2}{3}} \f$ and \f$ \bar{I}_2 = I_2 J^{-\frac{4}{3}} \f$ are the first
   * and second invariant of the isochoric right Cauchy-Green tensor, respectively, \f$ I_1 =
   * \text{tr}(\boldsymbol{C}) \f$ and \f$ I_2 = 0.5 (I_1^2 - \text{tr}(\boldsymbol{C}^2)) \f$ are the first and
   * second invariant of the right Cauchy-Green tensor \f$ \boldsymbol{C} = \boldsymbol{F}^T \boldsymbol{F} \f$, \f$ J
   * = \sqrt{\det(\boldsymbol{C})} = \det(\boldsymbol{F}) \f$ is the determinant of the deformation gradient, and \f$
   * C_1, C_2, D_1 \f$ are material parameters.
   *
   * It is the #DecoupledPotential of Isochoric::MooneyRivlinPotential and Volumetric::OgdenPotential with
   * \f$\kappa = 2/D_1\f$ and \f$\theta = -2\f$.
   *
   * @tparam T Scalar type, e.g. double, float, etc.
   * @param C Right Cauchy-Green tensor
   * @param C1 Mooney-Rivlin material parameter C1
   * @param C2 Mooney-Rivlin material parameter C2
   * @param D1 Mooney-Rivlin material parameter D1
   * @return Energy density
   */
  template < typename T >
  T MooneyRivlinPotential( const Math::Tensors::FastorTensors::Tensor33t< T >& C,
                           const double                                        C1,
                           const double                                        C2,
                           const double                                        D1 )
  {
    return Detail::decoupledPotential< 0 >( C, Isochoric::MooneyRivlin{ C1, C2 }, Volumetric::Ogden{ 2. / D1, -2. } );
  }

  /** @brief Yeoh hyperelastic energy density function.
   *
   *  The energy density function is given as
   *  \f[
   *    W = C_1 (\bar{I}_1 - 3) + C_2 (\bar{I}_1 - 3)^2 + C_3 (\bar{I}_1 - 3)^3 + \frac{1}{D_1}\left(
   *    \frac{J^2 - 1}{2} - \ln J \right)
   *  \f]
   *  where \f$\bar{I}_1 = I_1 J^{-2/3}\f$, \f$I_1 = \mathrm{tr}(\boldsymbol{C})\f$, and
   *  \f$J = \sqrt{\det(\boldsymbol{C})}\f$.
   *
   * It is the #DecoupledPotential of Isochoric::YeohPotential and Volumetric::OgdenPotential with
   * \f$\kappa = 2/D_1\f$ and \f$\theta = -2\f$.
   *
   * @tparam T Scalar type, e.g. double, float or an automatic differentiation number.
   * @param C Right Cauchy-Green tensor.
   * @param C1 Yeoh material parameter.
   * @param C2 Yeoh material parameter.
   * @param C3 Yeoh material parameter.
   * @param D1 Volumetric penalty parameter.
   * @return Energy density.
   */
  template < typename T >
  T YeohPotential( const Math::Tensors::FastorTensors::Tensor33t< T >& C,
                   const double                                        C1,
                   const double                                        C2,
                   const double                                        C3,
                   const double                                        D1 )
  {
    return Detail::decoupledPotential< 0 >( C, Isochoric::Yeoh{ C1, C2, C3 }, Volumetric::Ogden{ 2. / D1, -2. } );
  }

  /** @brief Compressible Arruda-Boyce 8-chain hyperelastic energy density function, with the same volumetric term
   * and \f$D_1\f$ convention as #MooneyRivlinPotential and #YeohPotential,
   * \f[
   *   W = W_{AB}(\bar I_1;\mu,\lambda_L) + \frac{1}{D_1}\left( \frac{J^2 - 1}{2} - \ln J \right),
   * \f]
   * where \f$W_{AB}\f$ is the isochoric Isochoric::ArrudaBoyce8ChainPotential.
   *
   * It is the #DecoupledPotential of Isochoric::ArrudaBoyce8ChainPotential and Volumetric::OgdenPotential with
   * \f$\kappa = 2/D_1\f$ and \f$\theta = -2\f$.
   *
   * @tparam T Scalar type, e.g. double, float or an automatic differentiation number.
   * @param C Right Cauchy-Green tensor.
   * @param mu Shear-modulus-like parameter.
   * @param lambdaL Locking stretch (must be > 1).
   * @param D1 Volumetric compliance parameter.
   * @return Energy density.
   */
  template < typename T >
  T ArrudaBoyce8ChainPotential( const Math::Tensors::FastorTensors::Tensor33t< T >& C,
                                const double                                        mu,
                                const double                                        lambdaL,
                                const double                                        D1 )
  {
    return Detail::decoupledPotential< 0 >( C,
                                            Isochoric::ArrudaBoyce8Chain{ mu, lambdaL },
                                            Volumetric::Ogden{ 2. / D1, -2. } );
  }

  /** @brief Hyperelastic potential given as a scalar function of the principal invariants
   * \f$(I_1, I_2, J)\f$ of \f$\boldsymbol C\f$ (PrincipalInvariantScalarFunction, e.g. Coupled::StandardNeoHooke),
   * for potentials that are not the sum of an isochoric and a volumetric part (see #DecoupledPotential for those).
   * The invariants are computed once; FirstOrder::CoupledPotential and SecondOrder::CoupledPotential also return
   * the derivatives w.r.t. the nine independent components of \f$\boldsymbol C\f$:
   * @code
   * const double psi = CoupledPotential( C, Coupled::StandardNeoHooke{ K, G } );
   * const auto [psi_, dPsi_dC, d2Psi_dCdC] = SecondOrder::CoupledPotential( C, Coupled::StandardNeoHooke{ K, G } );
   * @endcode
   *
   * @tparam T Scalar type, e.g. double, std::complex<double> or an automatic differentiation number.
   * @param C Right Cauchy-Green tensor.
   * @param potential Scalar function of the principal invariants, e.g. Coupled::StandardNeoHooke.
   * @return Energy density.
   */
  template < typename T, PrincipalInvariantScalarFunction Potential >
  T CoupledPotential( const Math::Tensors::FastorTensors::Tensor33t< T >& C, const Potential& potential )
  {
    return Detail::coupledPotential< 0 >( C, potential );
  }

  /** @cond */
  namespace Detail {
    /* PenceGouPotentialB as a function of (I1, J), with its scalar derivatives w.r.t. (I1, I2, J). */
    template < int Order, typename T >
    ScalarDerivatives< T, Order > penceGouB( const std::array< T, 3 >& x, const double K, const double G )
    {
      constexpr std::size_t iI1 = Invariants::RightCauchyGreen::indexI1;
      constexpr std::size_t iJ  = Invariants::RightCauchyGreen::indexJ;

      const T& I1     = x[iI1];
      const T& J      = x[iJ];
      const T  Jm23   = pow( J, -2. / 3. );
      const T  JmJInv = J - 1. / J;

      ScalarDerivatives< T, Order > f;
      f.f = K / 8. * JmJInv * JmJInv + G / 2. * ( I1 * Jm23 - 3. );
      if constexpr ( Order >= 1 ) {
        f.df[iI1] = G / 2. * Jm23;
        f.df[iJ]  = K / 4. * JmJInv * ( 1. + 1. / ( J * J ) ) - G / 3. * I1 * Jm23 / J;
      }
      if constexpr ( Order >= 2 ) {
        const T J2     = J * J;
        f.d2f[iJ][iJ]  = K / 4. * ( 1. + 3. / ( J2 * J2 ) ) + 5. / 9. * G * I1 * Jm23 / J2;
        f.d2f[iI1][iJ] = -G / 3. * Jm23 / J;
        f.d2f[iJ][iI1] = f.d2f[iI1][iJ];
      }
      return f;
    }
  } // namespace Detail
  /** @endcond */

  namespace FirstOrder {
    /** @brief Potential of the principal invariants (see the plain
     * Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::CoupledPotential) and its first derivative
     * w.r.t. C.
     *
     * @tparam T Scalar type, e.g. double, std::complex<double> or an automatic differentiation number.
     * @param C Right Cauchy-Green tensor.
     * @param potential Scalar function of the principal invariants, e.g. Coupled::StandardNeoHooke.
     * @return Tuple of energy density and its first derivative.
     */
    template < typename T, PrincipalInvariantScalarFunction Potential >
    std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T > > CoupledPotential(
      const Math::Tensors::FastorTensors::Tensor33t< T >& C,
      const Potential&                                    potential )
    {
      return Detail::coupledPotential< 1 >( C, potential );
    }

    /** @brief Single-pass decoupled potential of two scalar functions (see the plain
     * Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::DecoupledPotential) and its first derivative
     * w.r.t. C.
     *
     * @tparam T Scalar type, e.g. double, std::complex<double> or an automatic differentiation number.
     * @param C Right Cauchy-Green tensor.
     * @param isochoricPart Isochoric part, e.g. Isochoric::Yeoh.
     * @param volumetricPart Volumetric part, e.g. Volumetric::SimoTaylor.
     * @return Tuple of energy density and its first derivative.
     */
    template < typename T, IsochoricScalarFunction IsochoricPart, VolumetricScalarFunction VolumetricPart >
    std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T > > DecoupledPotential(
      const Math::Tensors::FastorTensors::Tensor33t< T >& C,
      const IsochoricPart&                                isochoricPart,
      const VolumetricPart&                               volumetricPart )
    {
      return Detail::decoupledPotential< 1 >( C, isochoricPart, volumetricPart );
    }

    /** @brief Hyperelastic Energy Density Function Wb acc. Pence & Gou (2015), Eq. (2.12) and its first derivative
     * w.r.t. C
     *
     *  The energy density function \f$W_b\f$ is given as
     *  \f[
     *    W_b = \frac{K}{8} \left(J - \frac{1}{J}\right)^2 + \frac{G}{2} \left(I_1 J^{-\frac{2}{3}} - 3\right)
     *  \f]
     *  where \f$ I_1 = \text{tr}(\boldsymbol{C}) \f$ is the first invariant of the right Cauchy-Green tensor
     *  \f$ \boldsymbol{C} = \boldsymbol{F}^T \boldsymbol{F} \f$, \f$ J = \sqrt{\det(\boldsymbol{C})} =
     * \det(\boldsymbol{F}) \f$ is the determinant of the deformation gradient, and \f$ K, G \f$ are the bulk and
     * shear modulus, respectively.
     *
     *  Additionally, the first derivative w.r.t. C is computed as
     *  \f[
     *    \frac{\partial W_b}{\partial \boldsymbol{C}} = \frac{\partial W_b}{\partial J} \frac{\partial J}{\partial
     * \boldsymbol{C}} + \frac{\partial W_b}{\partial I_1} \frac{\partial I_1}{\partial \boldsymbol{C}} \f] where \f[
     *    \frac{\partial J}{\partial \boldsymbol{C}} = \frac{1}{2} J \boldsymbol{C}^{-1}
     *  \f]
     *  and
     *  \f[
     *    \frac{\partial I_1}{\partial \boldsymbol{C}} = \boldsymbol{I}
     *  \f]
     *
     * @tparam T Scalar type, e.g. double, float, etc.
     * @param C Right Cauchy-Green tensor
     * @param K Bulk modulus
     * @param G Shear modulus
     * @return A tuple containing energy density and its first derivative w.r.t. C
     */
    template < typename T >
    std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T > > PenceGouPotentialB(
      const Math::Tensors::FastorTensors::Tensor33t< T >& C,
      const double                                        K,
      const double                                        G )
    {
      const auto inv = Invariants::RightCauchyGreen::computeInvariants< 1 >( C, false );
      return Invariants::RightCauchyGreen::composeWithInvariants( Detail::penceGouB< 1 >( inv.x, K, G ), inv );
    }

    /** @brief Mooney-Rivlin potential (see the plain
     * Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::MooneyRivlinPotential) and its first derivative
     * w.r.t. C.
     *
     * @tparam T Scalar type, e.g. double, std::complex<double>.
     * @param C Right Cauchy-Green tensor.
     * @param C1 Mooney-Rivlin material parameter C1.
     * @param C2 Mooney-Rivlin material parameter C2.
     * @param D1 Mooney-Rivlin material parameter D1.
     * @return A tuple containing the energy density and its first derivative w.r.t. C.
     */
    template < typename T >
    std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T > > MooneyRivlinPotential(
      const Math::Tensors::FastorTensors::Tensor33t< T >& C,
      const double                                        C1,
      const double                                        C2,
      const double                                        D1 )
    {
      return Detail::decoupledPotential< 1 >( C, Isochoric::MooneyRivlin{ C1, C2 }, Volumetric::Ogden{ 2. / D1, -2. } );
    }

    /** @brief Yeoh potential (see the plain
     * Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::YeohPotential) and its first derivative w.r.t.
     * C.
     *
     * @tparam T Scalar type, e.g. double, std::complex<double>.
     * @param C Right Cauchy-Green tensor.
     * @param C1 Yeoh material parameter.
     * @param C2 Yeoh material parameter.
     * @param C3 Yeoh material parameter.
     * @param D1 Volumetric penalty parameter.
     * @return A tuple containing the energy density and its first derivative w.r.t. C.
     */
    template < typename T >
    std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T > > YeohPotential(
      const Math::Tensors::FastorTensors::Tensor33t< T >& C,
      const double                                        C1,
      const double                                        C2,
      const double                                        C3,
      const double                                        D1 )
    {
      return Detail::decoupledPotential< 1 >( C, Isochoric::Yeoh{ C1, C2, C3 }, Volumetric::Ogden{ 2. / D1, -2. } );
    }

    /** @brief Compressible Arruda-Boyce 8-chain potential (see the plain
     * Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::ArrudaBoyce8ChainPotential) and its first
     * derivative w.r.t. C.
     *
     * @tparam T Scalar type, e.g. double, std::complex<double>.
     * @param C Right Cauchy-Green tensor.
     * @param mu Shear-modulus-like parameter.
     * @param lambdaL Locking stretch (must be > 1).
     * @param D1 Volumetric compliance parameter.
     * @return A tuple containing the energy density and its first derivative w.r.t. C.
     */
    template < typename T >
    std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T > > ArrudaBoyce8ChainPotential(
      const Math::Tensors::FastorTensors::Tensor33t< T >& C,
      const double                                        mu,
      const double                                        lambdaL,
      const double                                        D1 )
    {
      return Detail::decoupledPotential< 1 >( C,
                                              Isochoric::ArrudaBoyce8Chain{ mu, lambdaL },
                                              Volumetric::Ogden{ 2. / D1, -2. } );
    }
  } // namespace FirstOrder

  namespace SecondOrder {
    /** @brief Potential of the principal invariants (see the plain
     * Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::CoupledPotential) and its first and second
     * derivative w.r.t. C.
     *
     * @tparam T Scalar type, e.g. double, std::complex<double> or an automatic differentiation number.
     * @param C Right Cauchy-Green tensor.
     * @param potential Scalar function of the principal invariants, e.g. Coupled::StandardNeoHooke.
     * @return Tuple of energy density and its first and second derivative.
     */
    template < typename T, PrincipalInvariantScalarFunction Potential >
    std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T >, Math::Tensors::FastorTensors::Tensor3333t< T > > CoupledPotential(
      const Math::Tensors::FastorTensors::Tensor33t< T >& C,
      const Potential&                                    potential )
    {
      return Detail::coupledPotential< 2 >( C, potential );
    }

    /** @brief Single-pass decoupled potential of two scalar functions (see the plain
     * Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::DecoupledPotential) and its first and second
     * derivative w.r.t. C.
     *
     * @tparam T Scalar type, e.g. double, std::complex<double> or an automatic differentiation number.
     * @param C Right Cauchy-Green tensor.
     * @param isochoricPart Isochoric part, e.g. Isochoric::Yeoh.
     * @param volumetricPart Volumetric part, e.g. Volumetric::SimoTaylor.
     * @return Tuple of energy density and its first and second derivative.
     */
    template < typename T, IsochoricScalarFunction IsochoricPart, VolumetricScalarFunction VolumetricPart >
    std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T >, Math::Tensors::FastorTensors::Tensor3333t< T > > DecoupledPotential(
      const Math::Tensors::FastorTensors::Tensor33t< T >& C,
      const IsochoricPart&                                isochoricPart,
      const VolumetricPart&                               volumetricPart )
    {
      return Detail::decoupledPotential< 2 >( C, isochoricPart, volumetricPart );
    }

    /** @brief Mooney-Rivlin potential (see the plain
     * Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::MooneyRivlinPotential) and its first and
     * second derivative w.r.t. C.
     *
     * @tparam T Scalar type, e.g. double, std::complex<double>.
     * @param C Right Cauchy-Green tensor.
     * @param C1 Mooney-Rivlin material parameter C1.
     * @param C2 Mooney-Rivlin material parameter C2.
     * @param D1 Mooney-Rivlin material parameter D1.
     * @return A tuple containing the energy density and its first and second derivative w.r.t. C.
     */
    template < typename T >
    std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T >, Math::Tensors::FastorTensors::Tensor3333t< T > > MooneyRivlinPotential(
      const Math::Tensors::FastorTensors::Tensor33t< T >& C,
      const double                                        C1,
      const double                                        C2,
      const double                                        D1 )
    {
      return Detail::decoupledPotential< 2 >( C, Isochoric::MooneyRivlin{ C1, C2 }, Volumetric::Ogden{ 2. / D1, -2. } );
    }

    /** @brief Yeoh potential (see the plain
     * Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::YeohPotential) and its first and second
     * derivative w.r.t. C.
     *
     * @tparam T Scalar type, e.g. double, std::complex<double>.
     * @param C Right Cauchy-Green tensor.
     * @param C1 Yeoh material parameter.
     * @param C2 Yeoh material parameter.
     * @param C3 Yeoh material parameter.
     * @param D1 Volumetric penalty parameter.
     * @return A tuple containing the energy density and its first and second derivative w.r.t. C.
     */
    template < typename T >
    std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T >, Math::Tensors::FastorTensors::Tensor3333t< T > > YeohPotential(
      const Math::Tensors::FastorTensors::Tensor33t< T >& C,
      const double                                        C1,
      const double                                        C2,
      const double                                        C3,
      const double                                        D1 )
    {
      return Detail::decoupledPotential< 2 >( C, Isochoric::Yeoh{ C1, C2, C3 }, Volumetric::Ogden{ 2. / D1, -2. } );
    }

    /** @brief Compressible Arruda-Boyce 8-chain potential (see the plain
     * Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::ArrudaBoyce8ChainPotential) and its first
     * and second derivative w.r.t. C.
     *
     * @tparam T Scalar type, e.g. double, std::complex<double>.
     * @param C Right Cauchy-Green tensor.
     * @param mu Shear-modulus-like parameter.
     * @param lambdaL Locking stretch (must be > 1).
     * @param D1 Volumetric compliance parameter.
     * @return A tuple containing the energy density and its first and second derivative w.r.t. C.
     */
    template < typename T >
    std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T >, Math::Tensors::FastorTensors::Tensor3333t< T > > ArrudaBoyce8ChainPotential(
      const Math::Tensors::FastorTensors::Tensor33t< T >& C,
      const double                                        mu,
      const double                                        lambdaL,
      const double                                        D1 )
    {
      return Detail::decoupledPotential< 2 >( C,
                                              Isochoric::ArrudaBoyce8Chain{ mu, lambdaL },
                                              Volumetric::Ogden{ 2. / D1, -2. } );
    }

    /** @brief Hyperelastic Energy Density Function Wb acc. Pence & Gou (2015), Eq. (2.12) and its first and second
     * derivative w.r.t. C
     *
     * The energy density function \f$W_b\f$ is given as
     * \f[
     *   W_b = \frac{K}{8} \left(J - \frac{1}{J}\right)^2 + \frac{G}{2} \left(I_1 J^{-\frac{2}{3}} - 3\right)
     * \f]
     * where \f$ I_1 = \text{tr}(\boldsymbol{C}) \f$ is the first invariant of the right Cauchy-Green tensor
     * \f$ \boldsymbol{C} = \boldsymbol{F}^T \boldsymbol{F} \f$, \f$ J = \sqrt{\det(\boldsymbol{C})} =
     * \det(\boldsymbol{F}) \f$ is the determinant of the deformation gradient, and \f$ K, G \f$ are the bulk and
     * shear modulus, respectively.
     *
     * Additionally, the first and second derivative w.r.t. \f$\boldsymbol{C}\f$,
     * i.e.,
     * \f[
     * \frac{\partial W_b}{\partial \boldsymbol{C}} \quad \text{and} \quad
     * \frac{\partial^2 W_b}{\partial \boldsymbol{C} \partial \boldsymbol{C}}
     * \f]
     * are computed.
     *
     * @tparam T Scalar type, e.g. double, float, etc.
     * @param C Right Cauchy-Green tensor
     * @param K Bulk modulus
     * @param G Shear modulus
     * @return A tuple containing energy density, its first and second derivative w.r.t. C
     *
     */
    template < typename T >
    std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T >, Math::Tensors::FastorTensors::Tensor3333t< T > > PenceGouPotentialB(
      const Math::Tensors::FastorTensors::Tensor33t< T >& C,
      const double                                        K,
      const double                                        G )
    {
      const auto inv = Invariants::RightCauchyGreen::computeInvariants< 2 >( C, false );
      return Invariants::RightCauchyGreen::composeWithInvariants( Detail::penceGouB< 2 >( inv.x, K, G ), inv );
    }

    /** @brief Isotropic Biot-Neo-Hooke energy density and derivatives w.r.t. right stretch \f$\boldsymbol{U}\f$.
     *
     * @details Internally evaluates the standard Neo-Hooke potential in terms of
     * \f$\boldsymbol{C}=\boldsymbol{U}\boldsymbol{U}\f$ and applies chain-rule transformations.
     *
     * @tparam T Scalar type, e.g. double, float or an automatic differentiation number.
     * @param U Right stretch tensor.
     * @param K Bulk modulus.
     * @param G Shear modulus.
     * @return Tuple of energy density, first derivative and second derivative w.r.t. \f$\boldsymbol{U}\f$.
     */
    template < typename T >
    std::tuple< T, Marmot::Math::Tensors::FastorTensors::Tensor33t< T >, Marmot::Math::Tensors::FastorTensors::Tensor3333t< T > > BiotNeoHooke(
      const Marmot::Math::Tensors::FastorTensors::Tensor33t< T >& U,
      const double&                                               K,
      const double&                                               G )
    {
      using namespace Fastor;
      using namespace Marmot::Math::Tensors::FastorIndices;
      Math::Tensors::FastorTensors::Tensor3333t< T > I4 = Marmot::Math::Tensors::makeOtherScalarType< T >(
        Math::Tensors::FastorTensors::Spatial3D::I4 );
      Math::Tensors::FastorTensors::Tensor33t< T >   C        = U % U;
      Math::Tensors::FastorTensors::Tensor3333t< T > dC_dU    = 2 * einsum< iL, ijkl >( U, I4 );
      Math::Tensors::FastorTensors::Tensor3333t< T > d2C_dUdU = 2 * I4;

      const auto [psi, dPsi_dC, d2Psi_dCdCGeneral] = CoupledPotential( C, Coupled::StandardNeoHooke{ K, G } );

      // U and C are symmetric: the chain rule below needs the derivative w.r.t. a symmetric C, i.e. the
      // minor-symmetric part of the derivative w.r.t. the nine independent components of C
      Math::Tensors::FastorTensors::Tensor3333t< T > d2Psi_dCdC;
      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          for ( int k = 0; k < 3; k++ )
            for ( int l = 0; l < 3; l++ )
              d2Psi_dCdC( i, j, k, l ) = .25 * ( d2Psi_dCdCGeneral( i, j, k, l ) + d2Psi_dCdCGeneral( j, i, k, l ) +
                                                 d2Psi_dCdCGeneral( i, j, l, k ) + d2Psi_dCdCGeneral( j, i, l, k ) );

      Math::Tensors::FastorTensors::Tensor33t< T > dPsi_dU = einsum< kl, klmn >( dPsi_dC, dC_dU );
      Math::Tensors::FastorTensors::Tensor3333t< T >
        d2Psi_dUdU = einsum< ijkl, klmn >( einsum< ijkl, ijmn >( dC_dU, d2Psi_dCdC ), dC_dU ) +
                     einsum< jL, ijkl >( dPsi_dC, d2C_dUdU );

      return { psi, dPsi_dU, d2Psi_dUdU };
    }

  } // namespace SecondOrder
} // namespace Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions
