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

  /** @brief Standard compressible Neo-Hooke energy density function in terms of \f$\boldsymbol{C}\f$.
   *
   * @tparam T Scalar type, e.g. double, float, autodiff scalar.
   * @param C Right Cauchy-Green tensor.
   * @param K Bulk modulus.
   * @param G Shear modulus.
   * @return Energy density.
   */
  template < typename T >
  T standardNeoHooke( const Math::Tensors::FastorTensors::Tensor33t< T >& C, const double K, const double G )
  {
    using namespace Fastor;
    using namespace Marmot::Math::Tensors::FastorIndices;
    const double lambda = K - 2.0 / 3.0 * G;

    /*
     * we use the potential in terms of C
     * Psi = G/2 ( tr(C) - 3 - 2 ln(J) ) + lambda/2 ( 0.5 ( J^2 - 1 ) - ln(J) )
     * where J = sqrt( det(C) ) and ln(J) = 0.5 ln( det(C) )
     *
     * we can the rewrite as
     * Psi = G/2 ( tr(C) - 3 - ln(det(C)) ) + lambda/2 ( 0.5 ( det(C) - 1 ) - 0.5 ln(det(C)) )
     *
     */

    const T trC    = trace( C );
    const T detC   = determinant( C );
    const T lnDetC = log( detC );
    // energy density
    const T psi = G / 2 * ( trC - 3.0 - lnDetC ) + lambda / 4 * ( detC - 1.0 - lnDetC );

    return psi;
  }

  /** @cond */
  namespace Detail {
    /** @brief Shared closed-form energy and derivative for the Arruda-Boyce 8-chain
     * potential's isochoric part, as a function of the isochoric first invariant
     * \f$\bar I_1\f$ alone -- used by both the plain (energy-only) and
     * FirstOrder (energy + first derivative w.r.t. C) overloads of
     * Isochoric::ArrudaBoyce8ChainPotential below, so the two cannot independently drift out
     * of sync with each other.
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
     * @tparam T Scalar type, e.g. double, std::complex<double>, autodiff::dual3rd.
     * @param Ibar1 Isochoric first invariant, \f$\bar I_1 = I_1 J^{-2/3}\f$.
     * @param mu Shear-modulus-like parameter (small-strain limit -> standard
     * neo-Hookean shear modulus).
     * @param lambdaL Locking stretch (must be > 1).
     * @return Pair (Psi_iso, dPsi_iso/dIbar1).
     */
    template < typename T >
    std::pair< T, T > arrudaBoyce8ChainEnergyAndDerivative( const T& Ibar1, const double mu, const double lambdaL )
    {
      const double lambdaL2 = lambdaL * lambdaL;
      const double w0       = 1. - 1. / lambdaL2;
      const T      w        = w0 - ( Ibar1 - 3. ) / ( 3. * lambdaL2 );

      // w0 < 1 for lambdaL > 1, so the floor always lies below the reference state w = w0
      const double wFloor = 1e-2 * w0;

      if ( Math::Tensors::makeReal( w ) > wFloor ) {
        const T psi_iso     = mu / 6. * ( Ibar1 - 3. ) - mu * lambdaL2 * log( w / w0 );
        const T dPsi_dIbar1 = mu / 6. + mu / 3. / w;
        return { psi_iso, dPsi_dIbar1 };
      }

      // C1 linear continuation of ln(w) about wFloor: d/dIbar1 of the energy is the constant below
      const T psi_iso = mu / 6. * ( Ibar1 - 3. ) -
                        mu * lambdaL2 * ( std::log( wFloor / w0 ) + ( w - wFloor ) / wFloor );
      const T dPsi_dIbar1 = T( mu / 6. + mu / 3. / wFloor );
      return { psi_iso, dPsi_dIbar1 };
    }
    /** @brief Isochoric first invariant \f$\bar I_1 = I_1 J^{-2/3}\f$ and its derivative
     * \f$\partial \bar I_1/\partial \boldsymbol C = J^{-2/3}\boldsymbol I - \frac{\bar I_1}{3}\boldsymbol C^{-1}\f$.
     *
     * @tparam T Scalar type, e.g. double, std::complex<double>.
     * @param C Right Cauchy-Green tensor.
     * @return Pair (Ibar1, dIbar1/dC).
     */
    template < typename T >
    std::pair< T, Math::Tensors::FastorTensors::Tensor33t< T > > isochoricFirstInvariantAndDerivative(
      const Math::Tensors::FastorTensors::Tensor33t< T >& C )
    {
      using namespace Fastor;
      const T Jm23  = pow( determinant( C ), -1. / 3. );
      const T Ibar1 = trace( C ) * Jm23;

      const Math::Tensors::FastorTensors::Tensor33t< T > I = Marmot::Math::Tensors::fastorTensorFromDoubleTensor< T >(
        Math::Tensors::FastorTensors::Spatial3D::I );
      const Math::Tensors::FastorTensors::Tensor33t< T > CInv = inverse( C );
      return { Ibar1,
               Math::Tensors::multiplyFastorTensorWithScalar( I, Jm23 ) -
                 Math::Tensors::multiplyFastorTensorWithScalar( CInv, T( Ibar1 / 3. ) ) };
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
      using namespace Fastor;
      using namespace Marmot::Math::Tensors::FastorIndices;

      const T J   = sqrt( determinant( C ) );
      const T I1  = trace( C );
      const T I1_ = I1 * pow( J, -2. / 3. );
      const T I2_ = 0.5 * ( I1 * I1 - trace( C % C ) ) * pow( J, -4. / 3. );
      T       res = C10 * ( I1_ - 3. ) + C01 * ( I2_ - 3. );

      return res;
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
     * @tparam T Scalar type, e.g. double, float, autodiff scalar.
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
      const T J         = sqrt( determinant( C ) );
      const T I1        = trace( C );
      const T I1_minus3 = I1 * pow( J, -2. / 3. ) - 3.;
      T       res       = C10 * I1_minus3 + C20 * I1_minus3 * I1_minus3 + C30 * I1_minus3 * I1_minus3 * I1_minus3;
      return res;
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
     * @tparam T Scalar type, e.g. double, float, autodiff scalar.
     * @param C Right Cauchy-Green tensor.
     * @param mu Shear modulus.
     * @return Energy density (isochoric part only).
     */
    template < typename T >
    T NeoHookePotential( const Math::Tensors::FastorTensors::Tensor33t< T >& C, const double mu )
    {
      const T J     = sqrt( determinant( C ) );
      const T Ibar1 = trace( C ) * pow( J, -2. / 3. );
      return mu / 2. * ( Ibar1 - 3. );
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
     * part of #standardNeoHooke as \f$\lambda_L\to\infty\f$, and is exactly
     * stress-free (zero energy and zero gradient) at \f$\boldsymbol C=\boldsymbol I\f$
     * for any finite \f$\lambda_L\f$, since \f$\bar I_1\f$ is invariant under
     * \f$\boldsymbol C \to \lambda\boldsymbol C\f$ for any scalar \f$\lambda\f$.
     *
     * @tparam T Scalar type, e.g. double, float, autodiff scalar.
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
      using namespace Fastor;
      const T J     = sqrt( determinant( C ) );
      const T I1    = trace( C );
      const T Ibar1 = I1 * pow( J, -2. / 3 );

      return Detail::arrudaBoyce8ChainEnergyAndDerivative( Ibar1, mu, lambdaL ).first;
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
        const auto [Ibar1, dIbar1_dC] = Detail::isochoricFirstInvariantAndDerivative( C );

        const T psi = mu / 2. * ( Ibar1 - 3. );
        const Math::Tensors::FastorTensors::Tensor33t< T >
          dPsi_dC = Math::Tensors::multiplyFastorTensorWithScalar( dIbar1_dC, T( mu / 2. ) );
        return { psi, dPsi_dC };
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
        const auto [Ibar1, dIbar1_dC] = Detail::isochoricFirstInvariantAndDerivative( C );
        const T Ibar1m3               = Ibar1 - 3.;

        const T psi         = C10 * Ibar1m3 + C20 * Ibar1m3 * Ibar1m3 + C30 * Ibar1m3 * Ibar1m3 * Ibar1m3;
        const T dPsi_dIbar1 = C10 + 2. * C20 * Ibar1m3 + 3. * C30 * Ibar1m3 * Ibar1m3;
        const Math::Tensors::FastorTensors::Tensor33t< T >
          dPsi_dC = Math::Tensors::multiplyFastorTensorWithScalar( dIbar1_dC, dPsi_dIbar1 );
        return { psi, dPsi_dC };
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
        using namespace Fastor;
        const T J     = sqrt( determinant( C ) );
        const T Jm23  = pow( J, -2. / 3. );
        const T Jm43  = pow( J, -4. / 3. );
        const T I1    = trace( C );
        const T Ibar1 = I1 * Jm23;
        const T I2    = 0.5 * ( I1 * I1 - trace( C % C ) );
        const T Ibar2 = I2 * Jm43;

        const Math::Tensors::FastorTensors::Tensor33t< T > I = Marmot::Math::Tensors::fastorTensorFromDoubleTensor< T >(
          Math::Tensors::FastorTensors::Spatial3D::I );
        const Math::Tensors::FastorTensors::Tensor33t< T > CInv = inverse( C );

        const Math::Tensors::FastorTensors::Tensor33t< T >
          dIbar1_dC = Math::Tensors::multiplyFastorTensorWithScalar( I, Jm23 ) -
                      Math::Tensors::multiplyFastorTensorWithScalar( CInv, Ibar1 / 3. );
        // dIbar2/dC = J^(-4/3)*I1*I - J^(-4/3)*C - (2/3)*Ibar2*Cinv (standard
        // isochoric second-invariant derivative identity, I2's own derivative
        // dI2/dC = I1*I - C combined with the J^(-4/3) scaling's own C-derivative).
        const Math::Tensors::FastorTensors::Tensor33t< T >
          dIbar2_dC = Math::Tensors::multiplyFastorTensorWithScalar( I, T( Jm43 * I1 ) ) -
                      Math::Tensors::multiplyFastorTensorWithScalar( C, Jm43 ) -
                      Math::Tensors::multiplyFastorTensorWithScalar( CInv, T( 2. / 3. * Ibar2 ) );

        const T psi = C10 * ( Ibar1 - 3. ) + C01 * ( Ibar2 - 3. );
        const Math::Tensors::FastorTensors::Tensor33t< T >
          dPsi_dC = Math::Tensors::multiplyFastorTensorWithScalar( dIbar1_dC, T( C10 ) ) +
                    Math::Tensors::multiplyFastorTensorWithScalar( dIbar2_dC, T( C01 ) );
        return { psi, dPsi_dC };
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
        const auto [Ibar1, dIbar1_dC]     = Detail::isochoricFirstInvariantAndDerivative( C );
        const auto [psi_iso, dPsi_dIbar1] = Detail::arrudaBoyce8ChainEnergyAndDerivative( Ibar1, mu, lambdaL );

        const Math::Tensors::FastorTensors::Tensor33t< T >
          dPsi_dC = Math::Tensors::multiplyFastorTensorWithScalar( dIbar1_dC, dPsi_dIbar1 );

        return { psi_iso, dPsi_dC };
      }
    } // namespace FirstOrder
  }   // namespace Isochoric

  /** @cond */
  namespace Detail {
    /* Volumetric strain energy functions U(J) acc. Doll & Schweizerhof (2000), J. Appl. Mech. 67, 17-21,
     * https://doi.org/10.1115/1.321146, as scalar functions of J. Each returns the pair (U, dU/dJ); the plain
     * and FirstOrder C-overloads below are thin wrappers around these, so they cannot drift apart. */

    template < typename T >
    std::pair< T, T > sussmanBatheVolumetricEnergyAndDerivative( const T& J, const double kappa )
    {
      return { kappa / 2. * ( J - 1. ) * ( J - 1. ), kappa * ( J - 1. ) };
    }

    template < typename T >
    std::pair< T, T > simoTaylorVolumetricEnergyAndDerivative( const T& J, const double kappa )
    {
      const T lnJ = log( J );
      return { kappa / 4. * ( ( J - 1. ) * ( J - 1. ) + lnJ * lnJ ), kappa / 2. * ( ( J - 1. ) + lnJ / J ) };
    }

    template < typename T >
    std::pair< T, T > ogdenVolumetricEnergyAndDerivative( const T& J, const double kappa, const double theta )
    {
      const T JmTheta = pow( J, -theta );
      return { kappa / ( theta * theta ) * ( theta * log( J ) + JmTheta - 1. ), kappa / theta * ( 1. - JmTheta ) / J };
    }

    template < typename T >
    std::pair< T, T > liuHofstetterMangVolumetricEnergyAndDerivative( const T& J, const double kappa )
    {
      const T lnJ = log( J );
      return { kappa * ( J * lnJ - J + 1. ), kappa * lnJ };
    }

    template < typename T >
    std::pair< T, T > generalizedDollSchweizerhofVolumetricEnergyAndDerivative( const T&                        J,
                                                                                const double                    kappa,
                                                                                const std::span< const double > alphas,
                                                                                const std::span< const double > betas )
    {
      if ( alphas.size() != betas.size() || alphas.empty() )
        throw std::invalid_argument( "Generalized Doll-Schweizerhof volumetric potential: alphas and betas must be "
                                     "non-empty and of equal size." );

      double sumAlphaPlusBeta = 0.;
      T      U                = T( 0. );
      T      dU_dJ            = T( 0. );
      for ( std::size_t i = 0; i < alphas.size(); i++ ) {
        const double a = alphas[i];
        const double b = betas[i];
        sumAlphaPlusBeta += a + b;
        U += pow( J, a + 1. ) / ( a + 1. ) + pow( J, 1. - b ) / ( b - 1. ) - 1. / ( a + 1. ) - 1. / ( b - 1. );
        dU_dJ += pow( J, a ) - pow( J, -b );
      }
      return { kappa / sumAlphaPlusBeta * U, kappa / sumAlphaPlusBeta * dU_dJ };
    }

    template < typename T >
    std::pair< T, T > dollSchweizerhofExponentialVolumetricEnergyAndDerivative( const T& J, const double kappa )
    {
      const T expJm1 = exp( J - 1. );
      return { kappa / 2. * ( expJm1 - log( J ) - 1. ), kappa / 2. * ( expJm1 - 1. / J ) };
    }

    template < typename T >
    std::pair< T, T > dollSchweizerhofLogarithmicVolumetricEnergyAndDerivative( const T& J, const double kappa )
    {
      const T lnJ = log( J );
      return { kappa / 2. * ( J - 1. ) * lnJ, kappa / 2. * ( lnJ + ( J - 1. ) / J ) };
    }

    /* Chain rule dU/dC = dU/dJ dJ/dC with dJ/dC = J/2 C^{-T}. */
    template < typename T >
    std::tuple< T, Math::Tensors::FastorTensors::Tensor33t< T > > volumetricFirstOrderFromJ(
      const Math::Tensors::FastorTensors::Tensor33t< T >& C,
      const T&                                            J,
      const std::pair< T, T >&                            U_dUdJ )
    {
      using namespace Fastor;
      const Math::Tensors::FastorTensors::Tensor33t< T > CInvT = transpose( inverse( C ) );
      return { U_dUdJ.first, Math::Tensors::multiplyFastorTensorWithScalar( CInvT, T( U_dUdJ.second * J / 2. ) ) };
    }
  } // namespace Detail
  /** @endcond */

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
     * @tparam T Scalar type, e.g. double, float, autodiff scalar.
     * @param C Right Cauchy-Green tensor.
     * @param kappa Bulk modulus.
     * @return Energy density (volumetric part only).
     */
    template < typename T >
    T PenaltyPotential( const Math::Tensors::FastorTensors::Tensor33t< T >& C, const double kappa )
    {
      const T lnDetC = log( determinant( C ) );
      return kappa / 8. * lnDetC * lnDetC;
    }

    /** @brief Quadratic volumetric strain energy function, \f$U_1\f$ in Table 1 of Doll & Schweizerhof (2000),
     * https://doi.org/10.1115/1.321146 (Sussman & Bathe 1987, Simo 1988, among others),
     * \f[
     *   U_1(J) = \frac{\kappa}{2}(J-1)^2, \qquad J = \sqrt{\det\boldsymbol C}.
     * \f]
     * Note that \f$U_1\f$ stays bounded as \f$J\to 0\f$ (violating condition (9) of the paper), so it should not be
     * used for large compression.
     *
     * @tparam T Scalar type, e.g. double, std::complex<double>, autodiff scalar.
     * @param C Right Cauchy-Green tensor.
     * @param kappa Bulk modulus.
     * @return Energy density (volumetric part only).
     */
    template < typename T >
    T SussmanBathePotential( const Math::Tensors::FastorTensors::Tensor33t< T >& C, const double kappa )
    {
      return Detail::sussmanBatheVolumetricEnergyAndDerivative( T( sqrt( determinant( C ) ) ), kappa ).first;
    }

    /** @brief Volumetric strain energy function \f$U_3\f$ in Table 1 of Doll & Schweizerhof (2000),
     * https://doi.org/10.1115/1.321146 (Simo & Taylor 1982), the mean of \f$U_1\f$
     * (#SussmanBathePotential) and \f$U_2\f$ (#PenaltyPotential),
     * \f[
     *   U_3(J) = \frac{\kappa}{4}\left[(J-1)^2 + (\ln J)^2\right], \qquad J = \sqrt{\det\boldsymbol C}.
     * \f]
     *
     * @tparam T Scalar type, e.g. double, std::complex<double>, autodiff scalar.
     * @param C Right Cauchy-Green tensor.
     * @param kappa Bulk modulus.
     * @return Energy density (volumetric part only).
     */
    template < typename T >
    T SimoTaylorPotential( const Math::Tensors::FastorTensors::Tensor33t< T >& C, const double kappa )
    {
      return Detail::simoTaylorVolumetricEnergyAndDerivative( T( sqrt( determinant( C ) ) ), kappa ).first;
    }

    /** @brief Ogden (1972) volumetric strain energy function, \f$U_4\f$ in Table 1 of Doll & Schweizerhof (2000),
     * https://doi.org/10.1115/1.321146,
     * \f[
     *   U_4(J) = \frac{\kappa}{\theta^2}\left(\theta\ln J + J^{-\theta} - 1\right), \qquad J = \sqrt{\det\boldsymbol
     * C}. \f] All physical conditions of the paper hold only for \f$\theta < -1\f$; the popular choice \f$\theta=-2\f$
     * gives \f$\frac{\kappa}{4}(J^2 - 1 - 2\ln J)\f$ (Simo & Taylor 1991), and Ogden's rubber fit \f$\theta=9\f$ is
     * admissible for nearly incompressible applications. \f$\theta\neq0\f$ is required.
     *
     * @tparam T Scalar type, e.g. double, std::complex<double>, autodiff scalar.
     * @param C Right Cauchy-Green tensor.
     * @param kappa Bulk modulus.
     * @param theta Shape parameter \f$\theta\f$ (\f$\neq 0\f$).
     * @return Energy density (volumetric part only).
     */
    template < typename T >
    T OgdenPotential( const Math::Tensors::FastorTensors::Tensor33t< T >& C, const double kappa, const double theta )
    {
      return Detail::ogdenVolumetricEnergyAndDerivative( T( sqrt( determinant( C ) ) ), kappa, theta ).first;
    }

    /** @brief Volumetric strain energy function \f$U_5\f$ in Table 1 of Doll & Schweizerhof (2000),
     * https://doi.org/10.1115/1.321146 (Liu, Hofstetter & Mang 1992, 1994),
     * \f[
     *   U_5(J) = \kappa\left(J\ln J - J + 1\right), \qquad J = \sqrt{\det\boldsymbol C}.
     * \f]
     * Note that \f$U_5\f$ stays bounded as \f$J\to 0\f$ (violating condition (9) of the paper), so it should not be
     * used for large compression.
     *
     * @tparam T Scalar type, e.g. double, std::complex<double>, autodiff scalar.
     * @param C Right Cauchy-Green tensor.
     * @param kappa Bulk modulus.
     * @return Energy density (volumetric part only).
     */
    template < typename T >
    T LiuHofstetterMangPotential( const Math::Tensors::FastorTensors::Tensor33t< T >& C, const double kappa )
    {
      return Detail::liuHofstetterMangVolumetricEnergyAndDerivative( T( sqrt( determinant( C ) ) ), kappa ).first;
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
     * @tparam T Scalar type, e.g. double, std::complex<double>, autodiff scalar.
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
      return Detail::generalizedDollSchweizerhofVolumetricEnergyAndDerivative( T( sqrt( determinant( C ) ) ),
                                                                               kappa,
                                                                               alphas,
                                                                               betas )
        .first;
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
     * @tparam T Scalar type, e.g. double, std::complex<double>, autodiff scalar.
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
      return Detail::generalizedDollSchweizerhofVolumetricEnergyAndDerivative( T( sqrt( determinant( C ) ) ),
                                                                               kappa,
                                                                               alphas,
                                                                               betas )
        .first;
    }

    /** @brief Exponential volumetric strain energy function \f$U_7\f$ acc. Doll & Schweizerhof (2000),
     * https://doi.org/10.1115/1.321146, Eq. (16)\f$_1\f$,
     * \f[
     *   U_7(J) = \frac{\kappa}{2}\left(e^{J-1} - \ln J - 1\right), \qquad J = \sqrt{\det\boldsymbol C}.
     * \f]
     *
     * @tparam T Scalar type, e.g. double, std::complex<double>, autodiff scalar.
     * @param C Right Cauchy-Green tensor.
     * @param kappa Bulk modulus.
     * @return Energy density (volumetric part only).
     */
    template < typename T >
    T DollSchweizerhofExponentialPotential( const Math::Tensors::FastorTensors::Tensor33t< T >& C, const double kappa )
    {
      return Detail::dollSchweizerhofExponentialVolumetricEnergyAndDerivative( T( sqrt( determinant( C ) ) ), kappa )
        .first;
    }

    /** @brief Logarithmic volumetric strain energy function \f$U_8\f$ acc. Doll & Schweizerhof (2000),
     * https://doi.org/10.1115/1.321146, Eq. (16)\f$_2\f$,
     * \f[
     *   U_8(J) = \frac{\kappa}{2}(J-1)\ln J, \qquad J = \sqrt{\det\boldsymbol C}.
     * \f]
     *
     * @tparam T Scalar type, e.g. double, std::complex<double>, autodiff scalar.
     * @param C Right Cauchy-Green tensor.
     * @param kappa Bulk modulus.
     * @return Energy density (volumetric part only).
     */
    template < typename T >
    T DollSchweizerhofLogarithmicPotential( const Math::Tensors::FastorTensors::Tensor33t< T >& C, const double kappa )
    {
      return Detail::dollSchweizerhofLogarithmicVolumetricEnergyAndDerivative( T( sqrt( determinant( C ) ) ), kappa )
        .first;
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
        using namespace Fastor;
        const T                                            lnDetC = log( determinant( C ) );
        const Math::Tensors::FastorTensors::Tensor33t< T > CInv   = inverse( C );
        const T                                            psi    = kappa / 8. * lnDetC * lnDetC;
        const Math::Tensors::FastorTensors::Tensor33t< T >
          dPsi_dC = Math::Tensors::multiplyFastorTensorWithScalar( CInv, T( kappa / 4. * lnDetC ) );
        return { psi, dPsi_dC };
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
        const T J = sqrt( determinant( C ) );
        return Detail::volumetricFirstOrderFromJ( C, J, Detail::sussmanBatheVolumetricEnergyAndDerivative( J, kappa ) );
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
        const T J = sqrt( determinant( C ) );
        return Detail::volumetricFirstOrderFromJ( C, J, Detail::simoTaylorVolumetricEnergyAndDerivative( J, kappa ) );
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
        const T J = sqrt( determinant( C ) );
        return Detail::volumetricFirstOrderFromJ( C, J, Detail::ogdenVolumetricEnergyAndDerivative( J, kappa, theta ) );
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
        const T J = sqrt( determinant( C ) );
        return Detail::volumetricFirstOrderFromJ( C,
                                                  J,
                                                  Detail::liuHofstetterMangVolumetricEnergyAndDerivative( J, kappa ) );
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
        const T      J        = sqrt( determinant( C ) );
        return Detail::
          volumetricFirstOrderFromJ( C,
                                     J,
                                     Detail::generalizedDollSchweizerhofVolumetricEnergyAndDerivative( J,
                                                                                                       kappa,
                                                                                                       alphas,
                                                                                                       betas ) );
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
        const T J = sqrt( determinant( C ) );
        return Detail::
          volumetricFirstOrderFromJ( C,
                                     J,
                                     Detail::generalizedDollSchweizerhofVolumetricEnergyAndDerivative( J,
                                                                                                       kappa,
                                                                                                       alphas,
                                                                                                       betas ) );
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
        const T J = sqrt( determinant( C ) );
        return Detail::
          volumetricFirstOrderFromJ( C,
                                     J,
                                     Detail::dollSchweizerhofExponentialVolumetricEnergyAndDerivative( J, kappa ) );
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
        const T J = sqrt( determinant( C ) );
        return Detail::
          volumetricFirstOrderFromJ( C,
                                     J,
                                     Detail::dollSchweizerhofLogarithmicVolumetricEnergyAndDerivative( J, kappa ) );
      }
    } // namespace FirstOrder
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
   * @tparam T Scalar type, e.g. double, std::complex<double>, autodiff scalar.
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
    return DecoupledPotential(
      C,
      [C1, C2]( const auto& Ct ) { return Isochoric::MooneyRivlinPotential( Ct, C1, C2 ); },
      [D1]( const auto& Ct ) { return Volumetric::OgdenPotential( Ct, 2. / D1, -2. ); } );
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
   * @tparam T Scalar type, e.g. double, float, autodiff scalar.
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
    return DecoupledPotential(
      C,
      [C1, C2, C3]( const auto& Ct ) { return Isochoric::YeohPotential( Ct, C1, C2, C3 ); },
      [D1]( const auto& Ct ) { return Volumetric::OgdenPotential( Ct, 2. / D1, -2. ); } );
  }

  namespace FirstOrder {

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
      using namespace Fastor;
      using namespace Marmot::Math::Tensors::FastorIndices;

      const T J  = sqrt( determinant( C ) );
      const T I1 = trace( C );
      // energy density
      T psi = K / 8. * pow( J - 1. / J, 2. ) + G / 2. * ( I1 * pow( J, -2. / 3 ) - 3. );

      // first derivative w.r.t. C
      const T dPsi_dJ  = K / 4. * ( J - 1. / J ) * ( 1. + 1. / ( J * J ) ) - G / 3. * I1 * pow( J, -5. / 3. );
      const T dPsi_dI1 = G / 2. * pow( J, -2. / 3. );

      const Math::Tensors::FastorTensors::Tensor33t< T > CInv = inverse( C );
      const Math::Tensors::FastorTensors::Tensor33t< T >
        dJ_dC = Marmot::Math::Tensors::multiplyFastorTensorWithScalar( transpose( CInv ), T( J / 2. ) );
      const Math::Tensors::FastorTensors::Tensor33t< T >
        dI1_dC = Marmot::Math::Tensors::fastorTensorFromDoubleTensor< T >( Math::Tensors::FastorTensors::Spatial3D::I );

      Math::Tensors::FastorTensors::Tensor33t< T >
        dPsi_dC = Marmot::Math::Tensors::multiplyFastorTensorWithScalar( dJ_dC, dPsi_dJ ) +
                  Marmot::Math::Tensors::multiplyFastorTensorWithScalar( dI1_dC, dPsi_dI1 );

      return { psi, dPsi_dC };
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
      return DecoupledPotential(
        C,
        [C1, C2]( const auto& Ct ) { return Isochoric::FirstOrder::MooneyRivlinPotential( Ct, C1, C2 ); },
        [D1]( const auto& Ct ) { return Volumetric::FirstOrder::OgdenPotential( Ct, 2. / D1, -2. ); } );
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
      return DecoupledPotential(
        C,
        [C1, C2, C3]( const auto& Ct ) { return Isochoric::FirstOrder::YeohPotential( Ct, C1, C2, C3 ); },
        [D1]( const auto& Ct ) { return Volumetric::FirstOrder::OgdenPotential( Ct, 2. / D1, -2. ); } );
    }
  } // namespace FirstOrder

  namespace SecondOrder {

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
      using namespace Fastor;
      using namespace Marmot::Math::Tensors::FastorIndices;

      const T J  = sqrt( determinant( C ) );
      const T I1 = trace( C );
      // energy density
      T psi = K / 8. * pow( J - 1. / J, 2. ) + G / 2. * ( I1 * pow( J, -2. / 3 ) - 3. );

      // first derivative w.r.t. C
      const T dPsi_dJ  = K / 4. * ( J - 1. / J ) * ( 1 + 1. / ( J * J ) ) - G / 3. * I1 * pow( J, -5. / 3. );
      const T dPsi_dI1 = G / 2. * pow( J, -2. / 3. );

      const Math::Tensors::FastorTensors::Tensor33t< T > CInv   = inverse( C );
      const Math::Tensors::FastorTensors::Tensor33t< T > dJ_dC  = 0.5 * J * transpose( CInv );
      const Math::Tensors::FastorTensors::Tensor33t< T > dI1_dC = Math::Tensors::FastorTensors::Spatial3D::I;

      Math::Tensors::FastorTensors::Tensor33t< T > dPsi_dC = dPsi_dJ * dJ_dC + dPsi_dI1 * dI1_dC;

      // second derivative w.r.t. C
      const T d2Psi_dJdJ  = K / 4. * ( 1. + 3. / ( J * J * J * J ) ) + 5. / 9. * G * I1 * pow( J, -8. / 3. );
      const T d2Psi_dJdI1 = -G / 3. * pow( J, -5. / 3. );
      Math::Tensors::FastorTensors::Tensor3333t< T > d2J_dCdC = J / 4. * einsum< JI, LK, to_IJKL >( CInv, CInv ) -
                                                                J / 2. * einsum< JK, LI, to_IJKL >( CInv, CInv );

      Math::Tensors::FastorTensors::Tensor3333t< T > d2Psi_dCdC = d2Psi_dJdJ * einsum< IJ, KL >( dJ_dC, dJ_dC ) +
                                                                  dPsi_dJ * d2J_dCdC +
                                                                  d2Psi_dJdI1 * ( einsum< IJ, KL >( dJ_dC, dI1_dC ) +
                                                                                  einsum< IJ, KL >( dI1_dC, dJ_dC ) );

      return { psi, dPsi_dC, d2Psi_dCdC };
    }

    /** @brief Standard compressible Neo-Hooke energy density and first/second derivatives w.r.t.
     * \f$\boldsymbol{C}\f$.
     *
     * @tparam T Scalar type, e.g. double, float, autodiff scalar.
     * @param C Right Cauchy-Green tensor.
     * @param K Bulk modulus.
     * @param G Shear modulus.
     * @return Tuple of energy density, first derivative and second derivative w.r.t. \f$\boldsymbol{C}\f$.
     */
    template < typename T >
    std::tuple< T, Marmot::Math::Tensors::FastorTensors::Tensor33t< T >, Marmot::Math::Tensors::FastorTensors::Tensor3333t< T > > standardNeoHooke(
      const Marmot::Math::Tensors::FastorTensors::Tensor33t< T >& C,
      const double&                                               K,
      const double&                                               G )
    {
      using namespace Fastor;
      using namespace Marmot::Math::Tensors::FastorIndices;
      const double lambda = K - 2.0 / 3.0 * G;

      /*
       * we use the potential in terms of C
       * Psi = G/2 ( tr(C) - 3 - 2 ln(J) ) + lambda/2 ( 0.5 ( J^2 - 1 ) - ln(J) )
       * where J = sqrt( det(C) ) and ln(J) = 0.5 ln( det(C) )
       *
       * we can the rewrite as
       * Psi = G/2 ( tr(C) - 3 - ln(det(C)) ) + lambda/2 ( 0.5 ( det(C) - 1 ) - 0.5 ln(det(C)) )
       *
       */

      const T trC    = trace( C );
      const T detC   = determinant( C );
      const T lnDetC = log( detC );
      // energy density
      const T psi = G / 2 * ( trC - 3.0 - lnDetC ) + lambda / 4 * ( detC - 1.0 - lnDetC );

      // first derivative quantities
      const Math::Tensors::FastorTensors::Tensor33t< T >& dTrC_dC = Marmot::Math::Tensors::makeOtherScalarType< T >(
        Math::Tensors::FastorTensors::Spatial3D::I );
      const Math::Tensors::FastorTensors::Tensor33t< T > invCt = transpose( inverse( C ) );
      const Math::Tensors::FastorTensors::Tensor33t< T >
        dDetC_dC = Marmot::Math::Tensors::multiplyFastorTensorWithScalar( invCt, detC );
      const Math::Tensors::FastorTensors::Tensor33t< T >& dLnDetC_dC = invCt;

      // first derivative with respect to C
      const Math::Tensors::FastorTensors::Tensor33t< T > dPsi_dC = G / 2 * ( dTrC_dC - dLnDetC_dC ) +
                                                                   lambda / 4 * ( dDetC_dC - dLnDetC_dC );

      // second derivative quantities
      using lj                                                       = Index< l_, j_ >;
      using li                                                       = Index< l_, i_ >;
      const Math::Tensors::FastorTensors::Tensor3333t< T > dInvCt_dC = -0.5 *
                                                                       ( einsum< ik, lj, to_ijkl >( invCt, invCt ) +
                                                                         einsum< jk, li, to_ijkl >( invCt, invCt ) );
      const Math::Tensors::FastorTensors::Tensor3333t< T >
        d2DetC_dC2 = Marmot::Math::Tensors::multiplyFastorTensorWithScalar( dInvCt_dC, detC ) +
                     einsum< ij, kl, to_ijkl >( invCt, dDetC_dC );
      const Math::Tensors::FastorTensors::Tensor3333t< T >& d2LnDetC_dC2 = dInvCt_dC;

      // second derivative with respect to C
      const Math::Tensors::FastorTensors::Tensor3333t< T > d2Psi_dC2 = lambda / 4.0 * d2DetC_dC2 -
                                                                       ( G / 2.0 + lambda / 4.0 ) * d2LnDetC_dC2;
      return { psi, dPsi_dC, d2Psi_dC2 };
    }

    /** @brief Isotropic Biot-Neo-Hooke energy density and derivatives w.r.t. right stretch \f$\boldsymbol{U}\f$.
     *
     * @details Internally evaluates the standard Neo-Hooke potential in terms of
     * \f$\boldsymbol{C}=\boldsymbol{U}\boldsymbol{U}\f$ and applies chain-rule transformations.
     *
     * @tparam T Scalar type, e.g. double, float, autodiff scalar.
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

      auto [psi, dPsi_dC, d2Psi_dCdC] = standardNeoHooke< T >( C, K, G );

      Math::Tensors::FastorTensors::Tensor33t< T > dPsi_dU = einsum< kl, klmn >( dPsi_dC, dC_dU );
      Math::Tensors::FastorTensors::Tensor3333t< T >
        d2Psi_dUdU = einsum< ijkl, klmn >( einsum< ijkl, ijmn >( dC_dU, d2Psi_dCdC ), dC_dU ) +
                     einsum< jL, ijkl >( dPsi_dC, d2C_dUdU );

      return { psi, dPsi_dU, d2Psi_dUdU };
    }

  } // namespace SecondOrder

  namespace ThirdOrder {
    /** @brief Standard compressible Neo-Hooke energy density and first/second/third derivatives w.r.t.
     * \f$\boldsymbol{C}\f$.
     *
     * @param C Right Cauchy-Green tensor.
     * @param K Bulk modulus.
     * @param G Shear modulus.
     * @return Tuple of energy density, first derivative, second derivative and third derivative.
     */
    std::tuple< double,
                Marmot::Math::Tensors::FastorTensors::Tensor33d,
                Marmot::Math::Tensors::FastorTensors::Tensor3333d,
                Marmot::Math::Tensors::FastorTensors::Tensor333333d >
    standardNeoHooke( const Marmot::Math::Tensors::FastorTensors::Tensor33d& C, const double& K, const double& G );
  } // namespace ThirdOrder
} // namespace Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions
