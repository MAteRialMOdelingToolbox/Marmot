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
#include "Marmot/MarmotEigenSystems.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMath.h"

namespace Marmot::ContinuumMechanics {

  namespace EnergyDensityFunctions {

    using namespace Fastor;
    using namespace FastorStandardTensors;

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
    T PenceGouPotentialA( const Tensor33t< T >& C, const double K, const double G )
    {

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
    T PenceGouPotentialB( const Tensor33t< T >& C, const double K, const double G )
    {

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
    T PenceGouPotentialC( const Tensor33t< T >& C, const double K, const double G )
    {

      const T J  = sqrt( determinant( C ) );
      const T I1 = trace( C );

      T res = G / 2. * ( I1 - 3. ) + 3. * G * G / ( 3. * K - 2. * G ) * ( pow( J, 2. / 3 - K / G ) - 1 );

      return res;
    }

    /**
     * @brief Mooney-Rivlin Hyperelastic Energy Density Function
     * The energy density function \f$W_{MR}\f$ is given as
     * \f[
     *   W_{MR} = C_1 (\bar{I}_1 - 3) + C_2 (\bar{I}_2 - 3) + \frac{1}{D_1} (J - 1)^2
     * \f]
     * where \f$ \bar{I}_1 = I_1 J^{-\frac{2}{3}} \f$ and \f$ \bar{I}_2 = I_2 J^{-\frac{4}{3}} \f$ are the first
     * and second invariant of the isochoric right Cauchy-Green tensor, respectively, \f$ I_1 =
     * \text{tr}(\boldsymbol{C}) \f$ and \f$ I_2 = 0.5 (I_1^2 - \text{tr}(\boldsymbol{C}^2)) \f$ are the first and
     * second invariant of the right Cauchy-Green tensor \f$ \boldsymbol{C} = \boldsymbol{F}^T \boldsymbol{F} \f$, \f$ J
     * = \sqrt{\det(\boldsymbol{C})} = \det(\boldsymbol{F}) \f$ is the determinant of the deformation gradient, and \f$
     * C_1, C_2, D_1 \f$ are material parameters.
     *
     * @tparam T Scalar type, e.g. double, float, etc.
     * @param C Right Cauchy-Green tensor
     * @param C1 Mooney-Rivlin material parameter C1
     * @param C2 Mooney-Rivlin material parameter C2
     * @param D1 Mooney-Rivlin material parameter D1
     * @return Energy density
     */
    /** @brief Compressible Mooney-Rivlin hyperelastic energy density function
     * (isochoric part only -- no volumetric term; see #VolumetricPenaltyPotential
     * for the shared volumetric term every potential in this file now uses).
     *
     *  The isochoric energy density is
     *  \f[
     *    W_{\rm iso} = C_{10} (\bar{I}_1 - 3) + C_{01} (\bar{I}_2 - 3)
     *  \f]
     *  where \f$ \bar{I}_1 = I_1 J^{-2/3} \f$ and \f$ \bar{I}_2 = I_2 J^{-4/3} \f$ are the first
     *  and second invariant of the isochoric right Cauchy-Green tensor, respectively, \f$ I_1 =
     *  \text{tr}(\boldsymbol{C}) \f$ and \f$ I_2 = 0.5 (I_1^2 - \text{tr}(\boldsymbol{C}^2)) \f$ are the first and
     *  second invariant of the right Cauchy-Green tensor \f$ \boldsymbol{C} = \boldsymbol{F}^T \boldsymbol{F} \f$, \f$ J
     *  = \sqrt{\det(\boldsymbol{C})} = \det(\boldsymbol{F}) \f$ is the determinant of the deformation gradient.
     *
     * @tparam T Scalar type, e.g. double, float, etc.
     * @param C Right Cauchy-Green tensor
     * @param C10 Mooney-Rivlin material parameter C10
     * @param C01 Mooney-Rivlin material parameter C01
     * @return Energy density (isochoric part only).
     */
    template < typename T >
    T MooneyRivlinPotential( const Tensor33t< T >& C, const double C10, const double C01 )
    {

      const T J   = sqrt( determinant( C ) );
      const T I1  = trace( C );
      const T I1_ = I1 * pow( J, -2. / 3. );
      const T I2_ = 0.5 * ( I1 * I1 - trace( C % C ) ) * pow( J, -4. / 3. );
      T       res = C10 * ( I1_ - 3. ) + C01 * ( I2_ - 3. );

      return res;
    }

    /** @brief Yeoh hyperelastic energy density function (isochoric part only --
     * no volumetric term; see #VolumetricPenaltyPotential for the shared
     * volumetric term every potential in this file now uses).
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
    T YeohPotential( const Tensor33t< T >& C, const double C10, const double C20, const double C30 )
    {

      const T J         = sqrt( determinant( C ) );
      const T I1        = trace( C );
      const T I1_minus3 = I1 * pow( J, -2. / 3. ) - 3.;
      T       res = C10 * I1_minus3 + C20 * I1_minus3 * I1_minus3 + C30 * I1_minus3 * I1_minus3 * I1_minus3;
      return res;
    }

    /** @brief Compressible neo-Hookean hyperelastic energy density function
     * (isochoric part only -- no volumetric term; see #VolumetricPenaltyPotential
     * for the shared volumetric term every potential in this file now uses),
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
    T NeoHookePotential( const Tensor33t< T >& C, const double mu )
    {
      const T J     = sqrt( determinant( C ) );
      const T Ibar1 = trace( C ) * pow( J, -2. / 3. );
      return mu / 2. * ( Ibar1 - 3. );
    }

    /** @brief Shared volumetric penalty potential, used by every hyperelastic
     * base in this file (in place of the several different, mutually
     * inconsistent volumetric conventions historically used by individual
     * materials -- see the class-level docs of #Marmot::Materials::BergstromBoyce
     * and #Marmot::Materials::CompressibleFiniteStrainLinearViscoelasticity for
     * why this specific form was adopted).
     *
     * \f[
     *   W_{\rm vol} = \frac{\kappa}{8}\left(\ln I_3\right)^2, \qquad I_3 = J^2 = \det\boldsymbol C
     * \f]
     *
     * @tparam T Scalar type, e.g. double, float, autodiff scalar.
     * @param C Right Cauchy-Green tensor.
     * @param kappa Bulk modulus.
     * @return Energy density (volumetric part only).
     */
    template < typename T >
    T VolumetricPenaltyPotential( const Tensor33t< T >& C, const double kappa )
    {
      const T lnDetC = log( determinant( C ) );
      return kappa / 8. * lnDetC * lnDetC;
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
    T standardNeoHooke( const Tensor33t< T >& C, const double K, const double G )
    {
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
    namespace detail {
      /** @brief Shared closed-form energy and derivative for the Arruda-Boyce 8-chain
       * potential's isochoric part, as a function of the isochoric first invariant
       * \f$\bar I_1\f$ alone -- used by both the plain (energy-only) and
       * FirstOrderDerived (energy + first derivative w.r.t. C) overloads of
       * #ArrudaBoyce8ChainPotential below, so the two cannot independently drift out
       * of sync with each other.
       *
       * Closed-form Cohen (1991) Pade approximation to the inverse Langevin function,
       * integrated in \f$\bar I_1\f$. Regularized near the locking limit
       * (\f$\bar I_1 \to 3\lambda_L^2\f$, i.e. the average chain stretch approaching
       * the locking stretch) by shifting the real part of the argument \f$w\f$ up to a
       * small positive floor rather than overwriting it, so that any complex-step or
       * dual-number perturbation carried in \f$\bar I_1\f$ survives the regularization.
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
        const T      w        = 1. - Ibar1 / ( 3. * lambdaL2 );
        const double w0       = 1. - 1. / lambdaL2;

        constexpr double wFloor   = 1e-2;
        const double     wReal    = Math::makeReal( w );
        const T          wClamped = wReal > wFloor ? w : w + T( wFloor - wReal );

        const T psi_iso     = mu / 6. * ( Ibar1 - 3. ) - mu * lambdaL2 * log( wClamped / w0 );
        const T dPsi_dIbar1 = mu / 6. + mu / 3. / wClamped;

        return { psi_iso, dPsi_dIbar1 };
      }
    } // namespace detail
    /** @endcond */

    /** @brief Arruda-Boyce 8-chain hyperelastic energy density function (isochoric
     * part only -- no volumetric term; each consuming material adds its own
     * volumetric penalty, since the classical 8-chain derivation is purely isochoric
     * and different materials in this codebase use different volumetric conventions),
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
    T ArrudaBoyce8ChainPotential( const Tensor33t< T >& C, const double mu, const double lambdaL )
    {
      const T J     = sqrt( determinant( C ) );
      const T I1    = trace( C );
      const T Ibar1 = I1 * pow( J, -2. / 3 );

      return detail::arrudaBoyce8ChainEnergyAndDerivative( Ibar1, mu, lambdaL ).first;
    }

    /** @brief Classical 3-term Ogden isochoric hyperelastic energy density function
     * (isochoric part only -- no volumetric term, same convention as
     * #ArrudaBoyce8ChainPotential: each consuming material adds its own volumetric
     * penalty).
     *
     * The isochoric energy density is
     * \f[
     *   \Psi_\mathrm{iso}(\bar\lambda) = \sum_{p=1}^{3} \frac{\mu_p}{\alpha_p}\left(
     *   \bar\lambda_1^{\alpha_p} + \bar\lambda_2^{\alpha_p} + \bar\lambda_3^{\alpha_p} - 3\right)
     * \f]
     * with isochoric principal stretches \f$\bar\lambda_i = J^{-1/3}\lambda_i\f$,
     * \f$\lambda_i = \sqrt{\mathrm{eig}_i(\boldsymbol C)}\f$ the principal stretches of
     * \f$\boldsymbol C\f$, and \f$J=\sqrt{\det\boldsymbol C}\f$. Unlike every other
     * potential in this file, a general (non-quadratic) Ogden exponent is not
     * expressible purely in terms of \f$\boldsymbol C\f$'s invariants, so this requires
     * a spectral decomposition of \f$\boldsymbol C\f$ via
     * Marmot::Math::computeEigenSystemJacobi() -- generic over scalar type (so it
     * differentiates correctly through autodiff::dual3rd the same way the invariant-based
     * potentials above do) and numerically safe at repeated eigenvalues (e.g.
     * \f$\boldsymbol C=\boldsymbol I\f$ at the reference configuration, where the Jacobi
     * sweep converges immediately since all off-diagonal entries already vanish). The
     * decomposition is performed once and its principal stretches reused for all three
     * terms. Reduces exactly to the isochoric part of #standardNeoHooke if a single term
     * has \f$\alpha_p=2\f$ and the other two have \f$\mu_p=0\f$, and is exactly
     * stress-free at \f$\boldsymbol C=\boldsymbol I\f$ for any \f$\alpha_p\f$, since every
     * \f$\bar\lambda_i=1\f$ there.
     *
     * @tparam T Scalar type, e.g. double, autodiff::dual3rd.
     * @param C Right Cauchy-Green tensor.
     * @param mu1 First-term Ogden modulus.
     * @param alpha1 First-term Ogden exponent.
     * @param mu2 Second-term Ogden modulus.
     * @param alpha2 Second-term Ogden exponent.
     * @param mu3 Third-term Ogden modulus.
     * @param alpha3 Third-term Ogden exponent.
     * @return Energy density (isochoric part only).
     */
    template < typename T >
    T OgdenPotential( const Tensor33t< T >& C,
                      const double         mu1,
                      const double         alpha1,
                      const double         mu2,
                      const double         alpha2,
                      const double         mu3,
                      const double         alpha3 )
    {
      const T   J    = sqrt( determinant( C ) );
      const T   Jm13 = pow( J, -1. / 3. );
      const auto [eigC, Q] = Marmot::Math::computeEigenSystemJacobi( C );

      Tensor3t< T > lambdaBar;
      for ( int i = 0; i < 3; ++i )
        lambdaBar( i ) = Jm13 * sqrt( eigC( i ) );

      const double mu[3]    = { mu1, mu2, mu3 };
      const double alpha[3] = { alpha1, alpha2, alpha3 };

      T psi( 0. );
      for ( int p = 0; p < 3; ++p ) {
        T sumPow( 0. );
        for ( int i = 0; i < 3; ++i )
          sumPow += pow( lambdaBar( i ), alpha[p] );
        psi += mu[p] / alpha[p] * ( sumPow - 3.0 );
      }
      return psi;
    }

    namespace FirstOrderDerived {

      /** @brief Shared volumetric penalty potential (see the plain
       * #VolumetricPenaltyPotential) and its first derivative w.r.t. C,
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
      std::tuple< T, Tensor33t< T > > VolumetricPenaltyPotential( const Tensor33t< T >& C, const double kappa )
      {
        const T              lnDetC = log( determinant( C ) );
        const Tensor33t< T > CInv   = inverse( C );
        const T              psi    = kappa / 8. * lnDetC * lnDetC;
        const Tensor33t< T > dPsi_dC = multiplyFastorTensorWithScalar( CInv, T( kappa / 4. * lnDetC ) );
        return { psi, dPsi_dC };
      }

      /** @brief Compressible neo-Hookean potential (see the plain
       * #NeoHookePotential, isochoric part only) and its first derivative w.r.t. C,
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
      std::tuple< T, Tensor33t< T > > NeoHookePotential( const Tensor33t< T >& C, const double mu )
      {
        const T Jm23  = pow( determinant( C ), -1. / 3. );
        const T Ibar1 = trace( C ) * Jm23;

        const Tensor33t< T > I           = fastorTensorFromDoubleTensor< T >( Spatial3D::I );
        const Tensor33t< T > CInv         = inverse( C );
        const Tensor33t< T > dIbar1_dC    = multiplyFastorTensorWithScalar( I, Jm23 ) -
                                          multiplyFastorTensorWithScalar( CInv, Ibar1 / 3. );

        const T              psi     = mu / 2. * ( Ibar1 - 3. );
        const Tensor33t< T > dPsi_dC = multiplyFastorTensorWithScalar( dIbar1_dC, T( mu / 2. ) );
        return { psi, dPsi_dC };
      }

      /** @brief Yeoh potential (see the plain #YeohPotential, isochoric part
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
      std::tuple< T, Tensor33t< T > > YeohPotential( const Tensor33t< T >& C,
                                                      const double          C10,
                                                      const double          C20,
                                                      const double          C30 )
      {
        const T Jm23      = pow( determinant( C ), -1. / 3. );
        const T Ibar1     = trace( C ) * Jm23;
        const T Ibar1m3   = Ibar1 - 3.;

        const Tensor33t< T > I        = fastorTensorFromDoubleTensor< T >( Spatial3D::I );
        const Tensor33t< T > CInv     = inverse( C );
        const Tensor33t< T > dIbar1_dC = multiplyFastorTensorWithScalar( I, Jm23 ) -
                                          multiplyFastorTensorWithScalar( CInv, Ibar1 / 3. );

        const T psi          = C10 * Ibar1m3 + C20 * Ibar1m3 * Ibar1m3 + C30 * Ibar1m3 * Ibar1m3 * Ibar1m3;
        const T dPsi_dIbar1  = C10 + 2. * C20 * Ibar1m3 + 3. * C30 * Ibar1m3 * Ibar1m3;
        const Tensor33t< T > dPsi_dC = multiplyFastorTensorWithScalar( dIbar1_dC, dPsi_dIbar1 );
        return { psi, dPsi_dC };
      }

      /** @brief Mooney-Rivlin potential (see the plain #MooneyRivlinPotential,
       * isochoric part only) and its first derivative w.r.t. C.
       *
       * @tparam T Scalar type, e.g. double, std::complex<double>.
       * @param C Right Cauchy-Green tensor.
       * @param C10 Mooney-Rivlin material parameter.
       * @param C01 Mooney-Rivlin material parameter.
       * @return A tuple containing the (isochoric) energy density and its first derivative w.r.t. C.
       */
      template < typename T >
      std::tuple< T, Tensor33t< T > > MooneyRivlinPotential( const Tensor33t< T >& C,
                                                              const double          C10,
                                                              const double          C01 )
      {
        const T J     = sqrt( determinant( C ) );
        const T Jm23  = pow( J, -2. / 3. );
        const T Jm43  = pow( J, -4. / 3. );
        const T I1    = trace( C );
        const T Ibar1 = I1 * Jm23;
        const T I2    = 0.5 * ( I1 * I1 - trace( C % C ) );
        const T Ibar2 = I2 * Jm43;

        const Tensor33t< T > I    = fastorTensorFromDoubleTensor< T >( Spatial3D::I );
        const Tensor33t< T > CInv = inverse( C );

        const Tensor33t< T > dIbar1_dC = multiplyFastorTensorWithScalar( I, Jm23 ) -
                                          multiplyFastorTensorWithScalar( CInv, Ibar1 / 3. );
        // dIbar2/dC = J^(-4/3)*I1*I - J^(-4/3)*C - (2/3)*Ibar2*Cinv (standard
        // isochoric second-invariant derivative identity, I2's own derivative
        // dI2/dC = I1*I - C combined with the J^(-4/3) scaling's own C-derivative).
        const Tensor33t< T > dIbar2_dC = multiplyFastorTensorWithScalar( I, T( Jm43 * I1 ) ) -
                                          multiplyFastorTensorWithScalar( C, Jm43 ) -
                                          multiplyFastorTensorWithScalar( CInv, T( 2. / 3. * Ibar2 ) );

        const T              psi     = C10 * ( Ibar1 - 3. ) + C01 * ( Ibar2 - 3. );
        const Tensor33t< T > dPsi_dC = multiplyFastorTensorWithScalar( dIbar1_dC, T( C10 ) ) +
                                       multiplyFastorTensorWithScalar( dIbar2_dC, T( C01 ) );
        return { psi, dPsi_dC };
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
      std::tuple< T, Tensor33t< T > > PenceGouPotentialB( const Tensor33t< T >& C, const double K, const double G )
      {
        using namespace FastorIndices;

        const T J  = sqrt( determinant( C ) );
        const T I1 = trace( C );
        // energy density
        T psi = K / 8. * pow( J - 1. / J, 2. ) + G / 2. * ( I1 * pow( J, -2. / 3 ) - 3. );

        // first derivative w.r.t. C
        const T dPsi_dJ  = K / 4. * ( J - 1. / J ) * ( 1. + 1. / ( J * J ) ) - G / 3. * I1 * pow( J, -5. / 3. );
        const T dPsi_dI1 = G / 2. * pow( J, -2. / 3. );

        const Tensor33t< T > CInv   = inverse( C );
        const Tensor33t< T > dJ_dC  = multiplyFastorTensorWithScalar( transpose( CInv ), T( J / 2. ) );
        const Tensor33t< T > dI1_dC = fastorTensorFromDoubleTensor< T >( Spatial3D::I );

        Tensor33t< T > dPsi_dC = multiplyFastorTensorWithScalar( dJ_dC, dPsi_dJ ) +
                                 multiplyFastorTensorWithScalar( dI1_dC, dPsi_dI1 );

        return { psi, dPsi_dC };
      }

      /** @brief Arruda-Boyce 8-chain hyperelastic energy density function (isochoric
       * part only) and its first derivative w.r.t. C -- see the plain
       * #ArrudaBoyce8ChainPotential for the formula and its reduction/stress-free
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
      std::tuple< T, Tensor33t< T > > ArrudaBoyce8ChainPotential( const Tensor33t< T >& C,
                                                                    const double          mu,
                                                                    const double          lambdaL )
      {
        const T J     = sqrt( determinant( C ) );
        const T I1    = trace( C );
        const T Ibar1 = I1 * pow( J, -2. / 3 );

        const auto [ psi_iso,
                     dPsi_dIbar1 ] = EnergyDensityFunctions::detail::arrudaBoyce8ChainEnergyAndDerivative( Ibar1,
                                                                                                            mu,
                                                                                                            lambdaL );

        const Tensor33t< T > I      = fastorTensorFromDoubleTensor< T >( Spatial3D::I );
        const Tensor33t< T > CInv   = inverse( C );
        // dIbar1/dC = d(I1*J^(-2/3))/dC = J^(-2/3)*dI1/dC + I1*d(J^(-2/3))/dC
        //           = J^(-2/3)*I - (I1/3)*J^(-2/3)*Cinv
        // and (I1/3)*J^(-2/3) = Ibar1/3 exactly (since Ibar1 = I1*J^(-2/3) by
        // definition) -- so J^(-2/3) multiplies ONLY the I term, not the Cinv
        // term (an earlier version of this code mistakenly applied J^(-2/3)
        // to both, which is wrong away from J=1 and was caught by
        // TestBergstromBoyce's P-4/I-1/I-2/I-3 failing against the isochoric
        // formulation's own reduction check).
        const Tensor33t< T > dIbar1_dC = multiplyFastorTensorWithScalar( I, T( pow( J, -2. / 3 ) ) ) -
                                          multiplyFastorTensorWithScalar( CInv, Ibar1 / 3. );

        const Tensor33t< T > dPsi_dC = multiplyFastorTensorWithScalar( dIbar1_dC, dPsi_dIbar1 );

        return { psi_iso, dPsi_dC };
      }
    } // namespace FirstOrderDerived

    namespace SecondOrderDerived {

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
      std::tuple< T, Tensor33t< T >, Tensor3333t< T > > PenceGouPotentialB( const Tensor33t< T >& C,
                                                                            const double          K,
                                                                            const double          G )
      {
        using namespace FastorIndices;

        const T J  = sqrt( determinant( C ) );
        const T I1 = trace( C );
        // energy density
        T psi = K / 8. * pow( J - 1. / J, 2. ) + G / 2. * ( I1 * pow( J, -2. / 3 ) - 3. );

        // first derivative w.r.t. C
        const T dPsi_dJ  = K / 4. * ( J - 1. / J ) * ( 1 + 1. / ( J * J ) ) - G / 3. * I1 * pow( J, -5. / 3. );
        const T dPsi_dI1 = G / 2. * pow( J, -2. / 3. );

        const Tensor33t< T > CInv   = inverse( C );
        const Tensor33t< T > dJ_dC  = 0.5 * J * transpose( CInv );
        const Tensor33t< T > dI1_dC = Spatial3D::I;

        Tensor33t< T > dPsi_dC = dPsi_dJ * dJ_dC + dPsi_dI1 * dI1_dC;

        // second derivative w.r.t. C
        const T          d2Psi_dJdJ  = K / 4. * ( 1. + 3. / ( J * J * J * J ) ) + 5. / 9. * G * I1 * pow( J, -8. / 3. );
        const T          d2Psi_dJdI1 = -G / 3. * pow( J, -5. / 3. );
        Tensor3333t< T > d2J_dCdC    = J / 4. * einsum< JI, LK, to_IJKL >( CInv, CInv ) -
                                    J / 2. * einsum< JK, LI, to_IJKL >( CInv, CInv );

        Tensor3333t< T > d2Psi_dCdC = d2Psi_dJdJ * einsum< IJ, KL >( dJ_dC, dJ_dC ) + dPsi_dJ * d2J_dCdC +
                                      d2Psi_dJdI1 *
                                        ( einsum< IJ, KL >( dJ_dC, dI1_dC ) + einsum< IJ, KL >( dI1_dC, dJ_dC ) );

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
      std::tuple< T, FastorStandardTensors::Tensor33t< T >, FastorStandardTensors::Tensor3333t< T > > standardNeoHooke(
        const FastorStandardTensors::Tensor33t< T >& C,
        const double&                                K,
        const double&                                G )
      {
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
        const Tensor33t< T >& dTrC_dC    = makeOtherScalarType< T >( Spatial3D::I );
        const Tensor33t< T >  invCt      = transpose( inverse( C ) );
        const Tensor33t< T >  dDetC_dC   = multiplyFastorTensorWithScalar( invCt, detC );
        const Tensor33t< T >& dLnDetC_dC = invCt;

        // first derivative with respect to C
        const Tensor33t< T > dPsi_dC = G / 2 * ( dTrC_dC - dLnDetC_dC ) + lambda / 4 * ( dDetC_dC - dLnDetC_dC );

        // second derivative quantities
        using namespace FastorIndices;
        using lj                          = Index< l_, j_ >;
        using li                          = Index< l_, i_ >;
        const Tensor3333t< T > dInvCt_dC  = -0.5 * ( einsum< ik, lj, to_ijkl >( invCt, invCt ) +
                                                    einsum< jk, li, to_ijkl >( invCt, invCt ) );
        const Tensor3333t< T > d2DetC_dC2 = multiplyFastorTensorWithScalar( dInvCt_dC, detC ) +
                                            einsum< ij, kl, to_ijkl >( invCt, dDetC_dC );
        const Tensor3333t< T >& d2LnDetC_dC2 = dInvCt_dC;

        // second derivative with respect to C
        const Tensor3333t< T > d2Psi_dC2 = lambda / 4.0 * d2DetC_dC2 - ( G / 2.0 + lambda / 4.0 ) * d2LnDetC_dC2;
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
      std::tuple< T, FastorStandardTensors::Tensor33t< T >, FastorStandardTensors::Tensor3333t< T > > BiotNeoHooke(
        const FastorStandardTensors::Tensor33t< T >& U,
        const double&                                K,
        const double&                                G )
      {
        using namespace FastorIndices;
        Tensor3333t< T > I4       = makeOtherScalarType< T >( Spatial3D::I4 );
        Tensor33t< T >   C        = U % U;
        Tensor3333t< T > dC_dU    = 2 * einsum< iL, ijkl >( U, I4 );
        Tensor3333t< T > d2C_dUdU = 2 * I4;

        auto [psi, dPsi_dC, d2Psi_dCdC] = standardNeoHooke< T >( C, K, G );

        Tensor33t< T >   dPsi_dU    = einsum< kl, klmn >( dPsi_dC, dC_dU );
        Tensor3333t< T > d2Psi_dUdU = einsum< ijkl, klmn >( einsum< ijkl, ijmn >( dC_dU, d2Psi_dCdC ), dC_dU ) +
                                      einsum< jL, ijkl >( dPsi_dC, d2C_dUdU );

        return { psi, dPsi_dU, d2Psi_dUdU };
      }

    } // namespace SecondOrderDerived

    namespace ThirdOrderDerived {
      /** @brief Standard compressible Neo-Hooke energy density and first/second/third derivatives w.r.t.
       * \f$\boldsymbol{C}\f$.
       *
       * @param C Right Cauchy-Green tensor.
       * @param K Bulk modulus.
       * @param G Shear modulus.
       * @return Tuple of energy density, first derivative, second derivative and third derivative.
       */
      std::tuple< double, FastorStandardTensors::Tensor33d, FastorStandardTensors::Tensor3333d, FastorStandardTensors::Tensor333333d > standardNeoHooke(
        const FastorStandardTensors::Tensor33d& C,
        const double&                           K,
        const double&                           G );
    } // namespace ThirdOrderDerived
  }   // namespace EnergyDensityFunctions

} // namespace Marmot::ContinuumMechanics
