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
 * Matthias Neuner matthias.neuner@uibk.ac.at
 * Magdalena Schreter magdalena.schreter@uibk.ac.at
 * Alexander Dummer alexander.dummer@uibk.ac.at
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
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTypedefs.h"

#define VOIGTFROMDIM( x ) ( ( ( x * x ) + x ) >> 1 )

/**
 * \brief This file includes functions needed for calculations with stress and strain tensors written in voigt notation.
 */
namespace Marmot {
  namespace ContinuumMechanics::VoigtNotation {

    /* constexpr int VoigtSize = 6; */

    enum VoigtSize { OneD = 1, TwoD = 3, ThreeD = 6, Axial = 4 };

    constexpr VoigtSize voigtSizeFromDimension( int x )
    {
      return (VoigtSize)( ( ( x * x ) + x ) >> 1 );
    }

    extern const Marmot::Vector6d P;
    extern const Marmot::Vector6d PInv;

    extern const Marmot::Vector6d I;
    extern const Marmot::Vector6d IHyd;
    extern const Matrix6d         IDev;

    // Plane Stress handling

    /**
     * Converts a 3D voigt notated vector to a plane stress vector.
     *\f[
      \displaystyle  \begin{bmatrix}
        \sigma_{11}\\
        \sigma_{22}\\
        \sigma_{33} = 0\\
        \sigma_{12}\\
        \sigma_{13} = 0\\
        \sigma_{23} = 0
           \end{bmatrix} \hspace{.5cm} \Rightarrow \hspace{.5cm} \begin{bmatrix}
                    \sigma_{11}\\
                    \sigma_{22}\\
                    \sigma_{12}
                             \end{bmatrix}
      \f]
     */
    Eigen::Vector3d voigtToPlaneVoigt( const Marmot::Vector6d& voigt );
    Vector4d        voigtToAxisymmetricVoigt( const Vector6d& voigt );

    /**
     * Converts a voigt notated plane stress vector to a 3D vector.
     *\f[
      \displaystyle  \begin{bmatrix}
        \sigma_{11}\\
        \sigma_{22}\\
        \sigma_{12}
           \end{bmatrix} \hspace{.5cm} \Rightarrow \hspace{.5cm} \begin{bmatrix}
                    \sigma_{11}\\
                    \sigma_{22}\\
                    0\\
                    \sigma_{12}\\
                    0\\
                    0
                             \end{bmatrix}
      \f]
     */
    Marmot::Vector6d planeVoigtToVoigt( const Eigen::Vector3d& voigtPlane );
    Vector6d         axisymmetricVoigtToVoigt( const Vector4d& voigtAxisymmetric );

    /**
     * Reduces a 3D voigt notated vector to a lower dimension defined by the template parameter 'voigtSize'.
     *
     * Considered cases:
     * 	- voigtSize \f$ = 1 \f$:
     *\f[
      \displaystyle  \begin{bmatrix}
        \sigma_{11}\\
        \sigma_{22} = 0\\
        \sigma_{33} = 0\\
        \sigma_{12} = 0\\
        \sigma_{13} = 0\\
        \sigma_{23} = 0
           \end{bmatrix} \hspace{.5cm} \Rightarrow \hspace{.5cm} \sigma_{11}
      \f]
     * 	- voigtSize \f$ = 3 \f$: Calls function voigtToPlaneVoigt().
     * 	- voigtSize \f$ = 6 \f$: Returns the input vector
     */

    template < enum VoigtSize voigtSize >
    Eigen::Matrix< double, voigtSize, 1 > reduce3DVoigt( const Marmot::Vector6d& Voigt3D )
    {
      if constexpr ( voigtSize == OneD )
        return ( Eigen::Matrix< double, 1, 1 >() << Voigt3D( 0 ) ).finished();
      else if constexpr ( voigtSize == TwoD )
        return voigtToPlaneVoigt( Voigt3D );
      else if constexpr ( voigtSize == ThreeD )
        return Voigt3D;
      else
        throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid dimension specified" );
    }

    /**
     * Computes a 3D voigt notated vector from a vector of lower dimension defined by the template parameter
     'voigtSize'.
     *
     * Considered cases:
     * 	- voigtSize \f$ = 1 \f$:
     *\f[
      \displaystyle	\sigma_{11} \hspace{.5cm} \Rightarrow \hspace{.5cm} \begin{bmatrix}
                    \sigma_{11}\\
                    0\\
                    0\\
                    0\\
                    0\\
                    0
                            \end{bmatrix}
      \f]
     * 	- voigtSize \f$ = 3 \f$: Calls function planeVoigtToVoigt().
     * 	- voigtSize \f$ = 6 \f$: Returns the input vector
     */

    template < enum VoigtSize voigtSize >
    Marmot::Vector6d make3DVoigt( const Eigen::Matrix< double, voigtSize, 1 >& Voigt )
    {
      if constexpr ( voigtSize == OneD )
        return ( Marmot::Vector6d() << Voigt( 0 ), 0, 0, 0, 0, 0 ).finished();
      else if constexpr ( voigtSize == TwoD )
        return planeVoigtToVoigt( Voigt );
      else if constexpr ( voigtSize == ThreeD )
        return Voigt;
      else if constexpr ( voigtSize == Axial )
        return axisymmetricVoigtToVoigt( Voigt );
      else
        throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid dimension specified" );
    }

    // function prototypes for  Marmot::Vector6d handling

    /**
     * Converts a voigt notated strain vector to its corresponding strain tensor
     *\f[
      \displaystyle  \begin{bmatrix}
        \varepsilon_{11}\\
        \varepsilon_{22}\\
        \varepsilon_{33}\\
        \gamma_{12}\\
        \gamma_{13}\\
        \gamma_{23}
           \end{bmatrix} \hspace{.5cm} \Rightarrow \hspace{.5cm} \begin{bmatrix}
           \varepsilon_{11} & \gamma_{12}/2 & \gamma_{13}/2 \\
           \gamma_{12}/2 & \varepsilon_{22} & \gamma_{23}/2 \\
           \gamma_{13}/2 & \gamma_{23}/2 & \varepsilon_{33} \\
     \end{bmatrix}
      \f]
     */
    template < typename T >
    Eigen::Matrix< T, 3, 3 > voigtToStrain( const Eigen::Matrix< T, 6, 1 >& voigt )
    {
      Eigen::Matrix< T, 3, 3 > strain;
      // clang-format off
      strain << voigt[0],       voigt[3] * 0.5, voigt[4] * 0.5, 
                voigt[3] * 0.5, voigt[1],       voigt[5] * 0.5, 
                voigt[4] * 0.5, voigt[5] * 0.5, voigt[2];
      // clang-format on
      return strain;
    }

    /**
     * Converts a voigt notated stress vector to its corresponding stress tensor
     *\f[
      \displaystyle  \begin{bmatrix}
        \sigma_{11}\\
        \sigma_{22}\\
        \sigma_{33}\\
        \sigma_{12}\\
        \sigma_{13}\\
        \sigma_{23}
           \end{bmatrix} \hspace{.5cm} \Rightarrow \hspace{.5cm} \begin{bmatrix}
           \sigma_{11} & \sigma_{12} & \sigma_{13} \\
           \sigma_{12} & \sigma_{22} & \sigma_{23} \\
           \sigma_{13} & \sigma_{23} & \sigma_{33} \\
     \end{bmatrix}
      \f]
     */
    template < typename T >
    Eigen::Matrix< T, 3, 3 > voigtToStress( const Eigen::Matrix< T, 6, 1 >& voigt )
    {
      Eigen::Matrix< T, 3, 3 > stress;
      // clang-format off
      stress << voigt[0], voigt[3], voigt[4], 
                voigt[3], voigt[1], voigt[5], 
                voigt[4], voigt[5], voigt[2];
      // clang-format on
      return stress;
    }

    /**
     * Converts a strain tensor to its corresponding voigt notated strain vector
     *\f[
      \displaystyle  \begin{bmatrix}
              \varepsilon_{11} & \varepsilon_{12} & \varepsilon_{13}\\
              \varepsilon_{21} & \varepsilon_{22} & \varepsilon_{23} \\
              \varepsilon_{31} & \varepsilon_{32} & \varepsilon_{33} \\
           \end{bmatrix} \hspace{.5cm} \Rightarrow \hspace{.5cm} \begin{bmatrix}
                    \varepsilon_{11}\\
                    \varepsilon_{22}\\
                    \varepsilon_{33}\\
                    2\,\varepsilon_{12}\\
                    2\,\varepsilon_{13}\\
                    2\,\varepsilon_{23}
                             \end{bmatrix}
      \f]
     */
    Marmot::Vector6d strainToVoigt( const Eigen::Matrix3d& strainTensor );

    /**
     * Converts a stress tensor to its corresponding voigt notated stress vector
     *\f[
      \displaystyle  \begin{bmatrix}
              \sigma_{11} & \sigma_{12} & \sigma_{13}\\
              \sigma_{21} & \sigma_{22} & \sigma_{23} \\
              \sigma_{31} & \sigma_{32} & \sigma_{33} \\
           \end{bmatrix} \hspace{.5cm} \Rightarrow \hspace{.5cm} \begin{bmatrix}
                    \sigma_{11}\\
                    \sigma_{22}\\
                    \sigma_{33}\\
                    \sigma_{12}\\
                    \sigma_{13}\\
                    \sigma_{23}
                             \end{bmatrix}
      \f]
     */
    template < typename T = double >
    Eigen::Matrix< T, 6, 1 > stressToVoigt( const Eigen::Matrix< T, 3, 3 >& stressTensor )
    {
      Eigen::Matrix< T, 6, 1 > stress;
      // clang-format off
            stress << stressTensor( 0, 0 ), 
                      stressTensor( 1, 1 ), 
                      stressTensor( 2, 2 ), 
                      stressTensor( 0, 1 ),
                      stressTensor( 0, 2 ), 
                      stressTensor( 1, 2 );
      // clang-format on
      return stress;
    }

    Eigen::Matrix< double, 6, 6 > stiffnessToVoigt( const Eigen::Tensor< double, 4 >& C );

    Eigen::Tensor< double, 4 > voigtToStiffness( const Eigen::Matrix< double, 6, 6 >& voigtStiffness );

    template < int nDim >
    Eigen::Matrix< double, nDim, nDim > stressMatrixFromVoigt(
      const Eigen::Matrix< double, VOIGTFROMDIM( nDim ), 1 >& Voigt )
    {
      if constexpr ( nDim == 1 )
        return ( Eigen::Matrix< double, nDim, nDim >() << Voigt( 0 ) ).finished();
      else if constexpr ( nDim == 2 )
        return ( Eigen::Matrix< double, nDim, nDim >() << Voigt( 0 ), Voigt( 2 ), Voigt( 2 ), Voigt( 1 ) ).finished();
      else if constexpr ( nDim == 3 )
        return voigtToStress( Voigt );
      else
        throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid dimension specified" );
    }

    template < int nDim >
    Eigen::Matrix< double, VOIGTFROMDIM( nDim ), 1 > voigtFromStrainMatrix(
      const Eigen::Matrix< double, nDim, nDim >& strain )
    {
      if constexpr ( nDim == 1 )
        return ( Eigen::Matrix< double, VOIGTFROMDIM( nDim ), 1 >() << strain( 0, 0 ) ).finished();
      else if constexpr ( nDim == 2 )
        return ( Eigen::Matrix< double, VOIGTFROMDIM( nDim ), 1 >() << strain( 0, 0 ),
                 strain( 1, 1 ),
                 2 * strain( 0, 1 ) )
          .finished();
      else if constexpr ( nDim == 3 )
        return strainToVoigt( strain );
      else
        throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid dimension specified" );
    }

    namespace Invariants {

      /** Computes the principal strains by solving the eigenvalue problem.
       *\f[
           \displaystyle |\varepsilon_{ij} - \lambda\, \delta_{ij}| = 0 \hspace{.5cm} \Rightarrow \hspace{.5cm}
     \lambda^{(1)},\lambda^{(2)},\lambda^{(3)}\hspace{0.3cm} \widehat{=}\hspace{0.3cm} \varepsilon_1,\, \varepsilon_2,\,
     \varepsilon_3 \f]
       * The resulting principal strains are NOT sorted.
       */
      Eigen::Vector3d principalStrains( const Marmot::Vector6d& strain );

      /** Computes the principal stresses by solving the eigenvalue problem.
       *\f[
           \displaystyle |\sigma_{ij} - \lambda\, \delta_{ij}| = 0 \hspace{.5cm} \Rightarrow \hspace{.5cm}
     \lambda^{(1)},\lambda^{(2)},\lambda^{(3)}\hspace{0.3cm} \widehat{=}\hspace{0.3cm}\sigma_1,\, \sigma_2,\, \sigma_3
        \f]
       * The resulting principal stresses are NOT sorted.
       */
      Eigen::Vector3d principalStresses( const Marmot::Vector6d& stress );
      // principal strains calculated from haigh westergaard strains ( sorted --> e1 > e2 > e3 )

      /** Calculates the principal strains from its corresponding haigh westergaard coordinates.
       *\f[
          \displaystyle \begin{bmatrix}
          \varepsilon_{1}\\
        \varepsilon_{2}\\
        \varepsilon_{3}
      \end{bmatrix} = \frac{1}{\sqrt{3}} \begin{bmatrix}
                  \xi\\
                \xi\\
                \xi
                 \end{bmatrix} + \sqrt{\frac{2}{3}}\,\rho\,\begin{bmatrix}
                                \cos(\theta)\\
                          -\sin\left(\frac{\pi}{6} - \theta\right)\\
                          -\sin\left(\frac{\pi}{6} + \theta\right)
                             \end{bmatrix}
        \f]
       *The computation of \f$ \xi\f$ ,\ \f$ \rho\f$ and \f$\theta\f$ can be found in haighWestergaardFromStrain()
       */
      Eigen::Vector3d sortedPrincipalStrains( const Marmot::Vector6d& strain );
      // principal stressDirections calculated by solving eigenvalue problem ( !NOT sorted! )

      /** Computes the principal stress directions \f$\boldsymbol{x}^{(k)}\f$ of the eigenvalues \f$ \sigma_k \f$ by
       solving
       *\f[
           \displaystyle \left(\boldsymbol{\sigma} - \sigma_k \cdot \boldsymbol{I}\right) \cdot \boldsymbol{x}^{(k)}  =
       0 \f]
       * The resulting principal stress directions are NOT sorted.
       */
      Eigen::Matrix3d principalStressesDirections( const Marmot::Vector6d& stress );

      /** Computes the equivalent von Mises stress.
       *\f[
           \displaystyle \sigma^{(eq)} = \sqrt{3 \cdot J_2}
        \f]
       * Wherein \f$ J_2 \f$ denotes the second invariant of the deviator stress tensor (see J2()).
       */
      double vonMisesEquivalentStress( const Marmot::Vector6d& stress );

      /** Computes the equivalent von Mises strain from deviatoric part of the strain tensor \f$ e_{ij} \f$
       *\f[
           \displaystyle \varepsilon^{(eq)} = \sqrt{ \frac{2}{3} \cdot e_{ij}\,e_{ij}}
        \f]
       */
      double vonMisesEquivalentStrain( const Marmot::Vector6d& strain );

      /** Computes the euclidian norm of the strain tensor \f$ ||\boldsymbol{\varepsilon}|| \f$
       */
      template < typename T >
      T normStrain( const Eigen::Matrix< T, 6, 1 >& strain )
      {
        return ContinuumMechanics::VoigtNotation::voigtToStrain( strain ).norm();
      }

      /** Computes the euclidian norm of the stress tensor \f$ ||\boldsymbol{\sigma}|| \f$
       */
      double normStress( const Marmot::Vector6d& stress );
      // Trace of compressive strains

      /** Computes the volumetric plastic strains in compression
       *\f[
           \displaystyle \varepsilon^{vol}_{\ominus} = \sum^{3}_{i = 1} \left\langle -\varepsilon_i \right\rangle
        \f]
       * using the Macaulay brackets \f$ \left\langle \bullet \right\rangle \f$ and the principal values of the strain
       tensor \f$\varepsilon_i \f$
       */
      double StrainVolumetricNegative( const Marmot::Vector6d& strain );

      /** Computes the first invariant \f$ I_1 \f$ of the stress tensor \f$ \boldsymbol{\sigma} \f$.
       *\f[
           \displaystyle I_1 = tr(\boldsymbol{\sigma})
        \f]
       */
      template < typename T >
      T I1( const Eigen::Matrix< T, 6, 1 >& stress )
      {
        return stress( 0 ) + stress( 1 ) + stress( 2 );
      }

      /** Computes the first invariant \f$ I^{(\varepsilon)}_1 \f$ of the strain tensor \f$ \boldsymbol{\varepsilon}
       \f$.
       *\f[
           \displaystyle I^{(\varepsilon)}_1 = tr(\boldsymbol{\varepsilon})
        \f]
       */
      double I1Strain( const Marmot::Vector6d& strain ); //

      /** Computes the second invariant \f$ I_2 \f$ of the stress tensor \f$ \boldsymbol{\sigma} \f$.
       *\f[
           \displaystyle I_2 = \frac{1}{2} \left(tr(\boldsymbol{\sigma})^2 - tr(\boldsymbol{\sigma}^2)\right)
        \f]
       */
      template < typename T >
      T I2( const Eigen::Matrix< T, 6, 1 >& stress )
      {
        const Eigen::Matrix< T, 6, 1 >& s = stress;
        return s( 0 ) * s( 1 ) + s( 1 ) * s( 2 ) + s( 2 ) * s( 0 ) - s( 3 ) * s( 3 ) - s( 4 ) * s( 4 ) -
               s( 5 ) * s( 5 );
      }

      /** Computes the second invariant \f$ I^{(\varepsilon)}_2 \f$ from a voigt notated strain vector \f$
  \boldsymbol{\varepsilon} \f$.
       *\f[
           \displaystyle I^{(\varepsilon)}_2 = \varepsilon_{11}\,\varepsilon_{22} + \varepsilon_{22}\,\varepsilon_{33} +
  \varepsilon_{11}\,\varepsilon_{33} - \frac{1}{4}(\gamma^2_{12}  - \gamma^2_{13}  - \gamma^2_{23}) \f]
       */
      double I2Strain( const Marmot::Vector6d& strain );

      /** Computes the third invariant \f$ I_3 \f$ of the stress tensor \f$ \boldsymbol{\sigma} \f$.
       *\f[
           \displaystyle I_3 = det(\boldsymbol{\sigma})
        \f]
       */
      template < typename T >
      T I3( const Eigen::Matrix< T, 6, 1 >& stress )
      {
        const Eigen::Matrix< T, 6, 1 >& s = stress;
        return s( 0 ) * s( 1 ) * s( 2 ) + 2. * s( 3 ) * s( 4 ) * s( 5 ) - s( 0 ) * s( 5 ) * s( 5 ) -
               s( 1 ) * s( 4 ) * s( 4 ) - s( 2 ) * s( 3 ) * s( 3 );
      }

      /** Computes the third invariant \f$ I^{(\varepsilon)}_3 \f$ from a voigt notated strain vector \f$
       * \boldsymbol{\varepsilon} \f$ by calling voigtToStrain() and calculating the determinant.
       */
      double I3Strain( const Marmot::Vector6d& strain );

      /** Computes the second invariant \f$ J_2 \f$ of the deviatoric part of the stress tensor \f$ \boldsymbol{s} \f$.
       *\f[
           \displaystyle J_2 = \frac{1}{3} I^2_1 - I_2
        \f]
       */
      template < typename T >
      T J2( const Eigen::Matrix< T, 6, 1 >& stress )
      {
        const T I1_ = I1( stress );
        const T I2_ = I2( stress );
        const T res = ( 1. / 3 ) * I1_ * I1_ - I2_;
        return Marmot::Math::makeReal( res ) >= 0 ? res : T( 0.0 );
      }

      /** Computes the second invariant \f$ J^{(\varepsilon)}_2 \f$ of the deviatoric part of the strain tensor \f$
       \boldsymbol{\varepsilon} \f$.
       *\f[
           \displaystyle J^{(\varepsilon)}_2 = \frac{1}{3} I^{(\varepsilon) 2}_1 - I^{(\varepsilon)}_2
        \f]
       */

      double J2Strain( const Marmot::Vector6d& strain );

      /** Computes the third invariant \f$ J_3 \f$ of the deviatoric part of the stress tensor \f$ \boldsymbol{s} \f$.
       *\f[
           \displaystyle J_3 = \frac{2}{27} I^3_1 - \frac{1}{3} I_1\cdot I_2\cdot I_3
        \f]
       */
      template < typename T >
      T J3( const Eigen::Matrix< T, 6, 1 >& stress )
      {
        const T I1_ = I1( stress );
        const T I2_ = I2( stress );
        const T I3_ = I3( stress );

        return ( 2. / 27 ) * pow( I1_, 3 ) - ( 1. / 3 ) * I1_ * I2_ + I3_;
      }

      /** Computes the third invariant \f$ J^{(\varepsilon)}_3 \f$ of a voigt notated deviatoric strain vector \f$
       * \boldsymbol{e} \f$ by calling voigtToStrain() and calculating the determinant.
       */
      double J3Strain( const Marmot::Vector6d& strain );

      // principal values in voigt
      std::pair< Eigen::Vector3d, Eigen::Matrix< double, 3, 6 > > principalValuesAndDerivatives(
        const Eigen::Matrix< double, 6, 1 >& S );

    } // namespace Invariants

    namespace Derivatives {

      // derivatives of Haigh Westergaard stresses with respect to cauchy stress in eng. notation

      /**
       * Computes the derivative \f$ \frac{d\, \sigma_m}{d\, \boldsymbol{\sigma}} \f$ of the mean stress \f$ \sigma_m
       * \f$ with respect to the voigt notated stress vector \f$ \boldsymbol{\sigma} \f$
       */
      Vector6d dStressMean_dStress();
      /**
       * Computes the derivative \f$ \frac{d\, \rho}{d\, \boldsymbol{\sigma}}\f$ of the haigh westergaard coordinate \f$
       * \rho \f$ with respect to the voigt notated stress vector \f$ \boldsymbol{\sigma} \f$
       */
      template < typename T >
      Eigen::Matrix< T, 6, 1 > dRho_dStress( T rho, const Eigen::Matrix< T, 6, 1 >& stress )
      {
        if ( Marmot::Math::makeReal( rho ) <= 1e-16 )
          return Eigen::Matrix< T, 6, 1 >::Zero();

        Eigen::Matrix< T, 6, 1 > s = IDev * stress;
        // P array results from the derivative of the double contraction s:s in voigt notation
        return T( 1. / rho ) * P.array() * s.array();
      }
      /**
       * Computes the derivative \f$ \frac{d\, \varepsilon_\rho}{d\, \boldsymbol{\varepsilon}}\f$ of the haigh
       * westergaard coordinate \f$
       * \strain_\rho \f$ with respect to the voigt notated strain vector \f$ \boldsymbol{\varepsilon} \f$
       */
      template < typename T >
      Eigen::Matrix< T, 6, 1 > dRhoStrain_dStrain( T rhoStrain, const Eigen::Matrix< T, 6, 1 >& strain )
      {
        if ( Marmot::Math::makeReal( rhoStrain ) <= 1e-16 )
          return Eigen::Matrix< T, 6, 1 >::Zero();

        Eigen::Matrix< T, 6, 1 > e = IDev * strain;
        // P array results from the derivative of the double contraction e:e in voigt notation
        return T( 1. / rhoStrain ) * PInv.array() * e.array();
      }

      /**
       * Computes the derivative \f$ \frac{d\, \theta}{d\, \boldsymbol{\sigma}}\f$ of the haigh westergaard coordinate
       * \f$ \theta \f$ with respect to the voigt notated stress vector \f$ \boldsymbol{\sigma} \f$
       */
      Marmot::Vector6d dTheta_dStress( double theta, const Marmot::Vector6d& stress );

      /**
       * Computes the derivative \f$ \frac{d\, \theta}{d\, J_2}\f$ of the haigh westergaard coordinate \f$ \theta \f$
       * with respect to the second deviatoric invariant \f$ J_2 \f$
       */
      double dTheta_dJ2( const Marmot::Vector6d& stress );

      /**
       * Computes the derivative \f$ \frac{d\, \theta}{d\, J_3}\f$ of the haigh westergaard coordinate \f$ \theta \f$
       * with respect to the third deviatoric invariant \f$ J_3 \f$
       */
      double dTheta_dJ3( const Marmot::Vector6d& stress );

      /**
       * Computes the derivative \f$ \frac{d\, \theta^{(\varepsilon)}}{d\, J^{(\varepsilon)}_2}\f$ of the haigh
       * westergaard coordinate \f$ \theta^{(\varepsilon)} \f$ with respect to the second deviatoric invariant \f$
       * J^{(\varepsilon)}_2 \f$.
       */
      double dThetaStrain_dJ2Strain( const Marmot::Vector6d& strain );

      /**
       * Computes the derivative \f$ \frac{d\, \theta^{(\varepsilon)}}{d\, J^{(\varepsilon)}_3}\f$ of the haigh
       * westergaard coordinate \f$ \theta^{(\varepsilon)} \f$ with respect to the third deviatoric invariant \f$
       * J^{(\varepsilon)}_3 \f$.
       */
      double dThetaStrain_dJ3Strain( const Marmot::Vector6d& strain );

      /**
       * Computes the derivative \f$ \frac{d\, J_2}{d\, \boldsymbol{\sigma}}\f$ of the second deviatoric invariant \f$
       * J_2 \f$ with respect to the voigt notated stress vector \f$ \boldsymbol{\sigma} \f$.
       */
      Marmot::Vector6d dJ2_dStress( const Marmot::Vector6d& stress );

      /**
       * Computes the derivative \f$ \frac{d\, J_3}{d\, \boldsymbol{\sigma}}\f$ of the third deviatoric invariant \f$
       * J_3 \f$ with respect to the voigt notated stress vector \f$ \boldsymbol{\sigma} \f$.
       */
      Marmot::Vector6d dJ3_dStress( const Marmot::Vector6d& stress );

      /**
       * Computes the derivative \f$ \frac{d\, J^{(\varepsilon)}_2}{d\, \boldsymbol{\sigma}}\f$ of the second deviatoric
       * invariant \f$ J^{(\varepsilon)}_2 \f$ with respect to the voigt notated strain vector \f$
       * \boldsymbol{\varepsilon} \f$.
       */
      Marmot::Vector6d dJ2Strain_dStrain( const Marmot::Vector6d& strain );

      /**
       * Computes the derivative \f$ \frac{d\, J^{(\varepsilon)}_3}{d\, \boldsymbol{\sigma}}\f$ of the third deviatoric
       * invariant \f$ J^{(\varepsilon)}_3 \f$ with respect to the voigt notated strain vector \f$
       * \boldsymbol{\varepsilon} \f$.
       */
      Marmot::Vector6d dJ3Strain_dStrain( const Marmot::Vector6d& strain );

      /**
       * Computes the derivative \f$ \frac{d\, \theta^{(\varepsilon)}}{d\, \boldsymbol{\varepsilon}}\f$ of the haigh
       * westergaard coordinate \f$ \theta^{(\varepsilon)} \f$ with respect to the voigt notated strain vector \f$
       * \boldsymbol{\varepsilon} \f$
       */
      Marmot::Vector6d dThetaStrain_dStrain( const Marmot::Vector6d& strain );

      // derivatives of principalStess with respect to stress

      /**
       * Computes the derivative \f$ \frac{d\, \sigma_I}{d\, \boldsymbol{\sigma}}\f$ of the principal stresses  \f$
       * \sigma_I \f$ with respect to the voigt notated stress vector \f$ \boldsymbol{\sigma} \f$
       */
      Marmot::Matrix36 dStressPrincipals_dStress( const Marmot::Vector6d& stress );

      // derivatives of plastic strains with respect to strains

      /**
       * Computes the derivative \f$ \frac{d\, \varepsilon^{vol}_{\ominus}}{d\, \varepsilon_I}\f$ of the volumetric
       * strains in compression  \f$ \varepsilon^{vol}_{\ominus} \f$ with respect to the principal strains  \f$
       * \varepsilon_I \f$
       */
      Eigen::Vector3d dStrainVolumetricNegative_dStrainPrincipal( const Marmot::Vector6d& strain );

      /**
       * Computes the derivative \f$ \frac{d\, \boldsymbol{\varepsilon}^{p}}{d\, \boldsymbol{\varepsilon}}\f$ of the
       voigt notated plastic strain vector \f$ \boldsymbol{\varepsilon}^{p} \f$ with respect to the voigt notated strain
       vector  \f$ \boldsymbol{\varepsilon} \f$
       *
       *\f[
           \displaystyle \frac{d\, \boldsymbol{\varepsilon}^{p}}{d\, \boldsymbol{\varepsilon}} = \boldsymbol{I} -
       \mathbb{C}^{-1}\,\mathbb{C}^{(ep)} \f]

       *using the elastic compliance tensor \f$ \mathbb{C}^{-1} \f$ and the elastoplastic stiffness tensor \f$
       \mathbb{C}^{(ep)} \f$
       */
      Matrix6d dEp_dE( const Matrix6d& CelInv, const Matrix6d& Cep );

      /**
       * Computes the derivative \f$ \frac{d\, \Delta\, \varepsilon^{p, vol}}{d\, \boldsymbol{\varepsilon}}\f$ of the
       * volumetric plastic strain increment \f$ \Delta\, \varepsilon^{p, vol}\f$ with respect to the voigt notated
       * strain vector  \f$ \boldsymbol{\varepsilon} \f$
       */
      RowVector6d dDeltaEpv_dE( const Matrix6d& CelInv, const Matrix6d& Cep );

      /**
       * Computes the derivative \f$ \frac{d\, \varepsilon_I}{d\, \boldsymbol{\varepsilon}}\f$ of the principal strains
       * \f$ \varepsilon_I \f$ with respect to the voigt notated strain vector  \f$ \boldsymbol{\varepsilon} \f$
       */
      Marmot::Matrix36 dSortedStrainPrincipal_dStrain( const Marmot::Vector6d& dEp );

      /**
       * Computes the derivative \f$ \frac{d\, \Delta\, \varepsilon^{p, vol}_{\ominus}}{d\, \boldsymbol{\varepsilon}}\f$
       * of the volumetric plastic strain increment in compression \f$ \Delta\, \varepsilon^{p, vol}_{\ominus}\f$ with
       * respect to the voigt notated strain vector  \f$ \boldsymbol{\varepsilon} \f$
       */
      RowVector6d dDeltaEpvneg_dE( const Marmot::Vector6d& dEp, const Matrix6d& CelInv, const Matrix6d& Cep );

    } // namespace Derivatives

    namespace Transformations {

      /**
       * Computes the transformation matrix \f$ R_{\varepsilon} \f$ to transform a voigt notated strain vector \f$
       * \boldsymbol{\varepsilon} \f$ to another cartesian coordinate system
       */
      Matrix6d transformationMatrixStrainVoigt( const Matrix3d& transformedCoordinateSystem );

      /**
       * Computes the transformation matrix \f$ R_{\sigma} \f$ to transform a voigt notated stress vector \f$
       * \boldsymbol{\sigma} \f$  to another cartesian coordinate system
       */
      Matrix6d transformationMatrixStressVoigt( const Matrix3d& transformedCoordinateSystem );

      /**
       * Returns the projection matrix to calculate the stress vector \f$ \boldsymbol{t}^{(n)} \f$ effective on a plane
     orientated with the normal vector \f$ \boldsymbol{n} \f$ from a voigt notated stress vector following cauchy's
     formula.
       *\f[
     \displaystyle t^{(n)}_i = \sigma_{ij}\,n_j
        \f]
       */
      Matrix36d projectVoigtStressToPlane( const Vector3d& normalVector );

      /**
       * Returns the projection matrix to calculate the strain vector \f$ \boldsymbol{\varepsilon}^{(n)} \f$ effective
       * on a plane orientated with the normal vector \f$ \boldsymbol{n} \f$ from a voigt notated strain vector (see
       * projectVoigtStressToPlane())).
       */
      Matrix36d projectVoigtStrainToPlane( const Vector3d& normalVector );

      /**
       * Rotates a stress tensor \f$ \boldsymbol{\sigma} \f$ applying a rotation matrix \f$ \boldsymbol{Q} \f$ in voigt
     notation.
       *\f[
     \displaystyle \boldsymbol{\sigma}^{\prime} = \boldsymbol{Q} \cdot \boldsymbol{\sigma} \cdot \boldsymbol{Q}^{T}
        \f]
        */
      Marmot::Vector6d rotateVoigtStress( const Eigen::Matrix3d& Q, const Marmot::Vector6d& stress );

    } // namespace Transformations

  }   // namespace ContinuumMechanics::VoigtNotation
} // namespace Marmot
