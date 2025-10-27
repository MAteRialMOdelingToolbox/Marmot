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

namespace Marmot {
  /**
   * @namespace ContinuumMechanics::VoigtNotation
   * @brief This namespace includes functions needed for calculations with stress and strain tensors written in voigt
   * notation.
   */
  namespace ContinuumMechanics::VoigtNotation {

    /* constexpr int VoigtSize = 6; */

    enum VoigtSize { OneD = 1, TwoD = 3, ThreeD = 6, Axial = 4 };

    /**
     * @brief Computes the size of the Voigt notation array based on the spatial dimension.
     * @param x The spatial dimension (e.g., 2 for 2D, 3 for 3D).
     * @return The size of the Voigt notation array as a VoigtSize.
     * @note The calculation is based on the formula: \f$ \frac{x * x + x}{2} \f$.
     */
    constexpr VoigtSize voigtSizeFromDimension( int x )
    {
      return (VoigtSize)( ( ( x * x ) + x ) >> 1 );
    }

    /** @brief Predefined 6D vector with scaling factor (2x) for shear components in Voigt notation.
     * @details The vector contains scaling factors: {1, 1, 1, 2, 2, 2}.
     */
    extern const Marmot::Vector6d P;

    /** @brief Predefined 6D vector with inverse scaling factor (0.5x) for shear components in Voigt notation.
     * @details The vector contains inverse scaling factors: {1, 1, 1, 0.5, 0.5, 0.5}.
     */
    extern const Marmot::Vector6d PInv;

    /** @brief Predefined 6D Vector representing the identity tensor in Voigt notation.
     * @details The vector contains ones in all components: {1, 1, 1, 0, 0, 0}.
     */
    extern const Marmot::Vector6d I;

    /** @brief Predefined 6D Vector representing the hydrostatic projection tensor in Voigt notation.
     * @details The vector contains: {1/3, 1/3, 1/3, 0, 0, 0}.
     */
    extern const Marmot::Vector6d IHyd;

    /** @brief Deviatoric projection tensor in Voigt notation.
     * @details \f$ 6 \times 6 \f$ matrix.
     */
    extern const Matrix6d IDev;

    // Plane Stress handling

    /**
     * @brief Converts a 3D stress vector in Voigt notation to a plane stress vector.
     * @param voigt 6-component 3D stress vector in Voigt notation:
     * \f[
     * \begin{bmatrix}
     * \sigma_{11} & \sigma_{22} & \sigma_{33} & \sigma_{12} & \sigma_{13} & \sigma_{23}
     * \end{bmatrix}^T
     * \f]
     * where \f$\sigma_{33} = \sigma_{13} = \sigma_{23} = 0\f$ in plane stress conditions.
     * @return 3-component plane stress vector in Voigt notation:
     * \f[
     * \begin{bmatrix}
     * \sigma_{11} & \sigma_{22} & \sigma_{12}
     * \end{bmatrix}^T
     * \f]
     * @note This function extracts only the in-plane components of a 3D stress vector.
     */
    Eigen::Vector3d voigtToPlaneVoigt( const Marmot::Vector6d& voigt );
    Vector4d        voigtToAxisymmetricVoigt( const Vector6d& voigt );

    /**
     * @brief Converts a plane stress vector in Voigt notation back to a 3D Voigt vector.
     * @param voigtPlane 3-component plane stress vector in Voigt notation:
     * \f[
     * \begin{bmatrix}
     * \sigma_{11} & \sigma_{22} & \sigma_{12}
     * \end{bmatrix}^T
     * \f]
     * @return 6-component 3D stress vector in Voigt notation:
     * \f[
     * \begin{bmatrix}
     * \sigma_{11} & \sigma_{22} & 0 & \sigma_{12} & 0 & 0
     * \end{bmatrix}^T
     * \f]
     * @note The out-of-plane components are set to zero, suitable for plane stress conditions.
     */
    Marmot::Vector6d planeVoigtToVoigt( const Eigen::Vector3d& voigtPlane );
    Vector6d         axisymmetricVoigtToVoigt( const Vector4d& voigtAxisymmetric );

    /**
     * @brief Reduces a 3D Voigt vector to a lower-dimensional vector defined by `voigtSize`.
     * @tparam voigtSize Target Voigt dimension (1, 3, or 6).
     * @param Voigt3D Input 6-component 3D Voigt vector.
     * @return Reduced Voigt vector of size `voigtSize`.
     * @note Considered cases:
     * - `voigtSize = 1`: Returns the \f$\sigma_{11}\f$ component.
     * - `voigtSize = 3`: Returns the plane stress components using `voigtToPlaneVoigt()`.
     * - `voigtSize = 6`: Returns the input 3D vector unchanged.
     * @throws std::invalid_argument if `voigtSize` is not 1, 3, or 6.
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
     * @brief Constructs a 3D Voigt vector from a lower-dimensional Voigt vector.
     * @tparam voigtSize Dimension of the input Voigt vector (1, 3, or 6).
     * @param Voigt Input Voigt vector of size `voigtSize`.
     * @return 6-component 3D Voigt vector.
     * @note Considered cases:
     * - `voigtSize = 1`: Returns a 3D vector with only \f$\sigma_{11}\f$ non-zero:
     * \f[
     * \begin{bmatrix}
     * \sigma_{11} & 0 & 0 & 0 & 0 & 0
     * \end{bmatrix}^T
     * \f]
     * - `voigtSize = 3`: Returns a 3D Voigt vector by calling `planeVoigtToVoigt()`.
     * - `voigtSize = 6`: Returns the input vector unchanged.
     * @throws std::invalid_argument if `voigtSize` is not 1, 3, or 6.
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
     * @brief Converts a 6-component Voigt strain vector to a 3x3 strain tensor.
     * @tparam T Scalar type (e.g., double, float).
     * @param voigt 6-component strain vector in Voigt notation:
     * \f[
     * \begin{bmatrix}
     * \varepsilon_{11} & \varepsilon_{22} & \varepsilon_{33} & \gamma_{12} & \gamma_{13} & \gamma_{23}
     * \end{bmatrix}^T
     * \f]
     * @return 3x3 strain tensor:
     * \f[
     * \begin{bmatrix}
     * \varepsilon_{11} & \gamma_{12}/2 & \gamma_{13}/2 \\
     * \gamma_{12}/2 & \varepsilon_{22} & \gamma_{23}/2 \\
     * \gamma_{13}/2 & \gamma_{23}/2 & \varepsilon_{33} \\
     * \end{bmatrix}
     * \f]
     * @note Shear components in Voigt notation are engineering strains and are halved in the resulting tensor.
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
     * @brief Converts a 6-component Voigt stress vector to a 3x3 stress tensor.
     * @tparam T Scalar type (e.g., double, float).
     * @param voigt 6-component stress vector in Voigt notation:
     * \f[
     * \begin{bmatrix}
     * \sigma_{11} & \sigma_{22} & \sigma_{33} & \sigma_{12} & \sigma_{13} & \sigma_{23}
     * \end{bmatrix}^T
     * \f]
     * @return 3x3 stress tensor:
     * \f[
     * \begin{bmatrix}
     * \sigma_{11} & \sigma_{12} & \sigma_{13} \\
     * \sigma_{12} & \sigma_{22} & \sigma_{23} \\
     * \sigma_{13} & \sigma_{23} & \sigma_{33} \\
     * \end{bmatrix}
     * \f]
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
     * @brief Converts a 3x3 strain tensor to a 6-component Voigt strain vector.
     * @param strainTensor 3x3 strain tensor:
     * \f[
     * \boldsymbol{\varepsilon} =
     * \begin{bmatrix}
     * \varepsilon_{11} & \varepsilon_{12} & \varepsilon_{13} \\
     * \varepsilon_{21} & \varepsilon_{22} & \varepsilon_{23} \\
     * \varepsilon_{31} & \varepsilon_{32} & \varepsilon_{33} \\
     * \end{bmatrix}
     * \f]
     * @return 6-component strain vector in Voigt notation:
     * \f[
     * \begin{aligned}
     * \boldsymbol{\varepsilon}^{\text{Voigt}} &=
     * \begin{bmatrix}
     *   \varepsilon_{11} & \varepsilon_{22} & \varepsilon_{33} & \gamma_{12} & \gamma_{13} & \gamma_{23}
     * \end{bmatrix}^T \\
     * &=
     * \begin{bmatrix}
     *   \varepsilon_{11} & \varepsilon_{22} & \varepsilon_{33} & 2 \varepsilon_{12} & 2 \varepsilon_{13} & 2
     * \varepsilon_{23}
     * \end{bmatrix}^T
     * \end{aligned}
     * \f]
     * @note Shear components in the Voigt vector (\f$\gamma_{ij}\f$ ) are engineering strains,
     * which are twice the tensor shear components.
     */
    Marmot::Vector6d strainToVoigt( const Eigen::Matrix3d& strainTensor );

    /**
     * @brief Converts a 3x3 stress tensor to a 6-component Voigt stress vector.
     * @tparam T Scalar type (e.g., double, float).
     * @param stressTensor 3x3 stress tensor:
     * \f[
     * \boldsymbol{\sigma} =
     * \begin{bmatrix}
     * \sigma_{11} & \sigma_{12} & \sigma_{13} \\
     * \sigma_{21} & \sigma_{22} & \sigma_{23} \\
     * \sigma_{31} & \sigma_{32} & \sigma_{33} \\
     * \end{bmatrix}
     * \f]
     * @return 6-component stress vector in Voigt notation:
     * \f[
     * \boldsymbol{\sigma}^{\text{Voigt}} =
     * \begin{bmatrix}
     * \sigma_{11} & \sigma_{22} & \sigma_{33} & \sigma_{12} & \sigma_{13} & \sigma_{23}
     * \end{bmatrix}^T
     * \f]
     * @note Shear components are copied directly; no factor of 2 is applied.
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

    /**
     * @brief Converts a fourth-order stiffness tensor to its Voigt notation representation (\f$ 6 \times 6 \f$ matrix).
     * @param C The fourth-order stiffness tensor represented as an Eigen::Tensor<double, 4>.
     * @return An Eigen::Matrix<double, 6, 6> representing the stiffness tensor in Voigt notation.
     * @note The input tensor `C` must follow the symmetry properties of a stiffness tensor for the
     *       conversion to be valid (minor symmetry).
     */
    Eigen::Matrix< double, 6, 6 > stiffnessToVoigt( const Eigen::Tensor< double, 4 >& C );

    /**
     * @brief Converts a stiffness matrix in Voigt notation (\f$ 6 \times 6 \f$ matrix) to a 4th-order stiffness tensor
     * (\f$ 3 \times 3 \times 3 \times 3 \f$ tensor).
     * @param voigtStiffness The \f$ 6 \times 6 \f$ matrix representing the stiffness in Voigt notation.
     * @return An Eigen::Tensor of rank 4 (4th-order tensor) representing the stiffness tensor.
     *         The dimensions of the tensor are \f$ 3 \times 3 \times 3 \times 3 \f$.
     */
    Eigen::Tensor< double, 4 > voigtToStiffness( const Eigen::Matrix< double, 6, 6 >& voigtStiffness );

    /**
     * @brief Converts a stress vector in Voigt notation to its corresponding tensor form.
     *
     * @details This function handles 1D, 2D, and 3D stress vectors in Voigt notation and returns
     * the full symmetric stress tensor.
     * @tparam nDim Dimension of the stress tensor (1, 2, or 3).
     * @param Voigt Input stress vector in Voigt notation, size `VOIGTFROMDIM(nDim) x 1`:
     * - 1D: \f$ \boldsymbol{\sigma}^{\text{Voigt}} =
     *   \begin{bmatrix} \sigma_{11} \end{bmatrix}^T \f$
     * - 2D: \f$ \boldsymbol{\sigma}^{\text{Voigt}} =
     *   \begin{bmatrix} \sigma_{11} & \sigma_{22} & \sigma_{12} \end{bmatrix}^T \f$
     * - 3D: \f$ \boldsymbol{\sigma}^{\text{Voigt}} =
     *   \begin{bmatrix} \sigma_{11} & \sigma_{22} & \sigma_{33} & \sigma_{12} & \sigma_{13} & \sigma_{23}
     * \end{bmatrix}^T \f$
     * @return Symmetric stress tensor of size \f$ nDim \times nDim \f$.
     * @throws std::invalid_argument If `nDim` is not 1, 2, or 3.
     * @note For 2D and 3D, the off-diagonal shear components in the tensor are taken
     *       directly from the Voigt vector (no factor of 2 applied).
     */
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

    /**
     * @brief Converts a strain tensor to its corresponding Voigt notation vector.
     *
     * @details This function supports 1D, 2D, and 3D strain tensors and returns the corresponding
     * Voigt vector (column vector). Shear components are scaled appropriately for engineering strain.
     * @tparam nDim Dimension of the strain tensor (1, 2, or 3).
     * @param strain Input strain tensor of size `nDim x nDim`.
     * @return Column vector of size `VOIGTFROMDIM(nDim) x 1` representing the strain in Voigt notation.
     * - 1D: \f$ \boldsymbol{\varepsilon}^{\text{Voigt}} =
     *   \begin{bmatrix} \varepsilon_{11} \end{bmatrix}^T \f$
     * - 2D: \f$ \boldsymbol{\varepsilon}^{\text{Voigt}} =
     *   \begin{bmatrix} \varepsilon_{11} & \varepsilon_{22} & 2 \varepsilon_{12} \end{bmatrix}^T \f$
     * - 3D: \f$ \boldsymbol{\varepsilon}^{\text{Voigt}} =
     *   \begin{bmatrix} \varepsilon_{11} & \varepsilon_{22} & \varepsilon_{33} & 2 \varepsilon_{12} & 2
     * \varepsilon_{13} & 2 \varepsilon_{23} \end{bmatrix}^T \f$
     * @throws std::invalid_argument If `nDim` is not 1, 2, or 3.
     * @note Axisymmetric voigt vectors are not supported.
     */
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

      /**
       * @brief Computes the principal strains by solving the eigenvalue problem.
       *
       * @details This function computes the eigenvalues of the strain tensor reconstructed
       * from the given Voigt strain vector. The eigenvalues correspond to the
       * principal strains \f$ \varepsilon_1, \varepsilon_2, \varepsilon_3 \f$:
       *\f[
           \displaystyle |\varepsilon_{ij} - \lambda\, \delta_{ij}| = 0 \hspace{.5cm} \Rightarrow \hspace{.5cm}
       \lambda^{(1)},\lambda^{(2)},\lambda^{(3)}\hspace{0.3cm} \widehat{=}\hspace{0.3cm} \varepsilon_1,\,
       \varepsilon_2,\, \varepsilon_3 \f]
       * @param strain 6-component strain vector in Voigt notation.
       * @return A 3D vector containing the principal strains (unsorted).
       */
      Eigen::Vector3d principalStrains( const Marmot::Vector6d& strain );

      /**
       * @brief Computes the principal stresses by solving the eigenvalue problem.
       *
       * @details This function computes the eigenvalues of the stress tensor reconstructed
       * from the given Voigt stress vector. The eigenvalues correspond to the
       * principal stresses \f$ \sigma_1, \sigma_2, \sigma_3 \f$:
       * \f[
          \displaystyle |\sigma_{ij} - \lambda\, \delta_{ij}| = 0 \hspace{.5cm} \Rightarrow \hspace{.5cm}
      \lambda^{(1)},\lambda^{(2)},\lambda^{(3)}\hspace{0.3cm} \widehat{=}\hspace{0.3cm}\sigma_1,\, \sigma_2,\, \sigma_3
       \f]
       * @param stress 6-component stress vector in Voigt notation.
       * @return A 3D vector containing the principal stresses (unsorted).
       */
      Eigen::Vector3d principalStresses( const Marmot::Vector6d& stress );
      // principal strains calculated from haigh westergaard strains ( sorted --> e1 > e2 > e3 )

      /**
       * @brief Computes the principal strains from Haigh–Westergaard coordinates.
       *
       * @details This function transforms the Haigh–Westergaard coordinates
       * \f$ (\xi, \rho, \theta) \f$ into the corresponding principal strains
       * \f$ \varepsilon_1, \varepsilon_2, \varepsilon_3 \f$ using:
       * \f[
       *   \begin{bmatrix}
       *     \varepsilon_{1} \\
       *     \varepsilon_{2} \\
       *     \varepsilon_{3}
       *   \end{bmatrix}
       *   =
       *   \frac{1}{\sqrt{3}}
       *   \begin{bmatrix}
       *     \xi \\
       *     \xi \\
       *     \xi
       *   \end{bmatrix}
       *   + \sqrt{\tfrac{2}{3}}\, \rho \,
       *   \begin{bmatrix}
       *     \cos(\theta) \\
       *     -\sin\!\left(\tfrac{\pi}{6} - \theta\right) \\
       *     -\sin\!\left(\tfrac{\pi}{6} + \theta\right)
       *   \end{bmatrix}
       * \f]
       * @param strain 6-component strain vector in Voigt notation.
       * @return A 3D vector of the principal strains, sorted by magnitude.
       * @note The computation of \f$ \xi \f$, \f$ \rho \f$, and \f$ \theta \f$
       *       is performed in `haighWestergaardFromStrain()`.
       */
      Eigen::Vector3d sortedPrincipalStrains( const Marmot::Vector6d& strain );
      // principal stressDirections calculated by solving eigenvalue problem ( !NOT sorted! )

      /**
       * @brief Computes the principal stress directions corresponding to the principal stresses.
       *
       * @details This function solves the eigenvalue problem
       * \f[
       *   \left( \boldsymbol{\sigma} - \sigma_k \cdot \boldsymbol{I} \right) \cdot
       *   \boldsymbol{x}^{(k)} = 0
       * \f]
       * to obtain the eigenvectors \f$ \boldsymbol{x}^{(k)} \f$ associated with
       * the principal stresses \f$ \sigma_k \f$.
       * @param stress 6-component stress vector in Voigt notation.
       * @return A 3x3 matrix whose columns are the principal stress directions (unsorted).
       */
      Eigen::Matrix3d principalStressesDirections( const Marmot::Vector6d& stress );

      /**
       * @brief Computes the equivalent von Mises stress.
       *
       * @details Defined as:
       * \f[
       *   \sigma_\text{eq} = \sqrt{3 \cdot J_2}
       * \f]
       * where \f$ J_2 \f$ is the second invariant of the deviatoric stress tensor
       * (see `J2()`).
       * @param stress 6-component stress vector in Voigt notation.
       * @return The von Mises equivalent stress.
       */
      double vonMisesEquivalentStress( const Marmot::Vector6d& stress );

      /**
       * @brief Computes the equivalent von Mises strain.
       *
       * @details Based on the deviatoric strain tensor \f$ e_{ij} \f$:
       * \f[
       *   \varepsilon^\text{(eq)} = \sqrt{\tfrac{2}{3} \, e_{ij} e_{ij}}
       * \f]
       * @param strain 6-component strain vector in Voigt notation.
       * @return The von Mises equivalent strain.
       */
      double vonMisesEquivalentStrain( const Marmot::Vector6d& strain );

      /**
       * @brief Computes the Euclidean norm of the strain tensor.
       *
       * @details Defined as:
       * \f[
       *   ||\boldsymbol{\varepsilon}|| = \sqrt{ \varepsilon_{ij}\varepsilon_{ij}}
       * \f]
       * @tparam T Scalar type (e.g., double, float).
       * @param strain 6-component strain vector in Voigt notation.
       * @return The Euclidean norm of the strain tensor.
       */
      template < typename T >
      T normStrain( const Eigen::Matrix< T, 6, 1 >& strain )
      {
        return ContinuumMechanics::VoigtNotation::voigtToStrain( strain ).norm();
      }

      /**
       * @brief Computes the Euclidean norm of the stress tensor.
       *
       * @details Defined as:
       * \f[
       *   ||\boldsymbol{\sigma}|| = \sqrt{\sigma_{ij}\sigma_{ij}}
       * \f]
       * @param stress 6-component stress vector in Voigt notation.
       * @return The Euclidean norm of the stress tensor.
       */
      double normStress( const Marmot::Vector6d& stress );
      // Trace of compressive strains

      /**
       * @brief Computes the volumetric plastic strain in compression.
       *
       * @details Defined as:
       * \f[
       *   \varepsilon^{\text{vol}}_{\ominus}
       *   = \sum_{i=1}^3 \left\langle -\varepsilon_i \right\rangle
       * \f]
       * where \f$ \varepsilon_i \f$ are the principal strains and
       * \f$ \langle \bullet \rangle \f$ denotes the Macaulay brackets.
       * @param strain 6-component strain vector in Voigt notation.
       * @return The volumetric compressive strain.
       */
      double StrainVolumetricNegative( const Marmot::Vector6d& strain );

      /**
       * @brief Computes the first invariant \f$ I_1 \f$ of the stress tensor.
       *
       * @details The first invariant of the stress tensor \f$ \boldsymbol{\sigma} \f$ is the trace:
       * \f[
       *   I_1 = \operatorname{tr}(\boldsymbol{\sigma})
       *       = \sigma_{11} + \sigma_{22} + \sigma_{33}
       * \f]
       * @tparam T Scalar type (e.g., double, float).
       * @param stress 6-component stress vector in Voigt notation:
       *   \f$ \boldsymbol{\sigma}^{\text{Voigt}} =
       *     \begin{bmatrix}
       *       \sigma_{11} & \sigma_{22} & \sigma_{33} & \sigma_{12} & \sigma_{13} & \sigma_{23}
       *     \end{bmatrix}^T \f$
       * @return Value of the first invariant \f$ I_1 \f$.
       */
      template < typename T >
      T I1( const Eigen::Matrix< T, 6, 1 >& stress )
      {
        return stress( 0 ) + stress( 1 ) + stress( 2 );
      }

      /**
       * @brief Computes the first invariant \f$ I^{(\varepsilon)}_1 \f$ of the strain tensor.
       *
       * @details The first invariant of the strain tensor \f$ \boldsymbol{\varepsilon} \f$ is the trace:
       * \f[
       *   I^{(\varepsilon)}_1 = \operatorname{tr}(\boldsymbol{\varepsilon})
       *       = \varepsilon_{11} + \varepsilon_{22} + \varepsilon_{33}
       * \f]
       * @param strain 6-component strain vector in Voigt notation:
       *   \f$ \boldsymbol{\varepsilon}^{\text{Voigt}} =
       *     \begin{bmatrix}
       *       \varepsilon_{11} & \varepsilon_{22} & \varepsilon_{33} &
       *       \gamma_{12} & \gamma_{13} & \gamma_{23}
       *     \end{bmatrix}^T \f$
       * @return Value of the first strain invariant \f$ I^{(\varepsilon)}_1 \f$.
       */
      double I1Strain( const Marmot::Vector6d& strain ); //

      /**
       * @brief Computes the second invariant \f$ I_2 \f$ of the stress tensor.
       *
       * @details The second invariant of the stress tensor \f$ \boldsymbol{\sigma} \f$ is defined as:
       * \f[
       *   I_2 = \tfrac{1}{2}\left( \operatorname{tr}(\boldsymbol{\sigma})^2
       *         - \operatorname{tr}(\boldsymbol{\sigma}^2) \right)
       * \f]
       * In Voigt notation, this expands to:
       * \f[
       *   I_2 = \sigma_{11}\sigma_{22} + \sigma_{22}\sigma_{33} + \sigma_{33}\sigma_{11}
       *         - \sigma_{12}^2 - \sigma_{13}^2 - \sigma_{23}^2
       * \f]
       * @tparam T Scalar type (e.g., double, float).
       * @param stress 6-component stress vector in Voigt notation:
       *   \f$ \boldsymbol{\sigma}^{\text{Voigt}} =
       *     \begin{bmatrix}
       *       \sigma_{11} & \sigma_{22} & \sigma_{33} & \sigma_{12} & \sigma_{13} & \sigma_{23}
       *     \end{bmatrix}^T \f$
       *
       * @return Value of the second invariant \f$ I_2 \f$.
       */
      template < typename T >
      T I2( const Eigen::Matrix< T, 6, 1 >& stress )
      {
        const Eigen::Matrix< T, 6, 1 >& s = stress;
        return s( 0 ) * s( 1 ) + s( 1 ) * s( 2 ) + s( 2 ) * s( 0 ) - s( 3 ) * s( 3 ) - s( 4 ) * s( 4 ) -
               s( 5 ) * s( 5 );
      }

      /**
       * @brief Computes the second invariant \f$ I^{(\varepsilon)}_2 \f$ of the strain tensor.
       *
       * @details For a strain tensor \f$ \boldsymbol{\varepsilon} \f$ in Voigt notation,
       * the second invariant is:
       * \f[
       *   I^{(\varepsilon)}_2 =
       *     \varepsilon_{11}\varepsilon_{22} +
       *     \varepsilon_{22}\varepsilon_{33} +
       *     \varepsilon_{33}\varepsilon_{11}
       *     - \tfrac{1}{4}\left(\gamma_{12}^2 + \gamma_{13}^2 + \gamma_{23}^2\right)
       * \f]
       * @param strain 6-component strain vector in Voigt notation:
       *   \f$ \boldsymbol{\varepsilon}^{\text{Voigt}} =
       *     [ \varepsilon_{11}, \varepsilon_{22}, \varepsilon_{33},
       *       \gamma_{12}, \gamma_{13}, \gamma_{23} ]^T \f$
       * @return Value of the second invariant \f$ I^{(\varepsilon)}_2 \f$.
       */
      double I2Strain( const Marmot::Vector6d& strain );

      /**
       * @brief Computes the third invariant \f$ I_3 \f$ of the stress tensor.
       *
       * @details Defined as the determinant of the Cauchy stress tensor \f$ \boldsymbol{\sigma} \f$:
       * \f[
       *   I_3 = \det(\boldsymbol{\sigma})
       * \f]
       * Expanded in Voigt notation:
       * \f[
       *   I_3 =
       *     \sigma_{11}\sigma_{22}\sigma_{33}
       *     + 2\sigma_{12}\sigma_{13}\sigma_{23}
       *     - \sigma_{11}\sigma_{23}^2
       *     - \sigma_{22}\sigma_{13}^2
       *     - \sigma_{33}\sigma_{12}^2
       * \f]
       * @tparam T Scalar type (e.g., double, float).
       * @param stress 6-component stress vector in Voigt notation.
       * @return Value of the third invariant \f$ I_3 \f$.
       */
      template < typename T >
      T I3( const Eigen::Matrix< T, 6, 1 >& stress )
      {
        const Eigen::Matrix< T, 6, 1 >& s = stress;
        return s( 0 ) * s( 1 ) * s( 2 ) + 2. * s( 3 ) * s( 4 ) * s( 5 ) - s( 0 ) * s( 5 ) * s( 5 ) -
               s( 1 ) * s( 4 ) * s( 4 ) - s( 2 ) * s( 3 ) * s( 3 );
      }

      /**
       * @brief Computes the third invariant \f$ I^{(\varepsilon)}_3 \f$ of the strain tensor.
       *
       * @details The third invariant is given by the determinant of the strain tensor:
       * \f[
       *   I^{(\varepsilon)}_3 = \det(\boldsymbol{\varepsilon})
       * \f]
       * @param strain 6-component strain vector in Voigt notation.
       * @return Value of the third strain invariant \f$ I^{(\varepsilon)}_3 \f$.
       */
      double I3Strain( const Marmot::Vector6d& strain );

      /**
       * @brief Computes the second invariant \f$ J_2 \f$ of the deviatoric stress tensor.
       *
       * @details For the deviatoric stress tensor \f$ \boldsymbol{s} \f$:
       * \f[
       *   J_2 = \tfrac{1}{3} I_1^2 - I_2
       * \f]
       * @tparam T Scalar type (e.g., double, float).
       * @param stress 6-component stress vector in Voigt notation.
       * @return Value of the second deviatoric invariant \f$ J_2 \f$ (clamped to non-negative).
       */
      template < typename T >
      T J2( const Eigen::Matrix< T, 6, 1 >& stress )
      {
        const T I1_ = I1( stress );
        const T I2_ = I2( stress );
        const T res = ( 1. / 3 ) * I1_ * I1_ - I2_;
        return Marmot::Math::makeReal( res ) >= 0 ? res : T( 0.0 );
      }

      /**
       * @brief Computes the second invariant \f$ J^{(\varepsilon)}_2 \f$ of the deviatoric strain tensor.
       *
       * @details For the deviatoric strain tensor \f$ \boldsymbol{e} \f$:
       * \f[
       *   J^{(\varepsilon)}_2 =
       *     \tfrac{1}{3} \left(I^{(\varepsilon)}_1\right)^2
       *     - I^{(\varepsilon)}_2
       * \f]
       * @param strain 6-component strain vector in Voigt notation.
       * @return Value of the second deviatoric strain invariant \f$ J^{(\varepsilon)}_2 \f$.
       */
      double J2Strain( const Marmot::Vector6d& strain );

      /**
       * @brief Computes the third invariant \f$ J_3 \f$ of the deviatoric stress tensor.
       *
       * @details For the deviatoric stress tensor \f$ \boldsymbol{s} \f$:
       * \f[
       *   J_3 = \tfrac{2}{27} I_1^3 - \tfrac{1}{3} I_1 I_2 + I_3
       * \f]
       * @tparam T Scalar type (e.g., double, float).
       * @param stress 6-component stress vector in Voigt notation.
       * @return Value of the third deviatoric stress invariant \f$ J_3 \f$.
       */
      template < typename T >
      T J3( const Eigen::Matrix< T, 6, 1 >& stress )
      {
        const T I1_ = I1( stress );
        const T I2_ = I2( stress );
        const T I3_ = I3( stress );

        return ( 2. / 27 ) * pow( I1_, 3 ) - ( 1. / 3 ) * I1_ * I2_ + I3_;
      }

      /**
       * @brief Computes the third invariant \f$ J^{(\varepsilon)}_3 \f$ of the deviatoric strain tensor \f$
       * \boldsymbol{e} \f$.
       *
       * The third invariant is obtained by reconstructing the strain tensor
       * from Voigt notation and computing its determinant:
       * \f[
       *   J^{(\varepsilon)}_3 = \det(\boldsymbol{e})
       * \f]
       * @param strain 6-component strain vector in Voigt notation.
       * @return Value of the third deviatoric strain invariant \f$ J^{(\varepsilon)}_3 \f$.
       */
      double J3Strain( const Marmot::Vector6d& strain );

      /**
       * @brief Computes principal values and their derivatives for a symmetric 3×3 matrix in Voigt notation.
       *
       * @details This routine implements a fast algorithm to compute the principal values
       * (eigenvalues) of a symmetric second-order tensor in Voigt notation
       * (off-diagonal entries are expected without the factor of 2).
       * @param S 6-component symmetric matrix in Voigt notation.
       * @return A pair consisting of:
       *   - A vector containing the (unsorted) principal values.
       *   - An \f$ 3 \times 6 \f$ matrix containing the derivatives of the
       *     principal values with respect to the Voigt components of \f$ S \f$.
       */
      std::pair< Eigen::Vector3d, Eigen::Matrix< double, 3, 6 > > principalValuesAndDerivatives(
        const Eigen::Matrix< double, 6, 1 >& S );

    } // namespace Invariants

    namespace Derivatives {

      // derivatives of Haigh Westergaard stresses with respect to cauchy stress in eng. notation

      /**
       * @brief Computes the derivative of the mean stress with respect to the stress vector.
       *
       * \f[
       *   \frac{\partial \sigma_m}{\partial \boldsymbol{\sigma}}
       * \f]
       * @return 6-component vector representing \f$ \tfrac{\partial \sigma_m}{\partial \boldsymbol{\sigma}} \f$
       *         in Voigt notation.
       */
      Vector6d dStressMean_dStress();

      /**
       * @brief Computes the derivative of the Haigh–Westergaard deviatoric radius \f$ \rho \f$
       *        with respect to the stress vector.
       *
       * \f[
       *   \frac{\partial \rho}{\partial \boldsymbol{\sigma}}
       * \f]
       * @tparam T Scalar type (e.g., double, float).
       * @param rho Precomputed Haigh–Westergaard coordinate \f$ \rho \f$.
       * @param stress 6-component stress vector in Voigt notation.
       * @return 6-component vector of partial derivatives
       *         \f$ \tfrac{\partial \rho}{\partial \boldsymbol{\sigma}} \f$ in Voigt notation.
       * @note Returns zero if \f$ \rho \leq 10^{-16} \f$ to avoid division by zero.
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
       * @brief Computes the derivative of the lode angle \f$ \theta \f$
       *        with respect to the stress vector.
       *
       * @param theta Lode angle \f$ \theta \f$.
       * @param stress 6-component stress vector in Voigt notation.
       * @return 6-component vector of partial derivatives
       *         \f$ \tfrac{\partial \theta}{\partial \boldsymbol{\sigma}} \f$ in Voigt notation.
       * @note Returns zero if \f$ \theta \leq 10^{-15} \f$ or if \f$ \theta \geq \frac{\pi}{3} - 10^{-15} \f$.
       */
      Marmot::Vector6d dTheta_dStress( double theta, const Marmot::Vector6d& stress );

      /**
       * @brief Computes the derivative \f$ \frac{\partial \theta}{\partial J_2}\f$ of the lode angle \f$ \theta
       * \f$ with respect to the second deviatoric invariant \f$ J_2 \f$.
       * @param stress 6-component stress vector in Voigt notation.
       * @return Scalar derivative \f$ \tfrac{\partial \theta}{\partial J_2} \f$.
       */
      double dTheta_dJ2( const Marmot::Vector6d& stress );

      /**
       * @brief Computes the derivative \f$ \frac{\partial \theta}{\partial J_3}\f$ of the Lode angle
       * \f$ \theta \f$ with respect to the third deviatoric invariant \f$ J_3 \f$.
       *
       * @param stress Stress vector in Voigt notation (\f$ \boldsymbol{\sigma} \f$).
       * @return Value of the derivative \f$ \frac{\partial \theta}{\partial J_3} \f$.
       */
      double dTheta_dJ3( const Marmot::Vector6d& stress );

      /**
       * @brief Computes the derivative \f$ \frac{\partial \theta^{(\varepsilon)}}{\partial J^{(\varepsilon)}_2}\f$ of
       * the Lode angle in strain space \f$ \theta^{(\varepsilon)} \f$ with respect to the second deviatoric invariant
       * \f$ J^{(\varepsilon)}_2 \f$.
       *
       * @param strain Strain vector in Voigt notation (\f$ \boldsymbol{\varepsilon} \f$).
       * @return Value of the derivative \f$ \frac{\partial \theta^{(\varepsilon)}}{\partial J^{(\varepsilon)}_2} \f$.
       */
      double dThetaStrain_dJ2Strain( const Marmot::Vector6d& strain );

      /**
       * @brief Computes the derivative \f$ \frac{\partial \theta^{(\varepsilon)}}{\partial J^{(\varepsilon)}_3}\f$ of
       * the Lode angle in strain space \f$ \theta^{(\varepsilon)} \f$ with respect to the third deviatoric invariant
       * \f$ J^{(\varepsilon)}_3 \f$.
       *
       * @param strain Strain vector in Voigt notation (\f$ \boldsymbol{\varepsilon} \f$).
       * @return Value of the derivative \f$ \frac{\partial \theta^{(\varepsilon)}}{\partial J^{(\varepsilon)}_3} \f$.
       */
      double dThetaStrain_dJ3Strain( const Marmot::Vector6d& strain );

      /**
       * @brief Computes the derivative \f$ \frac{\partial J_2}{\partial \boldsymbol{\sigma}}\f$ of the second
       * deviatoric invariant \f$ J_2 \f$ with respect to the stress vector in Voigt notation.
       *
       * @param stress Stress vector in Voigt notation (\f$ \boldsymbol{\sigma} \f$).
       * @return Vector of partial derivatives \f$ \frac{\partial J_2}{\partial \boldsymbol{\sigma}} \f$ in Voigt
       * notation.
       */
      Marmot::Vector6d dJ2_dStress( const Marmot::Vector6d& stress );

      /**
       * @brief Computes the derivative \f$ \frac{\partial J_3}{\partial \boldsymbol{\sigma}}\f$ of the third deviatoric
       * invariant \f$ J_3 \f$ with respect to the stress vector in Voigt notation.
       *
       * @param stress Stress vector in Voigt notation (\f$ \boldsymbol{\sigma} \f$).
       * @return Vector of partial derivatives \f$ \frac{\partial J_3}{\partial \boldsymbol{\sigma}} \f$ in Voigt
       * notation.
       */
      Marmot::Vector6d dJ3_dStress( const Marmot::Vector6d& stress );

      /**
       * @brief Computes the derivative \f$ \frac{\partial J^{(\varepsilon)}_2}{\partial \boldsymbol{\varepsilon}} \f$
       * of the second deviatoric invariant \f$ J^{(\varepsilon)}_2 \f$ with respect to the strain vector in Voigt
       * notation.
       *
       * @param strain Strain vector in Voigt notation (\f$ \boldsymbol{\varepsilon} \f$).
       * @return Vector of partial derivatives \f$ \frac{\partial J^{(\varepsilon)}_2}{\partial
       * \boldsymbol{\varepsilon}} \f$ in Voigt notation.
       */
      Marmot::Vector6d dJ2Strain_dStrain( const Marmot::Vector6d& strain );

      /**
       * @brief Computes the derivative \f$ \frac{\partial J^{(\varepsilon)}_3}{\partial \boldsymbol{\varepsilon}} \f$
       * of the third deviatoric invariant \f$ J^{(\varepsilon)}_3 \f$ with respect to the strain vector in Voigt
       * notation.
       *
       * @param strain Strain vector in Voigt notation (\f$ \boldsymbol{\varepsilon} \f$).
       * @return Vector of partial derivatives \f$ \frac{\partial J^{(\varepsilon)}_3}{\partial
       * \boldsymbol{\varepsilon}} \f$ in Voigt notation.
       */
      Marmot::Vector6d dJ3Strain_dStrain( const Marmot::Vector6d& strain );

      /**
       * @brief Computes the derivative \f$ \frac{\partial \theta^{(\varepsilon)}}{\partial \boldsymbol{\varepsilon}}
       * \f$ of the Lode angle in strain space \f$ \theta^{(\varepsilon)} \f$ with respect to the strain vector in Voigt
       * notation.
       *
       * @param strain Strain vector in Voigt notation (\f$ \boldsymbol{\varepsilon} \f$).
       * @return Vector of partial derivatives \f$ \frac{\partial \theta^{(\varepsilon)}}{\partial
       * \boldsymbol{\varepsilon}} \f$ in Voigt notation.
       */
      Marmot::Vector6d dThetaStrain_dStrain( const Marmot::Vector6d& strain );

      // derivatives of principalStess with respect to stress

      /**
       * @brief Computes the derivative \f$ \frac{\partial \sigma_I}{\partial \boldsymbol{\sigma}} \f$ of the principal
       * stresses
       * \f$ \sigma_I \f$ with respect to the stress vector in Voigt notation.
       *
       * @param stress Stress vector in Voigt notation (\f$ \boldsymbol{\sigma} \f$).
       * @return Matrix of partial derivatives \f$ \frac{\partial \sigma_I}{\partial \boldsymbol{\sigma}} \f$,
       * where each row corresponds to the gradient of one principal stress with respect to \f$ \boldsymbol{\sigma} \f$.
       */
      Marmot::Matrix36 dStressPrincipals_dStress( const Marmot::Vector6d& stress );

      // derivatives of plastic strains with respect to strains

      /**
       * @brief Computes the derivative \f$ \frac{\partial \varepsilon^{vol}_{\ominus}}{\partial \varepsilon_I} \f$ of
       * the volumetric compressive strain \f$ \varepsilon^{vol}_{\ominus} \f$ with respect to the principal strains
       * \f$ \varepsilon_I \f$.
       *
       * @param strain Principal/total strain provided in Voigt notation (\f$ \boldsymbol{\varepsilon} \f$), 6
       * components.
       * @return A 3-component vector (\f$ \partial \varepsilon^{vol}_{\ominus} / \partial \varepsilon_I \f$), where
       * each entry is the partial derivative of the corresponding volumetric-compressive component with respect to the
       * principal strains \f$ \varepsilon_1, \varepsilon_2, \varepsilon_3 \f$.
       */
      Eigen::Vector3d dStrainVolumetricNegative_dStrainPrincipal( const Marmot::Vector6d& strain );

      /**
       * @brief Computes the derivative \f$ \frac{\partial \boldsymbol{\varepsilon}^{p}}{\partial
       * \boldsymbol{\varepsilon}} \f$ of the Voigt-notated plastic strain vector \f$ \boldsymbol{\varepsilon}^{p} \f$
       * with respect to the Voigt-notated total strain vector \f$ \boldsymbol{\varepsilon} \f$.
       *
       * @details The relation used is
       * \f[
       *   \frac{\partial \boldsymbol{\varepsilon}^{p}}{\partial \boldsymbol{\varepsilon}}
       *   \;=\; \boldsymbol{I} \;-\; \mathbb{C}^{-1}\,\mathbb{C}^{(ep)}
       * \f]
       * where \f$ \mathbb{C}^{-1} \f$ is the elastic compliance (6×6) and \f$ \mathbb{C}^{(ep)} \f$ the elastoplastic
       * stiffness (6×6).
       * @param CelInv Elastic compliance tensor \f$ \mathbb{C}^{-1} \f$ (6×6 matrix).
       * @param Cep Elastoplastic stiffness tensor \f$ \mathbb{C}^{(ep)} \f$ (6×6 matrix).
       * @return 6×6 matrix \f$ \partial \boldsymbol{\varepsilon}^{p} / \partial \boldsymbol{\varepsilon} \f$ in Voigt
       * form.
       */
      Matrix6d dEp_dE( const Matrix6d& CelInv, const Matrix6d& Cep );

      /**
       * @brief Computes the derivative \f$ \frac{\partial \Delta \varepsilon^{p,vol}}{\partial
       * \boldsymbol{\varepsilon}} \f$ of the volumetric plastic strain increment \f$ \Delta \varepsilon^{p,vol} \f$
       * with respect to the Voigt-notated total strain vector \f$ \boldsymbol{\varepsilon} \f$.
       *
       * @param CelInv Elastic compliance tensor \f$ \mathbb{C}^{-1} \f$ (6×6 matrix).
       * @param Cep Elastoplastic stiffness tensor \f$ \mathbb{C}^{(ep)} \f$ (6×6 matrix).
       * @return Row vector (1×6) containing \f$ \partial \Delta \varepsilon^{p,vol} / \partial \boldsymbol{\varepsilon}
       * \f$ in Voigt notation.
       */
      RowVector6d dDeltaEpv_dE( const Matrix6d& CelInv, const Matrix6d& Cep );

      /**
       * @brief Computes the derivative \f$ \frac{\partial \varepsilon_I}{\partial \boldsymbol{\varepsilon}} \f$ of the
       * principal strains \f$ \varepsilon_I \f$ with respect to the Voigt-notated strain vector \f$
       * \boldsymbol{\varepsilon} \f$.
       *
       * @param dEp Input strain vector in Voigt notation (6 components).
       * @return A 3×6 matrix where each row contains the partial derivatives of one principal strain
       * \f$ \varepsilon_1, \varepsilon_2, \varepsilon_3 \f$ with respect to the 6 Voigt strain components.
       */
      Marmot::Matrix36 dSortedStrainPrincipal_dStrain( const Marmot::Vector6d& dEp );

      /**
       * @brief Computes the derivative \f$ \frac{\partial \Delta \varepsilon^{p,vol}_{\ominus}}
       * {\partial \boldsymbol{\varepsilon}} \f$ of the volumetric plastic strain increment in compression
       * \f$ \Delta \varepsilon^{p,vol}_{\ominus} \f$ with respect to the strain vector in Voigt notation
       * \f$ \boldsymbol{\varepsilon} \f$.
       *
       * @param dEp Incremental plastic strain vector in Voigt notation (\f$\Delta\boldsymbol{\varepsilon^{p}} \f$), 6
       * components.
       * @param CelInv Elastic compliance tensor \f$ \mathbb{C}^{-1} \f$ (6×6 matrix).
       * @param Cep Elastoplastic stiffness tensor \f$ \mathbb{C}^{(ep)} \f$ (6×6 matrix).
       * @return Row vector (1×6) containing
       * \f$ \partial \Delta \varepsilon^{p,vol}_{\ominus} / \partial \boldsymbol{\varepsilon} \f$
       * in Voigt notation.
       */
      RowVector6d dDeltaEpvneg_dE( const Marmot::Vector6d& dEp, const Matrix6d& CelInv, const Matrix6d& Cep );

    } // namespace Derivatives

    namespace Transformations {

      /**
       * @brief Computes the transformation matrix \f$ T_{\varepsilon} \f$ to rotate a strain vector in Voigt notation
       * into another Cartesian coordinate system.
       *
       * @details This function modifies the stress transformation matrix \f$T_{\varepsilon}\f$ to account for the
       * engineering shear convention used for strains in Voigt notation. Specifically:
       * - The upper-right 3×3 block is scaled by 0.5
       * - The lower-left 3×3 block is scaled by 2
       *
       * This function builds the 6×6 transformation matrix from the direction cosines of the new basis:
       * \f[
       *   \boldsymbol{\varepsilon}'_{\text{voigt}} \;=\; T_{\varepsilon}\,\boldsymbol{\varepsilon}_{\text{voigt}}.
       * \f]
       * which is equivalent to the tensor transformation
       * \f[
       *   \boldsymbol{\varepsilon}' \;=\; \boldsymbol{Q}\,\boldsymbol{\varepsilon}\,\boldsymbol{Q}^{T}
       * \f]
       * @param transformedCoordinateSystem 3×3 rotation matrix \f$\boldsymbol{Q}\f$ that maps from
       * the original to the target Cartesian basis.
       * @return 6×6 transformation matrix \f$T_{\varepsilon}\f$ that transforms strain vectors in Voigt notation.
       */
      Matrix6d transformationMatrixStrainVoigt( const Matrix3d& transformedCoordinateSystem );

      /**
       * @brief Computes the transformation matrix \f$ T_{\sigma} \f$ to rotate a stress vector in Voigt notation
       * into another Cartesian coordinate system.
       * @details This function builds the 6×6 transformation matrix from the direction cosines of the new basis:
       * \f[
       *   \boldsymbol{\sigma}'_{\text{voigt}} \;=\; T_{\sigma}\,\boldsymbol{\sigma}_{\text{voigt}},
       * \f]
       * which is equivalent to the tensor transformation
       * \f[
       *   \boldsymbol{\sigma}' \;=\; \boldsymbol{Q}\,\boldsymbol{\sigma}\,\boldsymbol{Q}^{T}.
       * \f]
       * @param transformedCoordinateSystem 3×3 rotation matrix \f$\boldsymbol{Q}\f$ that maps from
       * the original to the target Cartesian basis.
       * @return 6×6 transformation matrix \f$T_{\sigma}\f$ that transforms stress vectors in Voigt notation.
       */
      Matrix6d transformationMatrixStressVoigt( const Matrix3d& transformedCoordinateSystem );

      /**
       * @brief Calculates the projection matrix for mapping Voigt stress to a traction vector.
       * @details This function returns the \f$3 \times 6\f$ projection matrix \f$\mathbf{P}^{(\mathbf{n})}\f$
       * that maps a **Voigt-notated stress vector** \f$\boldsymbol{\sigma}_{\text{voigt}} \in \mathbb{R}^6\f$
       * to the **traction vector** \f$\mathbf{t}^{(\mathbf{n})} \in \mathbb{R}^3\f$ on a plane with
       * an outward unit normal vector \f$\mathbf{n}\f$.
       * The relationship is given by:
       * \f[
       * \mathbf{t}^{(\mathbf{n})} = \mathbf{P}^{(\mathbf{n})} \, \boldsymbol{\sigma}_{\text{voigt}}.
       * \f]
       * This expression is mathematically equivalent to Cauchy's formula \f$ \boldsymbol{t}^{(n)}_i =
       * \sigma_{ji} \, n_j \f$.
       * @param normalVector The \f$3 \times 1\f$ **unit vector** \f$\mathbf{n}\f$ representing the **normal**
       * of the plane.
       * @return The \f$3 \times 6\f$ projection matrix \f$\mathbf{P}^{(\mathbf{n})}\f$.
       */
      Matrix36d projectVoigtStressToPlane( const Vector3d& normalVector );

      /**
       * @brief Calculates the projection matrix for mapping a Voigt strain vector to the strain on a plane.
       * @details This function is similar to `projectVoigtStressToPlane` but for **strains** in voigt notation. For
       * more details, see the documentation of `projectVoigtStressToPlane`.
       * @param normalVector The \f$3 \times 1\f$ **unit vector** \f$\mathbf{n}\f$ representing the **normal**
       * of the plane.
       * @return The \f$3 \times 6\f$ projection matrix \f$\mathbf{P}^{(\mathbf{n})}\f$.
       */
      Matrix36d projectVoigtStrainToPlane( const Vector3d& normalVector );

      /**
       * @brief Rotates a stress tensor in Voigt notation.
       * @details This function rotates a stress tensor \f$ \boldsymbol{\sigma} \f$ using a rotation matrix
       * \f$ \boldsymbol{Q} \f$, while keeping the tensor in **Voigt notation**.
       * The rotated stress \f$ \boldsymbol{\sigma}^{\prime} \f$ is computed as:
       * \f[
       * \boldsymbol{\sigma}^{\prime} = \boldsymbol{Q} \, \boldsymbol{\sigma} \, \boldsymbol{Q}^{T}.
       * \f]
       * @param Q A \f$3 \times 3\f$ rotation matrix \f$ \boldsymbol{Q} \f$.
       * @param stress A Voigt-notated stress vector \f$ \boldsymbol{\sigma} \in \mathbb{R}^6 \f$.
       * @return The rotated Voigt-notated stress vector \f$ \boldsymbol{\sigma}^{\prime} \in \mathbb{R}^6 \f$.
       */
      Marmot::Vector6d rotateVoigtStress( const Eigen::Matrix3d& Q, const Marmot::Vector6d& stress );

      /**
       * @brief Transforms a stress tensor in Voigt notation from the global to a local coordinate system.
       * @param stress The stress tensor in Voigt notation (\f$ \mathbb{R}^6 \f$) in the global system.
       * @param transformedCoordinateSystem A \f$3 \times 3\f$ orthonormal matrix whose columns are the axes of the
       * local coordinate system expressed in the global system. Typically, this is a rotation matrix.
       * @return The stress tensor in Voigt notation (\f$ \mathbb{R}^6 \f$) in the local system.
       * @details The transformation convention used is: \f$ \boldsymbol{\sigma}_{\text{local}} = Q^T \,
       * \boldsymbol{\sigma}_{\text{global}} \, Q \f$, where \f$ Q \f$ is the transformation matrix.
       */
      Marmot::Vector6d transformStressToLocalSystem( const Marmot::Vector6d& stress,
                                                     const Matrix3d&         transformedCoordinateSystem );

      /**
       * @brief Transforms a strain tensor in Voigt notation from the global to a local coordinate system.
       * @param stress The strain tensor in Voigt notation (\f$ \mathbb{R}^6 \f$) in the global system.
       * @param transformedCoordinateSystem A \f$3 \times 3\f$ orthonormal matrix whose columns are the axes of the
       * local coordinate system expressed in the global system. Typically, this is a rotation matrix.
       * @return The strain tensor in Voigt notation (\f$ \mathbb{R}^6 \f$) in the local system.
       * @details The transformation convention used is: \f$ \boldsymbol{\varepsilon}_{\text{local}} = Q^T \,
       * \boldsymbol{\varepsilon}_{\text{global}} \, Q \f$, where \f$ Q \f$ is the transformation matrix.
       */
      Marmot::Vector6d transformStrainToLocalSystem( const Marmot::Vector6d& stress,
                                                     const Matrix3d&         transformedCoordinateSystem );

      /**
       * @brief Transforms a stress tensor in Voigt notation from a local to the global coordinate system.
       * @param stress The stress tensor in Voigt notation (\f$ \mathbb{R}^6 \f$) in the local system.
       * @param transformedCoordinateSystem A \f$3 \times 3\f$ orthonormal matrix whose columns are the axes of the
       * local coordinate system expressed in the global system. Typically, this is a rotation matrix.
       * @return The stress tensor in Voigt notation (\f$ \mathbb{R}^6 \f$) in the global system.
       * @details The transformation convention used is: \f$ \boldsymbol{\sigma}_{\text{global}} = Q \,
       * \boldsymbol{\sigma}_{\text{local}} \, Q^T \f$, where \f$ Q \f$ is the transformation matrix.
       */
      Marmot::Vector6d transformStressToGlobalSystem( const Marmot::Vector6d& stress,
                                                      const Matrix3d&         transformedCoordinateSystem );

      /**
       * @brief Transforms a strain tensor in Voigt notation from a local to the global coordinate system.
       * @param stress The strain tensor in Voigt notation (\f$ \mathbb{R}^6 \f$) in the local system.
       * @param transformedCoordinateSystem A \f$3 \times 3\f$ orthonormal matrix whose columns are the axes of the
       * local coordinate system expressed in the global system. Typically, this is a rotation matrix.
       * @return The strain tensor in Voigt notation (\f$ \mathbb{R}^6 \f$) in the global system.
       * @details The transformation convention used is: \f$ \boldsymbol{\varepsilon}_{\text{global}} = Q \,
       * \boldsymbol{\varepsilon}_{\text{local}} \, Q^T \f$, where \f$ Q \f$ is the transformation matrix.
       */
      Marmot::Vector6d transformStrainToGlobalSystem( const Marmot::Vector6d& stress,
                                                      const Matrix3d&         transformedCoordinateSystem );

      /**
       * @brief Transforms a stiffness tensor in Voigt notation from a local to the global coordinate system.
       * @param stiffness The stiffness tensor in Voigt notation (\f$ \mathbb{R}^{6 \times 6} \f$) in the local system.
       * @param transformedCoordinateSystem A \f$3 \times 3\f$ orthonormal matrix whose columns are the axes of the
       * local coordinate system expressed in the global system. Typically, this is a rotation matrix.
       * @return The stiffness tensor in Voigt notation (\f$ \mathbb{R}^{6 \times 6} \f$) in the global system.
       * @details The transformation convention used is: \f$ \mathbf{C}_{\text{global}} = T \, \mathbf{C}_{\text{local}}
       * \, T^T \f$, where \f$ T \f$ is the transformation matrix derived from \f$ Q \f$.
       */
      Marmot::Matrix6d transformStiffnessToGlobalSystem( const Marmot::Matrix6d& stiffness,
                                                         const Matrix3d&         transformedCoordinateSystem );

    } // namespace Transformations

  }   // namespace ContinuumMechanics::VoigtNotation
} // namespace Marmot
