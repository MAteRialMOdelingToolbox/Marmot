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
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include <utility>

namespace Marmot {
  namespace ContinuumMechanics::CommonTensors {
    /**
     * @brief Initializes the fourth-order tensor \f$I_{ijkl} = \delta_{ij}\delta_{kl}\f$.
     */
    EigenTensors::Tensor3333d Initialize_I2xI2();

    /**
     * @brief Fourth-order tensor \f$I_{ijkl} = \delta_{ij}\delta_{kl}\f$.
     */
    inline const EigenTensors::Tensor3333d I2xI2 = Initialize_I2xI2();

    /**
     * @brief Initializes the symmetric fourth-order identity tensor \f$
     * I_{ijkl}^{sym}=\frac{1}{2}(\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk}) \f$.
     */
    EigenTensors::Tensor3333d Initialize_Isym();

    /**
     * @brief Symmetric fourth-order identity tensor
     * \f$ I_{ijkl}^{sym}=\frac{1}{2}(\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk}) \f$.
     */
    inline const EigenTensors::Tensor3333d Isym = Initialize_Isym();

    /**
     * @brief Initializes the skew-symmetric fourth-order identity tensor
     * \f$ I_{ijkl}^{skew}=\frac{1}{2}(\delta_{ik}\delta_{jl}-\delta_{il}\delta_{jk}) \f$.
     */
    EigenTensors::Tensor3333d Initialize_Iskew();

    /** @brief Skew-symmetric part of the fourth-order identity tensor \f$
     * I_{ijkl}^{skew}=\frac{1}{2}(\delta_{ik}\delta_{jl}-\delta_{il}\delta_{jk}) \f$.
     */
    inline const EigenTensors::Tensor3333d Iskew = Initialize_Iskew();

    /**
     * @brief Initializes the fourth-order identity tensor \f$ I_{ijkl} = \delta_{ik}\delta_{jl} \f$.
     */
    EigenTensors::Tensor3333d Initialize_IFourthOrder();

    /**
     * @brief Fourth-order identity tensor \f$ I_{ijkl} = \delta_{ik}\delta_{jl} \f$.
     */
    inline const EigenTensors::Tensor3333d IFourthOrder = Initialize_IFourthOrder();

    /**
     * @brief Initializes the transposed fourth-order identity tensor \f$ I_{ijkl}^{T} = \delta_{il}\delta_{jk} \f$.
     */
    EigenTensors::Tensor3333d Initialize_IFourthOrderTranspose();

    /// @brief Transposed fourth-order identity tensor \f$ I_{ijkl}^{T} = \delta_{il}\delta_{jk} \f$.
    inline const EigenTensors::Tensor3333d IFourthOrderTranspose = Initialize_IFourthOrderTranspose();

    /**
     * @brief Initializes the derivative tensor of deviatoric stress w.r.t. stress  \f$ \frac{\partial
     * s_{ij}}{\partial\sigma_{kl}} = \delta_{ik}\delta_{jl} - \frac{1}{3} \delta_{ij}\delta_{kl} \f$.
     */
    EigenTensors::Tensor3333d Initialize_dDeviatoricStress_dStress();

    /** @brief Derivative of the deviatoric stress with respect to stress \f$ \frac{\partial
     * s_{ij}}{\partial\sigma_{kl}} =
     * \delta_{ik}\delta_{jl} - \frac{1}{3} \delta_{ij}\delta_{kl} \f$.
     */
    inline const EigenTensors::Tensor3333d dDeviatoricStress_dStress = Initialize_dDeviatoricStress_dStress();

    /**
     * @brief Initializes the 3D Levi-Civita permutation tensor \f$E_{ijk}\f$.
     * @details A fully antisymmetric third-order tensor defined as
     * \f[
     *   E_{ijk} =
     *   \begin{cases}
     *     +1 & \text{if } (i,j,k) \text{ is an even permutation of } (1,2,3), \\
     *     -1 & \text{if } (i,j,k) \text{ is an odd permutation of } (1,2,3), \\
     *      0 & \text{if any two indices are equal.}
     *   \end{cases}
     * \f]
     */
    EigenTensors::Tensor333d Initialize_LeviCivita3D();

    /**
     * @brief 3D Levi-Civita permutation tensor \f$E_{ijk}\f$.
     * @copydetails Initialize_LeviCivita3D
     */
    inline const EigenTensors::Tensor333d LeviCivita3D = Initialize_LeviCivita3D();

    /**
     * @brief 2D Levi-Civita permutation tensor \f$E_{ij}\f$.
     * @details A fully antisymmetric second-order tensor defined as
     * \f[
     *   E_{ij} =
     *   \begin{cases}
     *     +1 & \text{if } (i,j) = (1,2), \\
     *     -1 & \text{if } (i,j) = (2,1), \\
     *      0 & \text{if } i=j.
     *   \end{cases}
     * \f]
     * Commonly used to represent 2D cross products and rotations in tensor notation.
     */
    EigenTensors::Tensor122d Initialize_LeviCivita2D();

    /**
     * @brief 2D Levi-Civita permutation tensor \f$\varepsilon_{ij}\f$.
     * @copydetails Initialize_LeviCivita2D
     */
    inline const EigenTensors::Tensor122d LeviCivita2D = Initialize_LeviCivita2D();

    /**
     * @brief Initializes the second-order identity tensor.
     * @details \f$ I_{ij} = \delta_{ij} \f$.
     */
    EigenTensors::Tensor33d Initialize_I2();

    /** @brief Second-order identity tensor.
     * @copydetails Initialize_I2
     */
    inline const EigenTensors::Tensor33d I2 = Initialize_I2();

    /**
     * @brief Returns the number of rotational DOFs for the given dimension.
     * @param nDim Problem dimension (2 or 3).
     * @return 1 in 2D, 3 in 3D.
     */
    constexpr int getNumberOfDofForRotation( int nDim )
    {
      if ( nDim == 2 )
        return 1;
      else
        return 3;
    }

    /**
     * @brief Constructs an aux. matrix, which helps to swap indices in Eigen::Matrices abused as higher order Tensors
     * by multiplication
     *
     * The transformation is defined as:
     * \f[
     *   T_{(ij)(kl)} \cdot P_{(kl)(lk)} = T_{(ij)(lk)}
     * \f]
     *
     * @tparam sizeI Number of rows in the tensor index block.
     * @tparam sizeJ Number of columns in the tensor index block.
     * @return Eigen::Matrix<double, sizeI * sizeJ, sizeI * sizeJ>
     *         The index swap tensor.
     *
     * @code
     * Eigen::Matrix3d t;
     * t << 1,2,3,
     *      4,5,6,
     *      7,8,9;
     * const auto P = makeIndexSwapTensor<3,3>();
     * std::cout << t.reshaped().transpose() << std::endl;
     * std::cout << t.reshaped().transpose() * P << std::endl;
     *
     * // Output:
     * // 1,4,7,2,5,8,3,6,9
     * // 1,2,3,4,5,6,7,8,9
     * @endcode
     */
    template < int sizeI, int sizeJ >
    Eigen::Matrix< double, sizeI * sizeJ, sizeI * sizeJ > makeIndexSwapTensor()
    {
      Eigen::Matrix< double, sizeI * sizeJ, sizeI * sizeJ > P;

      auto d = []( int a, int b ) -> double { return a == b ? 1.0 : 0.0; };

      for ( int i = 0; i < sizeI; i++ )
        for ( int j = 0; j < sizeJ; j++ )
          for ( int l = 0; l < sizeI; l++ )
            for ( int k = 0; k < sizeJ; k++ )
              P( i + j * sizeI, l * sizeJ + k ) = d( i, l ) * d( k, j );

      return P;
    }

    /**
     * @brief Provides the reference Levi-Civita permutation tensor for the given spatial dimension.
     *
     * @tparam nDim Spatial dimension (2 or 3).
     * @return The Levi-Civita tensor:
     * - In 2D: \f$E_{ij}\f$, a second-order antisymmetric tensor with components
     * @copydetails Initialize_LeviCivita2D
     * - In 3D: \f$E_{ijk}\f$, a third-order antisymmetric tensor with components
     * @copydetails Initialize_LeviCivita3D
     */
    template < int nDim >
    constexpr Eigen::TensorFixedSize< double, Eigen::Sizes< getNumberOfDofForRotation( nDim ), nDim, nDim > > getReferenceToCorrectLeviCivita()
    // template <int nDim>
    // template <int nDim=2>
    // constexpr Eigen::TensorFixedSize< double, Eigen::Sizes<getNumberOfDofForRotation (nDim), nDim, nDim> >
    // getReferenceToCorrectLeviCivita()
    {
      if constexpr ( nDim == 2 )
        return LeviCivita2D;
      else
        return LeviCivita3D;
    }

  } // namespace ContinuumMechanics::CommonTensors

  namespace ContinuumMechanics::TensorUtility {
    /** @brief Kronecker delta function \f$ \delta_{ab} \f$.
     * @return 1 if a == b, otherwise 0.
     */
    constexpr int d( int a, int b )
    {
      return a == b ? 1 : 0;
    }

    /**
     * @brief Map an object's raw data as a fixed-size Eigen matrix.
     *
     * @details Creates an `Eigen::Map` view of @p t reinterpreted as an `x`-by-`y` matrix. No data is copied; the map
     * directly references the storage returned by `t.data()`.
     *
     * @tparam x Number of rows (compile-time).
     * @tparam y Number of columns (compile-time).
     * @tparam T Object type, must define `Scalar` and provide `data()`.
     * @param t Input object.
     * @return `Eigen::Map<Eigen::Matrix<typename T::Scalar, x, y>>`.
     *
     * @note Caller must ensure `t.data()` has at least `x*y` elements
     *       and matches Eigen’s default (column-major) layout.
     */
    template < int x,
               int y,
               typename T,
               typename = std::enable_if< !std::is_const< std::remove_reference< T > >::value > >
    auto as( T& t )
    {
      return Eigen::Map< Eigen::Matrix< typename T::Scalar, x, y > >( t.data() );
    }

    /**
     * @brief Map an object's raw data as a fixed-size Eigen matrix (const version).
     * @copydetails as
     */
    template < int x, int y, typename T, typename = void >
    auto as( const T& t )
    {
      return Eigen::Map< const Eigen::Matrix< typename T::Scalar, x, y > >( t.data() );
    }

    /**
     * @brief Flattens an Eigen object into a 1D column vector map.
     *
     * Creates an `Eigen::Map` that reinterprets the data of @p t as a
     * column vector of size `RowsAtCompileTime * ColsAtCompileTime`.
     * No copy is made; modifications through the map affect @p t.
     *
     * @tparam Derived Eigen matrix/array type (non-const).
     * @param t Eigen object to be flattened.
     * @return `Eigen::Map<Eigen::Matrix<Scalar, Rows*Cols, 1>>`
     *         referencing the data of @p t.
     */
    template < typename Derived,
               typename = std::enable_if< !std::is_const< std::remove_reference< Derived > >::value > >
    auto flatten( Derived& t )
    {
      return Eigen::Map<
        Eigen::Matrix< typename Derived::Scalar, Derived::RowsAtCompileTime * Derived::ColsAtCompileTime, 1 > >(
        t.data() );
    }

    /**
     * @brief Flattens an Eigen object into a 1D column vector map (const version).
     */
    template < typename Derived, typename = void >
    auto flatten( const Derived& t )
    {
      return Eigen::Map<
        const Eigen::Matrix< typename Derived::Scalar, Derived::RowsAtCompileTime * Derived::ColsAtCompileTime, 1 > >(
        t.data() );
    }

    /**
     * @brief Build contraction dimension pairs for Eigen tensor operations.
     *
     * Generates an `Eigen::array<Eigen::IndexPair<int>, N>` from the
     * compile-time parameter pack @p Pairs, where each consecutive pair
     * defines a contraction index mapping. Evaluated entirely at compile time.
     *
     * @tparam Pairs Sequence of integers; must contain an even number of values
     *         forming index pairs.
     * @return Array of index pairs usable in Eigen tensor contraction.
     */
    template < int... Pairs >
    constexpr auto contractionDims() // should be evaluated at compile time
    {
      static_assert( sizeof...( Pairs ) % 2 == 0, "Pairs must contain an even number of elements." );

      constexpr int numPairs = sizeof...( Pairs ) / 2;

      Eigen::array< Eigen::IndexPair< int >, numPairs > result{};

      if constexpr ( numPairs > 0 ) {          // return empty array if no pairs are given
        constexpr int values[] = { Pairs... }; // Pairs... expands the parameter pack

        // Fill the result array at compile-time
        for ( int i = 0; i < numPairs; ++i ) {
          result[i] = { values[2 * i], values[2 * i + 1] };
        }
      }
      return result;
    }

    /**
     * @brief Compute the dyadic (outer) product of two 3D vectors.
     *
     * @details Forms a 3x3 matrix where each entry is given by \f$ a_{ij} = b_i c_j \f$.
     *
     * @param vector1 First 3D vector.
     * @param vector2 Second 3D vector.
     * @return 3x3 matrix representing the dyadic product.
     */
    Eigen::Matrix3d dyadicProduct( const Eigen::Vector3d& vector1, const Eigen::Vector3d& vector2 );

    namespace IndexNotation {
      /**
       * @brief Convert a Voigt index to tensor indices.
       *
       * Maps a Voigt notation index @p ij to the corresponding
       * `(i, j)` tensor indices for a given dimension @p nDim.
       *
       * @tparam nDim Problem dimension (1, 2, or 3).
       * @param ij Voigt index.
       * @return Pair of tensor indices (i, j).
       * @throws std::invalid_argument if @p nDim or @p ij is invalid.
       */
      template < int nDim >
      constexpr std::pair< int, int > fromVoigt( int ij )
      {
        if constexpr ( nDim == 1 )
          return std::pair< int, int >( 0, 0 );
        else if ( nDim == 2 )
          switch ( ij ) {
          case 0: return std::pair< int, int >( 0, 0 );
          case 1: return std::pair< int, int >( 1, 1 );
          case 2: return std::pair< int, int >( 0, 1 );
          }

        else if ( nDim == 3 ) {
          switch ( ij ) {
          case 0: return std::pair< int, int >( 0, 0 );
          case 1: return std::pair< int, int >( 1, 1 );
          case 2: return std::pair< int, int >( 2, 2 );
          case 3: return std::pair< int, int >( 0, 1 );
          case 4: return std::pair< int, int >( 0, 2 );
          case 5: return std::pair< int, int >( 1, 2 );
          }
        }

        throw std::invalid_argument( MakeString()
                                     << __PRETTY_FUNCTION__ << ": invalid dimension / voigt index specified" );
      }

      /**
       * @brief Maps tensor indices (i, j) to the corresponding Voigt
       * notation index for a given dimension @p nDim.
       *
       * @tparam nDim Problem dimension (1, 2, or 3).
       * @param i Row index of the tensor.
       * @param j Column index of the tensor.
       * @return Voigt index corresponding to (i, j).
       * @throws std::invalid_argument if @p nDim is invalid.
       */
      template < int nDim >
      constexpr int toVoigt( int i, int j )
      {
        if constexpr ( nDim == 1 )
          return 0;
        else if ( nDim == 2 )
          return ( i == j ) ? ( i == 0 ? 0 : 1 ) : 2;

        else if ( nDim == 3 ) {
          constexpr int tensor2VoigtNotationIndicesMapping[3][3] = { { 0, 3, 4 }, { 3, 1, 5 }, { 4, 5, 2 } };
          return tensor2VoigtNotationIndicesMapping[i][j];
        }

        throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid dimension specified" );
      }

      /**
       * @brief Construct the Voigt mapping tensor.
       *
       * Creates a 3rd-order tensor that maps tensor indices (i, j)
       * to their corresponding Voigt index. Each entry is 1 at
       * `(toVoigt<nDim>(i, j), i, j)` and 0 elsewhere.
       *
       * @tparam nDim Problem dimension (1, 2, or 3).
       * @return A tensor of shape (VoigtSize, nDim, nDim) encoding
       *         the Voigt mapping.
       */
      template < int nDim >
      Eigen::TensorFixedSize< double, Eigen::Sizes< VOIGTFROMDIM( nDim ), nDim, nDim > > voigtMap()
      {
        using namespace Eigen;
        Eigen::TensorFixedSize< double, Eigen::Sizes< VOIGTFROMDIM( nDim ), nDim, nDim > > result;
        result.setZero();
        for ( int i = 0; i < nDim; i++ )
          for ( int j = 0; j < nDim; j++ )
            result( toVoigt< nDim >( i, j ), i, j ) = 1;
        return result;
      }

    } // namespace IndexNotation

    // namespace ContinuumMechanics::VoigtNotation

  } // namespace ContinuumMechanics::TensorUtility

} // namespace Marmot
