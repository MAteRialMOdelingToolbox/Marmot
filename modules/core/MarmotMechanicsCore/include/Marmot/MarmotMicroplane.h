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
 * Konstantin Basche konstantin.basche@uibk.ac.at
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

#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotTypedefs.h"

using namespace Marmot::ContinuumMechanics::CommonTensors;
using namespace Marmot::ContinuumMechanics::TensorUtility;
using namespace Eigen;

namespace Marmot {
  /**
   * @namespace Marmot::Microplane
   * @brief This namespace provides the hard-coded microplane quadrature schemes (orientations and
   * integration weights on the unit hemisphere) and the two kinematic constraints built on top of
   * them: the volumetric-deviatoric (V-D) split (VDSplit) and the normal-tangential (N-T) split
   * (NTSplit). Both classes construct their projection tensors from the same discretization
   * (Discretization::get()) and provide the corresponding virtual-work homogenization back to the
   * macroscopic stress/tangent.
   */
  namespace Microplane {

    /**
     * @brief The available hard-coded microplane quadrature schemes (integration point
     * orientations and weights on the unit hemisphere) according to Bažant (1986).
     * @details #N21Symmetric and #N21Asymmetric both discretize the hemisphere into 21
     * microplanes but differ in whether the orientations exhibit orthogonal symmetries (Bažant
     * 1986, Table 1 vs. Table 4); #N28, #N37 and #N61 use 28, 37 and 61 microplanes, respectively.
     * See #nqpFromScheme for the corresponding number of microplanes and Discretization::get() for
     * the hard-coded orientations/weights.
     */
    enum class QuadratureScheme { N21Symmetric, N21Asymmetric, N28, N37, N61 };

    /**
     * @brief Returns the number of microplanes \f$ n_{qp} \f$ (integration points on the unit
     * hemisphere) for a given #QuadratureScheme.
     * @param scheme The quadrature scheme.
     * @return Number of microplanes \f$ n_{qp} \in \{21, 28, 37, 61\} \f$.
     */
    constexpr size_t nqpFromScheme( QuadratureScheme scheme )
    {
      return scheme == QuadratureScheme::N28   ? 28
             : scheme == QuadratureScheme::N37 ? 37
             : scheme == QuadratureScheme::N61 ? 61
                                               : 21; // N21Symmetric / N21Asymmetric
    }

    namespace Discretization {
      /**
       * @brief Hard-coded microplane orientations and integration weights for one
       * #QuadratureScheme, as returned by #get(). Dynamically sized since it is shared by
       * #VDSplit and #NTSplit instantiations for different quadrature schemes; the number of rows
       * is given by #nqpFromScheme.
       */
      struct Data {
        /// @brief Microplane unit normal vectors \f$ n^{(m)} \f$, one per row (\f$ n_{qp} \times 3 \f$).
        Eigen::MatrixXd normals;
        /// @brief Integration weights \f$ w^{(m)} \f$, one per microplane (\f$ n_{qp} \times 1 \f$), normalized
        /// such that \f$ \sum_m w^{(m)} = 3 \f$.
        Eigen::VectorXd weights;
      };

      /**
       * @brief Returns the hard-coded microplane orientations and integration weights for @p
       * scheme.
       * @details Single source for the discretization data shared by #VDSplit and
       * #NTSplit, so both splits sample exactly the same set of microplane orientations
       * for a given #QuadratureScheme.
       * @param scheme The quadrature scheme.
       * @return The microplane normal vectors and integration weights; see #Data.
       */
      Data get( QuadratureScheme scheme );
    } // namespace Discretization

    /**
     * @brief Normal-tangential (N-T) split and projection tensors of the microplane
     * theory, as used by the M7-family models (MicroplaneM7, GM7).
     * @details On each microplane, the macroscopic strain is projected onto the normal direction
     * \f$ n^{(m)} \f$ and the two in-plane orthogonal shear directions \f$ l^{(m)} \f$, \f$ m^{(m)} \f$
     * (spanning the plane orthogonal to \f$ n^{(m)} \f$), giving a normal strain \f$ \varepsilon_N \f$
     * and two shear strains \f$ \varepsilon_M, \varepsilon_L \f$. #pTensors holds the corresponding
     * Voigt-notation projection tensors \f$ N_{ij} = n_i n_j \f$,
     * \f$ M_{ij} = \tfrac{1}{2}(m_i n_j + m_j n_i) \f$, \f$ L_{ij} = \tfrac{1}{2}(l_i n_j + l_j n_i) \f$,
     * built from the hard-coded discretization (Discretization::get()), and #homogenizeStress /
     * #homogenizeTangent perform the virtual-work homogenization of per-microplane normal/shear
     * stresses and tangents back to the macroscopic stress and tangent.
     * @tparam scheme The #QuadratureScheme (and thus the number of microplanes #nqp) this
     * instantiation is built for.
     */
    template < QuadratureScheme scheme >
    class NTSplit {
    public:
      /// @brief Number of microplanes for scheme; see #nqpFromScheme.
      static constexpr size_t nqp = nqpFromScheme( scheme );

      /// @brief One scalar per microplane, e.g. #projectionTensors::w.
      typedef Eigen::Matrix< double, 1, nqp > nmp_sized;
      /// @brief One Voigt-notation vector per microplane, e.g. #projectionTensors::N.
      typedef Eigen::Matrix< double, 6, nqp > nmp6_sized;

      NTSplit() : pTensors( initiateSystem() ){};

      /**
       * @brief Kinematic projection tensors (normal \f$ N \f$, and in-plane shear \f$ M \f$/\f$ L \f$
       * directions) and integration weights \f$ w \f$ for all #nqp microplanes, in Voigt notation.
       */
      struct projectionTensors {
        /**
         * @brief Normal projection tensor \f$ N_{ij}^{(m)} = n_i^{(m)} n_j^{(m)} \f$ in Voigt notation, one
         * column per microplane.
         * @details Projects the macroscopic strain tensor \f$ \eps \f$ onto the microplane normal strain
         * \f$ \varepsilon_N^{(m)} \f$ of microplane @p m:
         * \f[
         *   \varepsilon_N^{(m)} = N_{ij}^{(m)} \varepsilon_{ij}
         * \f]
         */
        nmp6_sized N;
        /**
         * @brief In-plane shear projection tensor \f$ M_{ij}^{(m)} = \tfrac{1}{2}(m_i^{(m)} n_j^{(m)} +
         * m_j^{(m)} n_i^{(m)}) \f$ in Voigt notation, one column per microplane.
         * @details Projects the macroscopic strain tensor \f$ \eps \f$ onto the microplane shear strain
         * \f$ \varepsilon_M^{(m)} \f$ (direction \f$ m^{(m)} \f$) of microplane @p m:
         * \f[
         *   \varepsilon_M^{(m)} = M_{ij}^{(m)} \varepsilon_{ij}
         * \f]
         */
        nmp6_sized M;
        /**
         * @brief In-plane shear projection tensor \f$ L_{ij}^{(m)} = \tfrac{1}{2}(l_i^{(m)} n_j^{(m)} +
         * l_j^{(m)} n_i^{(m)}) \f$ in Voigt notation, one column per microplane.
         * @details Projects the macroscopic strain tensor \f$ \eps \f$ onto the microplane shear strain
         * \f$ \varepsilon_L^{(m)} \f$ (direction \f$ l^{(m)} \f$) of microplane @p m:
         * \f[
         *   \varepsilon_L^{(m)} = L_{ij}^{(m)} \varepsilon_{ij}
         * \f]
         */
        nmp6_sized L;
        /// @brief Integration weights \f$ w^{(m)} \f$ for the chosen scheme; see Discretization::Data::weights.
        nmp_sized w;
      };
      /// @brief The N/M/L projection tensors and integration weights for this instance; see #projectionTensors.
      projectionTensors pTensors;

      /**
       * @brief Homogenizes per-microplane normal and in-plane shear stresses to the macroscopic
       * stress, in Voigt notation.
       * @details Approximates the continuum microplane virtual-work homogenization over the unit
       * hemisphere:
       * \f[
       *   \sigma_{ij} = \frac{3}{2\pi} \int_\Omega \left( \sigma_N N_{ij} + \sigma_L L_{ij} + \sigma_M M_{ij} \right)
       * \, d\Omega
       * \f]
       * by the quadrature sum actually evaluated here, in which the hemisphere-integral
       * normalization constant is already absorbed into the integration weights \f$ w^{(m)} \f$
       * (see Discretization::get(), where \f$ \sum_m w^{(m)} = 3 \f$):
       * \f[
       *   \sigma_{ij} \approx \sum_{m=1}^{n_{qp}} w^{(m)} \left( \sigma_N^{(m)} N_{ij}^{(m)}
       *     + \sigma_M^{(m)} M_{ij}^{(m)} + \sigma_L^{(m)} L_{ij}^{(m)} \right)
       * \f]
       * @tparam T Scalar type of the per-microplane stresses, e.g. `double` or
       * `autodiff::dual`; the projection tensors #pTensors stay `double` (they are
       * geometric constants, never differentiated), so `T` may differ from `double`.
       * @param stressN Per-microplane normal stress \f$ \sigma_N^{(m)} \f$, size \f$ 1 \times n_{qp} \f$.
       * @param stressM Per-microplane shear stress \f$ \sigma_M^{(m)} \f$ (direction \f$ M \f$), size \f$ 1 \times
       * n_{qp} \f$.
       * @param stressL Per-microplane shear stress \f$ \sigma_L^{(m)} \f$ (direction \f$ L \f$), size \f$ 1 \times
       * n_{qp} \f$.
       * @return Homogenized macroscopic stress in Voigt notation.
       */
      template < typename T >
      Eigen::Matrix< T, 6, 1 > homogenizeStress( const Eigen::Matrix< T, 1, nqp >& stressN,
                                                 const Eigen::Matrix< T, 1, nqp >& stressM,
                                                 const Eigen::Matrix< T, 1, nqp >& stressL ) const;

      /**
       * @brief Homogenizes per-microplane normal and in-plane shear tangents to the macroscopic
       * tangent, in Voigt notation.
       * @details Approximates the continuum microplane virtual-work homogenization over the unit
       * sphere \f$ \Omega \f$ of the macroscopic tangent, in tensor (index) notation:
       * \f[
       *   \frac{\partial \sigma_{ij}}{\partial \varepsilon_{kl}} = \int_\Omega \left(
       *       N_{ij} \frac{\partial \sigma_N}{\partial \varepsilon_{kl}}
       *     + M_{ij} \frac{\partial \sigma_M}{\partial \varepsilon_{kl}}
       *     + L_{ij} \frac{\partial \sigma_L}{\partial \varepsilon_{kl}} \right) \, d\Omega
       * \f]
       * approximated by the quadrature sum actually evaluated here (see #homogenizeStress for how
       * the hemisphere/sphere-integral normalization is absorbed into the weights \f$ w^{(m)} \f$).
       * @tparam T Scalar type of the per-microplane tangents (see #homogenizeStress).
       * @param dSigNddStrain Per-microplane normal-direction tangent \f$ \partial \sigma_N^{(m)} /
       * \partial \boldsymbol{\varepsilon} \f$, size \f$ 6 \times n_{qp} \f$; taken by value since
       * each row is scaled in place by the integration weight before contraction.
       * @param dSigMddStrain Per-microplane \f$ M \f$-direction tangent \f$ \partial \sigma_M^{(m)}
       * / \partial \boldsymbol{\varepsilon} \f$, size \f$ 6 \times n_{qp} \f$.
       * @param dSigLddStrain Per-microplane \f$ L \f$-direction tangent \f$ \partial \sigma_L^{(m)}
       * / \partial \boldsymbol{\varepsilon} \f$, size \f$ 6 \times n_{qp} \f$.
       * @return Homogenized macroscopic tangent \f$ \partial \sigma_{ij} / \partial
       * \varepsilon_{kl} \f$ in Voigt notation.
       */
      template < typename T >
      Eigen::Matrix< T, 6, 6 > homogenizeTangent( Eigen::Matrix< T, 6, nqp > dSigNddStrain,
                                                  Eigen::Matrix< T, 6, nqp > dSigMddStrain,
                                                  Eigen::Matrix< T, 6, nqp > dSigLddStrain ) const;

    private:
      /**
       * @brief Builds the microplane normal vectors, the fixed in-plane shear basis, and the N/M/L
       * projection tensors and integration weights (#projectionTensors) from
       * Discretization::get(scheme).
       */
      projectionTensors initiateSystem();
    };

    template < QuadratureScheme scheme >
    template < typename T >
    Eigen::Matrix< T, 6, 1 > NTSplit< scheme >::homogenizeStress( const Eigen::Matrix< T, 1, nqp >& stressN,
                                                                  const Eigen::Matrix< T, 1, nqp >& stressM,
                                                                  const Eigen::Matrix< T, 1, nqp >& stressL ) const
    {
      return pTensors.N * ( stressN.array() * pTensors.w.array() ).matrix().transpose() +
             pTensors.M * ( stressM.array() * pTensors.w.array() ).matrix().transpose() +
             pTensors.L * ( stressL.array() * pTensors.w.array() ).matrix().transpose();
    }

    template < QuadratureScheme scheme >
    template < typename T >
    Eigen::Matrix< T, 6, 6 > NTSplit< scheme >::homogenizeTangent( Eigen::Matrix< T, 6, nqp > dSigNddStrain,
                                                                   Eigen::Matrix< T, 6, nqp > dSigMddStrain,
                                                                   Eigen::Matrix< T, 6, nqp > dSigLddStrain ) const
    {
      for ( size_t t = 0; t < 6; t++ ) {
        dSigNddStrain.row( t ) = dSigNddStrain.row( t ).array() * pTensors.w.array();
        dSigMddStrain.row( t ) = dSigMddStrain.row( t ).array() * pTensors.w.array();
        dSigLddStrain.row( t ) = dSigLddStrain.row( t ).array() * pTensors.w.array();
      }
      return pTensors.N * dSigNddStrain.transpose() + pTensors.M * dSigMddStrain.transpose() +
             pTensors.L * dSigLddStrain.transpose();
    }

    template < QuadratureScheme scheme >
    NTSplit< scheme >::projectionTensors NTSplit< scheme >::initiateSystem()
    {
      const Discretization::Data data = Discretization::get( scheme );

      projectionTensors pTensors;
      pTensors.w = data.weights;

      Matrix< int, 6, 2 > ij;
      ij << 1, 1, 2, 2, 3, 3, 1, 2, 1, 3, 2, 3;

      // Fixed in-plane basis vector spanning the M/L directions together with each microplane
      // normal - not actually "random", just a vector not parallel to any of the hard-coded normals.
      const Vector3d rand_vec( -0.46, 0.29, 0.82 );
      Vector3d       n_vec, m_vec, l_vec;

      for ( size_t qp_i = 0; qp_i < nqp; qp_i++ ) {
        n_vec = data.normals.row( qp_i );

        m_vec                = rand_vec - n_vec.dot( rand_vec ) * n_vec;
        const double lengthm = m_vec.norm();
        m_vec /= lengthm;

        // calculate l_vec as cross product of n_vec x m_vec
        l_vec = n_vec.cross( m_vec );
        l_vec /= l_vec.norm();

        // calculate projection Matrices N_ij, M_ij, L_ij for kinematic constraint
        for ( size_t k = 0; k < 6; k++ ) {
          int i                 = ij( k, 0 );
          int j                 = ij( k, 1 );
          pTensors.N( k, qp_i ) = n_vec( i - 1 ) * n_vec( j - 1 ); // 6 x nqp matrix
          pTensors.L( k, qp_i ) = 0.5 * ( l_vec( i - 1 ) * n_vec( j - 1 ) + l_vec( j - 1 ) * n_vec( i - 1 ) );
          pTensors.M( k, qp_i ) = 0.5 * ( m_vec( i - 1 ) * n_vec( j - 1 ) + m_vec( j - 1 ) * n_vec( i - 1 ) );
        }
      }
      return pTensors;
    }

    /**
     * @brief Volumetric-deviatoric (V-D) split and projection tensors of the microplane
     * theory, as used e.g. by damage-plasticity constitutive microplane models.
     * @details The macroscopic strain is decomposed on each microplane into a volumetric part
     * \f$ \varepsilon_V = \tfrac{1}{3}\operatorname{tr}(\eps) \f$, identical
     * for every microplane, and a deviatoric microplane strain vector \f$ \eps_D \f$. #pTensors holds the corresponding
     * projection tensors (#projectionTensors::Vol, #projectionTensors::Dev) built from the hard-coded discretization
     * (Discretization::get()), and the `homogenize*` methods perform the virtual-work homogenization of per-microplane
     * volumetric/deviatoric quantities back to macroscopic 3x3 tensors.
     * @tparam scheme The #QuadratureScheme (and thus the number of microplanes #nqp) this
     * instantiation is built for.
     */
    template < QuadratureScheme scheme >
    class VDSplit {
    public:
      /// @brief Number of microplanes for scheme; see #nqpFromScheme.
      static constexpr size_t nqp = nqpFromScheme( scheme );

      VDSplit() : pTensors( initiateSystem() ){};

      /**
       * @brief Volumetric/deviatoric projection tensors and integration weights for all #nqp
       * microplanes.
       */
      struct projectionTensors {
        /**
         * @brief Holds the deviatoric projection tensors for each microplane.
         * @details Projects the macroscopic strain tensor \f$ \eps \f$ onto the microplane deviatoric strain vector \f$
         * \eps_D^{(m)} \f$ of microplane @p m:
         * \f[
         *   (\varepsilon_D^{(m)})_i = Dev_{ijk}^{(m)} \varepsilon_{jk}
         * \f]
         * Built from the microplane normal vectors #nMatrix in #initiateSystem(). For all \f$ n_{qp} \f$ microplanes,
         * these are stored together in a single 4th-order tensor of size \f$ n_{qp} \times 3 \times 3 \times 3 \f$,
         * where the first index selects the microplane and \f$ Dev^{(m)} \f$ above denotes the corresponding
         * 3rd-order slice.
         */
        TensorFixedSize< double, Sizes< nqp, 3, 3, 3 > > Dev;
        /**
         * @brief Volumetric projection tensor \f$ V_{ij} = \delta_{ij}/3 \f$, tensor form of #Vol.
         * @details Projects the macroscopic strain tensor \f$ \eps \f$ onto the (microplane-independent)
         * volumetric strain \f$ \varepsilon_V \f$:
         * \f[
         *   \varepsilon_V = V_{ij} \varepsilon_{ij} = \tfrac{1}{3} \operatorname{tr}( \eps )
         * \f]
         */
        TensorFixedSize< double, Sizes< 3, 3 > > VolT;
        /**
         * @brief Volumetric projection tensor \f$ V_{ij} = \delta_{ij}/3 \f$, matrix form; see #VolT for the
         * kinematic projection \f$ \varepsilon_V = V_{ij} \varepsilon_{ij} \f$.
         */
        Matrix< double, 3, 3 > Vol;
        /// @brief Integration weights \f$ w^{(m)} \f$ for the chosen scheme; see Discretization::Data::weights.
        Matrix< double, nqp, 1 > w;
        /// @brief Microplane normal vectors \f$ n^{(m)} \f$, one per row; see Discretization::Data::normals.
        Matrix< double, nqp, 3 > nMatrix;
      };
      /// @brief The volumetric/deviatoric projection tensors and integration weights for this instance; see
      /// #projectionTensors.
      projectionTensors pTensors;

      /**
       * @brief Homogenizes per-microplane volumetric and deviatoric values to the macroscopic
       * tensor.
       * @details Evaluates the discrete virtual-work homogenization:
       * \f[
       *   X_{jk} = \left( \sum_{m=1}^{n_{qp}} w^{(m)} X_V^{(m)} \right) \frac{\delta_{jk}}{3}
       *     + \sum_{m=1}^{n_{qp}} w^{(m)}\, Dev_{ijk}^{(m)}\, \big( X_D^{(m)} \big)_i
       * \f]
       * Denoting the \f$ i \f$-th component of \f$ X_D^{(m)} \f$ by \f$ \big( X_D^{(m)} \big)_i \f$ (same index
       * convention as the kinematic projection \f$ (\varepsilon_D^{(m)})_i = Dev_{ijk}^{(m)} \varepsilon_{jk} \f$
       * documented at #projectionTensors::Dev, contracted here on \f$ i \f$ instead of \f$ j,k \f$). \f$ Dev^{(m)}
       * \f$ is #projectionTensors::Dev for microplane \f$ m \f$. This is used both to homogenize stress (\f$
       * \boldsymbol{X} = \sig \f$) and to homogenize strain (\f$ \boldsymbol{X} = \eps \f$) - the method itself is
       * agnostic to the physical meaning of @p volValue / @p devValue.
       * @param volValue Per-microplane volumetric value \f$ X_V^{(m)} \f$, size \f$ n_{qp} \times 1 \f$.
       * @param devValue Per-microplane deviatoric value \f$ X_D^{(m)} \f$, size \f$ n_{qp} \times 3 \f$.
       * @return Homogenized macroscopic tensor \f$ X_{jk} \f$ (3x3).
       */
      Matrix3d homogenize( const auto& volValue, const auto& devValue ) const;

      /**
       * @brief Applies the projection tensors of a single microplane @p m to a volumetric/deviatoric
       * value pair, without homogenizing over all microplanes.
       * @details Computes, the (unweighted) macroscopic-space contribution of microplane @p m alone:
       * \f[
       *   X_{jk} = X_V^{(m)} \, V_{jk} + Dev_{ijk}^{(m)}\, \big( X_D^{(m)} \big)_i
       * \f]
       * denoting the \f$ i \f$-th component of \f$ X_D^{(m)} \f$ by \f$ \big( X_D^{(m)} \big)_i \f$ (same index
       * convention as #homogenize and the kinematic projection documented at #projectionTensors::Dev).
       * @param volValue Volumetric value \f$ X_V^{(m)} \f$ of microplane @p m.
       * @param devValue Deviatoric value \f$ X_D^{(m)} \f$ of microplane @p m.
       * @param m Index of the microplane.
       * @return The 3x3 tensor contribution \f$ X_{jk} \f$ of microplane @p m.
       */
      Matrix3d applyProjectionTensors( const double volValue, const Vector3d& devValue, const int m ) const;

    private:
      /**
       * @brief Builds the microplane normal vectors, integration weights, and volumetric/deviatoric
       * projection tensors (#projectionTensors) from Discretization::get(scheme).
       */
      projectionTensors initiateSystem();
    };

    template < QuadratureScheme scheme >
    Matrix3d VDSplit< scheme >::homogenize( const auto& volValue, const auto& devValue ) const
    {
      Matrix3d vol_3x3 = pTensors.Vol * pTensors.w.dot( volValue );

      Matrix< double, nqp, 3 > devValue_weighted = ( devValue.array().colwise() * pTensors.w.array() ).matrix();
      TensorFixedSize< double, Sizes< 3, 3 > > dev_3x3T;
      TensorMap< Tensor< double, 2 > >         devValue_weightedT( devValue_weighted.data(), nqp, 3 );
      dev_3x3T = pTensors.Dev.contract( devValue_weightedT, contractionDims< 0, 0, 1, 1 >() );
      Map< Matrix< double, 3, 3 > > dev_3x3( dev_3x3T.data() );

      return vol_3x3 + dev_3x3;
    }

    template < QuadratureScheme scheme > // V * sigv + Dev^T * sigd
    Matrix3d VDSplit< scheme >::applyProjectionTensors( const double    volValue,
                                                        const Vector3d& devValue,
                                                        const int       m ) const
    { // without damage
      TensorMap< const Tensor< double, 1 > > devValueT( devValue.data(), 3 );
      TensorFixedSize< double, Sizes< 3, 3 > >
        deviatoricPart = pTensors.Dev.chip( m, 0 ).eval().contract( devValueT, contractionDims< 0, 0 >() );

      // Map Tensor back to Matrix
      Map< const Matrix< double, 3, 3 > > deviatoricPartM( deviatoricPart.data(), 3, 3 );

      return pTensors.Vol * volValue + deviatoricPartM;
    }

    template < QuadratureScheme scheme >
    VDSplit< scheme >::projectionTensors VDSplit< scheme >::initiateSystem()
    {
      const Discretization::Data data = Discretization::get( scheme );

      projectionTensors pTensors;
      pTensors.nMatrix = data.normals;
      pTensors.w       = data.weights;

      // Deviatoric projection tensor
      TensorFixedSize< double, Sizes< 3, 3 > > kronDelta_33;
      for ( int i = 0; i < 3; ++i ) {
        for ( int j = 0; j < 3; ++j ) {
          kronDelta_33( i, j ) = ( i == j );
        }
      }

      TensorFixedSize< double, Sizes< 3, 3, 3, 3 > > delta_ijDelta_kl;
      delta_ijDelta_kl = kronDelta_33.contract( kronDelta_33,
                                                contractionDims<>() ); // outer product: delta_ij delta_kl

      // map the tensor on matrix to use Eigen's contraction
      TensorMap< Tensor< double, 2 > > normVec_tensor( pTensors.nMatrix.data(), nqp, 3 );

      // Deviatoric projection tensor nqp x 3 x 3 x 3
      pTensors.Dev = normVec_tensor.contract( Isym, contractionDims< 1, 0 >() ) -
                     1. / 3. * normVec_tensor.contract( delta_ijDelta_kl, contractionDims< 1, 0 >() );

      // Volumetric projection tensor
      pTensors.Vol = 1. / 3. * Matrix3d::Identity();
      pTensors.VolT.setZero();           // same as Vol but of type Tensor
      for ( int i = 0; i < 3; ++i ) {
        pTensors.VolT( i, i ) = 1. / 3.; // Set diagonal elements to 1
      }

      return pTensors;
    }

  } // namespace Microplane
} // namespace Marmot
