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
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include <Eigen/QR>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <utility>

namespace Marmot::ContinuumMechanics {

  /**
   * @brief Reduced stress states by static condensation: some stress components vanish, the strains at them follow.
   *
   * @details In Marmot's Voigt order \f$[11, 22, 33, 12, 13, 23]\f$ (engineering shear strains), a reduced stress
   * state is defined by the prescribed strain components \f$p\f$ and the stress-free components \f$f\f$: the strains
   * \f$\boldsymbol{\varepsilon}_p\f$ are given, \f$\boldsymbol{\sigma}_f = \mathbf{0}\f$ determines the strains
   * \f$\boldsymbol{\varepsilon}_f\f$, and the strains at the components in neither set keep their given values
   * (usually zero). With the 3D tangent partitioned into the blocks \f$\mathbf{C}_{pp}, \mathbf{C}_{pf},
   * \mathbf{C}_{fp}, \mathbf{C}_{ff}\f$, linearizing \f$\boldsymbol{\sigma}_f = \mathbf{0}\f$ gives the sensitivity
   * of the free strains and the condensed tangent:
   * \f[
   *   \frac{\partial \boldsymbol{\varepsilon}_f}{\partial \boldsymbol{\varepsilon}_p}
   *     = -\mathbf{C}_{ff}^{-1}\, \mathbf{C}_{fp}, \qquad
   *   \frac{\partial \boldsymbol{\sigma}_p}{\partial \boldsymbol{\varepsilon}_p}
   *     = \mathbf{C}_{pp} - \mathbf{C}_{pf}\, \mathbf{C}_{ff}^{-1}\, \mathbf{C}_{fp}.
   * \f]
   * Any other derivative is condensed the same way: e.g. a quantity \f$\mathbf{y}(\boldsymbol{\varepsilon},
   * \mathbf{K})\f$ with an additional argument \f$\mathbf{K}\f$ has \f$\mathrm{d}\boldsymbol{\varepsilon}_f /
   * \mathrm{d}\mathbf{K} = -\mathbf{C}_{ff}^{-1}\, \partial\boldsymbol{\sigma}_f / \partial \mathbf{K}\f$ and
   * \f$\mathrm{d}\mathbf{y} / \mathrm{d}\boldsymbol{\varepsilon}_p = \partial\mathbf{y} /
   * \partial\boldsymbol{\varepsilon}_p + \partial\mathbf{y} / \partial\boldsymbol{\varepsilon}_f\, \partial
   * \boldsymbol{\varepsilon}_f / \partial \boldsymbol{\varepsilon}_p\f$.
   *
   * A singular \f$\mathbf{C}_{ff}\f$ is admitted if the reduced problem stays consistent, e.g. a fully softened
   * smeared crack without out-of-plane shear stiffness: then \f$\mathbf{C}_{ff}^{-1}\f$ above is the pseudo-inverse
   * \f$\mathbf{C}_{ff}^{+}\f$, so the strains along the null space of \f$\mathbf{C}_{ff}\f$ are not changed (minimum
   * norm corrections), which is exact as long as neither \f$\boldsymbol{\sigma}_p\f$ nor \f$\boldsymbol{\sigma}_f\f$
   * depend on them. Otherwise, the reduced stress state is not determined, and StressUpdateFailed is thrown.
   *
   * computeReducedStress() is the iteration of a 3D stress update to such a state, condense() the algebra of a given
   * tangent. The standard states are defined below: plane stress (planeStressPrescribed, planeStressFree) and uniaxial
   * stress (uniaxialStressPrescribed, uniaxialStressFree).
   */
  namespace ReducedStress {

    /// @brief Plane stress, prescribed components: \f$\varepsilon_{11}, \varepsilon_{22}, \gamma_{12}\f$ (in the
    /// order of the plane Voigt notation).
    inline constexpr std::array< int, 3 > planeStressPrescribed = { 0, 1, 3 };
    /// @brief Plane stress, stress-free components: \f$\sigma_{33} = \sigma_{13} = \sigma_{23} = 0\f$.
    inline constexpr std::array< int, 3 > planeStressFree = { 2, 4, 5 };
    /// @brief Uniaxial stress, prescribed component: \f$\varepsilon_{11}\f$.
    inline constexpr std::array< int, 1 > uniaxialStressPrescribed = { 0 };
    /// @brief Uniaxial stress, stress-free components: all but \f$\sigma_{11}\f$.
    inline constexpr std::array< int, 5 > uniaxialStressFree = { 1, 2, 3, 4, 5 };

    /// @brief Convergence tolerance of computeReducedStress(): \f$\max_i |\sigma_{f_i}|\f$ relative to the stress
    /// scale.
    inline constexpr double convergenceTolerance = 1e-12;
    /// @brief Relaxed convergence tolerance of computeReducedStress(), from the 9th evaluation on.
    inline constexpr double relaxedConvergenceTolerance = 1e-9;
    /// @brief Rank threshold of \f$\mathbf{C}_{ff}\f$: a pivot of its complete orthogonal decomposition below this
    /// fraction of the largest pivot counts as zero.
    inline constexpr double rankTolerance = 1e-12;
    /// @brief Consistency of a singular \f$\mathbf{C}_{ff}\f$ (condense()): the largest component of
    /// \f$\mathbf{C}_{pf}\, \mathbf{N}\f$ and of \f$(\mathbf{I} - \mathbf{C}_{ff} \mathbf{C}_{ff}^{+})\,
    /// \mathbf{C}_{fp}\f$, with \f$\mathbf{N}\f$ the projector on the null space of \f$\mathbf{C}_{ff}\f$, relative to
    /// \f$\max_{ij} |C_{ij}|\f$; and of the stress residual outside the range of \f$\mathbf{C}_{ff}\f$
    /// (computeReducedStress()) relative to the stress scale.
    inline constexpr double consistencyTolerance = 1e-9;

    /**
     * @brief The static condensation of a 3D tangent to a reduced stress state.
     * @tparam nP number of prescribed components
     * @tparam nF number of stress-free components
     */
    template < int nP, int nF >
    struct Condensation {
      /// @brief The condensed tangent \f$\mathbf{C}_{pp} - \mathbf{C}_{pf} \mathbf{C}_{ff}^{+} \mathbf{C}_{fp}\f$
      /// (\f$\mathbf{C}_{ff}^{+} = \mathbf{C}_{ff}^{-1}\f$ if it is regular).
      Eigen::Matrix< double, nP, nP > dStressP_dStrainP;
      /// @brief The sensitivity of the free strains \f$\partial \boldsymbol{\varepsilon}_f / \partial
      /// \boldsymbol{\varepsilon}_p = -\mathbf{C}_{ff}^{+} \mathbf{C}_{fp}\f$.
      Eigen::Matrix< double, nF, nP > dStrainF_dStrainP;
      /// @brief The complete orthogonal decomposition of \f$\mathbf{C}_{ff}\f$ (threshold rankTolerance), to condense
      /// further derivatives with the same pseudo-inverse, see freeStrainSensitivity().
      Eigen::CompleteOrthogonalDecomposition< Eigen::Matrix< double, nF, nF > > CffCOD;
      /// @brief computeReducedStress(): the number of evaluations of the 3D stress update (1 if the initial guess
      /// satisfies the stress-free conditions); 0 from condense().
      int nIterations = 0;
      /// @brief The tangent of the stress-free components \f$\mathbf{C}_{ff}\f$, decomposed in CffCOD.
      Eigen::Matrix< double, nF, nF > Cff;
    };

    /**
     * @brief Checks the index sets of a reduced stress state.
     * @tparam nP_ number of prescribed components
     * @tparam nF_ number of stress-free components
     * @param p the prescribed components
     * @param f the stress-free components
     * @throws std::invalid_argument unless @p p and @p f are distinct Voigt indices in [0, 6)
     */
    template < std::size_t nP_, std::size_t nF_ >
    void checkIndexSets( const std::array< int, nP_ >& p, const std::array< int, nF_ >& f )
    {
      static_assert( nP_ > 0 && nF_ > 0 && nP_ + nF_ <= 6, "ReducedStress: invalid number of components" );
      std::array< bool, 6 > used{};
      auto                  use = [&]( int i ) {
        if ( i < 0 || i > 5 || used[i] )
          throw std::invalid_argument( Marmot::Utilities::MakeString()
                                       << __PRETTY_FUNCTION__
                                       << ": the prescribed and the free components must be distinct Voigt indices "
                                                           "in [0, 6)" );
        used[i] = true;
      };
      std::for_each( p.begin(), p.end(), use );
      std::for_each( f.begin(), f.end(), use );
    }

    /**
     * @brief The static condensation of a 3D tangent: the condensed tangent, the sensitivity of the free strains and
     * the complete orthogonal decomposition of \f$\mathbf{C}_{ff}\f$ (see the namespace ReducedStress).
     * @tparam nP_ number of prescribed components
     * @tparam nF_ number of stress-free components
     * @param C the 3D tangent \f$\partial\boldsymbol{\sigma} / \partial\boldsymbol{\varepsilon}\f$ (Voigt, engineering
     * shear strains)
     * @param p the prescribed components, distinct Voigt indices in [0, 6)
     * @param f the stress-free components, distinct Voigt indices in [0, 6), disjoint from @p p
     * @return the condensation, a Condensation<nP, nF>, with Condensation::nIterations = 0
     * @throws std::invalid_argument if @p p and @p f are no valid index sets
     * @throws Marmot::Utilities::StressUpdateFailed if @p C is not finite, or if \f$\mathbf{C}_{ff}\f$ is singular
     * (rankTolerance) and the condensation is inconsistent (consistencyTolerance): \f$\boldsymbol{\sigma}_p\f$ depends
     * on strains along the null space of \f$\mathbf{C}_{ff}\f$, or \f$\mathbf{C}_{fp}\f$ is not in its range. The
     * reduced stress state is then not determined by the material; this happens within a stress update, so the
     * caller can cut back. A consistent singular \f$\mathbf{C}_{ff}\f$ is condensed with its pseudo-inverse.
     */
    template < std::size_t nP_, std::size_t nF_ >
    auto condense( const Marmot::Math::Matrix6d& C, const std::array< int, nP_ >& p, const std::array< int, nF_ >& f )
    {
      constexpr int nP = static_cast< int >( nP_ );
      constexpr int nF = static_cast< int >( nF_ );
      checkIndexSets( p, f );

      Eigen::Matrix< double, nF, nF > Cff;
      Eigen::Matrix< double, nF, nP > Cfp;
      Eigen::Matrix< double, nP, nF > Cpf;
      Eigen::Matrix< double, nP, nP > Cpp;
      for ( int i = 0; i < nF; i++ ) {
        for ( int j = 0; j < nF; j++ )
          Cff( i, j ) = C( f[i], f[j] );
        for ( int j = 0; j < nP; j++ )
          Cfp( i, j ) = C( f[i], p[j] );
      }
      for ( int i = 0; i < nP; i++ ) {
        for ( int j = 0; j < nF; j++ )
          Cpf( i, j ) = C( p[i], f[j] );
        for ( int j = 0; j < nP; j++ )
          Cpp( i, j ) = C( p[i], p[j] );
      }

      if ( !C.allFinite() )
        throw Marmot::Utilities::StressUpdateFailed( Marmot::Utilities::MakeString()
                                                     << __PRETTY_FUNCTION__ << ": the tangent is not finite" );

      Condensation< nP, nF > condensation;
      condensation.CffCOD.setThreshold( rankTolerance );
      condensation.CffCOD.compute( Cff );
      condensation.Cff = Cff;

      if ( condensation.CffCOD.rank() < nF ) {
        // consistency of the singular block: sigma_p independent of the null space of C_ff, C_fp in its range
        using MatrixFF                    = Eigen::Matrix< double, nF, nF >;
        const MatrixFF nullSpaceProjector = MatrixFF::Identity() - condensation.CffCOD.solve( Cff );
        const MatrixFF complementOfRange  = MatrixFF::Identity() - Cff * condensation.CffCOD.pseudoInverse();
        const double   scale              = C.cwiseAbs().maxCoeff();
        const double   pDependsOnNull     = ( Cpf * nullSpaceProjector ).cwiseAbs().maxCoeff();
        const double   fpOutsideRange     = ( complementOfRange * Cfp ).cwiseAbs().maxCoeff();
        if ( pDependsOnNull > consistencyTolerance * scale || fpOutsideRange > consistencyTolerance * scale )
          throw Marmot::Utilities::StressUpdateFailed(
            Marmot::Utilities::MakeString()
            << __PRETTY_FUNCTION__ << ": the tangent of the stress-free components is singular (rank "
            << condensation.CffCOD.rank() << " of " << nF << ") and the condensation is inconsistent" );
      }

      condensation.dStrainF_dStrainP = -condensation.CffCOD.solve( Cfp );
      condensation.dStressP_dStrainP = Cpp + Cpf * condensation.dStrainF_dStrainP;
      return condensation;
    }

    /**
     * @brief A reduced stress state of a 3D stress update: Newton's method on the strains at the stress-free
     * components (see the namespace ReducedStress).
     *
     * @details Each iteration evaluates @p computeStress3D at the current strain increment and corrects the free
     * strains by \f$\Delta\boldsymbol{\varepsilon}_f \leftarrow \Delta\boldsymbol{\varepsilon}_f -
     * \mathbf{C}_{ff}^{+} \boldsymbol{\sigma}_f\f$ (minimum norm, see condense()), until \f$\max_i |\sigma_{f_i}|
     * \leq 10^{-12}\, s\f$ (convergenceTolerance), or \f$10^{-9}\, s\f$ (relaxedConvergenceTolerance) from the 9th
     * evaluation on. The stress scale \f$s\f$ is the largest stress component
     * \f$|\sigma_i|\f$ (all six) of all evaluations of this call, so the tolerance is relative and unit free, and it is
     * met by a stress state of round-off size; if all stresses vanish, the residual does too. After 15 evaluations
     * without convergence, it throws.
     *
     * The last evaluation of @p computeStress3D is the one at the converged strain, so additional outputs of the
     * callable (state variables, energies, other derivatives) refer to the converged 3D state.
     *
     * @tparam nP_ number of prescribed components
     * @tparam nF_ number of stress-free components
     * @tparam StressUpdate3D callable `void( const Marmot::Math::Vector6d& dStrain, Marmot::Math::Vector6d& stress,
     * Marmot::Math::Matrix6d& dStress_dStrain )`: the 3D stress update for the strain increment `dStrain` from the
     * start of the increment, i.e., every call must reset what the previous call updated (state variables, the stress
     * at the start of the increment, energies). `stress` and `dStress_dStrain` are zeroed before each call.
     * @param computeStress3D the 3D stress update
     * @param[in,out] dStrain in: the prescribed strain increments at @p p, an initial guess at @p f, and the given
     * increments (usually zero) at the other components; out: the converged increment
     * @param p the prescribed components
     * @param f the stress-free components
     * @return the condensation of the tangent of the converged 3D state, a Condensation<nP, nF>, with the number of
     * evaluations
     * @throws std::invalid_argument if @p p and @p f are no valid index sets (checkIndexSets(), before the first
     * evaluation)
     * @throws Marmot::Utilities::StressUpdateFailed if the iteration does not converge, a stress is not finite, the
     * condensation is inconsistent (condense()), or \f$\mathbf{C}_{ff}\f$ is singular and the residual
     * \f$\boldsymbol{\sigma}_f\f$ has a component outside its range larger than consistencyTolerance \f$s\f$, which no
     * correction of the free strains can remove; the stress update can then be repeated with a smaller increment.
     */
    template < std::size_t nP_, std::size_t nF_, typename StressUpdate3D >
    auto computeReducedStress( StressUpdate3D&&              computeStress3D,
                               Marmot::Math::Vector6d&       dStrain,
                               const std::array< int, nP_ >& p,
                               const std::array< int, nF_ >& f )
    {
      constexpr int    nF               = static_cast< int >( nF_ );
      constexpr double tolerance        = convergenceTolerance;
      constexpr double relaxedTolerance = relaxedConvergenceTolerance;
      constexpr int    nRelaxed         = 8;
      constexpr int    nMaxIterations   = 15;
      checkIndexSets( p, f ); // before the first access at f

      Marmot::Math::Vector6d         stress;
      Marmot::Math::Matrix6d         dStress_dStrain;
      Eigen::Matrix< double, nF, 1 > stressF;
      double                         stressScale = 0;
      for ( int iteration = 1;; iteration++ ) {
        stress.setZero();
        dStress_dStrain.setZero();
        computeStress3D( std::as_const( dStrain ), stress, dStress_dStrain );

        if ( !stress.allFinite() )
          throw Marmot::Utilities::StressUpdateFailed( Marmot::Utilities::MakeString()
                                                       << __PRETTY_FUNCTION__ << ": the stress is not finite" );

        stressScale = std::max( stressScale, stress.cwiseAbs().maxCoeff() );
        for ( int i = 0; i < nF; i++ )
          stressF( i ) = stress( f[i] );
        const double residual = stressF.cwiseAbs().maxCoeff();

        auto condensation = condense( dStress_dStrain, p, f );

        if ( residual <= tolerance * stressScale ||
             ( iteration > nRelaxed && residual <= relaxedTolerance * stressScale ) ) {
          condensation.nIterations = iteration;
          return condensation;
        }

        if ( iteration >= nMaxIterations ) {
          throw Marmot::Utilities::StressUpdateFailed(
            Marmot::Utilities::MakeString()
            << __PRETTY_FUNCTION__ << ": no convergence of the " << nF << " stress-free components after "
            << nMaxIterations << " iterations, max |sigma_f| = " << residual << ", stress scale " << stressScale );
        }

        const Eigen::Matrix< double, nF, 1 > correction = condensation.CffCOD.solve( stressF );
        if ( condensation.CffCOD.rank() < nF ) {
          const double outsideRange = ( stressF - condensation.Cff * correction ).cwiseAbs().maxCoeff();
          if ( outsideRange > consistencyTolerance * stressScale )
            throw Marmot::Utilities::StressUpdateFailed(
              Marmot::Utilities::MakeString()
              << __PRETTY_FUNCTION__ << ": the tangent of the stress-free components is singular (rank "
              << condensation.CffCOD.rank() << " of " << nF << ") and the stress residual " << outsideRange
              << " is outside its range" );
        }
        for ( int i = 0; i < nF; i++ )
          dStrain( f[i] ) -= correction( i );
      }
    }

    /**
     * @brief The sensitivity of the free strains to a further variable \f$\mathbf{x}\f$ of the stress update (e.g.,
     * nonlocal variables), with the same pseudo-inverse as the condensation: \f$\mathrm{d}\boldsymbol{\varepsilon}_f /
     * \mathrm{d}\mathbf{x} = -\mathbf{C}_{ff}^{+}\, \partial\boldsymbol{\sigma}_f / \partial\mathbf{x}\f$, such that
     * the free stresses stay zero.
     * @tparam nP number of prescribed components
     * @tparam nF number of stress-free components
     * @tparam nX number of further variables
     * @param condensation the condensation of the tangent (condense(), computeReducedStress())
     * @param dStressF_dX \f$\partial\boldsymbol{\sigma}_f / \partial\mathbf{x}\f$
     * @return \f$\mathrm{d}\boldsymbol{\varepsilon}_f / \mathrm{d}\mathbf{x}\f$
     * @throws Marmot::Utilities::StressUpdateFailed if \f$\mathbf{C}_{ff}\f$ is singular and a column of
     * @p dStressF_dX has a component outside its range larger than consistencyTolerance times the column's scale
     * (\f$\max(|\partial\sigma_{f_i}/\partial x_j|, \max|\mathbf{C}_{ff}|\, \max_i|\mathrm{d}\varepsilon_{f_i} /
     * \mathrm{d}x_j|)\f$): no change of the free strains then keeps the free stresses zero when \f$\mathbf{x}\f$
     * changes.
     */
    template < int nP, int nF, int nX >
    Eigen::Matrix< double, nF, nX > freeStrainSensitivity( const Condensation< nP, nF >&          condensation,
                                                           const Eigen::Matrix< double, nF, nX >& dStressF_dX )
    {
      const Eigen::Matrix< double, nF, nX > dStrainF_dX = -condensation.CffCOD.solve( dStressF_dX );
      if ( condensation.CffCOD.rank() < nF ) {
        const Eigen::Matrix< double, nF, nX > residual = dStressF_dX + condensation.Cff * dStrainF_dX;
        const double                          CffMax   = condensation.Cff.cwiseAbs().maxCoeff();
        for ( int j = 0; j < nX; j++ ) {
          const double scale = std::max( dStressF_dX.col( j ).cwiseAbs().maxCoeff(),
                                         CffMax * dStrainF_dX.col( j ).cwiseAbs().maxCoeff() );
          if ( residual.col( j ).cwiseAbs().maxCoeff() > consistencyTolerance * scale )
            throw Marmot::Utilities::StressUpdateFailed(
              Marmot::Utilities::MakeString()
              << __PRETTY_FUNCTION__ << ": the tangent of the stress-free components is singular (rank "
              << condensation.CffCOD.rank() << " of " << nF << ") and dStressF_dX (column " << j
              << ") is outside its range" );
        }
      }
      return dStrainF_dX;
    }

  } // namespace ReducedStress

  namespace PlaneStrain {

    /**
     * @brief Extracts the plane strain stiffness matrix from a 3D stiffness matrix.
     *
     * @details Reduces the full 3D stiffness matrix \f$\mathbb{C}\f$ to the plane strain
     * stiffness matrix corresponding to the in-plane components.
     *
     * @param C 3D stiffness matrix (\f$6\times6\f$) in Voigt notation.
     * @return Plane strain stiffness matrix (\f$3\times3\f$) in Voigt notation.
     */
    Eigen::Matrix3d getPlaneStrainTangent( const Marmot::Math::Matrix6d& C );

    /**
     * @brief Extracts the plane strain derivative of the stress with respect to the deformation gradient.
     *
     * @details Reduces the derivative of the 3D stress tensor in Voigt notation
     * with respect to the deformation gradient \f$F_{ij}\f$ to the plane strain setting.
     *
     * @param dStressdDeformationGradient3D Derivative of 3D stress with respect to the deformation gradient
     *                                      (\f$6\times3\times3\f$ tensor).
     * @return Reduced derivative of stress with respect to deformation gradient (\f$3\times2\times2\f$ tensor)
     *         in the plane strain setting.
     */
    Marmot::Math::Tensors::Tensor322d reduce3D_dStress_dDeformationGradient(
      const Marmot::Math::Tensors::Tensor633d& dStressdDeformationGradient3D );

    /**
     * @brief Computes the derivative of the 3D strain tensor with respect to the plane strain tensor.
     *
     * @return Matrix (\f$6\times3\f$) representing \f$\frac{\partial \varepsilon_{3D}}{\partial
     * \varepsilon_{plane}}\f$.
     */
    Eigen::Matrix< double, 6, 3 > dStrainDStrainPlaneStrain();
  } // namespace PlaneStrain

  namespace PlaneStress {

    /**
     * @brief Computes the plane stress derivative of the stress with respect to the deformation gradient.
     *
     * @details Reduces the 3D derivative of the stress tensor with respect to the deformation gradient
     * \f$ F_{ij} \f$ to a plane stress setting.
     *
     * @param dStressdDeformationGradient3D Derivative of 3D stress w.r.t. deformation gradient
     *                                      (\f$6 \times 3 \times 3\f$ tensor).
     * @return Reduced derivative of stress w.r.t. deformation gradient (\f$3 \times 2 \times 2\f$ tensor)
     *         for plane stress.
     */
    Marmot::Math::Tensors::Tensor322d compute_dStress_dDeformationGradient(
      const Marmot::Math::Tensors::Tensor633d& dStressdDeformationGradient3D );

    /**
     * @brief Computes the out-of-plane strain correction component for plane stress conditions.
     *
     * @details Determines the correction \f$\varepsilon_{33}^{\text{corr}}\f$ such that the resulting
     * stress satisfies plane stress conditions.
     * \f$ \boldsymbol{\varepsilon}^{\text{el}} \f$, so that the resulting stress satisfies plane stress:
     * \f[
     *     \sigma^\text{plane} = \mathbb{C} : (\boldsymbol{\varepsilon}^{\text{el}} +
     * \boldsymbol{\varepsilon}^{\text{corr}})
     * \f]
     *
     * @param elasticStrain The in-plane elastic strain vector \f$ \boldsymbol{\varepsilon}^{\text{el}} \f$ (6
     * components, Voigt notation).
     * @param nu Poisson's ratio of the material.
     * @return A 6-component strain vector \f$ \boldsymbol{\varepsilon}^{\text{corr}} \f$ containing the computed
     * out-of-plane correction component \f$\varepsilon_{33}^{\text{corr}}\f$ (Voigt notation).
     */
    Marmot::Math::Vector6d planeStressCompensationStrain( const Marmot::Math::Vector6d& elasticStrain, double nu );

    /**
     * @brief Computes the transformation matrix for the plane stress tangent.
     *
     * @details Returns the 6x6 transformation matrix \f$\mathbf{T}\f$ relating arbitrary strain increments
     * to plane stress increments using the material stiffness.
     *
     * @param tangent 6x6 material tangent stiffness matrix (\f$\mathbf{C}\f$).
     * @return 6x6 transformation matrix (\f$\mathbf{T}\f$) for plane stress correction.
     */
    Marmot::Math::Matrix6d planeStressTangentTransformationMatrix( const Marmot::Math::Matrix6d& tangent );

    /**
     * @brief Computes the derivative of the 3D strain tensor with respect to the plane stress strain tensor.
     *
     * @param tangent 6x6 material tangent stiffness matrix (\f$\mathbf{C}\f$).
     * @return 6x3 matrix representing \f$\frac{\partial \varepsilon_{3D}}{\partial \varepsilon_{\sigma_{plane}}}\f$.
     */
    Eigen::Matrix< double, 6, 3 > dStrainDStrainPlaneStress( const Marmot::Math::Matrix6d& tangent );

    /**
     * @brief Computes the derivative of the plane stress tensor with respect to the 3D stress tensor.
     *
     * @return 3x6 matrix representing \f$\frac{\partial \sigma_{plane}}{\partial \sigma_{3D}}\f$.
     */
    Eigen::Matrix< double, 3, 6 > dStressPlaneStressDStress();
  } // namespace PlaneStress
} // namespace Marmot::ContinuumMechanics
