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
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <stdexcept>

/**
 * @brief The critical time step of the central difference method for a continuum element, from the largest
 * eigenfrequency of the element.
 *
 * @details The central difference method with a lumped mass \f$\mathbf{M}\f$ is stable for \f$\Delta t \leq 2 /
 * \omega_{\max}\f$, with the largest eigenfrequency \f$\omega_{\max}\f$ of the assembled system, and the largest
 * eigenfrequency of the system is bounded by the largest of its elements: \f$\omega_{\max}^2 \leq \max_e
 * \lambda_{\max}(\mathbf{M}_e^{-1} \mathbf{K}_e)\f$. An element's own eigenvalue is therefore a safe estimate, for
 * every polynomial order, every quadrature rule and every distortion, since it is computed from the element's own
 * strain operator and its own lumped mass -- where a characteristic length over a wave speed is the exact limit of a
 * linear bar only.
 *
 * The stiffness of an isotropic elastic element splits along the two Lame constants,
 * \f[ \mathbf{K}_e = \lambda\, \mathbf{K}_\lambda + \mu\, \mathbf{K}_\mu, \qquad
 *     \mathbf{K}_\lambda = \int \mathbf{b}\, \mathbf{b}^T \,\mathrm{d}V, \qquad
 *     \mathbf{K}_\mu = \int \mathbf{B}^T \mathbf{D}\, \mathbf{B} \,\mathrm{d}V, \f]
 * with the strain operator \f$\mathbf{B}\f$ in Voigt notation, the divergence \f$\mathbf{b}^T\f$ (the sum of its
 * rows of the normal strains) and \f$\mathbf{D}\f$ = 2 on the normal and 1 on the engineering shear strains, so
 * that \f$\mathbf{u}^T \mathbf{K}_e \mathbf{u} = \int \lambda (\nabla \cdot \mathbf{u})^2 + 2 \mu\,
 * \boldsymbol{\varepsilon} : \boldsymbol{\varepsilon} \,\mathrm{d}V\f$. (\f$\mathbf{K}_\mu\f$ is not the deviatoric
 * operator: \f$2\,\boldsymbol{\varepsilon} : \boldsymbol{\varepsilon}\f$ contains the volumetric strain too.) A
 * material tells its largest wave speed \f$c\f$, i.e. the P-wave modulus \f$\rho c^2 = \lambda + 2\mu\f$, but not how
 * it divides into \f$\lambda\f$ and \f$\mu\f$. With \f$\lambda = \rho c^2 - 2\mu\f$, the stiffness is affine in
 * \f$\mu \in [0, \rho c^2 / 2]\f$ (\f$\nu \in [0, 1/2]\f$), and the largest eigenvalue of a matrix is convex in it,
 * so it is largest at one of the two ends: the incompressible end \f$\rho c^2\, \mathbf{K}_\lambda\f$ or the end
 * without lateral contraction \f$\rho c^2 \mathbf{K}_\mu / 2\f$. Hence
 * \f[ \omega_{\max}^2 \leq \rho c^2\, \Lambda, \qquad
 *     \Lambda = \max\Big( \lambda_{\max}(\mathbf{M}_e^{-1} \mathbf{K}_\lambda),\;
 *                         \tfrac{1}{2} \lambda_{\max}(\mathbf{M}_e^{-1} \mathbf{K}_\mu) \Big), \qquad
 *     \Delta t = \frac{2}{\sqrt{\rho c^2\, \Lambda}} . \f]
 * The bound is exact for a material at the end whose term is the larger one and conservative everywhere else, by as
 * much as a compact linear element at \f$\nu = 0\f$, where the incompressible end governs (its time increment is then
 * about 0.6 of the element's true one). It assumes an isotropic material with the same \f$\nu\f$ at all quadrature
 * points of the element, \f$\nu \geq 0\f$, and takes the largest P-wave modulus of its quadrature points.
 *
 * A plane stress section is covered by the plane operators: its in-plane modulus \f$2\lambda\mu / (\lambda +
 * 2\mu)\f$ is smaller than \f$\lambda\f$.
 */
namespace Marmot::Discretizations::Elements::CriticalTimeStep {

  /**
   * @brief The P-wave modulus \f$\rho c^2\f$ of a quadrature point, checked before it is combined with those of the
   * other points.
   * @details Every check is written so that a NaN fails it: a quadrature point whose wave speed or density is not a
   * number must stop the analysis, not drop out of the largest P-wave modulus of the element (a NaN loses every
   * comparison, so `std::max` would pass over it, and the time increment would come out too large).
   * @param density The density \f$\rho\f$ of the material at the point.
   * @param waveSpeed The largest wave speed \f$c\f$ of the material at the point.
   * @param elementLabel The label of the element, for the message.
   * @param quadraturePoint The number of the quadrature point, for the message.
   * @return The P-wave modulus \f$\rho c^2\f$.
   * @throws std::runtime_error if the wave speed or the density is not a positive number.
   */
  inline double pWaveModulusOfQuadraturePoint( double density, double waveSpeed, int elementLabel, int quadraturePoint )
  {
    if ( !( waveSpeed > 0.0 ) || !( density > 0.0 ) )
      throw std::runtime_error( Marmot::Utilities::MakeString()
                                << "CriticalTimeStep: element " << elementLabel << ", quadrature point "
                                << quadraturePoint << ": the wave speed (" << waveSpeed << ") and the density ("
                                << density << ") must be positive numbers." );
    return density * waveSpeed * waveSpeed;
  }

  /**
   * @brief The operators \f$\mathbf{K}_\lambda\f$ and \f$\mathbf{K}_\mu\f$ of an element, accumulated
   * quadrature point by quadrature point, and the critical time step they give with the element's lumped mass;
   * see the namespace documentation.
   */
  class IsotropicElasticStiffnessBound {

  public:
    /**
     * @brief Start with zero operators.
     * @param nDof The number of displacement degrees of freedom of the element, the columns of its strain operator.
     */
    explicit IsotropicElasticStiffnessBound( int nDof )
      : lambdaOperator( Eigen::MatrixXd::Zero( nDof, nDof ) ), muOperator( Eigen::MatrixXd::Zero( nDof, nDof ) )
    {
    }

    /**
     * @brief Add the contribution of a quadrature point.
     * @tparam Derived The Eigen expression type of the strain operator.
     * @param B The strain operator in Voigt notation, the normal strains first, then the engineering shear strains;
     * its columns are the displacement degrees of freedom of the element, in the order of the lumped mass.
     * @param nNormalStrains The number of rows of the normal strains: the spatial dimension, or 3 for an
     * axisymmetric section (rr, zz, \f$\theta\theta\f$).
     * @param dV The integration weight of the quadrature point as the element integrates its internal force with it.
     */
    template < typename Derived >
    void addQuadraturePoint( const Eigen::MatrixBase< Derived >& B, int nNormalStrains, double dV )
    {
      const Eigen::RowVectorXd divergence = B.topRows( nNormalStrains ).colwise().sum();
      lambdaOperator += divergence.transpose() * divergence * dV;
      muOperator += 2.0 * B.topRows( nNormalStrains ).transpose() * B.topRows( nNormalStrains ) * dV;
      const int nShearStrains = static_cast< int >( B.rows() ) - nNormalStrains;
      muOperator += B.bottomRows( nShearStrains ).transpose() * B.bottomRows( nShearStrains ) * dV;
    }

    /**
     * @brief The critical time step \f$2 / \sqrt{\rho c^2 \Lambda}\f$.
     * @param lumpedMass The lumped mass of the displacement degrees of freedom, as the element assembles it.
     * @param pWaveModulus The largest P-wave modulus \f$\rho c^2\f$ of the quadrature points.
     * @return The critical time step.
     * @throws std::invalid_argument for a lumped mass or a P-wave modulus that is not positive, NaN included: no time
     * increment is meaningful then.
     * @throws std::runtime_error if the eigenvalue solver fails.
     */
    double criticalTimeStep( const Eigen::VectorXd& lumpedMass, double pWaveModulus ) const
    {
      // every entry: minCoeff() does not promise to return a NaN entry
      if ( !( lumpedMass.array() > 0.0 ).all() )
        throw std::invalid_argument( "CriticalTimeStep: the element has a lumped mass that is not positive." );
      if ( !( pWaveModulus > 0.0 ) )
        throw std::invalid_argument( "CriticalTimeStep: non-positive wave speed or density." );

      const Eigen::VectorXd scaling = lumpedMass.cwiseInverse().cwiseSqrt();

      // Lambda, the bound of omega_max^2 per unit P-wave modulus: the larger of the two ends
      const double largestFrequencySquaredPerPWaveModulus = std::max( largestEigenvalue( scaling, lambdaOperator ),
                                                                      0.5 * largestEigenvalue( scaling, muOperator ) );

      return 2.0 / std::sqrt( pWaveModulus * largestFrequencySquaredPerPWaveModulus );
    }

  private:
    /// The operator \f$\mathbf{K}_\lambda\f$ of the Lame constant \f$\lambda\f$, \f$\int \mathbf{b}\,\mathbf{b}^T\f$.
    Eigen::MatrixXd lambdaOperator;
    /// The operator \f$\mathbf{K}_\mu\f$ of the shear modulus, \f$\int \mathbf{B}^T \mathbf{D}\,\mathbf{B}\f$.
    Eigen::MatrixXd muOperator;

    /**
     * @brief The largest eigenvalue of \f$\mathbf{M}^{-1}\mathbf{K}\f$, as that of the symmetric
     * \f$\mathbf{M}^{-1/2}\mathbf{K}\mathbf{M}^{-1/2}\f$.
     * @param scaling The diagonal \f$\mathbf{M}^{-1/2}\f$.
     * @param K The symmetric operator.
     * @return The largest eigenvalue.
     */
    static double largestEigenvalue( const Eigen::VectorXd& scaling, const Eigen::MatrixXd& K )
    {
      const Eigen::MatrixXd                                  scaled = scaling.asDiagonal() * K * scaling.asDiagonal();
      const Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > solver( scaled, Eigen::EigenvaluesOnly );
      if ( solver.info() != Eigen::Success )
        throw std::runtime_error( "CriticalTimeStep: the eigenvalues of the element's operators did not converge." );
      return solver.eigenvalues().maxCoeff();
    }
  };

} // namespace Marmot::Discretizations::Elements::CriticalTimeStep
