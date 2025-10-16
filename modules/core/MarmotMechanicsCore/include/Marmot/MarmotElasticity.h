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
 * Thomas Mader thomas.mader@uibk.ac.at
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
#include "Marmot/MarmotVoigt.h"

namespace Marmot {

  namespace ContinuumMechanics {

    /**
     * @brief Functions for the description of isotropic elastic behavior
     */
    namespace Elasticity::Isotropic {

      /**
       * @brief Computes the isotropic Young's modulus \f$ E \f$ from the compression modulus \f$ K \f$
       * and shear modulus \f$ G \f$.
       *
       *\f[
          \displaystyle E = \frac{9\,K\,G}{3\,K + G}
        \f]
       *
       * @param K Compression modulus \f$ K \f$.
       * @param G Shear modulus \f$ G \f$.
       * @return Young's modulus \f$ E \f$.
       */
      double constexpr E( const double K, const double G )
      {
        return 9. * K * G / ( 3. * K + G );
      }

      /**
       * @brief Computes the isotropic Poisson's ratio \f$ \nu \f$ from the compression modulus \f$ K \f$
       * and shear modulus \f$ G \f$.
       *
       *\f[
          \displaystyle \nu = \frac{3\,K - 2\,G}{6\,K + 2\,G}
         \f]
       *
       * @param K Compression modulus \f$ K \f$.
       * @param G Shear modulus \f$ G \f$.
       * @return Poisson's ratio \f$ \nu \f$.
       */
      double constexpr nu( const double K, const double G )
      {
        return ( 3 * K - 2 * G ) / ( 6 * K + 2 * G );
      }

      /**
       * @brief Computes the isotropic shear modulus \f$ G \f$ from the Young's modulus \f$ E \f$
       * and Poisson's ratio \f$ \nu \f$.
       *
       *\f[
         \displaystyle G = \frac{E}{2\,(1 + \nu)}
        \f]
       * @param E Young's modulus \f$ E \f$.
       * @param nu Poisson's ratio \f$ \nu \f$.
       * @return Shear modulus \f$ G \f$.
      */
      double constexpr shearModulus( const double E, const double nu )
      {
        return E / ( 2 * ( 1 + nu ) );
      }

      /**
       * @brief Computes the isotropic Lamé parameter \f$ \lambda \f$ from the Young's modulus \f$ E \f$ and Poisson's
       * ratio \f$ \nu \f$.
       *
       *\f[
         \displaystyle \lambda = \frac{E\,\nu}{(1 + \nu)(1 - 2\,\nu)}
        \f]
       *
       * @param E Young's modulus \f$ E \f$.
       * @param nu Poisson's ratio \f$ \nu \f$.
       * @return Lamé parameter \f$ \lambda \f$.
       */
      double constexpr lameParameter( const double E, const double nu )
      {
        return E * nu / ( ( 1 + nu ) * ( 1 - 2 * nu ) );
      }

      /**
       * @brief Computes the isotropic compliance tensor \f$\mathbb{C}^{-1}\f$.
       *\f[
         \displaystyle \mathbb{ C }^{-1} = \begin{bmatrix}
                        \frac{1}{E} & \frac{-\nu}{E} & \frac{-\nu}{E} & 0 & 0 & 0 \\
                        \frac{-\nu}{E} & \frac{1}{E} & \frac{-\nu}{E} & 0 & 0 & 0 \\
                        \frac{-\nu}{E} & \frac{-\nu}{E} & \frac{1}{E} & 0 & 0 & 0 \\
                  0 & 0 & 0 & \frac{1}{G} & 0 & 0 \\
                  0 & 0 & 0 & 0 & \frac{1}{G} & 0 \\
                  0 & 0 & 0 & 0 & 0 & \frac{1}{G}
                   \end{bmatrix}
        \f]
       *
       * with
       *
       * \f[
       *   G = \frac{E}{2\,(1 + \nu)}.
       * \f]
       *
       * @param E Young's modulus \f$E\f$.
       * @param nu Poisson's ratio \f$\nu\f$.
       * @return Compliance tensor \f$\mathbb{C}^{-1}\f$ as a 6x6 matrix.
       */
      Matrix6d complianceTensor( const double E, const double nu );

      /**
       * @brief Computes the isotropic stiffness tensor \f$\mathbb{C}\f$ from
       * Young's modulus \f$E\f$ and Poisson's ratio \f$\nu\f$.
       *
       *\f[
       \displaystyle \mathbb{ C } = \frac{E\,(1-\nu)}{(1+\nu)(1-2\,\nu)} \begin{bmatrix}
                                       1 & \frac{\nu}{1-\nu} & \frac{\nu}{1-\nu} & 0 & 0 & 0 \\
                                                   \frac{\nu}{1-\nu} & 1 & \frac{\nu}{1-\nu} & 0 & 0 & 0 \\
                                             \frac{\nu}{1-\nu} & \frac{\nu}{1-\nu} & 1 & 0 & 0 & 0 \\
                                         0 & 0 & 0 & \frac{1-2\,\nu}{2\,(1-\nu)} & 0 & 0 \\
                                             0 & 0 & 0 & 0 & \frac{1-2\,\nu}{2\,(1-\nu)} & 0 \\
                                             0 & 0 & 0 & 0 & 0 & \frac{1-2\,\nu}{2\,(1-\nu)}
                                        \end{bmatrix}
        \f]

       * @param E Young's modulus \f$E\f$.
       * @param nu Poisson's ratio \f$\nu\f$.
       * @return Stiffness tensor \f$\mathbb{C}\f$ as a 6x6 matrix.
       */
      Matrix6d stiffnessTensor( const double E, const double nu );

      /**
       * @brief Computes the isotropic stiffness tensor \f$\mathbb{C}\f$ from
       * the bulk modulus \f$K\f$ and shear modulus \f$G\f$.
       *
       * @param K Bulk modulus \f$K\f$.
       * @param G Shear modulus \f$G\f$.
       * @return Stiffness tensor \f$\mathbb{C}\f$ as a 6x6 matrix.
       */
      Matrix6d stiffnessTensorKG( const double K, const double G );

    } // namespace Elasticity::Isotropic

    /**
     * @brief Functions for the description of transversely isotropic elastic behavior
     */
    namespace Elasticity::TransverseIsotropic {

      /**
       * @brief Computes the transversely isotropic compliance tensor \f$\mathbb{C}^{-1}\f$.
       *\f[
          \displaystyle \mathbb{ C }^{-1} = \begin{bmatrix}
                       \frac{1}{E_1} & \frac{-\nu_{12}}{E_1} & \frac{-\nu_{12}}{E_1} & 0 & 0 & 0 \\
                             \frac{-\nu_{12}}{E_1} & \frac{1}{E_2} & \frac{-\nu_{23}}{E_2} & 0 & 0 & 0 \\
                             \frac{-\nu_{12}}{E_1} & \frac{-\nu_{23}}{E_2} & \frac{1}{E_2} & 0 & 0 & 0 \\
                         0 & 0 & 0 & \frac{1}{G_{12}} & 0 & 0 \\
                         0 & 0 & 0 & 0 & \frac{1}{G_{12}} & 0 \\
                         0 & 0 & 0 & 0 & 0 & \frac{1}{G_{23}}
                                          \end{bmatrix}
       \f]
       *
       * The isotropic plane is defined with respect to the \f$x_2\f$–\f$x_3\f$ axes of a local coordinate system.
       * The tensor is expressed in terms of the out-of-plane Young's modulus \f$E_1\f$, shear modulus \f$G_{12}\f$,
       * and Poisson's ratio \f$\nu_{12}\f$, together with the in-plane Young's modulus \f$E_2\f$ and
       * Poisson's ratio \f$\nu_{23}\f$.
       *
       *The in-plane shear modulus \f$ G_{23} \f$ is given by
       *
       \f[
         \displaystyle G_{23} = \frac{E_2}{2\,(1 + \nu_{23})}
       \f]
       *
       * @param E1 Out-of-plane Young's modulus \f$E_1\f$.
       * @param E2 In-plane Young's modulus \f$E_2\f$.
       * @param nu12 Poisson's ratio \f$\nu_{12}\f$.
       * @param nu23 Poisson's ratio \f$\nu_{23}\f$.
       * @param G12 shear modulus \f$G_{12}\f$.
       * @return Compliance tensor \f$\mathbb{C}^{-1}\f$ as a 6x6 matrix.
       */
      Matrix6d complianceTensor( const double E1,
                                 const double E2,
                                 const double nu12,
                                 const double nu23,
                                 const double G12 );
      /**
       * @brief Computes the transversely isotropic stiffness tensor \f$ \mathbb{C} \f$ as inverse of the transversely
       * isotropic compliance tensor \f$ \mathbb{C}^{-1} \f$.
       *
       * @param E1 Out-of-plane Young's modulus \f$E_1\f$.
       * @param E2 In-plane Young's modulus \f$E_2\f$.
       * @param nu12 Poisson's ratio \f$\nu_{12}\f$.
       * @param nu23 Poisson's ratio \f$\nu_{23}\f$.
       * @param G12 shear modulus \f$G_{12}\f$.
       * @return Stiffness tensor \f$\mathbb{C}\f$ as a 6x6 matrix.
       */
      Matrix6d stiffnessTensor( const double E1,
                                const double E2,
                                const double nu12,
                                const double nu23,
                                const double G12 );
    } // namespace Elasticity::TransverseIsotropic

    /**
     * @brief Functions for the description of orthotropic elastic behavior
     */
    namespace Elasticity::Orthotropic {
      /**
       * @brief Computes the orthotropic compliance tensor \f$\mathbb{C}^{-1}\f$,
       * defined in the principal material directions \f$x_1\f$, \f$x_2\f$, and \f$x_3\f$.
       *
       *\f[
          \displaystyle \mathbb{ C }^{-1} = \begin{bmatrix}
                         \frac{1}{E_1} & \frac{-\nu_{12}}{E_2} & \frac{-\nu_{13}}{E_3} & 0 & 0 & 0 \\
                         \frac{-\nu_{12}}{E_2} & \frac{1}{E_2} & \frac{-\nu_{23}}{E_3} & 0 & 0 & 0 \\
                         \frac{-\nu_{13}}{E_3} & \frac{-\nu_{23}}{E_3} & \frac{1}{E_3} & 0 & 0 & 0 \\
                       0 & 0 & 0 & \frac{1}{G_{12}} & 0 & 0 \\
                       0 & 0 & 0 & 0 & \frac{1}{G_{13}} & 0 \\
                       0 & 0 & 0 & 0 & 0 & \frac{1}{G_{23}}
                         \end{bmatrix}
        \f]
       *
       * @param E1 Young's modulus \f$E_1\f$ in \f$ x_1 \f$ - direction.
       * @param E2 Young's modulus \f$E_2\f$ in \f$ x_2 \f$ - direction.
       * @param E3 Young's modulus \f$E_3\f$ in \f$ x_3 \f$ - direction.
       * @param nu12 Poisson's ratio \f$\nu_{12}\f$ effective between \f$ x_1 \f$ and \f$ x_2 \f$ - direction.
       * @param nu23 Poisson's ratio \f$\nu_{23}\f$ effective between \f$ x_2 \f$ and \f$ x_3 \f$ - direction.
       * @param nu13 Poisson's ratio \f$\nu_{13}\f$ effective between \f$ x_1 \f$ and \f$ x_3 \f$ - direction.
       * @param G12 Shear modulus \f$G_{12}\f$ effective between \f$ x_1 \f$ and \f$ x_2 \f$ - direction.
       * @param G23 Shear modulus \f$G_{23}\f$ effective between \f$ x_2 \f$ and \f$ x_3 \f$ - direction.
       * @param G31 Shear modulus \f$G_{31}\f$ effective between \f$ x_1 \f$ and \f$ x_3 \f$ - direction.
       * @return Compliance tensor \f$\mathbb{C}^{-1}\f$ as a 6x6 matrix.
       */
      Matrix6d complianceTensor( const double E1,
                                 const double E2,
                                 const double E3,
                                 const double nu12,
                                 const double nu23,
                                 const double nu13,
                                 const double G12,
                                 const double G23,
                                 const double G31 );
      /**
       * @brief Computes the orthotropic stiffness tensor \f$ \mathbb{C} \f$ as inverse of the orthotropic compliance
       * tensor
       * \f$ \mathbb{C}^{-1} \f$.
       *
       * @param E1 Young's modulus \f$E_1\f$ in \f$ x_1 \f$ - direction.
       * @param E2 Young's modulus \f$E_2\f$ in \f$ x_2 \f$ - direction.
       * @param E3 Young's modulus \f$E_3\f$ in \f$ x_3 \f$ - direction.
       * @param nu12 Poisson's ratio \f$\nu_{12}\f$ effective between \f$ x_1 \f$ and \f$ x_2 \f$ - direction.
       * @param nu23 Poisson's ratio \f$\nu_{23}\f$ effective between \f$ x_2 \f$ and \f$ x_3 \f$ - direction.
       * @param nu13 Poisson's ratio \f$\nu_{13}\f$ effective between \f$ x_1 \f$ and \f$ x_3 \f$ - direction.
       * @param G12 Shear modulus \f$G_{12}\f$ effective between \f$ x_1 \f$ and \f$ x_2 \f$ - direction.
       * @param G23 Shear modulus \f$G_{23}\f$ effective between \f$ x_2 \f$ and \f$ x_3 \f$ - direction.
       * @param G31 Shear modulus \f$G_{31}\f$ effective between \f$ x_1 \f$ and \f$ x_3 \f$ - direction.
       * @return Stiffness tensor \f$\mathbb{C}\f$ as a 6x6 matrix.
       */
      Matrix6d stiffnessTensor( const double E1,
                                const double E2,
                                const double E3,
                                const double nu12,
                                const double nu23,
                                const double nu13,
                                const double G12,
                                const double G23,
                                const double G31 );

    } // namespace Elasticity::Orthotropic
  }   // namespace ContinuumMechanics
} // namespace Marmot
