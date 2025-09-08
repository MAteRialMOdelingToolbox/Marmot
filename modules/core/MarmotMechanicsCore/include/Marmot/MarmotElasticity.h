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
     * \brief Functions for the description of isotropic elastic behavior
     */
    namespace Elasticity::Isotropic {

      /**
       * Computes the isotropic young's modulus E from the compression modulus K and shear modulus G
       *\f[
          \displaystyle E = \frac{9\,K\,G}{3\,K + G}
        \f]
       */
      double constexpr E( const double K, const double G )
      {
        return 9. * K * G / ( 3. * K + G );
      }

      /**
       * Computes the isotropic poisson's ratio \f$ \nu \f$ from the compression modulus K and shear modulus G
       *\f[
         \displaystyle \nu = \frac{3\,K - 2\,G}{6\,K + 2\,G}
        \f]
       */
      double constexpr nu( const double K, const double G )
      {
        return ( 3 * K - 2 * G ) / ( 6 * K + 2 * G );
      }

      /**
       * Computes the isotropic shear modulus G from the young's modulus E and poisson's ratio \f$ \nu \f$.
       *\f[
         \displaystyle G = \frac{E}{2\,(1 + \nu)}
        \f]
       */
      double constexpr shearModulus( const double E, const double nu )
      {
        return E / ( 2 * ( 1 + nu ) );
      }

      /**
       * Computes the isotropic lame parameter \f$ \lambda \f$ from the young's modulus E and poisson's ratio \f$ \nu
       \f$
       *\f[
         \displaystyle \lambda = \frac{E\,\nu}{(1 + \nu)(1 - 2\,\nu)}
        \f]
       */
      double constexpr lameParameter( const double E, const double nu )
      {
        return E * nu / ( ( 1 + nu ) * ( 1 - 2 * nu ) );
      }

      /**
       *Computes the isotropic compliance tensor:
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

        from the young's modulus E and the poisson's ratio \f$ \nu \f$ with

        \f[
         \displaystyle G = \frac{E}{2\,(1 + \nu)}
        \f]
       */
      Matrix6d complianceTensor( const double E, const double nu );

      /**
       *Computes the isotropic stiffness tensor:
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

       from the young's modulus E and the poisson's ratio \f$ \nu \f$.
      */
      Matrix6d stiffnessTensor( const double E, const double nu );

      /**
       *Computes the isotropic stiffness tensor \f$\mathbb{ C }\f$ from the bulk modulus K and the shear modulus G.
       */
      Matrix6d stiffnessTensorKG( const double K, const double G );

    } // namespace Elasticity::Isotropic

    /**
     * \brief Functions for the description of transversely isotropic elastic behavior
     */
    namespace Elasticity::TransverseIsotropic {

      /**
       *Computes the transversely isotropic compliance tensor:
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
       *The isotropic plane is defined with respect to the  \f$ x_2 \f$ and \f$ x_3 \f$ axes of a local coordinate
       system. It is computed from the young's modulus \f$ E_1 \f$, the shear modulus \f$ G_{12} \f$ and poisson's ratio
       \f$ \nu_{12} \f$ effective out of the isotropic plane and the in - plane young's modulus \f$ E_2 \f$ and
       poisson's ratio \f$ \nu_{23} \f$.
       *
       *The in - plane shear modulus \f$ G_{23} \f$ can be expressed by

       \f[
         \displaystyle G_{23} = \frac{E_2}{2\,(1 + \nu_{23})}
       \f]
       */
      Matrix6d complianceTensor( const double E1,
                                 const double E2,
                                 const double nu12,
                                 const double nu23,
                                 const double G12 );
      /**
       * Computes the transversely isotropic stiffness tensor \f$ \mathbb{C} \f$ as inverse of the transversely
       * isotropic compliance tensor \f$ \mathbb{C}^{-1} \f$.
       */
      Matrix6d stiffnessTensor( const double E1,
                                const double E2,
                                const double nu12,
                                const double nu23,
                                const double G12 );
    } // namespace Elasticity::TransverseIsotropic

    /**
     * \brief Functions for the description of orthotropic elastic behavior
     */
    namespace Elasticity::Orthotropic {
      /**
       *Computes the orthotropic compliance tensor defined in the principal directions of the material \f$ x_1 \f$, \f$
       x_2 \f$ and \f$ x_3 \f$ :
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
       *from the following independent parameters:
       *<table>
       *<tr><th>Parameter        <th>Description
       *<tr><td>\f$ E_1 \f$      <td> Young's modulus effective in \f$ x_1 \f$ - direction
       *<tr><td>\f$ E_2 \f$      <td> Young's modulus effective in \f$ x_2 \f$ - direction
       *<tr><td>\f$ E_3 \f$ 	 <td> Young's modulus effective in \f$ x_3 \f$ - direction
       *<tr><td>\f$ \nu_{12} \f$ <td> Poisson's ratio effective between \f$ x_1 \f$ and \f$ x_2 \f$ - direction
       *<tr><td>\f$ \nu_{13} \f$ <td> Poisson's ratio effective between \f$ x_1 \f$ and \f$ x_3 \f$ - direction
       *<tr><td>\f$ \nu_{23} \f$ <td> Poisson's ratio effective between \f$ x_2 \f$ and \f$ x_3 \f$ - direction
       *<tr><td>\f$ G_{12} \f$ <td> Shear modulus effective between \f$ x_1 \f$ and \f$ x_2 \f$ - direction
       *<tr><td>\f$ G_{13} \f$ <td> Shear modulus effective between \f$ x_1 \f$ and \f$ x_3 \f$ - direction
       *<tr><td>\f$ G_{23} \f$ <td> Shear modulus effective between \f$ x_2 \f$ and \f$ x_3 \f$ - direction
       *</table>
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
       * Computes the orthotropic stiffness tensor \f$ \mathbb{C} \f$ as inverse of the orthotropic compliance tensor
       * \f$ \mathbb{C}^{-1} \f$.
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
