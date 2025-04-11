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
#include "Marmot/HaighWestergaard.h"
#include <utility>

namespace Marmot {
  namespace ContinuumMechanics::CommonConstitutiveModels {

    /**
     * \brief Generalized failure criterion proposed by Menétrey and Willam
     *
     * Implementation of the generalized failure criterion proposed by Menétrey and Willam \cite{MenetreyWillam1995}
     *
     * \f[
     *  f(\xi,\rho,\theta) = \left(A_f\,\rho\right)^2+m\left(B_f\,\rho\,r(\theta,e)+C_f\,\xi\right) - 1
     * \f]
     *
     * with the Haigh-Westergaard coordinates \f$\xi,\,\rho,\,\theta\f$ and the parameters \f$A_f,\,B_f,\,C_f,\,m,\,e\f$
     * which are automatically evaluated dependent on the desired formulation (supported types see \ref
     * MenetreyWillamType).
     *
     * \image html yieldSurfaceDeviatoric_MenetreyWillam.png width=75%
     * \image latex yieldSurfaceDeviatoric_MenetreyWillam.png width=75%
     *
     * # Example
     *
     * For a given pair of principal stress components in the Haigh-Westergaard stress space, evaluate the Mohr-Coulomb
     * failure criterion:
     */
    // clang-format off
     /**
     * @code
     * #include <Marmot/MenetreyWillam.h>
     * #include <Marmot/HaighWestergaard.h>
     *
     * using namespace Marmot;
     *
     * void main()
     * {
     *   const double ft = 1.0;
     *   const double fc = 10.0;
     *   auto mw = MenetreyWillam( ft , MenetreyWillamType::MohrCoulomb, fc );
     *   const HaighWestergaard::HaighWestergaardCoordinates hw = { 1.0, -10, 0 };
     *   if ( mw.yieldCriterion( hw ) >= 0 ) 
     *     std::cout << "Material is yielding!" << std::endl;
     *   else
     *     std::cout << "Material is elastic!" << std::endl;
     * }
     * @endcode
     */
    // clang-format on
    class MenetreyWillam {
    public:
      /**
       * Aggregate of the parameters for the generalized strength criterion
       * proposed by Menétrey and Willam.
       */
      struct MenetreyWillamParameters {
        double Af; /**< additional parameter to capture the influence of the
                      deviatoric stress invariant \f$\rho\f$ */
        double Bf; /**< addditional parameter to capture the influence of the
                      deviatoric stress invariant \f$\rho\f$ and the lode angle */
        double Cf; /**< addditional parameter to capture the influence of the
                      hydrostatic stress invariant \f$\xi\f$ */
        double m;  /**< Friction parameter \f$m\f$ */
        double e;  /**< Eccentricity parameter \f$e\f$; to obtain a smooth and
                      convex surface $e$ has to be in the range of \f$0.5\leq e
                      \leq 1\f$ */
      } param;

      // Reduction of the generalized failure criterion to a specific type.
      enum class MenetreyWillamType {
        Mises,         /**< von-Mises failure criterion; only the tensile strength @ref
                          ft has to be specified. */
        Rankine,       /**< Rankine failure criterion; only the tensile strength @ref
                          ft has to be specified. \note In deviatoric sections, compared to the original Rankine criterion the
                          vertices are slightly rounded (specified by the eccentricity parameter \f$e\f$). */
        DruckerPrager, /**< Drucker-Prager failure criterion; both the tensile
                          strength @ref ft as well as the compressive strength
                          @ref fc have to be specified. */
        MohrCoulomb    /**< Drucker-Prager failure criterion; both the tensile
                          strength @ref ft as well as the compressive strength @ref
                          fc have to be specified. \note In deviatoric sections, compared to the original Mohr-Coulomb
                          criterion the vertices are slightly rounded (specified by the eccentricity parameter \f$e\f$). */
      };

      /**
       * Constructor that takes the uniaxial tensile strength \f$f_t\f$ and two
       * optional arguments consisting of the specific type of failure criterion \ref MenetreyWillamType
       * and the uniaxial compressive strength \f$f_c\f$. The call of the
       * constructor automatically fills the corresponding Menetrey-Willam \ref param.
       */
      MenetreyWillam( const double              ft,
                      const MenetreyWillamType& type = MenetreyWillamType::Rankine,
                      const double              fc   = 0 );

      /**
       * This function can be used to reset the type of the specified failure
       * criterion by entering the uniaxial compressive strength \ref fc, the uniaxial tensile strength \ref ft
       * and the \ref MenetreyWillamType.
       */
      void setParameters( const double ft, const double fc, const MenetreyWillamType& type );

      /**
       * Compute the polar radius \f$r\f$ from the Lode angle \f$\theta\f$. The
       * eccentricity parameter will be used from the chosen Menetrey-Willam parameters \ref param.
       */
      template < typename T >
      T polarRadius( const double& theta ) const
      {
        return polarRadius( theta, param.e );
      }

      /**
       * Static version for computing the polar radius \f$r\f$ from the Lode angle
       * \f$\theta\f$ and a specified value for the eccentricity parameter
       * \f$e\f$.
       */
      template < typename T >
      static T polarRadius( const T& theta, const double& e )
      {
        // computes the deviatoric shape roundness for a given eccentricity (e) at a certain
        // position (theta) the numerator and denominator are stored for performance reasons, as
        // they are also needed for the derivative dRdTheta
        if ( e >= 1.0 )
          return 1;

        const double e2        = e * e;
        const T      cosTheta  = cos( theta );
        const T      cos2Theta = cosTheta * cosTheta;

        const T numerator   = ( 4. * ( 1. - e2 ) * cos2Theta + ( 2. * e - 1. ) * ( 2. * e - 1. ) );
        const T denominator = ( 2. * ( 1. - e2 ) * cosTheta +
                                ( 2. * e - 1. ) * sqrt( 4. * ( 1. - e2 ) * cos2Theta + 5. * e2 - 4. * e ) );

        const T r = numerator / denominator;

        return r;
      }
      /**
       * Compute the polar radius \f$r\f$ and its derivative
       * \f$\frac{dr}{d\theta}\f$ from the Lode angle \f$\theta\f$. The
       * eccentricity parameter will be used from the @ref param struct.
       */
      template < typename T >
      std::pair< T, T > dPolarRadius_dTheta( const T& theta ) const
      {
        return dPolarRadius_dTheta( theta, param.e );
      }

      /**
       * Static version for computing the polar radius \f$r$\f$  and its
       * derivative \f$\frac{dr}{d\theta}\f$ from the Lode angle \f$\theta\f$
       * and a specified value for the eccentricity parameter \f$e\f$.
       */
      template < typename T >
      static std::pair< T, T > dPolarRadius_dTheta( const T& theta, const double& e )
      {
        //
        // computes the deviatoric shape roundness for a given eccentricity (e) at a certain
        // position (theta) the numerator and denominator are stored for performance reasons, as
        // they are also needed for the derivative dRdTheta
        T r, dRdTheta;

        if ( e >= 1.0 ) {
          r        = 1;
          dRdTheta = 0;
          return std::make_pair( r, dRdTheta );
        }

        const double e2        = e * e;
        const T      cosTheta  = cos( theta );
        const T      cos2Theta = cosTheta * cosTheta;

        T numerator   = ( 4. * ( 1. - e2 ) * cos2Theta + ( 2. * e - 1. ) * ( 2. * e - 1. ) );
        T denominator = ( 2. * ( 1. - e2 ) * cosTheta +
                          ( 2. * e - 1. ) * sqrt( 4. * ( 1. - e2 ) * cos2Theta + 5. * e2 - 4. * e ) );

        r = numerator / denominator;

        // compute the derivative
        const T sinTheta = sin( theta );

        const T a          = numerator;
        const T b          = 1. / denominator;
        const T aux        = 4. * ( 1. - e2 ) * cos2Theta + 5. * e2 - 4. * e;
        const T dAuxdTheta = 4. * ( 1. - e2 ) * 2. * cosTheta * -sinTheta;
        const T dadTheta   = 2. * cosTheta * ( -sinTheta ) * 4. * ( 1. - e2 );
        const T dbdTheta   = -pow( denominator, -2. ) * ( 2. * ( 1. - e2 ) * -sinTheta +
                                                        ( 2. * e - 1 ) * 1. / 2 * pow( aux, -1. / 2 ) * dAuxdTheta );
        dRdTheta           = b * dadTheta + a * dbdTheta;

        return { r, dRdTheta };
      }

      template < typename T >
      static std::tuple< T, T, T > d2PolarRadius_dTheta2( const T& theta, const double& e )
      {

        const auto [r, dr_dTheta] = dPolarRadius_dTheta< T >( theta, e );
        if ( e >= 1.0 ) {
          return { r, dr_dTheta, 0 };
        }

        const double e2   = pow( e, 2 );
        const double e2m1 = e2 - 1;
        const double tem1 = 2 * e - 1;

        const T cos_theta   = cos( theta );
        const T sin_theta   = sin( theta );
        const T cos_2_theta = cos_theta * cos_theta;
        const T sin_2_theta = sin_theta * sin_theta;

        const T x      = 5. * e2 - 4. * e - 4. * e2m1 * cos_2_theta;
        const T sqrt_x = sqrt( x );

        const T one_over_sqrt_x = 1. / sqrt_x;

        const T d2r_dTheta2 = e2m1 *
                              ( -16. * e2m1 * ( 4.0 * tem1 * one_over_sqrt_x * cos_theta + 2. ) * sin_2_theta *
                                  cos_theta / ( tem1 * sqrt_x - 2 * e2m1 * cos_theta ) -
                                8. * sin_2_theta + 8. * cos_2_theta +
                                ( pow( 2. * e - 1., 2 ) - 4 * e2m1 * cos_2_theta ) *
                                  ( 16.0 * tem1 * e2m1 * pow( x, -1.5 ) * sin_2_theta * cos_2_theta +
                                    4.0 * tem1 * one_over_sqrt_x * sin_2_theta -
                                    4.0 * tem1 * one_over_sqrt_x * cos_2_theta +
                                    e2m1 * ( 4.0 * tem1 * one_over_sqrt_x * cos_theta + 2. ) *
                                      ( 8.0 * tem1 * one_over_sqrt_x * cos_theta + 4. ) * sin_2_theta /
                                      ( tem1 * sqrt_x - 2. * e2m1 * cos_theta ) -
                                    2. * cos_theta ) /
                                  ( tem1 * sqrt_x - 2. * e2m1 * cos_theta ) ) /
                              ( tem1 * sqrt_x - 2. * e2m1 * cos_theta );

        return { r, dr_dTheta, d2r_dTheta2 };
      }

      /**
       * Evaluate the yield function \f$f\f$ depending on the Haigh-Westergaard stress
       * coordinates @ref hw. \f$f<0\f$ means no yielding while \f$f\geq0\f$ means
       * yielding. When the optional fillet parameter \ref varEps is given, a potential vertex
       * along the hydrostatic axis is rounded and thus a smooth failure criterion is obtained.
       *
       * \note The yield function can be also used as plastic potential function if needed.
       */
      template < typename T >
      T yieldFunction( const ContinuumMechanics::HaighWestergaard::HaighWestergaardCoordinates< T >& hw,
                       const double varEps = 0.0 ) const
      {
        const T r_ = polarRadius( hw.theta, param.e );
        if ( varEps == 0 )
          return ( param.Af * hw.rho ) * ( param.Af * hw.rho ) +
                 param.m * ( param.Bf * hw.rho * r_ + param.Cf * hw.xi ) - 1.;
        else
          return param.Af * param.Af * hw.rho * hw.rho +
                 param.m *
                   ( sqrt( param.Bf * hw.rho * r_ * param.Bf * hw.rho * r_ + varEps * varEps ) + param.Cf * hw.xi ) -
                 1.;
      }

      /**
       * Evaluate the derivatives of the yield function with respect to the
       * Haigh-Westergaard stress coordinates, i.e.
       * \f$df/d\xi,\,df/d\rho,\,df/d\theta\f$. When the optional fillet parameter \ref varEps is given, a potential
       * vertex along the hydrostatic axis is rounded. Thus, a smooth failure criterion is obtained, where derivatives
       * can be calculated uniquely.
       */
      template < typename T >
      std::tuple< T, T, T > dYieldFunction_dHaighWestergaard(
        const ContinuumMechanics::HaighWestergaard::HaighWestergaardCoordinates< T >& hw,
        const double                                                                  varEps = 0.0 ) const
      {
        const auto [r_, dRdTheta_] = dPolarRadius_dTheta( hw.theta, param.e );

        T dFdXi, dFdRho, dFdTheta;
        dFdXi = param.m * param.Cf;

        if ( varEps == 0.0 ) {
          dFdRho   = 2. * param.Af * hw.rho * param.Af + param.m * param.Bf * r_;
          dFdTheta = param.m * param.Bf * hw.rho * dRdTheta_;
        }
        else {
          const T auxTerm1 = param.m * 0.5 *
                             pow( param.Bf * param.Bf * hw.rho * hw.rho * r_ * r_ + varEps * varEps, -0.5 ) * 2 *
                             param.Bf * hw.rho * r_ * param.Bf;

          dFdRho   = param.Af * param.Af * 2. * hw.rho + auxTerm1 * r_;
          dFdTheta = auxTerm1 * hw.rho * dRdTheta_;
        }

        return { dFdXi, dFdRho, dFdTheta };
      }

      /**
       * Compute a fillet parameter for the vertex of the yield surface along the hydrostatic axis in the same way as
       * Abaqus does. This parameter is only relevant in the case of the Drucker-Prager or the Mohr-Coulomb criterion.
       * The calculated smoothing value can then be used as the optional input argument \ref varEps in \ref
       * yieldFunction and \ref dYieldFunction_dHaighWestergaard.
       */
      static inline double abaqusMohrCoulombPotentialVarEpsToMenetreyWillam( const double varEps, const double psi )
      {
        return varEps * 2 * std::sin( psi );
      }

      /**
       * Compute the eccentricity parameter based on a given compressive strength \ref fc and
       * tensile strength \ref ft.
       */
      static inline double e( const double fc, const double ft ) { return ( fc + 2 * ft ) / ( 2 * fc + ft ); }

      /**
       * Compute the cohesion based on a given compressive strength \ref fc and
       * tensile strength \ref ft.
       */
      static inline double c( const double fc, const double ft )
      {
        const double phi_ = phi( fc, ft );
        return ft / 2. * ( 1 + std::sin( phi_ ) ) / ( std::cos( phi_ ) );
      }

      /**
       * Compute the friction angle in radian based on a given compressive strength \ref fc and
       * tensile strength \ref ft.
       */
      static inline double phi( const double fc, const double ft ) { return std::asin( ( fc - ft ) / ( fc + ft ) ); }

      /**
       * Compute the tensile strength based on a given cohesion \ref c and a friction angle \ref phi.
       */
      static inline double ft( const double c, const double phi )
      {
        return 2 * c * std::cos( phi ) / ( 1 + std::sin( phi ) );
      }

      /**
       * Compute the compressive strength based on a given cohesion \ref c and a friction angle \ref phi.
       */
      static inline double fc( const double c, const double phi )
      {
        return 2 * c * std::cos( phi ) / ( 1 - std::sin( phi ) );
      }
    };

  } // namespace ContinuumMechanics::CommonConstitutiveModels
} // namespace Marmot
