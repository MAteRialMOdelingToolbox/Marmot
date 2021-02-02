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

namespace Marmot {
    namespace ContinuumMechanics::CommonConstitutiveModels::MenetreyWillam {

        double r( const double theta, const double e, double& numerator, double& denominator );

        double dRdTheta( const double theta, const double e, const double rNumerator, const double rDenominator );

        double e( const double fc, const double ft );

        double c( const double fc, const double ft );

        double phi( const double fc, const double ft );

        double fc( const double c, const double phi );

        double ft( const double c, const double phi );

        double f( const double Af, const double Bf, const double Cf, const double m, const double e, const double xi, const double rho, const double theta );

        void dFdHaighWestergaard( double& dFdXi,
                                  double& dFdRho,
                                  double& dFdTheta,
                                  const double  Af,
                                  const double  Bf,
                                  const double  Cf,
                                  const double  m,
                                  const double  e,
                                  const double  xi,
                                  const double  rho,
                                  const double  theta );

        double fRounded( const double Af,
                         const double Bf,
                         const double Cf,
                         const double m,
                         const double e,
                         const double xi,
                         const double rho,
                         const double theta,
                         const double varEps );

        void dFRoundeddHaighWestergaard( double& dFdXi,
                                         double& dFdRho,
                                         double& dFdTheta,
                                         const double  Af,
                                         const double  Bf,
                                         const double  Cf,
                                         const double  m,
                                         const double  e,
                                         const double  xi,
                                         const double  rho,
                                         const double  theta,
                                         const double  varEps );

        double abaqusMohrCoulombPotentialVarEpsToMenetreyWillam( const double varEps, const double psi );

        void RankineParameters( double& Af, double& Bf, double& Cf, double& m, double& e, const double ft, const double fc = 0 );

        void MisesParameters( double& Af, double& Bf, double& Cf, double& m, double& e, const double ft, const double fc = 0 );

        void DruckerPragerParameters( double& Af, double& Bf, double& Cf, double& m, double& e, const double ft, const double fc );

        void MohrCoulombParameters( double& Af, double& Bf, double& Cf, double& m, double& e, const double ft, const double fc );

    } // namespace ContinuumMechanics::CommonConstitutiveModels::MenetreyWillam
} // namespace Marmot
