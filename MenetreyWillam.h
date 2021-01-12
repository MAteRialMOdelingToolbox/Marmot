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

        double r( double theta, double e, double& numerator, double& denominator );

        double dRdTheta( double theta, double e, double rNumerator, double rDenominator );

        double e( double fc, double ft );

        double c( double fc, double ft );

        double phi( double fc, double ft );

        double fc( double c, double phi );

        double ft( double c, double phi );

        double f( double Af, double Bf, double Cf, double m, double e, double xi, double rho, double theta );

        void dFdHaighWestergaard( double& dFdXi,
                                  double& dFdRho,
                                  double& dFdTheta,
                                  double  Af,
                                  double  Bf,
                                  double  Cf,
                                  double  m,
                                  double  e,
                                  double  xi,
                                  double  rho,
                                  double  theta );

        double fRounded( double Af,
                         double Bf,
                         double Cf,
                         double m,
                         double e,
                         double xi,
                         double rho,
                         double theta,
                         double varEps );

        void dFRoundeddHaighWestergaard( double& dFdXi,
                                         double& dFdRho,
                                         double& dFdTheta,
                                         double  Af,
                                         double  Bf,
                                         double  Cf,
                                         double  m,
                                         double  e,
                                         double  xi,
                                         double  rho,
                                         double  theta,
                                         double  varEps );

        double abaqusMohrCoulombPotentialVarEpsToMenetreyWillam( double varEps, double psi );

        void RankineParameters( double& Af, double& Bf, double& Cf, double& m, double& e, double ft, double fc = 0 );

        void MisesParameters( double& Af, double& Bf, double& Cf, double& m, double& e, double ft, double fc = 0 );

        void DruckerPragerParameters( double& Af, double& Bf, double& Cf, double& m, double& e, double ft, double fc );

        void MohrCoulombParameters( double& Af, double& Bf, double& Cf, double& m, double& e, double ft, double fc );

    } // namespace ContinuumMechanics::CommonConstitutiveModels::MenetreyWillam
} // namespace Marmot
