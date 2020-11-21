#pragma once

namespace Marmot {
    namespace MW {

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

    } // namespace MW
} // namespace Marmot
