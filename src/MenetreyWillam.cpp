#include "MenetreyWillam.h"
#include <cmath>
#include "bftConstants.h"

namespace bft{
    namespace MW
    {
        using namespace Constants;
        double r(double theta, double e, double& numerator, double &denominator)
        {
            // computes the deviatoric shape roundness for a given eccentricity (e) at a certain position (theta)
            // the numerator and denominator are stored for performance reasons, as they are also needed for the derivative dRdTheta
            if(e>=1.0){
                numerator=1;
                denominator = 1;
                return 1;
            }

            const double e2 =               e*e;
            const double cosTheta =         std::cos(theta);
            const double cos2Theta =        cosTheta * cosTheta;

            numerator =                     (4*(1-e2)*cos2Theta + (2*e -1)*(2*e -1));
            denominator =                   (2*(1-e2)*cosTheta+(2*e-1)*std::sqrt(4*(1-e2)*cos2Theta+ 5*e2 -4*e));

            return numerator / denominator;
        }

        double dRdTheta(double theta, double e,  double rNumerator, double rDenominator)
        {
                // computes the derivate for a given r (defined by its numerator, denominator, calculated by r()), 

            if(e>=1.0)
                    return 0;

            const double cosTheta =         std::cos(theta);
            const double cos2Theta =        cosTheta * cosTheta;
            const double sinTheta =         std::sin(theta);
            const double e2 = e*e;
            const double a = rNumerator;
            const double b = 1./rDenominator;
            const double aux = 4*(1-e2)*cos2Theta+ 5*e2 -4*e;
            const double dAuxdTheta = 4 * (1-e2)*2*cosTheta*-sinTheta;
            const double dadTheta =         2 * cosTheta * (-sinTheta) * 4 * (1-e2);
            const double dbdTheta = - pow(rDenominator,-2) * ( 2*(1-e2)*-sinTheta    + (2*e-1)*1./2 * pow(aux, -1./2) * dAuxdTheta);
             
            // r = a * b 
            // = num * den^-1
            return b * dadTheta + a * dbdTheta;
        }
        double e(double fc, double ft)
        {
            return (fc+2*ft)/(2*fc + ft);
        }

        double c(double fc, double ft)
        {
            const double phi_ = phi(fc, ft);
            return ft * ( 1 + std::sin(phi_)) / (2 + std::cos(phi_));
        }

        double phi(double fc, double ft)
        {
            return std::asin ( (fc - ft) / (fc + ft));
        }

        double ft(double c, double phi)
        {
            return 2*c*std::cos(phi) / (1 + std::sin(phi));
        }

        double fc(double c, double phi)
        {
            return 2*c*std::cos(phi) / (1 - std::sin(phi));
        }


        double f(double Af, double Bf, double Cf, double m, double e, double xi, double rho, double theta)
        {
            double num = 0;
            double den = 0;
            return (Af * rho)*(Af*rho) + m * (Bf * rho * r(theta, e, num, den) + Cf * xi) - 1; 
        }

        void dFdHaighWestergaard(double&dFdXi, double&dFdRho, double&dFdTheta, double Af, double Bf, 
                double Cf, double m, double e, double xi, double rho, double theta)
        {
            double rNum=0, rDen=0;
            const double r_ = r(theta, e, rNum, rDen);
            const double dRdTheta_ = dRdTheta(theta, e, rNum, rDen);

            dFdXi =         m * Cf;
            dFdRho =        2*Af*rho*Af +  m*Bf*r_;
            dFdTheta =      m*Bf*rho*dRdTheta_;
        }

        double fRounded(double Af, double Bf, double Cf, double m, double e, double xi, double rho, double theta, double varEps)
        {
            double rNum=0, rDen=0;
            const double r_ = r(theta, e, rNum, rDen);
            return Af*Af*rho*rho + m* (std::sqrt( Bf*rho*r_*Bf*rho*r_ + varEps*varEps) + Cf*xi) -1;
        }

        void dFRoundeddHaighWestergaard(double&dFdXi, double&dFdRho, double&dFdTheta, double Af, double Bf, 
                double Cf, double m, double e, double xi, double rho, double theta, double varEps)
        {
            double rNum=0, rDen=0;
            const double r_ = r(theta, e, rNum, rDen);
            const double dRdTheta_ = dRdTheta(theta, e, rNum, rDen);

            const double auxTerm1 = m * 1./2 * std::pow( Bf*Bf*rho*rho*r_*r_ + varEps*varEps, -1./2) * 2*Bf*rho*r_ * Bf;

            dFdXi =         m * Cf;
            dFdRho =        Af*Af*rho*rho + auxTerm1 * r_;
            dFdTheta =      auxTerm1 * rho * dRdTheta_;
        }

        void RankineParameters(double& Af, double&Bf, double& Cf, double&m, double& e, double ft, double fc)
        {
            Af =    0;
            Bf =    1./ ( sqrt6 * ft);
            Cf =    1./ ( sqrt3 * ft);
            m =     1;
            e =     0.51;
        }

        void MisesParameters(double& Af, double&Bf, double& Cf, double&m, double& e, double ft, double fc)
        {
            Af =    0;
            Bf =    sqrt3_2 / ft;
            Cf =    0;
            m =     1;
            e =     1;
        }

        void DruckerPragerParameters(double& Af, double&Bf, double& Cf, double&m, double& e, double ft, double fc)
        {
            Af =    0;
            Bf =    sqrt3_8 * (fc + ft)/(fc*ft);
            Cf =    3./2 * (fc-ft)/(fc*ft);
            m =     1;
            e =     1;
        }
        void MohrCoulombParameters(double& Af, double&Bf, double& Cf, double&m, double& e, double ft, double fc)
        {
            Af =    0;
            Bf =    1./sqrt6 * (fc +2*ft)/(fc*ft);
            Cf =    1./sqrt3 * (fc-ft)/(fc*ft);
            m =     1;
            e =     (fc+2*ft)/(2*fc +ft);
        }
    }
}
