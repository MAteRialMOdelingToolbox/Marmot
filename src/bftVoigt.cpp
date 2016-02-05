#include "bftVoigt.h"
#include "bftMechanics.h"
#include "bftConstants.h"

namespace bft{
    namespace Vgt{

        using namespace Constants;

        const Matrix6& P()
        {
            const static Matrix6 P(P_data);
            return P;
        }

        const Matrix6& PInv()
        {
            const static Matrix6 P(PInv_data);
            return P;
        }


        const Vector6& I()
        {
            const static Vector6 I(I_data);
            return I;
        }
        const Vector6& Ihyd()
        {
            const static Vector6 Ihyd(Ihyd_data);
            return Ihyd;
        }


        const Matrix6& Idev()
        {
            const static Matrix6 Idev(Idev_data);
            return Idev;
        }

        Matrix3d voigtToStrain(const Vector6& voigt)
        {
            Matrix3d strain;
            strain <<       
                    voigt[0],   voigt[3]/2, voigt[4]/2,
                    voigt[3]/2, voigt[1],   voigt[5]/2,
                    voigt[4]/2, voigt[5]/2, voigt[2];
            return strain;
        }

        Matrix3d voigtToStress(const Vector6& voigt)
        {
            Matrix3d stress;
            stress <<       
                    voigt[0],   voigt[3], voigt[4],
                    voigt[3],   voigt[1], voigt[5],
                    voigt[4],   voigt[5], voigt[2];
            return stress;
        }

        Vector6 strainToVoigt(const Matrix3d& strainTensor)
        {
            Vector6 strain;
            strain <<       
                    strainTensor(0,0), strainTensor(1,1), strainTensor(2,2), 
                    2*strainTensor(0,1), 2*strainTensor(0,2), 2*strainTensor(1,2);
            return strain;
        }

        Vector6 stressToVoigt(const Matrix3d& stressTensor)
        {
            Vector6 stress;
            stress <<       
                    stressTensor(0,0), stressTensor(1,1), stressTensor(2,2), 
                    stressTensor(0,1), stressTensor(0,2), stressTensor(1,2);
            return stress;
        }
        Vector3d principalStrains(const Vector6& voigtStrain)
        {
            SelfAdjointEigenSolver<Matrix3d> es(voigtToStrain(voigtStrain));
            return es.eigenvalues();
        }

        Vector3d principalStresses(const Vector6& voigtStress)
        {
            SelfAdjointEigenSolver<Matrix3d> es(voigtToStress(voigtStress));
            return es.eigenvalues();
        }
        Vector3d haighWestergaard(const Vector6& stress)
        {
            Vector3d hw;
            const double J2_ = J2(stress);
            const double J3_ = J3(stress) ;
            hw(0) = I1(stress) / sqrt3;
            hw(1) = sqrt(2 * J2_);

            const double x =  3*(sqrt3/2) * J3_ / (pow(J2_, 3./2));
            if(x<= -1)
                hw(2) = 1./3 * Pi ;
            else if(x>= 1)
                hw(2) = 0;
            else if(x!=x)
                hw(2) = 1./3 * Pi;
            else
                hw(2) = 1./3 * acos(x);
            return hw;
        }


        double I1(const Vector6& stress)
        {
            return stress.head(3).sum();
        }

        double I2(const Vector6& stress)
        {
            const Vector6& s = stress;

            return s(0)*s(1) + s(1)*s(2) + s(2)*s(0)
                    -s(3)*s(3) -s(4)*s(4) -s(5)*s(5);
        }

        double I3(const Vector6& stress)
        {
            const Vector6& s = stress;
            return  s(0)*s(1)*s(2)
                    +2*s(3)*s(4)*s(5)
                    -s(0)*s(5)*s(5)
                    -s(1)*s(4)*s(4)
                    -s(2)*s(3)*s(3);
        }

        double J2(const Vector6& stress)
        {
            double I1_ = I1(stress);
            double I2_ = I2(stress);
            double res = (1./3)* I1_*I1_ - I2_;
            return res >= 0 ? res : 0.0;
        }

        double J3(const Vector6& stress)
        {
            double I1_ = I1(stress);
            double I2_ = I2(stress);
            double I3_ = I3(stress);

            return (2./27) * pow(I1_,3)     - (1./3)*I1_*I2_        + I3_;
        }

        Vector6 dSigmaMdSigma()
        {
            return 1./3 * I();
        }

        Vector6 dRhodSigma(double rho, const Vector6& stress)
        {

            if(rho <= 1e-10)
                return Vector6::Zero();

            Vector6 s = Idev() * stress;

            return 1./rho  * P() * s;
        }

        Vector6 dThetadSigma(double theta, const Vector6& stress)
        {
            if(theta <= 0 || theta >= Pi/3)
                return Vector6::Zero();

            const double J2_ = J2(stress);
            const double J3_ = J3(stress);

            const double cos2_3theta = std::cos(3*theta) * std::cos(3*theta);
            const double dThetadJ2 = 3*sqrt3/4 * J3_/(std::pow(J2_, 2.5) * std::sqrt(1.0 - cos2_3theta)); 
            const double dThetadJ3 = - sqrt3/2 * 1./ (std::pow(J2_, 1.5) * std::sqrt(1.0 - cos2_3theta));

            if(isNaN(dThetadJ2) || isNaN(dThetadJ3))
                return Vector6::Zero(); 


            Vector6 s_ = Idev() * stress;
            const Matrix3d s = voigtToStress (s_);

            return dThetadJ2*P()*s_ + dThetadJ3*(P() * stressToVoigt(s*s) - 2./3 * J2_ *I());

        }

    }
}
