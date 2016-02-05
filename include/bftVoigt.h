#pragma once
#include "bftTypedefs.h"

namespace bft
{
    namespace Vgt{

        const int VoigtSize = 6;
        
        const double P_data[]= 
        {1,0,0,0,0,0,
        0,1,0,0,0,0,
        0,0,1,0,0,0,
        0,0,0,2,0,0,
        0,0,0,0,2,0,
        0,0,0,0,0,2};

        const double PInv_data[]= 
        {1,0,0,0,0,0,
        0,1,0,0,0,0,
        0,0,1,0,0,0,
        0,0,0,0.5,0,0,
        0,0,0,0,0.5,0,
        0,0,0,0,0,0.5};

        const Matrix6& P();
        const Matrix6& PInv();
        // I Vector
        const double I_data[]= {1,1,1,0,0,0};

        const Vector6& I();
        // hydrostatic Vector
        const double Ihyd_data[] = {1./3, 1./3, 1./3, 0, 0, 0};
        const Vector6& Ihyd();

        // Deviatoric Operator
        const double Idev_data[]= {
                2./3,    -1./3,   -1./3,    0,  0,  0,
                -1./3,   2./3,    -1./3,    0,  0,  0,
                -1./3,   -1./3,   2./3,     0,  0,  0,
                0,          0,      0,      1,  0,  0,
                0,          0,      0,      0,  1,  0,
                0,          0,      0,      0,  0,  1};
        const Matrix6& Idev();
        // function prototypes for Vector6 handling

        Matrix3d voigtToStrain(const Vector6& strainVector);
        Matrix3d voigtToStress(const Vector6& stressVector);
        Vector6 strainToVoigt(const Matrix3d& strainTensor);
        Vector6 stressToVoigt(const Matrix3d& stressTensor);
        Vector3d principalStrains(const Vector6& strain);
        Vector3d principalStresses(const Vector6& stress);
        Vector3d haighWestergaard(const Vector6& stress);

        double I1(const Vector6& stress);
        double I2(const Vector6& stress);
        double I3(const Vector6& stress);

        double J2(const Vector6& stress);
        double J3(const Vector6& stress);

        Vector6 dSigmaMdSigma();
        Vector6 dRhodSigma(double rho, const Vector6& stress);
        Vector6 dThetadSigma(double theta, const Vector6& stress);

    };

}
