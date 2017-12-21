#pragma once
#include "bftTypedefs.h"

#define isNaN(x) (x!=x)

namespace bft
{
	//****************************************************
	namespace mechanics
	{
        template <typename T> 
            int sgn(T val) {
                return (T(0) < val) - (val < T(0));
            }

        Matrix6 Cel(double E, double nu);
        Matrix6 CelInverse(double E, double nu);
        double macauly(double scalar);
        int heaviside(double scalar);
        Matrix3d getPlaneStressTangent(const Matrix6& C);
        Matrix3d getPlaneStrainTangent(const Matrix6& C);
    }
	//****************************************************
    namespace Vgt{

        const int VoigtSize = 6;
       
        extern const Vector6 P;
        extern const Vector6 PInv;

        extern const Vector6 I;
        extern const Vector6 IHyd;
        extern const Matrix6 IDev;

        // Plane Stress handling
        Vector3d voigtToPlaneVoigt(const Vector6& voigt);
        Vector6 planeVoigtToVoigt(const Vector3d& voigtPlane);

        /*compute E33 for a given elastic strain, to compute the compensation for 
         * planeStress = Cel : (elasticStrain + compensationStrain) */
        Vector6 planeStressCompensationStrain(const Vector6& elasticStrain, double nu);
        /* Returns the transformation Matrix T which fullfills
         * planeStressIncrement = C : (T : arbitraryStrainIncrement) */
        Matrix6 planeStressTangentTransformationMatrix(const Matrix6& tangent);
        Matrix<double, 6, 3> dStrainDStrainPlaneStress(const Matrix6& tangent);
        Matrix<double, 3, 6> dStressPlaneStressDStress();
        Matrix<double, 6, 3> dStrainDStrainPlaneStrain();
        
        // function prototypes for Vector6 handling
        Matrix3d voigtToStrain(const Vector6& strainVector);
        Matrix3d voigtToStress(const Vector6& stressVector);
        Vector6 strainToVoigt(const Matrix3d& strainTensor);
        Vector6 stressToVoigt(const Matrix3d& stressTensor);
		Vector3d haighWestergaard(const Vector6& stress);
		Vector3d haighWestergaardStrain(const Vector6& strain);

		// principal strains calculated by solving eigenvalue problem ( !NOT sorted! )
        Vector3d principalStrains(const Vector6& strain);
		// principal strains calculated from haigh westergaard strains ( sorted --> e1 > e2 > e3 )
		Vector3d principalStrainsHW(const Vector6& strain);
		// principal stresses calculated by solving eigenvalue problem ( !NOT sorted! )
        Vector3d principalStresses(const Vector6& stress);

		// Euclidian norm of strain
		double normStrain(const Vector6& strain);
		// Trace of compressive strains
		double Evolneg(const Vector6& strain);

		// Invariants - keep atention: different for stress/strain tensor
        double I1(const Vector6& stress);
        double I2(const Vector6& stress);
		double I2strain(const Vector6& strain);
        double I3(const Vector6& stress);
		double I3strain(const Vector6& strain);

		// Invariants of the deviatoric part of the stress/strain tensor
        double J2(const Vector6& stress);
		double J2strain(const Vector6& strain);
        double J3(const Vector6& stress);
		double J3strain(const Vector6& strain);

		// derivatives of Haigh Westergaard stresses with respect to cauchy stress in eng. notation
        Vector6 dSigmaMdSigma();
        Vector6 dRhodSigma(double rho, const Vector6& stress);
        Vector6 dThetadSigma(double theta, const Vector6& stress);
		// derivatives of Haigh Westergaard stresses with respect to deviatoric invariants
		double dTheta_dJ2(const Vector6& stress);
		double dTheta_dJ3(const Vector6& stress);
		// derivatives of Haigh Westergaard strains with respect to deviatoric invariants
		double dThetaE_dJ2E(const Vector6& strain);
		double dThetaE_dJ3E(const Vector6& strain);
		
        // derivatives of deviatoric invariants with respect to eng. stresses
		Vector6 dJ2_dStress(const Vector6& stress);
		Vector6 dJ3_dStress(const Vector6& stress);
		// derivatives of deviatoric invariants with respect to eng. strains
		Vector6 dJ2E_dE(const Vector6& strain);
		Vector6 dJ3E_dE(const Vector6& strain);
		// derivatives of Haigh Westergaard strains with respect to eng. strains
		Vector6 dThetaE_dE(const Vector6& strain);

		// derivatives of plastic strains with respect to strains
		Vector3d dDeltaEpvneg_dDeltaEpPrincipals(const Vector6& strain);
		Matrix6 dEp_dE(const Matrix6& CelInv, const Matrix6& Cep);
		RowVector6d dDeltaEpv_dE(const Matrix6& CelInv, const Matrix6& Cep);
		Matrix36 dDeltaEpPrincipals_dDeltaEp(const Vector6& dEp);
		RowVector6d dDeltaEpvneg_dE(const Vector6& dEp, const Matrix6& CelInv, const Matrix6& Cep);
    }

}
