#pragma once
#include "MarmotFunctions.h"
#include "MarmotTypedefs.h"

#define isNaN( x ) ( x != x )
#define VOIGTFROMDIM( x ) ( ( ( x * x ) + x ) >> 1 )
//#define DIMFROMVOIGT( x ) ( x<<1  )

namespace marmot {
    namespace Vgt {

        // TODO: Remove, only valid for 3D!
        const int VoigtSize = 6;

        extern const marmot::Vector6 P;
        extern const marmot::Vector6 PInv;

        extern const marmot::Vector6 I;
        extern const marmot::Vector6 IHyd;
        extern const Matrix6      IDev;

        // Plane Stress handling
        Eigen::Vector3d voigtToPlaneVoigt( const marmot::Vector6& voigt );
        marmot::Vector6    planeVoigtToVoigt( const Eigen::Vector3d& voigtPlane );

        template <int voigtSize>
        Eigen::Matrix<double, voigtSize, 1> reduce3DVoigt( const marmot::Vector6& Voigt3D )
        {
            if constexpr ( voigtSize == 1 )
                return ( Eigen::Matrix<double, 1, 1>() << Voigt3D( 0 ) ).finished();
            else if constexpr ( voigtSize == 3 )
                return voigtToPlaneVoigt( Voigt3D );
            else if constexpr ( voigtSize == 6 )
                return Voigt3D;
            else
                throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid dimension specified" );
        }

        template <int voigtSize>
        marmot::Vector6 make3DVoigt( const Eigen::Matrix<double, voigtSize, 1>& Voigt )
        {
            if constexpr ( voigtSize == 1 )
                return ( marmot::Vector6() << Voigt( 0 ), 0, 0, 0, 0, 0 ).finished();
            else if constexpr ( voigtSize == 3 )
                return planeVoigtToVoigt( Voigt );
            else if constexpr ( voigtSize == 6 )
                return Voigt;
            else
                throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid dimension specified" );
        }

        /*compute E33 for a given elastic strain, to compute the compensation for
         * planeStress = Cel : (elasticStrain + compensationStrain) */
        marmot::Vector6 planeStressCompensationStrain( const marmot::Vector6& elasticStrain, double nu );
        /* Returns the transformation Matrix T which fullfills
         * planeStressIncrement = C : (T : arbitraryStrainIncrement) */
        Matrix6                     planeStressTangentTransformationMatrix( const Matrix6& tangent );
        Eigen::Matrix<double, 6, 3> dStrainDStrainPlaneStress( const Matrix6& tangent );
        Eigen::Matrix<double, 3, 6> dStressPlaneStressDStress();
        Eigen::Matrix<double, 6, 3> dStrainDStrainPlaneStrain();

        // function prototypes for  marmot::Vector6 handling
        Eigen::Matrix3d voigtToStrain( const marmot::Vector6& strainVector );
        Eigen::Matrix3d voigtToStress( const marmot::Vector6& stressVector );
        marmot::Vector6    strainToVoigt( const Eigen::Matrix3d& strainTensor );
        marmot::Vector6    stressToVoigt( const Eigen::Matrix3d& stressTensor );
        Eigen::Vector3d haighWestergaard( const marmot::Vector6& stress );
        Eigen::Vector3d haighWestergaardStrain( const marmot::Vector6& strain );

        template <int nDim>
        Eigen::Matrix<double, nDim, nDim> StressMatrixFromVoigt(
            const Eigen::Matrix<double, VOIGTFROMDIM( nDim ), 1>& Voigt )
        {
            if constexpr ( nDim == 1 )
                return ( Eigen::Matrix<double, nDim, nDim>() << Voigt( 0 ) ).finished();
            else if constexpr ( nDim == 2 )
                return ( Eigen::Matrix<double, nDim, nDim>() << Voigt( 0 ), Voigt( 2 ), Voigt( 2 ), Voigt( 1 ) )
                    .finished();
            else if constexpr ( nDim == 3 )
                return voigtToStress( Voigt );
            else
                throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid dimension specified" );
        }

        template <int nDim>
        Eigen::Matrix<double, VOIGTFROMDIM( nDim ), 1> VoigtFromStrainMatrix(
            const Eigen::Matrix<double, nDim, nDim>& strain )
        {
            if constexpr ( nDim == 1 )
                return ( Eigen::Matrix<double, VOIGTFROMDIM( nDim ), 1>() << strain( 0, 0 ) ).finished();
            else if constexpr ( nDim == 2 )
                return ( Eigen::Matrix<double, VOIGTFROMDIM( nDim ), 1>() << strain( 0, 0 ),
                         strain( 1, 1 ),
                         2 * strain( 0, 1 ) )
                    .finished();
            else if constexpr ( nDim == 3 )
                return strainToVoigt( strain );
            else
                throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid dimension specified" );
        }

        // principal strains calculated by solving eigenvalue problem ( !NOT sorted! )
        Eigen::Vector3d principalStrains( const marmot::Vector6& strain );
        // principal strains calculated from haigh westergaard strains ( sorted --> e1 > e2 > e3 )
        Eigen::Vector3d principalStrainsHW( const marmot::Vector6& strain );
        // principal stresses calculated by solving eigenvalue problem ( !NOT sorted! )
        Eigen::Vector3d principalStresses( const marmot::Vector6& stress );
        // principal stressDirections calculated by solving eigenvalue problem ( !NOT sorted! )
        Eigen::Matrix3d principalStressesDirections( const marmot::Vector6& stress );
        // rotate a 2nd order tensor T in voigt notation by : T' = Q * T * Q^T
        marmot::Vector6 rotateVoigtStress( const Eigen::Matrix3d Q, const marmot::Vector6& stress );

        // equivalent von Mises stress
        double vonMisesEquivalentStress( const marmot::Vector6& stress );
        // equivalent von Mises strain
        double vonMisesEquivalentStrain( const marmot::Vector6& strain );
        // Euclidian norm of strain
        double normStrain( const marmot::Vector6& strain );
        // Euclidian norm of stress
        double normStress( const marmot::Vector6& stress );
        // Trace of compressive strains
        double Evolneg( const marmot::Vector6& strain );

        // Invariants - keep atention: different for stress/strain tensor
        double I1( const marmot::Vector6& stress );
        double I2( const marmot::Vector6& stress );
        double I2strain( const marmot::Vector6& strain );
        double I3( const marmot::Vector6& stress );
        double I3strain( const marmot::Vector6& strain );

        // Invariants of the deviatoric part of the stress/strain tensor
        double J2( const marmot::Vector6& stress );
        double J2strain( const marmot::Vector6& strain );
        double J3( const marmot::Vector6& stress );
        double J3strain( const marmot::Vector6& strain );

        // derivatives of Haigh Westergaard stresses with respect to cauchy stress in eng. notation
        marmot::Vector6 dSigmaMdSigma();
        marmot::Vector6 dRhodSigma( double rho, const marmot::Vector6& stress );
        marmot::Vector6 dThetadSigma( double theta, const marmot::Vector6& stress );
        // derivatives of Haigh Westergaard stresses with respect to deviatoric invariants
        double dTheta_dJ2( const marmot::Vector6& stress );
        double dTheta_dJ3( const marmot::Vector6& stress );
        // derivatives of Haigh Westergaard strains with respect to deviatoric invariants
        double dThetaE_dJ2E( const marmot::Vector6& strain );
        double dThetaE_dJ3E( const marmot::Vector6& strain );
        // derivatives of deviatoric invariants with respect to eng. stresses
        marmot::Vector6 dJ2_dStress( const marmot::Vector6& stress );
        marmot::Vector6 dJ3_dStress( const marmot::Vector6& stress );
        // derivatives of deviatoric invariants with respect to eng. strains
        marmot::Vector6 dJ2E_dE( const marmot::Vector6& strain );
        marmot::Vector6 dJ3E_dE( const marmot::Vector6& strain );
        // derivatives of Haigh Westergaard strains with respect to eng. strains
        marmot::Vector6 dThetaE_dE( const marmot::Vector6& strain );

        // derivatives of principalStess with respect to stress
        marmot::Matrix36 dStressPrincipals_dStress( const marmot::Vector6& stress );

        // derivatives of plastic strains with respect to strains
        Eigen::Vector3d dDeltaEpvneg_dDeltaEpPrincipals( const marmot::Vector6& strain );
        Matrix6         dEp_dE( const Matrix6& CelInv, const Matrix6& Cep );
        RowVector6d     dDeltaEpv_dE( const Matrix6& CelInv, const Matrix6& Cep );
        marmot::Matrix36   dDeltaEpPrincipals_dDeltaEp( const marmot::Vector6& dEp );
        RowVector6d     dDeltaEpvneg_dE( const marmot::Vector6& dEp, const Matrix6& CelInv, const Matrix6& Cep );

    } // namespace Vgt

    namespace mechanics {
        template <typename T>
        int sgn( T val )
        {
            return ( T( 0 ) < val ) - ( val < T( 0 ) );
        }

        Matrix6         Cel( double E, double nu );
        Matrix6         CelInverse( double E, double nu );
        Eigen::Matrix3d getPlaneStressTangent( const Matrix6& C );

        Eigen::Matrix3d getPlaneStrainTangent( const Matrix6& C );
        double          getUniaxialStressTangent( const Eigen::Ref<const Matrix6>& C );
        double          E( const double K, const double G );
        double          nu( const double K, const double G );

        namespace PlaneStrain {
            Tensor322d dStressdDeformationGradient( const Tensor633d& dStressdDeformationGradient3D );
        }
        namespace PlaneStress {
            Tensor322d dStressdDeformationGradient( const Tensor633d& dStressdDeformationGradient3D );
        }
    } // namespace mechanics
} // namespace marmot
