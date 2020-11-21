#include "MarmotVoigt.h"
#include "MarmotConstants.h"
#include "MarmotMath.h"
#include "MarmotTensor.h"
#include <iostream>

using namespace Eigen;

namespace Marmot {
    namespace mechanics {

        Matrix6 Cel( double E, double nu )
        {
            Matrix6 Cel;
            // clang-format off
            Cel <<  (1-nu), nu, nu, 0, 0, 0,
                nu, (1-nu), nu, 0, 0, 0,
                nu, nu, (1-nu), 0, 0, 0,
                0, 0, 0, (1-2*nu)/2, 0, 0,
                0, 0, 0, 0, (1-2*nu)/2, 0,
                0, 0, 0, 0, 0, (1-2*nu)/2;
            // clang-format on
            Cel *= E / ( ( 1 + nu ) * ( 1 - 2 * nu ) );
            return Cel;
        }

        Matrix6 CelInverse( double E, double nu )
        {
            Matrix6      CelInv;
            const double G = E / ( 2 * ( 1 + nu ) );
            // clang-format off
            CelInv <<   1./E,   -nu/E,  -nu/E,  0,      0,      0,
                   -nu/E,  1./E,   -nu/E,  0,      0,      0,
                   -nu/E,  -nu/E,  1./E,   0,      0,      0,
                   0,      0,      0,      1./G,   0,      0,
                   0,      0,      0,      0,      1./G,   0,
                   0,      0,      0,      0,      0,      1./G;
            // clang-format on
            return CelInv;
        }

        Matrix3d getPlaneStressTangent( const Matrix6& C )
        {
            return Vgt::dStressPlaneStressDStress() * C * Vgt::dStrainDStrainPlaneStress( C );
        }

        double getUniaxialStressTangent( const Ref<const Matrix6>& C )
        {
            Vector6 b;
            b << 1, 0, 0, 0, 0, 0;
            Matrix6 A = C;
            A.row( 0 ) << 1, 0, 0, 0, 0, 0;

            Vector6 dEdEUniaxial = A.colPivHouseholderQr().solve( b );

            return C.row( 0 ) * dEdEUniaxial;
        }

        Matrix3d getPlaneStrainTangent( const Matrix6& C )
        {
            Matrix3d CPlaneStrain              = Matrix3d::Zero();
            CPlaneStrain.topLeftCorner( 2, 2 ) = C.topLeftCorner( 2, 2 );
            CPlaneStrain( 2, 2 )               = C( 3, 3 );
            CPlaneStrain.block<1, 2>( 2, 0 )   = C.block<1, 2>( 3, 0 );
            CPlaneStrain.block<2, 1>( 0, 2 )   = C.block<2, 1>( 0, 3 );

            return CPlaneStrain;
        }

        namespace PlaneStrain {

            Tensor322d dStressdDeformationGradient( const Tensor633d& dStressdDeformationGradient3D )
            {
                static constexpr int planeVoigtIndices[] = {0, 1, 3};
                Tensor322d           tangent2D;
                for ( int i = 0; i < 3; i++ )
                    for ( int j = 0; j < 2; j++ )
                        for ( int k = 0; k < 2; k++ )
                            tangent2D( i, j, k ) = dStressdDeformationGradient3D( planeVoigtIndices[i], j, k );

                return tangent2D;
            }
        } // namespace PlaneStrain

        namespace PlaneStress {

            Tensor322d dStressdDeformationGradient( const Tensor633d& dS_dF_3D )
            {
                /*  dS^PS    dS^PS    dS    dF^Comp
                 *  ----- == ----- * ---- * -------
                 *  dF^PS    dS       dF    dF^PS
                 *
                 *  dF^Comp_ij = dF_ij^PS       if  ij =/= 33
                 *             = dF_^Comp_33    else
                 *
                 *             with dS_33 = 0 =  dS33/dFij * dFij^PS + dS33/dF33 * dF^Comp_33
                 *
                 *             --> dF^Comp_33 = -dS33/dFij * dFij^PS * dF33/dS33
                 *
                 *             --> dF^Comp
                 *                 -------
                 *                  dF^PS
                 * */

                // projection to plane
                Tensor322d dS_dF = PlaneStrain::dStressdDeformationGradient( dS_dF_3D );

                using namespace TensorUtility::IndexNotation;

                // clang-format off
                for ( int m = 0; m < 2; m ++ )
                    for ( int n = 0; n < 2; n ++ )
                        for ( int k = 0; k < 2; k ++ )
                            for ( int l = 0; l < 2; l ++ )
                                dS_dF(              toVoigt<2> (m,n), k, l ) 
                                    =  -  dS_dF_3D( toVoigt<3> (m,n), 2, 2 ) 
                                    * 1./ dS_dF_3D( toVoigt<3> (2,2), 2, 2 )
                                    *     dS_dF_3D( toVoigt<3> (2,2), k, l );
                // clang-format on
                return dS_dF;
            }

        } // namespace PlaneStress

        double E( const double K, const double G ) { return 9. * K * G / ( 3. * K + G ); }

        double nu( const double K, const double G ) { return ( 3 * K - 2 * G ) / ( 6 * K + 2 * G ); }
    } // namespace mechanics
    namespace Vgt {

        using namespace Constants;

        const Vector6 P    = ( Vector6() << 1, 1, 1, 2, 2, 2 ).finished();
        const Vector6 PInv = ( Vector6() << 1, 1, 1, .5, .5, .5 ).finished();

        const Vector6 I    = ( Vector6() << 1, 1, 1, 0, 0, 0 ).finished();
        const Vector6 IHyd = ( Vector6() << 1. / 3, 1. / 3, 1. / 3, 0, 0, 0 ).finished();

        const Matrix6 IDev = ( Matrix6() <<
                                   // clang-format off
                2./3,    -1./3,   -1./3,    0,  0,  0,
                -1./3,   2./3,    -1./3,    0,  0,  0,
                -1./3,   -1./3,   2./3,     0,  0,  0,
                0,          0,      0,      1,  0,  0,
                0,          0,      0,      0,  1,  0,
                0,          0,      0,      0,  0,  1).finished();
        // clang-format on

        Matrix3d voigtToStrain( const Vector6& voigt )
        {
            Matrix3d strain;
            strain << voigt[0], voigt[3] / 2, voigt[4] / 2, voigt[3] / 2, voigt[1], voigt[5] / 2, voigt[4] / 2,
                voigt[5] / 2, voigt[2];
            return strain;
        }

        Matrix3d voigtToStress( const Vector6& voigt )
        {
            Matrix3d stress;
            // clang-format off
            stress << voigt[0], voigt[3], voigt[4], 
                      voigt[3], voigt[1], voigt[5], 
                      voigt[4], voigt[5], voigt[2];
            // clang-format on
            return stress;
        }

        Vector6 strainToVoigt( const Matrix3d& strainTensor )
        {
            Vector6 strain;
            strain << strainTensor( 0, 0 ), strainTensor( 1, 1 ), strainTensor( 2, 2 ), 2 * strainTensor( 0, 1 ),
                2 * strainTensor( 0, 2 ), 2 * strainTensor( 1, 2 );
            return strain;
        }

        Vector6 stressToVoigt( const Matrix3d& stressTensor )
        {
            Vector6 stress;
            stress << stressTensor( 0, 0 ), stressTensor( 1, 1 ), stressTensor( 2, 2 ), stressTensor( 0, 1 ),
                stressTensor( 0, 2 ), stressTensor( 1, 2 );
            return stress;
        }

        Vector3d voigtToPlaneVoigt( const Vector6& voigt )
        {
            /* converts a 6d voigt Vector with Abaqus notation
               S11, S22, S33, S12, S13, S23
               to a 3d Vector S11, S22, S12 */
            Vector3d voigtPlane;
            voigtPlane << voigt[0], voigt[1], voigt[3];
            return voigtPlane;
        }

        Vector6 planeVoigtToVoigt( const Vector3d& voigtPlane )
        {
            /* converts a 3d voigt Vector with notation
               S11, S22, S12 to a Vector6 with
               S11, S22, S33, S12, S13, S23
               !!! Don't use if 3rd component is NOT ZERO !!!*/
            Vector6 voigt;
            voigt << voigtPlane[0], voigtPlane[1], 0, voigtPlane[2], 0, 0;
            return voigt;
        }

        Vector6 planeStressCompensationStrain( const Vector6& strain, double nu )
        {
            /*compute E33 for a given strain, to compute the compensation for
             * planeStress = Cel : (strain + compensationStrain) */
            const Vector6& e                = strain;
            const double   strainCorrComp33 = -nu / ( 1 - nu ) * ( e( 0 ) + e( 1 ) ) - e( 2 );
            Vector6        result;
            result << 0, 0, strainCorrComp33, 0, 0, 0;
            return result;
        }

        Matrix6 planeStressTangentTransformationMatrix( const Matrix6& tangent )
        {
            /* Returns the transformation Matrix T which fullfills
             * planeStressIncrement = C : T * strainIncrement
             * for isotropic material behavior only!
             * */
            Matrix6      T   = Matrix6::Identity();
            const double t31 = -tangent( 2, 0 ) / tangent( 2, 2 );
            const double t32 = -tangent( 2, 1 ) / tangent( 2, 2 );
            T.row( 2 ).head( 3 ) << t31, t32, 0.0;
            return T;
        }

        Matrix<double, 6, 3> dStrainDStrainPlaneStress( const Matrix6& tangent )
        {
            Matrix<double, 6, 3> T = Matrix<double, 6, 3>::Zero();
            T( 0, 0 )              = 1;
            T( 1, 1 )              = 1;
            T( 3, 2 )              = 1;
            T( 2, 0 )              = -tangent( 2, 0 ) / tangent( 2, 2 );
            T( 2, 1 )              = -tangent( 2, 1 ) / tangent( 2, 2 );
            return T;
        }

        Matrix<double, 6, 3> dStrainDStrainPlaneStrain()
        {
            Matrix<double, 6, 3> T = Matrix<double, 6, 3>::Zero();
            T( 0, 0 )              = 1;
            T( 1, 1 )              = 1;
            T( 3, 2 )              = 1;
            return T;
        }

        Matrix<double, 3, 6> dStressPlaneStressDStress()
        {
            Matrix<double, 3, 6> T = Matrix<double, 3, 6>::Zero();
            T( 0, 0 )              = 1;
            T( 1, 1 )              = 1;
            T( 2, 3 )              = 1;
            return T;
        }

        Vector3d haighWestergaard( const Vector6& stress )
        {
            Vector3d     hw;
            const double J2_ = J2( stress );
            hw( 0 )          = I1( stress ) / sqrt3;
            hw( 1 )          = sqrt( 2 * J2_ );

            if ( hw( 1 ) != 0 ) {
                const double J3_ = J3( stress );
                const double x   = 3 * ( sqrt3 / 2 ) * J3_ / ( pow( J2_, 3. / 2 ) );
                if ( x <= -1 )
                    hw( 2 ) = 1. / 3 * Pi;
                else if ( x >= 1 )
                    hw( 2 ) = 0;
                else if ( x != x )
                    hw( 2 ) = 1. / 3 * Pi;
                else
                    hw( 2 ) = 1. / 3 * acos( x );
            }
            else
                hw( 2 ) = 0;

            return hw;
        }

        Vector3d haighWestergaardStrain( const Vector6& strain )
        {
            Vector3d     hw;
            const double J2_ = J2strain( strain );

            hw( 0 ) = I1( strain ) / sqrt3; // will be changed to eM if used
            hw( 1 ) = sqrt( 2 * J2_ );

            if ( hw( 1 ) != 0 ) {
                const double J3_ = J3strain( strain );
                const double x   = 3. * ( sqrt3 / 2. ) * J3_ / ( pow( J2_, 3. / 2 ) );
                if ( x <= -1 )
                    hw( 2 ) = 1. / 3 * Constants::Pi;
                else if ( x >= 1 )
                    hw( 2 ) = 0;
                else if ( x != x )
                    hw( 2 ) = 1. / 3 * Constants::Pi;
                else
                    hw( 2 ) = 1. / 3 * acos( x );
            }
            else
                hw( 2 ) = 0;

            return hw;
        }

        Vector3d principalStrains( const Vector6& voigtStrain )
        {
            SelfAdjointEigenSolver<Matrix3d> es( voigtToStrain( voigtStrain ) );
            return es.eigenvalues();
        }

        Vector3d principalStrainsHW( const Vector6& voigtStrain )
        {
            // if you wanna sort your eigenvalues after size
            Vector3d hw = haighWestergaardStrain( voigtStrain );
            Vector3d strainPrinc;
            strainPrinc << hw( 0 ) / sqrt3 + sqrt2_3 * hw( 1 ) * std::cos( hw( 2 ) ),
                hw( 0 ) / sqrt3 + sqrt2_3 * hw( 1 ) * ( -std::sin( Constants::Pi / 6. - hw( 2 ) ) ),
                hw( 0 ) / sqrt3 + sqrt2_3 * hw( 1 ) * ( -std::sin( Constants::Pi / 6. + hw( 2 ) ) );
            return strainPrinc;
        }

        Vector3d principalStresses( const Vector6& voigtStress )
        {
            SelfAdjointEigenSolver<Matrix3d> es( voigtToStress( voigtStress ) );
            return es.eigenvalues();
        }
        Eigen::Matrix3d principalStressesDirections( const Marmot::Vector6& voigtStress )
        {
            SelfAdjointEigenSolver<Matrix3d> es( voigtToStress( voigtStress ) );
            Matrix3d Q = es.eigenvectors();  Q.col(2) = Q.col(0).cross(Q.col(1)); // for a clockwise coordinate system
            return Q; 
        }
        Marmot::Vector6 rotateVoigtStress( const Eigen::Matrix3d Q, const Marmot::Vector6& voigtStress )
        {
            const Matrix3d& T  = voigtToStress( voigtStress );
            const Matrix3d& TR = Q * T * Q.transpose();
            return stressToVoigt( TR );
        }

        double vonMisesEquivalentStress( const Vector6& stress ) { return sqrt( 3. * J2( stress ) ); }

        double vonMisesEquivalentStrain( const Vector6& strain )
        {
            // e_eq = sqrt( 2/3 * e_ij * e_ij )
            const Vector6& e = strain;
            return std::sqrt( 2. / 3. * ( e( 0 ) * e( 0 ) + e( 1 ) * e( 1 ) + e( 2 ) * e( 2 ) ) +
                              1. / 3. * ( e( 3 ) * e( 3 ) + e( 4 ) * e( 4 ) + e( 5 ) * e( 5 ) ) );
        }

        double normStrain( const Vector6& strain ) { return Vgt::voigtToStrain( strain ).norm(); }

        double normStress( const Vector6& stress ) { return Vgt::voigtToStress( stress ).norm(); }

        double Evolneg( const Vector6& strain )
        {
            Vector3d dEpPrincipal = principalStrains( strain );

            return Math::macauly( -dEpPrincipal( 0 ) ) + Math::macauly( -dEpPrincipal( 1 ) ) +
                   Math::macauly( -dEpPrincipal( 2 ) );
        }

        double I1( const Vector6& stress ) { return stress.head( 3 ).sum(); }

        double I2( const Vector6& stress )
        {
            const Vector6& s = stress;

            return s( 0 ) * s( 1 ) + s( 1 ) * s( 2 ) + s( 2 ) * s( 0 ) - s( 3 ) * s( 3 ) - s( 4 ) * s( 4 ) -
                   s( 5 ) * s( 5 );
        }

        double I2strain( const Vector6& strain ) // you could also use normal I2, but with epsilon12
                                                 // instead of 2*epsilon12
        {
            const Vector6& e = strain;

            return e( 0 ) * e( 1 ) + e( 1 ) * e( 2 ) + e( 2 ) * e( 0 ) - e( 3 ) / 2. * e( 3 ) / 2. -
                   e( 4 ) / 2. * e( 4 ) / 2. - e( 5 ) / 2. * e( 5 ) / 2.;
        }

        double I3( const Vector6& stress )
        {
            const Vector6& s = stress;
            return s( 0 ) * s( 1 ) * s( 2 ) + 2 * s( 3 ) * s( 4 ) * s( 5 ) - s( 0 ) * s( 5 ) * s( 5 ) -
                   s( 1 ) * s( 4 ) * s( 4 ) - s( 2 ) * s( 3 ) * s( 3 );
        }

        double I3strain( const Vector6& strain ) // you could also use normal I3, but with epsilon12
                                                 // instead of 2*epsilon12
        {
            return voigtToStrain( strain ).determinant();
        }

        double J2( const Vector6& stress )
        {
            double I1_ = I1( stress );
            double I2_ = I2( stress );
            double res = ( 1. / 3 ) * I1_ * I1_ - I2_;
            return res >= 0 ? res : 0.0;
        }

        double J2strain( const Vector6& strain )
        {
            const double res = 1. / 3. * std::pow( I1( strain ), 2. ) - I2strain( strain );
            return res > 0 ? res : 0;
        }

        double J3( const Vector6& stress )
        {
            double I1_ = I1( stress );
            double I2_ = I2( stress );
            double I3_ = I3( stress );

            return ( 2. / 27 ) * pow( I1_, 3 ) - ( 1. / 3 ) * I1_ * I2_ + I3_;
        }

        double J3strain( const Vector6& strain ) // determinant of the deviatoric strain tensor
        {
            return voigtToStrain( IDev * strain ).determinant();
        }

        Vector6 dSigmaMdSigma() { return 1. / 3 * I; }

        Vector6 dRhodSigma( double rho, const Vector6& stress )
        {

            if ( rho <= 1e-16 )
                return Vector6::Zero();

            Vector6 s = IDev * stress;

            return 1. / rho * P.array() * s.array();
        }

        Vector6 dThetadSigma( double theta, const Vector6& stress )
        {
            if ( theta <= 1e-15 || theta >= Pi / 3 - 1e-15 )
                return Vector6::Zero();

            // const double J2_ = J2(stress);
            // const double J3_ = J3(stress);

            const double dThetadJ2 = dTheta_dJ2( stress );
            const double dThetadJ3 = dTheta_dJ3( stress );

            if ( isNaN( dThetadJ2 ) || isNaN( dThetadJ3 ) )
                return Vector6::Zero();

            return dThetadJ2 * dJ2_dStress( stress ) + dThetadJ3 * dJ3_dStress( stress );
        }

        double dTheta_dJ2( const Vector6& stress )
        {
            const Vector3d hw    = haighWestergaard( stress );
            const double&  theta = hw( 2 );

            if ( theta <= 1e-14 || theta >= Pi / 3 - 1e-14 )
                return 1e16;

            const double J2_ = J2( stress );
            const double J3_ = J3( stress );

            const double cos2_3theta = std::cos( 3 * theta ) * std::cos( 3 * theta );
            const double dThetadJ2   = 3 * sqrt3 / 4 * J3_ / ( std::pow( J2_, 2.5 ) * std::sqrt( 1.0 - cos2_3theta ) );
            return dThetadJ2;
        }

        double dTheta_dJ3( const Vector6& stress )
        {
            const Vector3d hw    = haighWestergaard( stress );
            const double&  theta = hw( 2 );

            if ( theta <= 1e-14 || theta >= Pi / 3 - 1e-14 )
                return -1e16;

            const double J2_ = J2( stress );
            // const double J3_ = J3(stress);

            const double cos2_3theta = std::cos( 3 * theta ) * std::cos( 3 * theta );
            const double dThetadJ3   = -sqrt3 / 2. * 1. / ( std::pow( J2_, 1.5 ) * std::sqrt( 1.0 - cos2_3theta ) );
            return dThetadJ3;
        }

        double dThetaE_dJ2E( const Vector6& strain )
        {
            const Vector3d hw    = haighWestergaardStrain( strain );
            const double   theta = hw( 2 );

            if ( theta <= 1e-15 || theta >= Pi / 3 - 1e-15 )
                return 1e16;
            else
                return 3. * std::sqrt( 3. ) / 4. * J3strain( strain ) /
                       ( std::pow( J2strain( strain ), 5. / 2 ) *
                         std::sqrt( 1. - std::pow( std::cos( 3. * theta ), 2. ) ) );
        }

        double dThetaE_dJ3E( const Vector6& strain )
        {
            const Vector3d hw    = haighWestergaardStrain( strain );
            const double&  theta = hw( 2 );

            if ( theta <= 1e-15 || theta >= Pi / 3 - 1e-15 )
                return -1e16;
            else
                return -std::sqrt( 3. ) / 2. * 1. /
                       ( std::pow( J2strain( strain ), 3. / 2 ) *
                         std::sqrt( 1. - std::pow( std::cos( 3. * theta ), 2. ) ) );
        }

        Vector6 dJ2_dStress( const Vector6& stress ) { return P.array() * ( IDev * stress ).array(); }

        Vector6 dJ3_dStress( const Vector6& stress )
        {
            Vector6 s = IDev * stress;
            return ( P.array() * stressToVoigt( voigtToStress( s ) * voigtToStress( s ) ).array() ).matrix() -
                   2. / 3. * J2( stress ) * I;
        }

        Vector6 dJ2E_dE( const Vector6& strain ) { return PInv.array() * ( IDev * strain ).array(); }

        Vector6 dJ3E_dE( const Vector6& strain )
        {
            Vector6 e = IDev * strain;
            return ( PInv.array() * strainToVoigt( voigtToStrain( e ) * voigtToStrain( e ) ).array() ).matrix() -
                   2. / 3. * J2strain( strain ) * I;
        }

        Vector6 dThetaE_dE( const Vector6& strain )

        {
            return dThetaE_dJ2E( strain ) * dJ2E_dE( strain ) + dThetaE_dJ3E( strain ) * dJ3E_dE( strain );
        }

        Matrix36 dStressPrincipals_dStress(
            const Vector6& stress ) // derivative when principal stresses are computed from solving Eigenvalue-Problem
        {
            MatrixXd J( 3, 6 );

            Vector6 leftX;
            Vector6 rightX;

            for ( size_t i = 0; i < 6; i++ ) {
                double volatile h = std::max( 1.0, std::abs( stress( i ) ) ) * Constants::cubicRootEps();
                leftX             = stress;
                leftX( i ) -= h;
                rightX = stress;
                rightX( i ) += h;

                J.col( i ) = 1. / ( 2 * h ) * ( principalStresses( rightX ) - principalStresses( leftX ) );
            }
            return J;
        }

        Vector3d dDeltaEpvneg_dDeltaEpPrincipals( const Vector6& strain )
        {
            Vector3d       dEvdEpPrinc  = Vector3d::Zero();
            const Vector3d deltaEpPrinc = principalStrainsHW( strain );

            for ( int i = 0; i < dEvdEpPrinc.size(); i++ )
                dEvdEpPrinc( i ) = -Math::heaviside( -deltaEpPrinc( i ) );

            return dEvdEpPrinc;
        }

        Matrix6 dEp_dE( const Matrix6& CelInv, const Matrix6& Cep ) { return Matrix6::Identity() - CelInv * Cep; }

        RowVector6d dDeltaEpv_dE( const Matrix6& CelInv, const Matrix6& Cep )
        {
            return I.transpose() * ( Matrix6::Identity() - CelInv * Cep );
        }

        Matrix36 dDeltaEpPrincipals_dDeltaEp(
            const Vector6& dEp ) // equations from page 218-219 PhD Thesis David Unteregger
        {
            Vector3d dEpPrinc_dEpvol   = Vector3d::Zero();
            Vector3d dEpPrinc_dEprho   = Vector3d::Zero();
            Vector3d dEPprinc_dEptheta = Vector3d::Zero();

            const double   sqrt2_3 = std::sqrt( 2. / 3. );
            const Vector3d hw      = haighWestergaardStrain( dEp );
            // const double& epsM =		hw(0);
            const double& rhoE = hw( 1 );
            // const double& thetaE =		hw(2);

            dEpPrinc_dEpvol = 1. / 3. * Vector3d::Ones();
            dEpPrinc_dEprho << sqrt2_3 * std::cos( hw( 2 ) ), sqrt2_3 * std::cos( hw( 2 ) - 2. * Constants::Pi / 3. ),
                sqrt2_3 * std::cos( hw( 2 ) + 2. * Constants::Pi / 3. );

            dEPprinc_dEptheta << -sqrt2_3 * hw( 1 ) * std::sin( hw( 2 ) ),
                -sqrt2_3 * hw( 1 ) * std::sin( hw( 2 ) - 2. * Constants::Pi / 3. ),
                -sqrt2_3 * hw( 1 ) * std::sin( hw( 2 ) + 2. * Constants::Pi / 3. );

            RowVector6d dEpvol_dEp   = RowVector6d::Zero();
            dEpvol_dEp               = I;
            RowVector6d dEprho_dEp   = RowVector6d::Zero();
            RowVector6d dEptheta_dEp = RowVector6d::Zero();

            if ( std::abs( rhoE ) > 1e-16 ) {
                dEprho_dEp   = 1. / rhoE * dJ2E_dE( dEp ).transpose();
                dEptheta_dEp = ( dThetaE_dJ2E( dEp ) * dJ2E_dE( dEp ).transpose() ) +
                               ( dThetaE_dJ3E( dEp ) * dJ3E_dE( dEp ).transpose() );
            }
            else {
                dEprho_dEp << 1.e16, 1.e16, 1.e16, 1.e16, 1.e16,
                    1.e16; // 1e16 from Code David (Line 67, D_2_Umatsub_damage3_derivatives)
                dEptheta_dEp << 0., 0., 0., 0., 0., 0.;
            }

            Matrix36 dEpPrincdEpAna;
            dEpPrincdEpAna = ( dEpPrinc_dEpvol * dEpvol_dEp ) + ( dEpPrinc_dEprho * dEprho_dEp ) +
                             ( dEPprinc_dEptheta * dEptheta_dEp );

            return dEpPrincdEpAna;
        }

        RowVector6d dDeltaEpvneg_dE( const Vector6& dEp, const Matrix6& CelInv, const Matrix6& Cep )
        {
            return dDeltaEpvneg_dDeltaEpPrincipals( dEp ).transpose() * dDeltaEpPrincipals_dDeltaEp( dEp ) *
                   dEp_dE( CelInv, Cep );
        }

    } // namespace Vgt
} // namespace Marmot
