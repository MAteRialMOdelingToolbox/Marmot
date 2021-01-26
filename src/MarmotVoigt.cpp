#include "Marmot/MarmotVoigt.h"
#include "Marmot/HaighWestergaard.h"
#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTensor.h"
#include <iostream>

using namespace Eigen;

namespace Marmot {
    namespace ContinuumMechanics::VoigtNotation {
        using namespace Constants;
        using namespace ContinuumMechanics::HaighWestergaard;

        const Vector6d P    = ( Vector6d() << 1, 1, 1, 2, 2, 2 ).finished();
        const Vector6d PInv = ( Vector6d() << 1, 1, 1, .5, .5, .5 ).finished();

        const Vector6d I    = ( Vector6d() << 1, 1, 1, 0, 0, 0 ).finished();
        const Vector6d IHyd = ( Vector6d() << 1. / 3, 1. / 3, 1. / 3, 0, 0, 0 ).finished();

        const Matrix6d IDev = ( Matrix6d() <<
                                    // clang-format off
        2./3,    -1./3,   -1./3,    0,  0,  0,
        -1./3,   2./3,    -1./3,    0,  0,  0,
        -1./3,   -1./3,   2./3,     0,  0,  0,
        0,          0,      0,      1,  0,  0,
        0,          0,      0,      0,  1,  0,
        0,          0,      0,      0,  0,  1).finished();
        // clang-format on

        Matrix3d voigtToStrain( const Vector6d& voigt )
        {
            Matrix3d strain;
            strain << voigt[0], voigt[3] / 2, voigt[4] / 2, voigt[3] / 2, voigt[1], voigt[5] / 2, voigt[4] / 2,
                voigt[5] / 2, voigt[2];
            return strain;
        }

        Matrix3d voigtToStress( const Vector6d& voigt )
        {
            Matrix3d stress;
            // clang-format off
            stress << voigt[0], voigt[3], voigt[4], 
                      voigt[3], voigt[1], voigt[5], 
                      voigt[4], voigt[5], voigt[2];
            // clang-format on
            return stress;
        }

        Vector6d strainToVoigt( const Matrix3d& strainTensor )
        {
            Vector6d strain;
            // clang-format off
            strain << strainTensor( 0, 0 ), 
                      strainTensor( 1, 1 ), 
                      strainTensor( 2, 2 ), 
                  2 * strainTensor( 0, 1 ),
                  2 * strainTensor( 0, 2 ), 
                  2 * strainTensor( 1, 2 );
            // clang-format on
            return strain;
        }

        Vector6d stressToVoigt( const Matrix3d& stressTensor )
        {
            Vector6d stress;
            // clang-format off
            stress << stressTensor( 0, 0 ), 
                      stressTensor( 1, 1 ), 
                      stressTensor( 2, 2 ), 
                      stressTensor( 0, 1 ),
                      stressTensor( 0, 2 ), 
                      stressTensor( 1, 2 );
            // clang-format on
            return stress;
        }

        Vector3d voigtToPlaneVoigt( const Vector6d& voigt )
        {
            /* converts a 6d voigt Vector with Abaqus notation
               S11, S22, S33, S12, S13, S23
               to a 3d Vector S11, S22, S12 */
            Vector3d voigtPlane;
            voigtPlane << voigt[0], voigt[1], voigt[3];
            return voigtPlane;
        }

        Vector6d planeVoigtToVoigt( const Vector3d& voigtPlane )
        {
            /* converts a 3d voigt Vector with notation
               S11, S22, S12 to a Vector6d with
               S11, S22, S33, S12, S13, S23
               !!! Don't use if 3rd component is NOT ZERO !!!*/
            Vector6d voigt;
            voigt << voigtPlane[0], voigtPlane[1], 0, voigtPlane[2], 0, 0;
            return voigt;
        }

        namespace Invariants {

            Vector3d principalStrains( const Vector6d& voigtStrain )
            {
                SelfAdjointEigenSolver< Matrix3d > es( voigtToStrain( voigtStrain ) );
                return es.eigenvalues();
            }

            Vector3d principalStresses( const Vector6d& voigtStress )
            {
                SelfAdjointEigenSolver< Matrix3d > es( voigtToStress( voigtStress ) );
                return es.eigenvalues();
            }
            Eigen::Matrix3d principalStressesDirections( const Marmot::Vector6d& voigtStress )
            {
                SelfAdjointEigenSolver< Matrix3d > es( voigtToStress( voigtStress ) );
                Matrix3d                           Q = es.eigenvectors();
                Q.col( 2 ) = Q.col( 0 ).cross( Q.col( 1 ) ); // for a clockwise coordinate system
                return Q;
            }

            std::pair< Eigen::Vector3d, Eigen::Matrix< double, 3, 6 > > principalsOfVoigtAndDerivatives(
                const Eigen::Matrix< double, 6, 1 >& S )
            {
                // This is a fast implementation of the classical algorithm for determining
                // the principal components of a symmetric 3x3 Matrix in Voigt notation (off diagonals expected with
                // factor 1) as well as its respective derivatives

                using namespace Eigen;

                using Vector6d = Matrix< double, 6, 1 >;
                using Matrix36 = Matrix< double, 3, 6 >;
                using Matrix66 = Matrix< double, 6, 6 >;

                const static auto dS0_dS = ( Vector6d() << 1, 0, 0, 0, 0, 0 ).finished();
                const static auto dS1_dS = ( Vector6d() << 0, 1, 0, 0, 0, 0 ).finished();
                const static auto dS2_dS = ( Vector6d() << 0, 0, 1, 0, 0, 0 ).finished();

                const double p1 = S( 3 ) * S( 3 ) + S( 4 ) * S( 4 ) + S( 5 ) * S( 5 );

                if ( p1 <= 1e-16 ) { // matrix is already diagonal
                    // clang-format off
                    const static auto dE_dS = ( Matrix36() << 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 ).finished();
                    return { S.head( 3 ), dE_dS };
                    // clang-format on
                }

                const Vector6d dP1_dS = { 0, 0, 0, 2 * S( 3 ), 2 * S( 4 ), 2 * S( 5 ) };

                const double   q     = S.head( 3 ).sum() / 3;
                const Vector6d dQ_dS = Marmot::ContinuumMechanics::VoigtNotation::I / 3;

                const double p2 = std::pow( S( 0 ) - q, 2 ) + std::pow( S( 1 ) - q, 2 ) + std::pow( S( 2 ) - q, 2 ) +
                                  2 * p1;
                const Vector6d dP2_dS = 2 * ( S( 0 ) - q ) * ( dS0_dS - dQ_dS ) +
                                        2 * ( S( 1 ) - q ) * ( dS1_dS - dQ_dS ) +
                                        2 * ( S( 2 ) - q ) * ( dS2_dS - dQ_dS ) + 2 * dP1_dS;

                const double   p     = std::sqrt( p2 / 6 );
                const Vector6d dP_dS = 1. / 12 * 1. / p * dP2_dS;

                const Vector6d B     = 1. / p * ( S - q * Marmot::ContinuumMechanics::VoigtNotation::I );
                const Matrix66 dB_dS = -B * 1. / p * dP_dS.transpose() +
                                       1. / p *
                                           ( Matrix66::Identity() -
                                             Marmot::ContinuumMechanics::VoigtNotation::I * dQ_dS.transpose() );

                const double detB = B( 0 ) * B( 1 ) * B( 2 ) + B( 3 ) * B( 4 ) * B( 5 ) * 2 - B( 2 ) * B( 3 ) * B( 3 ) -
                                    B( 1 ) * B( 4 ) * B( 4 ) - B( 0 ) * B( 5 ) * B( 5 );

                Vector6d dDetB_dB;
                dDetB_dB( 0 ) = B( 1 ) * B( 2 ) - B( 5 ) * B( 5 );
                dDetB_dB( 1 ) = B( 0 ) * B( 2 ) - B( 4 ) * B( 4 );
                dDetB_dB( 2 ) = B( 0 ) * B( 1 ) - B( 3 ) * B( 3 );
                dDetB_dB( 3 ) = B( 4 ) * B( 5 ) * 2 - B( 2 ) * B( 3 ) * 2;
                dDetB_dB( 4 ) = B( 3 ) * B( 5 ) * 2 - B( 1 ) * B( 4 ) * 2;
                dDetB_dB( 5 ) = B( 3 ) * B( 4 ) * 2 - B( 0 ) * B( 5 ) * 2;

                const double   r     = detB * 1. / 2;
                const Vector6d dR_dS = 1. / 2 * dDetB_dB.transpose() * dB_dS;

                double phi;
                double dPhi_dR;
                if ( r <= -1 ) {
                    phi     = Constants::Pi / 3;
                    dPhi_dR = 0.0;
                }
                else if ( r >= 1 ) {
                    phi     = 1.0;
                    dPhi_dR = 0.0;
                }
                else {
                    phi     = std::acos( r ) / 3;
                    dPhi_dR = 1. / 3 * -1. / ( std::sqrt( 1 - r * r ) );
                }
                const Vector6d dPhi_dS = dPhi_dR * dR_dS;

                Vector3d e;
                Matrix36 dE_dS;

                e( 0 )         = q + 2 * p * std::cos( phi );
                dE_dS.row( 0 ) = dQ_dS + 2 * ( dP_dS * std::cos( phi ) - p * std::sin( phi ) * dPhi_dS );

                const double phiShifted = phi + ( 2. / 3 * Constants::Pi );
                e( 2 )                  = q + 2 * p * std::cos( phiShifted );
                dE_dS.row( 2 ) = dQ_dS + 2 * ( dP_dS * std::cos( phiShifted ) - p * std::sin( phiShifted ) * dPhi_dS );

                e( 1 )         = 3 * q - e( 0 ) - e( 2 );
                dE_dS.row( 1 ) = 3 * dQ_dS - dE_dS.row( 0 ).transpose() - dE_dS.row( 2 ).transpose();

                return { e, dE_dS };
            }

            double vonMisesEquivalentStress( const Vector6d& stress ) { return sqrt( 3. * J2( stress ) ); }

            double vonMisesEquivalentStrain( const Vector6d& strain )
            {
                // e_eq = sqrt( 2/3 * e_ij * e_ij )
                const Vector6d& e = strain;
                return std::sqrt( 2. / 3. * ( e( 0 ) * e( 0 ) + e( 1 ) * e( 1 ) + e( 2 ) * e( 2 ) ) +
                                  1. / 3. * ( e( 3 ) * e( 3 ) + e( 4 ) * e( 4 ) + e( 5 ) * e( 5 ) ) );
            }

            double normStrain( const Vector6d& strain )
            {
                return ContinuumMechanics::VoigtNotation::voigtToStrain( strain ).norm();
            }

            double normStress( const Vector6d& stress )
            {
                return ContinuumMechanics::VoigtNotation::voigtToStress( stress ).norm();
            }

            double Evolneg( const Vector6d& strain )
            {
                Vector3d dEpPrincipal = principalStrains( strain );

                return Math::macauly( -dEpPrincipal( 0 ) ) + Math::macauly( -dEpPrincipal( 1 ) ) +
                       Math::macauly( -dEpPrincipal( 2 ) );
            }

            double I1( const Vector6d& stress ) { return stress.head( 3 ).sum(); }

            double I2( const Vector6d& stress )
            {
                const Vector6d& s = stress;

                return s( 0 ) * s( 1 ) + s( 1 ) * s( 2 ) + s( 2 ) * s( 0 ) - s( 3 ) * s( 3 ) - s( 4 ) * s( 4 ) -
                       s( 5 ) * s( 5 );
            }

            double I2strain( const Vector6d& strain ) // you could also use normal I2, but with epsilon12
                                                      // instead of 2*epsilon12
            {
                const Vector6d& e = strain;

                return e( 0 ) * e( 1 ) + e( 1 ) * e( 2 ) + e( 2 ) * e( 0 ) - e( 3 ) / 2. * e( 3 ) / 2. -
                       e( 4 ) / 2. * e( 4 ) / 2. - e( 5 ) / 2. * e( 5 ) / 2.;
            }

            double I3( const Vector6d& stress )
            {
                const Vector6d& s = stress;
                return s( 0 ) * s( 1 ) * s( 2 ) + 2 * s( 3 ) * s( 4 ) * s( 5 ) - s( 0 ) * s( 5 ) * s( 5 ) -
                       s( 1 ) * s( 4 ) * s( 4 ) - s( 2 ) * s( 3 ) * s( 3 );
            }

            double I3strain( const Vector6d& strain ) // you could also use normal I3, but with epsilon12
                                                      // instead of 2*epsilon12
            {
                return voigtToStrain( strain ).determinant();
            }

            double J2( const Vector6d& stress )
            {
                double I1_ = I1( stress );
                double I2_ = I2( stress );
                double res = ( 1. / 3 ) * I1_ * I1_ - I2_;
                return res >= 0 ? res : 0.0;
            }

            double J2strain( const Vector6d& strain )
            {
                const double res = 1. / 3. * std::pow( I1( strain ), 2. ) - I2strain( strain );
                return res > 0 ? res : 0;
            }

            double J3( const Vector6d& stress )
            {
                double I1_ = I1( stress );
                double I2_ = I2( stress );
                double I3_ = I3( stress );

                return ( 2. / 27 ) * pow( I1_, 3 ) - ( 1. / 3 ) * I1_ * I2_ + I3_;
            }

            double J3strain( const Vector6d& strain ) // determinant of the deviatoric strain tensor
            {
                return voigtToStrain( IDev * strain ).determinant();
            }
        } // namespace Invariants

        namespace Derivatives {
            using namespace Invariants;

            Vector6d dSigmaMdSigma() { return 1. / 3 * I; }

            Vector6d dRhodSigma( double rho, const Vector6d& stress )
            {

                if ( rho <= 1e-16 )
                    return Vector6d::Zero();

                Vector6d s = IDev * stress;

                return 1. / rho * P.array() * s.array();
            }

            Vector6d dThetadSigma( double theta, const Vector6d& stress )
            {
                if ( theta <= 1e-15 || theta >= Pi / 3 - 1e-15 )
                    return Vector6d::Zero();

                // const double J2_ = J2(stress);
                // const double J3_ = J3(stress);

                const double dThetadJ2 = dTheta_dJ2( stress );
                const double dThetadJ3 = dTheta_dJ3( stress );

                if ( isNaN( dThetadJ2 ) || isNaN( dThetadJ3 ) )
                    return Vector6d::Zero();

                return dThetadJ2 * dJ2_dStress( stress ) + dThetadJ3 * dJ3_dStress( stress );
            }

            double dTheta_dJ2( const Vector6d& stress )
            {
                const Vector3d hw    = haighWestergaard( stress );
                const double&  theta = hw( 2 );

                if ( theta <= 1e-14 || theta >= Pi / 3 - 1e-14 )
                    return 1e16;

                const double J2_ = J2( stress );
                const double J3_ = J3( stress );

                const double cos2_3theta = std::cos( 3 * theta ) * std::cos( 3 * theta );
                const double dThetadJ2   = 3 * sqrt3 / 4 * J3_ /
                                         ( std::pow( J2_, 2.5 ) * std::sqrt( 1.0 - cos2_3theta ) );
                return dThetadJ2;
            }

            double dTheta_dJ3( const Vector6d& stress )
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

            double dThetaE_dJ2E( const Vector6d& strain )
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

            double dThetaE_dJ3E( const Vector6d& strain )
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

            Vector6d dJ2_dStress( const Vector6d& stress ) { return P.array() * ( IDev * stress ).array(); }

            Vector6d dJ3_dStress( const Vector6d& stress )
            {
                Vector6d s = IDev * stress;
                return ( P.array() * stressToVoigt( voigtToStress( s ) * voigtToStress( s ) ).array() ).matrix() -
                       2. / 3. * J2( stress ) * I;
            }

            Vector6d dJ2E_dE( const Vector6d& strain ) { return PInv.array() * ( IDev * strain ).array(); }

            Vector6d dJ3E_dE( const Vector6d& strain )
            {
                Vector6d e = IDev * strain;
                return ( PInv.array() * strainToVoigt( voigtToStrain( e ) * voigtToStrain( e ) ).array() ).matrix() -
                       2. / 3. * J2strain( strain ) * I;
            }

            Vector6d dThetaE_dE( const Vector6d& strain )

            {
                return dThetaE_dJ2E( strain ) * dJ2E_dE( strain ) + dThetaE_dJ3E( strain ) * dJ3E_dE( strain );
            }

            Matrix36 dStressPrincipals_dStress( const Vector6d& stress ) // derivative when principal stresses are
                                                                         // computed from solving Eigenvalue-Problem
            {
                MatrixXd J( 3, 6 );

                Vector6d leftX;
                Vector6d rightX;

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

            Vector3d dDeltaEpvneg_dDeltaEpPrincipals( const Vector6d& strain )
            {
                Vector3d       dEvdEpPrinc  = Vector3d::Zero();
                const Vector3d deltaEpPrinc = principalStrainsHW( strain );

                for ( int i = 0; i < dEvdEpPrinc.size(); i++ )
                    dEvdEpPrinc( i ) = -Math::heaviside( -deltaEpPrinc( i ) );

                return dEvdEpPrinc;
            }

            Matrix6d dEp_dE( const Matrix6d& CelInv, const Matrix6d& Cep )
            {
                return Matrix6d::Identity() - CelInv * Cep;
            }

            RowVector6d dDeltaEpv_dE( const Matrix6d& CelInv, const Matrix6d& Cep )
            {
                return I.transpose() * ( Matrix6d::Identity() - CelInv * Cep );
            }

            Matrix36 dDeltaEpPrincipals_dDeltaEp(
                const Vector6d& dEp ) // equations from page 218-219 PhD Thesis David Unteregger
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
                dEpPrinc_dEprho << sqrt2_3 * std::cos( hw( 2 ) ),
                    sqrt2_3 * std::cos( hw( 2 ) - 2. * Constants::Pi / 3. ),
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

            RowVector6d dDeltaEpvneg_dE( const Vector6d& dEp, const Matrix6d& CelInv, const Matrix6d& Cep )
            {
                return dDeltaEpvneg_dDeltaEpPrincipals( dEp ).transpose() * dDeltaEpPrincipals_dDeltaEp( dEp ) *
                       dEp_dE( CelInv, Cep );
            }

        } // namespace Derivatives
        namespace Transformations {
            Matrix6d transformationMatrixStrainVoigt( const Matrix3d& transformedCoordinateSystem )
            {
                Matrix6d transformationMatrix = transformationMatrixStressVoigt( transformedCoordinateSystem );
                transformationMatrix.topRightCorner( 3, 3 ) *= 0.5;
                transformationMatrix.bottomLeftCorner( 3, 3 ) *= 2;

                return transformationMatrix;
            }

            Matrix6d transformationMatrixStressVoigt( const Matrix3d& transformedCoordinateSystem )
            {
                const Matrix3d N = Math::directionCosines( transformedCoordinateSystem );

                Matrix6d transformationMatrix;

                // clang-format off
                transformationMatrix << 
                    pow(N(0,0),2), pow(N(0,1),2), pow(N(0,2),2), 2*N(0,0)*N(0,1), 2*N(0,2)*N(0,1), 2*N(0,2)*N(0,0),
                    pow(N(1,0),2), pow(N(1,1),2), pow(N(1,2),2), 2*N(1,0)*N(1,1), 2*N(1,2)*N(1,1), 2*N(1,0)*N(1,2),
                    pow(N(2,0),2), pow(N(2,1),2), pow(N(2,2),2), 2*N(2,0)*N(2,1), 2*N(2,2)*N(2,1), 2*N(2,0)*N(2,2),
                    N(0,0)*N(1,0), N(0,1)*N(1,1), N(0,2)*N(1,2), N(0,0)*N(1,1)+N(0,1)*N(1,0), N(0,1)*N(1,2)+N(0,2)*N(1,1), N(0,0)*N(1,2)+N(0,2)*N(1,0),
                    N(1,0)*N(2,0), N(1,1)*N(2,1), N(1,2)*N(2,2), N(1,0)*N(2,1)+N(1,1)*N(2,0), N(1,1)*N(2,2)+N(1,2)*N(2,1), N(1,0)*N(2,2)+N(1,2)*N(2,0),
                    N(2,0)*N(0,0), N(2,1)*N(0,1), N(2,2)*N(0,2), N(2,0)*N(0,1)+N(2,1)*N(0,0), N(2,1)*N(0,2)+N(2,2)*N(0,1), N(2,0)*N(0,2)+N(2,2)*N(0,0);
                // clang-format on

                return transformationMatrix;
            }

            Matrix36d projectVoigtStressToPlane( const Vector3d& normalVector )
            {
                const Vector3d& n = normalVector;

                // clang-format off
                Matrix36d projectMatrix;
                projectMatrix << n( 0 ),      0,      0, n( 1 ), 0,      n( 2 ), 
                                      0, n( 1 ),      0, n( 0 ), n( 2 ),      0, 
                                      0,      0, n( 2 ),      0, n( 1 ), n( 0 );
                // clang-format on
                return projectMatrix;
            }

            Matrix36d projectVoigtStrainToPlane( const Vector3d& normalVector )
            {
                Matrix36d projectMatrix = projectVoigtStressToPlane( normalVector );
                projectMatrix.topRightCorner( 3, 3 ) *= 0.5;

                return projectMatrix;
            }

            Marmot::Vector6d rotateVoigtStress( const Eigen::Matrix3d Q, const Marmot::Vector6d& voigtStress )
            {
                const Matrix3d& T  = voigtToStress( voigtStress );
                const Matrix3d& TR = Q * T * Q.transpose();
                return stressToVoigt( TR );
            }
        } // namespace Transformations
    }     // namespace ContinuumMechanics::VoigtNotation
} // namespace Marmot
