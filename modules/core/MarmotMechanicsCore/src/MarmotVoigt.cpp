#include "Marmot/MarmotVoigt.h"
#include "Marmot/HaighWestergaard.h"
#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotTypedefs.h"

using namespace Eigen;

namespace Marmot::ContinuumMechanics::Voigt {
  using namespace Marmot::Math::Constants;
  using namespace ContinuumMechanics::Invariants::HaighWestergaard;

  Marmot::Math::Vector6d strainToVoigt( const Marmot::Math::Matrix3d& strainTensor )
  {
    Marmot::Math::Vector6d strain;
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

  Marmot::Math::Vector3d voigtToPlaneVoigt( const Marmot::Math::Vector6d& voigt )
  {
    /* converts a 6d voigt Vector with Abaqus notation
       S11, S22, S33, S12, S13, S23
       to a 3d Vector S11, S22, S12 */
    Marmot::Math::Vector3d voigtPlane;
    voigtPlane << voigt[0], voigt[1], voigt[3];
    return voigtPlane;
  }

  Marmot::Math::Vector4d voigtToAxisymmetricVoigt( const Marmot::Math::Vector6d& voigt )
  {
    // Converts full 3D Voigt vector to axisymmetric Voigt vector
    /* Assumes Voigt order:
       [0] ε_rr
       [1] ε_zz
       [2] ε_θθ
       [3] γ_rz
       [4] γ_rθ (ignored)
       [5] γ_zθ (ignored)

       Axisymmetric Voigt order:
       [0] ε_rr
       [1] ε_zz
       [2] ε_θθ
       [3] γ_rz
    */
    Marmot::Math::Vector4d voigtAxisym;
    voigtAxisym << voigt[0], voigt[1], voigt[2], voigt[3];
    return voigtAxisym;
  }

  Marmot::Math::Vector6d planeVoigtToVoigt( const Marmot::Math::Vector3d& voigtPlane )
  {
    /* converts a 3d voigt Vector with notation
       S11, S22, S12 to a Marmot::Math::Vector6d with
       S11, S22, S33, S12, S13, S23
       !!! Don't use if 3rd component is NOT ZERO !!!*/
    Marmot::Math::Vector6d voigt;
    voigt << voigtPlane[0], voigtPlane[1], 0, voigtPlane[2], 0, 0;
    return voigt;
  }

  Marmot::Math::Vector6d axisymmetricVoigtToVoigt( const Marmot::Math::Vector4d& voigtAxisymmetric )
  {
    // Input: ε_rr, ε_zz, ε_tt (θθ), γ_rz
    // Output: ε_rr, ε_zz, ε_tt, γ_rz, 0, 0
    Marmot::Math::Vector6d voigt;
    voigt << voigtAxisymmetric[0], // ε_rr
      voigtAxisymmetric[1],        // ε_zz
      voigtAxisymmetric[2],        // ε_θθ
      voigtAxisymmetric[3],        // γ_rz
      0.0,                         // γ_rθ
      0.0;                         // γ_zθ
    return voigt;
  }

} // namespace Marmot::ContinuumMechanics::Voigt

namespace Marmot::ContinuumMechanics::Invariants {
  using namespace Marmot::Math::Constants;
  using namespace ContinuumMechanics::Voigt;
  using namespace ContinuumMechanics::Invariants::HaighWestergaard;

  Marmot::Math::Vector3d principalStrains( const Marmot::Math::Vector6d& voigtStrain )
  {
    SelfAdjointEigenSolver< Marmot::Math::Matrix3d > es( voigtToStrain( voigtStrain ) );
    return es.eigenvalues();
  }

  Marmot::Math::Vector3d principalStresses( const Marmot::Math::Vector6d& voigtStress )
  {
    SelfAdjointEigenSolver< Marmot::Math::Matrix3d > es( voigtToStress( voigtStress ) );
    return es.eigenvalues();
  }

  Marmot::Math::Vector3d sortedPrincipalStrains( const Marmot::Math::Vector6d& voigtStrain )
  {
    HaighWestergaardCoordinates hw = haighWestergaardFromStrain( voigtStrain );
    Marmot::Math::Vector3d      strainPrinc;

    strainPrinc << hw.xi / sqrt3 + sqrt2_3 * hw.rho * std::cos( hw.theta ),
      hw.xi / sqrt3 + sqrt2_3 * hw.rho * ( -std::sin( Marmot::Math::Constants::Pi / 6. - hw.theta ) ),
      hw.xi / sqrt3 + sqrt2_3 * hw.rho * ( -std::sin( Marmot::Math::Constants::Pi / 6. + hw.theta ) );

    return strainPrinc;
  }

  Eigen::Matrix3d principalStressesDirections( const Marmot::Math::Vector6d& voigtStress )
  {
    SelfAdjointEigenSolver< Marmot::Math::Matrix3d > es( voigtToStress( voigtStress ) );
    Marmot::Math::Matrix3d                           Q = es.eigenvectors();
    Q.col( 2 ) = Q.col( 0 ).cross( Q.col( 1 ) ); // for a clockwise coordinate system
    return Q;
  }

  std::pair< Eigen::Vector3d, Eigen::Matrix< double, 3, 6 > > principalValuesAndDerivatives(
    const Eigen::Matrix< double, 6, 1 >& S )
  {
    // This is a fast implementation of the classical algorithm for determining
    // the principal components of a symmetric 3x3 Matrix in Voigt notation (off
    // diagonals expected with factor 1) as well as its respective derivatives

    using namespace Eigen;

    using Matrix66 = Matrix< double, 6, 6 >;

    const static auto dS0_dS = ( Marmot::Math::Vector6d() << 1, 0, 0, 0, 0, 0 ).finished();
    const static auto dS1_dS = ( Marmot::Math::Vector6d() << 0, 1, 0, 0, 0, 0 ).finished();
    const static auto dS2_dS = ( Marmot::Math::Vector6d() << 0, 0, 1, 0, 0, 0 ).finished();

    const double p1 = S( 3 ) * S( 3 ) + S( 4 ) * S( 4 ) + S( 5 ) * S( 5 );

    if ( p1 <= 1e-16 ) { // matrix is already diagonal
                         // clang-format off
                    const static auto dE_dS = ( Marmot::Math::Matrix36() << 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 ).finished();
                    return { S.head( 3 ), dE_dS };
                         // clang-format on
    }

    const Marmot::Math::Vector6d dP1_dS = { 0, 0, 0, 2 * S( 3 ), 2 * S( 4 ), 2 * S( 5 ) };

    const double                 q     = S.head( 3 ).sum() / 3;
    const Marmot::Math::Vector6d dQ_dS = Marmot::ContinuumMechanics::Voigt::I / 3;

    const double p2 = std::pow( S( 0 ) - q, 2 ) + std::pow( S( 1 ) - q, 2 ) + std::pow( S( 2 ) - q, 2 ) + 2 * p1;
    const Marmot::Math::Vector6d dP2_dS = 2 * ( S( 0 ) - q ) * ( dS0_dS - dQ_dS ) +
                                          2 * ( S( 1 ) - q ) * ( dS1_dS - dQ_dS ) +
                                          2 * ( S( 2 ) - q ) * ( dS2_dS - dQ_dS ) + 2 * dP1_dS;

    const double                 p     = std::sqrt( p2 / 6 );
    const Marmot::Math::Vector6d dP_dS = 1. / 12 * 1. / p * dP2_dS;

    const Marmot::Math::Vector6d B     = 1. / p * ( S - q * Marmot::ContinuumMechanics::Voigt::I );
    const Matrix66               dB_dS = -B * 1. / p * dP_dS.transpose() +
                           1. / p * ( Matrix66::Identity() - Marmot::ContinuumMechanics::Voigt::I * dQ_dS.transpose() );

    const double detB = B( 0 ) * B( 1 ) * B( 2 ) + B( 3 ) * B( 4 ) * B( 5 ) * 2 - B( 2 ) * B( 3 ) * B( 3 ) -
                        B( 1 ) * B( 4 ) * B( 4 ) - B( 0 ) * B( 5 ) * B( 5 );

    Marmot::Math::Vector6d dDetB_dB;
    dDetB_dB( 0 ) = B( 1 ) * B( 2 ) - B( 5 ) * B( 5 );
    dDetB_dB( 1 ) = B( 0 ) * B( 2 ) - B( 4 ) * B( 4 );
    dDetB_dB( 2 ) = B( 0 ) * B( 1 ) - B( 3 ) * B( 3 );
    dDetB_dB( 3 ) = B( 4 ) * B( 5 ) * 2 - B( 2 ) * B( 3 ) * 2;
    dDetB_dB( 4 ) = B( 3 ) * B( 5 ) * 2 - B( 1 ) * B( 4 ) * 2;
    dDetB_dB( 5 ) = B( 3 ) * B( 4 ) * 2 - B( 0 ) * B( 5 ) * 2;

    const double                 r     = detB * 1. / 2;
    const Marmot::Math::Vector6d dR_dS = 1. / 2 * dDetB_dB.transpose() * dB_dS;

    double phi;
    double dPhi_dR;
    if ( r <= -1 ) {
      phi     = Marmot::Math::Constants::Pi / 3;
      dPhi_dR = 0.0;
    }
    else if ( r >= 1 ) {
      phi     = 0.0;
      dPhi_dR = 0.0;
    }
    else {
      phi     = std::acos( r ) / 3;
      dPhi_dR = 1. / 3 * -1. / ( std::sqrt( 1 - r * r ) );
    }
    const Marmot::Math::Vector6d dPhi_dS = dPhi_dR * dR_dS;

    Marmot::Math::Vector3d e;
    Marmot::Math::Matrix36 dE_dS;

    e( 0 )         = q + 2 * p * std::cos( phi );
    dE_dS.row( 0 ) = dQ_dS + 2 * ( dP_dS * std::cos( phi ) - p * std::sin( phi ) * dPhi_dS );

    const double phiShifted = phi + ( 2. / 3 * Marmot::Math::Constants::Pi );
    e( 2 )                  = q + 2 * p * std::cos( phiShifted );
    dE_dS.row( 2 )          = dQ_dS + 2 * ( dP_dS * std::cos( phiShifted ) - p * std::sin( phiShifted ) * dPhi_dS );

    e( 1 )         = 3 * q - e( 0 ) - e( 2 );
    dE_dS.row( 1 ) = 3 * dQ_dS - dE_dS.row( 0 ).transpose() - dE_dS.row( 2 ).transpose();

    return { e, dE_dS };
  }

  double vonMisesEquivalentStress( const Marmot::Math::Vector6d& stress )
  {
    return sqrt( 3. * J2( stress ) );
  }

  double vonMisesEquivalentStrain( const Marmot::Math::Vector6d& strain )
  {
    // e_eq = sqrt( 2/3 * e_ij * e_ij )
    const Marmot::Math::Vector6d& e = strain;
    return std::sqrt( 2. / 3. * ( e( 0 ) * e( 0 ) + e( 1 ) * e( 1 ) + e( 2 ) * e( 2 ) ) +
                      1. / 3. * ( e( 3 ) * e( 3 ) + e( 4 ) * e( 4 ) + e( 5 ) * e( 5 ) ) );
  }

  double normStress( const Marmot::Math::Vector6d& stress )
  {
    return ContinuumMechanics::Voigt::voigtToStress( stress ).norm();
  }

  double StrainVolumetricNegative( const Marmot::Math::Vector6d& strain )
  {
    Marmot::Math::Vector3d dEpPrincipal = principalStrains( strain );

    return Math::macauly( -dEpPrincipal( 0 ) ) + Math::macauly( -dEpPrincipal( 1 ) ) +
           Math::macauly( -dEpPrincipal( 2 ) );
  }

  double I1Strain( const Marmot::Math::Vector6d& strain )
  {
    return strain.head( 3 ).sum();
  }

  double I2Strain( const Marmot::Math::Vector6d& strain ) // you could also use normal I2, but
                                                          // with epsilon12 instead of 2*epsilon12
  {
    const Marmot::Math::Vector6d& e = strain;

    return e( 0 ) * e( 1 ) + e( 1 ) * e( 2 ) + e( 2 ) * e( 0 ) - e( 3 ) / 2. * e( 3 ) / 2. - e( 4 ) / 2. * e( 4 ) / 2. -
           e( 5 ) / 2. * e( 5 ) / 2.;
  }

  double I3Strain( const Marmot::Math::Vector6d& strain ) // you could also use normal I3, but
                                                          // with epsilon12 instead of 2*epsilon12
  {
    return voigtToStrain( strain ).determinant();
  }

  double J2Strain( const Marmot::Math::Vector6d& strain )
  {
    const double res = 1. / 3. * std::pow( I1( strain ), 2. ) - I2Strain( strain );
    return res > 0 ? res : 0;
  }

  double J3Strain( const Marmot::Math::Vector6d& strain ) // determinant of the deviatoric strain tensor
  {
    return voigtToStrain< double >( IDev * strain ).determinant();
  }

} // namespace Marmot::ContinuumMechanics::Invariants

namespace Marmot::ContinuumMechanics::Invariants::Derivatives {
  using namespace Marmot::Math::Constants;
  using namespace ContinuumMechanics::Voigt;
  using namespace ContinuumMechanics::Invariants::HaighWestergaard;
  using namespace Invariants;

  Marmot::Math::Vector6d dStressMean_dStress()
  {
    return 1. / 3 * I;
  }

  double dThetaStrain_dJ2Strain( const Marmot::Math::Vector6d& strain )
  {
    const HaighWestergaardCoordinates hw    = haighWestergaardFromStrain( strain );
    const double                      theta = hw.theta;

    if ( theta <= 1e-15 || theta >= Pi / 3 - 1e-15 )
      return 1e16;
    else
      return 3. * std::sqrt( 3. ) / 4. * J3Strain( strain ) /
             ( std::pow( J2Strain( strain ), 5. / 2 ) * std::sqrt( 1. - std::pow( std::cos( 3. * theta ), 2. ) ) );
  }

  double dThetaStrain_dJ3Strain( const Marmot::Math::Vector6d& strain )
  {
    const HaighWestergaardCoordinates hw    = haighWestergaardFromStrain( strain );
    const double&                     theta = hw.theta;

    if ( theta <= 1e-15 || theta >= Pi / 3 - 1e-15 )
      return -1e16;
    else
      return -std::sqrt( 3. ) / 2. * 1. /
             ( std::pow( J2Strain( strain ), 3. / 2 ) * std::sqrt( 1. - std::pow( std::cos( 3. * theta ), 2. ) ) );
  }

  Marmot::Math::Vector6d dJ2Strain_dStrain( const Marmot::Math::Vector6d& strain )
  {
    return PInv.array() * ( IDev * strain ).array();
  }

  Marmot::Math::Vector6d dJ3Strain_dStrain( const Marmot::Math::Vector6d& strain )
  {
    Marmot::Math::Vector6d e = IDev * strain;
    return ( PInv.array() * strainToVoigt( voigtToStrain( e ) * voigtToStrain( e ) ).array() ).matrix() -
           2. / 3. * J2Strain( strain ) * I;
  }

  Marmot::Math::Vector6d dThetaStrain_dStrain( const Marmot::Math::Vector6d& strain )

  {
    return dThetaStrain_dJ2Strain( strain ) * dJ2Strain_dStrain( strain ) +
           dThetaStrain_dJ3Strain( strain ) * dJ3Strain_dStrain( strain );
  }

  Marmot::Math::Matrix36 dStressPrincipals_dStress(
    const Marmot::Math::Vector6d& stress ) // derivative when principal stresses are
                                           // computed from solving Eigenvalue-Problem
  {
    MatrixXd J( 3, 6 );

    Marmot::Math::Vector6d leftX;
    Marmot::Math::Vector6d rightX;

    for ( size_t i = 0; i < 6; i++ ) {
      double volatile h = std::max( 1.0, std::abs( stress( i ) ) ) * Marmot::Math::Constants::cubicRootEps();
      leftX             = stress;
      leftX( i ) -= h;
      rightX = stress;
      rightX( i ) += h;

      J.col( i ) = 1. / ( 2 * h ) * ( principalStresses( rightX ) - principalStresses( leftX ) );
    }
    return J;
  }

  Marmot::Math::Vector3d dStrainVolumetricNegative_dStrainPrincipal( const Marmot::Math::Vector6d& strain )
  {
    Marmot::Math::Vector3d       dEvdEpPrinc  = Marmot::Math::Vector3d::Zero();
    const Marmot::Math::Vector3d deltaEpPrinc = Invariants::sortedPrincipalStrains( strain );

    for ( int i = 0; i < dEvdEpPrinc.size(); i++ )
      dEvdEpPrinc( i ) = -Math::heaviside( -deltaEpPrinc( i ) );

    return dEvdEpPrinc;
  }

  Marmot::Math::Matrix6d dEp_dE( const Marmot::Math::Matrix6d& CelInv, const Marmot::Math::Matrix6d& Cep )
  {
    return Marmot::Math::Matrix6d::Identity() - CelInv * Cep;
  }

  Marmot::Math::RowVector6d dDeltaEpv_dE( const Marmot::Math::Matrix6d& CelInv, const Marmot::Math::Matrix6d& Cep )
  {
    return I.transpose() * ( Marmot::Math::Matrix6d::Identity() - CelInv * Cep );
  }

  Marmot::Math::Matrix36 dSortedStrainPrincipal_dStrain(
    const Marmot::Math::Vector6d& dEp ) // equations from page 218-219 PhD Thesis David Unteregger
  {
    Marmot::Math::Vector3d dEpPrinc_dEpvol   = Marmot::Math::Vector3d::Zero();
    Marmot::Math::Vector3d dEpPrinc_dEprho   = Marmot::Math::Vector3d::Zero();
    Marmot::Math::Vector3d dEPprinc_dEptheta = Marmot::Math::Vector3d::Zero();

    const double                      sqrt23 = std::sqrt( 2. / 3. );
    const HaighWestergaardCoordinates hw     = haighWestergaardFromStrain( dEp );
    // const double& epsM =		hw(0);
    const double& rhoE = hw.rho;
    // const double& thetaE =		hw(2);

    dEpPrinc_dEpvol = 1. / 3. * Marmot::Math::Vector3d::Ones();
    dEpPrinc_dEprho << sqrt23 * std::cos( hw.theta ),
      sqrt23 * std::cos( hw.theta - 2. * Marmot::Math::Constants::Pi / 3. ),
      sqrt23 * std::cos( hw.theta + 2. * Marmot::Math::Constants::Pi / 3. );

    dEPprinc_dEptheta << -sqrt23 * hw.rho * std::sin( hw.theta ),
      -sqrt23 * hw.rho * std::sin( hw.theta - 2. * Marmot::Math::Constants::Pi / 3. ),
      -sqrt23 * hw.rho * std::sin( hw.theta + 2. * Marmot::Math::Constants::Pi / 3. );

    Marmot::Math::RowVector6d dEpvol_dEp   = Marmot::Math::RowVector6d::Zero();
    dEpvol_dEp                             = I;
    Marmot::Math::RowVector6d dEprho_dEp   = Marmot::Math::RowVector6d::Zero();
    Marmot::Math::RowVector6d dEptheta_dEp = Marmot::Math::RowVector6d::Zero();

    if ( std::abs( rhoE ) > 1e-16 ) {
      dEprho_dEp   = 1. / rhoE * dJ2Strain_dStrain( dEp ).transpose();
      dEptheta_dEp = ( dThetaStrain_dJ2Strain( dEp ) * dJ2Strain_dStrain( dEp ).transpose() ) +
                     ( dThetaStrain_dJ3Strain( dEp ) * dJ3Strain_dStrain( dEp ).transpose() );
    }
    else {
      dEprho_dEp << 1.e16, 1.e16, 1.e16, 1.e16, 1.e16,
        1.e16; // 1e16 from Code David (Line 67,
               // D_2_Umatsub_damage3_derivatives)
      dEptheta_dEp << 0., 0., 0., 0., 0., 0.;
    }

    Marmot::Math::Matrix36 dEpPrincdEpAna;
    dEpPrincdEpAna = ( dEpPrinc_dEpvol * dEpvol_dEp ) + ( dEpPrinc_dEprho * dEprho_dEp ) +
                     ( dEPprinc_dEptheta * dEptheta_dEp );

    return dEpPrincdEpAna;
  }

  Marmot::Math::RowVector6d dDeltaEpvneg_dE( const Marmot::Math::Vector6d& dEp,
                                             const Marmot::Math::Matrix6d& CelInv,
                                             const Marmot::Math::Matrix6d& Cep )
  {
    return dStrainVolumetricNegative_dStrainPrincipal( dEp ).transpose() * dSortedStrainPrincipal_dStrain( dEp ) *
           dEp_dE( CelInv, Cep );
  }

} // namespace Marmot::ContinuumMechanics::Invariants::Derivatives
namespace Marmot::ContinuumMechanics::Voigt::Transformations {
  using namespace Marmot::Math::Constants;
  Marmot::Math::Matrix6d transformationMatrixStrainVoigt( const Marmot::Math::Matrix3d& transformedCoordinateSystem )
  {
    Marmot::Math::Matrix6d transformationMatrix = transformationMatrixStressVoigt( transformedCoordinateSystem );
    transformationMatrix.topRightCorner( 3, 3 ) *= 0.5;
    transformationMatrix.bottomLeftCorner( 3, 3 ) *= 2;

    return transformationMatrix;
  }

  Marmot::Math::Matrix6d transformationMatrixStressVoigt( const Marmot::Math::Matrix3d& transformedCoordinateSystem )
  {
    const Marmot::Math::Matrix3d N = Math::directionCosines( transformedCoordinateSystem );

    Marmot::Math::Matrix6d transformationMatrix;

    // clang-format off
        transformationMatrix <<
                    pow(N(0,0),2), pow(N(0,1),2), pow(N(0,2),2), 2*N(0,0)*N(0,1), 2*N(0,0)*N(0,2), 2*N(0,2)*N(0,1),
                    pow(N(1,0),2), pow(N(1,1),2), pow(N(1,2),2), 2*N(1,0)*N(1,1), 2*N(1,0)*N(1,2), 2*N(1,2)*N(1,1),
                    pow(N(2,0),2), pow(N(2,1),2), pow(N(2,2),2), 2*N(2,0)*N(2,1), 2*N(2,0)*N(2,2), 2*N(2,2)*N(2,1),
                    N(0,0)*N(1,0), N(0,1)*N(1,1), N(0,2)*N(1,2), N(0,0)*N(1,1)+N(0,1)*N(1,0), N(0,0)*N(1,2)+N(0,2)*N(1,0), N(0,1)*N(1,2)+N(0,2)*N(1,1),
                    N(2,0)*N(0,0), N(2,1)*N(0,1), N(2,2)*N(0,2), N(2,0)*N(0,1)+N(2,1)*N(0,0), N(2,0)*N(0,2)+N(2,2)*N(0,0), N(2,1)*N(0,2)+N(2,2)*N(0,1),
                    N(1,0)*N(2,0), N(1,1)*N(2,1), N(1,2)*N(2,2), N(1,0)*N(2,1)+N(1,1)*N(2,0), N(1,0)*N(2,2)+N(1,2)*N(2,0), N(1,1)*N(2,2)+N(1,2)*N(2,1);
    // clang-format on

    return transformationMatrix;
  }

  Marmot::Math::Matrix36d projectVoigtStressToPlane( const Marmot::Math::Vector3d& normalVector )
  {
    const Marmot::Math::Vector3d& n = normalVector;

    // clang-format off
                Marmot::Math::Matrix36d projectMatrix;
                projectMatrix << n( 0 ),      0,      0, n( 1 ), n( 2 ),      0,
                                      0, n( 1 ),      0, n( 0 ),      0, n( 2 ),
                                      0,      0, n( 2 ),      0, n( 0 ), n( 1 );
    // clang-format on
    return projectMatrix;
  }

  Marmot::Math::Matrix36d projectVoigtStrainToPlane( const Marmot::Math::Vector3d& normalVector )
  {
    Marmot::Math::Matrix36d projectMatrix = projectVoigtStressToPlane( normalVector );
    projectMatrix.topRightCorner( 3, 3 ) *= 0.5;

    return projectMatrix;
  }

  Marmot::Math::Vector6d rotateVoigtStress( const Eigen::Matrix3d& Q, const Marmot::Math::Vector6d& voigtStress )
  {
    const Marmot::Math::Matrix3d& T  = voigtToStress( voigtStress );
    const Marmot::Math::Matrix3d& TR = Q * T * Q.transpose();
    return stressToVoigt( TR );
  }

  Marmot::Math::Vector6d transformStressToLocalSystem( const Marmot::Math::Vector6d& stress,
                                                       const Marmot::Math::Matrix3d& transformedCoordinateSystem )
  {
    const Marmot::Math::Matrix3d s_prime = Math::transformToLocalSystem( voigtToStress( stress ),
                                                                         transformedCoordinateSystem );

    return stressToVoigt( s_prime );
  }

  Marmot::Math::Vector6d transformStrainToLocalSystem( const Marmot::Math::Vector6d& strain,
                                                       const Marmot::Math::Matrix3d& transformedCoordinateSystem )
  {
    const Marmot::Math::Matrix3d e_prime = Math::transformToLocalSystem( voigtToStrain( strain ),
                                                                         transformedCoordinateSystem );

    return strainToVoigt( e_prime );
  }

  Marmot::Math::Vector6d transformStressToGlobalSystem( const Marmot::Math::Vector6d& stress,
                                                        const Marmot::Math::Matrix3d& transformedCoordinateSystem )
  {
    const Marmot::Math::Matrix3d s_prime = Math::transformToGlobalSystem( voigtToStress( stress ),
                                                                          transformedCoordinateSystem );

    return stressToVoigt( s_prime );
  }
  Marmot::Math::Vector6d transformStrainToGlobalSystem( const Marmot::Math::Vector6d& strain,
                                                        const Marmot::Math::Matrix3d& transformedCoordinateSystem )
  {
    const Marmot::Math::Matrix3d e_prime = Math::transformToGlobalSystem( voigtToStrain( strain ),
                                                                          transformedCoordinateSystem );

    return strainToVoigt( e_prime );
  }

  Marmot::Math::Matrix6d transformStiffnessToGlobalSystem( const Marmot::Math::Matrix6d& stiffness,
                                                           const Marmot::Math::Matrix3d& transformedCoordinateSystem )
  {
    const Marmot::Math::Tensors::Tensor3333d
      stiffnessTensorLocal = voigtToStiffness< Marmot::Math::Tensors::Tensor3333d >( stiffness );
    Marmot::Math::Tensors::Tensor3333d stiffnessTensorGlobal;
    stiffnessTensorGlobal.setZero();
    Marmot::Math::Matrix3d N = transformedCoordinateSystem.transpose();

    for ( int i = 0; i < 3; i++ )
      for ( int j = 0; j < 3; j++ )
        for ( int k = 0; k < 3; k++ )
          for ( int l = 0; l < 3; l++ )
            for ( int m = 0; m < 3; m++ )
              for ( int n = 0; n < 3; n++ )
                for ( int o = 0; o < 3; o++ )
                  for ( int p = 0; p < 3; p++ )
                    stiffnessTensorGlobal( i, j, k, l ) += N( i, m ) * N( j, n ) * N( k, o ) * N( l, p ) *
                                                           stiffnessTensorLocal( m, n, o, p );

    Marmot::Math::Matrix6d res = stiffnessToVoigt( stiffnessTensorGlobal );
    return res;
  }

} // namespace Marmot::ContinuumMechanics::Voigt::Transformations
