#include "Marmot/MarmotVoigt.h"
#include "Marmot/HaighWestergaard.h"
#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotTypedefs.h"

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

    Eigen::Matrix< double, 6, 6 > stiffnessToVoigt( const Eigen::Tensor< double, 4 >& C )
    {
      // Ordering for Voigt notation (0->xx, 1->yy, 2->zz, 3->xy, 4->yz, 5->xz)
      std::array< std::pair< int, int >, 6 > ordering = {
        { { 0, 0 }, { 1, 1 }, { 2, 2 }, { 0, 1 }, { 2, 0 }, { 1, 2 } } };

      Eigen::Matrix< double, 6, 6 > voigtStiffness;
      voigtStiffness.setZero(); // Initialize with zeros

      for ( int a = 0; a < 6; ++a ) {
        int i = ordering[a].first;
        int j = ordering[a].second;
        for ( int b = 0; b < 6; ++b ) {
          int k = ordering[b].first;
          int l = ordering[b].second;

          // Populate the Voigt stiffness matrix
          voigtStiffness( a, b ) = C( i, j, k, l );
          voigtStiffness( a, b ) += C( j, i, k, l );
          voigtStiffness( a, b ) += C( j, i, l, k );
          voigtStiffness( a, b ) += C( i, j, l, k );
          voigtStiffness( a, b ) /= 4.0;
        }
      }

      return voigtStiffness;
    }

    Eigen::Tensor< double, 4 > voigtToStiffness( const Eigen::Matrix< double, 6, 6 >& voigtStiffness )
    {
      using namespace TensorUtility::IndexNotation;

      EigenTensors::Tensor3333d stiffness;
      stiffness.setZero(); // Initialize with zeros

      int row;
      int col;
      for ( int i = 0; i < 3; i++ ) {
        for ( int j = 0; j < 3; j++ ) {
          row = toVoigt< 3 >( i, j );
          for ( int k = 0; k < 3; k++ ) {
            for ( int l = 0; l < 3; l++ ) {
              col = toVoigt< 3 >( k, l );
              stiffness( i, j, k, l ) += voigtStiffness( row, col );
            };
          };
        };
      };

      return stiffness;
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

    Vector4d voigtToAxisymmetricVoigt( const Vector6d& voigt )
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
      Vector4d voigtAxisym;
      voigtAxisym << voigt[0], voigt[1], voigt[2], voigt[3];
      return voigtAxisym;
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

    Vector6d axisymmetricVoigtToVoigt( const Vector4d& voigtAxisymmetric )
    {
      // Input: ε_rr, ε_zz, ε_tt (θθ), γ_rz
      // Output: ε_rr, ε_zz, ε_tt, γ_rz, 0, 0
      Vector6d voigt;
      voigt << voigtAxisymmetric[0], // ε_rr
        voigtAxisymmetric[1],        // ε_zz
        voigtAxisymmetric[2],        // ε_θθ
        voigtAxisymmetric[3],        // γ_rz
        0.0,                         // γ_rθ
        0.0;                         // γ_zθ
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

      Vector3d sortedPrincipalStrains( const Vector6d& voigtStrain )
      {
        HaighWestergaardCoordinates hw = haighWestergaardFromStrain( voigtStrain );
        Vector3d                    strainPrinc;

        strainPrinc << hw.xi / sqrt3 + sqrt2_3 * hw.rho * std::cos( hw.theta ),
          hw.xi / sqrt3 + sqrt2_3 * hw.rho * ( -std::sin( Constants::Pi / 6. - hw.theta ) ),
          hw.xi / sqrt3 + sqrt2_3 * hw.rho * ( -std::sin( Constants::Pi / 6. + hw.theta ) );

        return strainPrinc;
      }

      Eigen::Matrix3d principalStressesDirections( const Marmot::Vector6d& voigtStress )
      {
        SelfAdjointEigenSolver< Matrix3d > es( voigtToStress( voigtStress ) );
        Matrix3d                           Q = es.eigenvectors();
        Q.col( 2 )                           = Q.col( 0 ).cross( Q.col( 1 ) ); // for a clockwise coordinate system
        return Q;
      }

      std::pair< Eigen::Vector3d, Eigen::Matrix< double, 3, 6 > > principalValuesAndDerivatives(
        const Eigen::Matrix< double, 6, 1 >& S )
      {
        // This is a fast implementation of the classical algorithm for determining
        // the principal components of a symmetric 3x3 Matrix in Voigt notation (off
        // diagonals expected with factor 1) as well as its respective derivatives

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

        const double   p2 = std::pow( S( 0 ) - q, 2 ) + std::pow( S( 1 ) - q, 2 ) + std::pow( S( 2 ) - q, 2 ) + 2 * p1;
        const Vector6d dP2_dS = 2 * ( S( 0 ) - q ) * ( dS0_dS - dQ_dS ) + 2 * ( S( 1 ) - q ) * ( dS1_dS - dQ_dS ) +
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
        dE_dS.row( 2 )          = dQ_dS + 2 * ( dP_dS * std::cos( phiShifted ) - p * std::sin( phiShifted ) * dPhi_dS );

        e( 1 )         = 3 * q - e( 0 ) - e( 2 );
        dE_dS.row( 1 ) = 3 * dQ_dS - dE_dS.row( 0 ).transpose() - dE_dS.row( 2 ).transpose();

        return { e, dE_dS };
      }

      double vonMisesEquivalentStress( const Vector6d& stress )
      {
        return sqrt( 3. * J2( stress ) );
      }

      double vonMisesEquivalentStrain( const Vector6d& strain )
      {
        // e_eq = sqrt( 2/3 * e_ij * e_ij )
        const Vector6d& e = strain;
        return std::sqrt( 2. / 3. * ( e( 0 ) * e( 0 ) + e( 1 ) * e( 1 ) + e( 2 ) * e( 2 ) ) +
                          1. / 3. * ( e( 3 ) * e( 3 ) + e( 4 ) * e( 4 ) + e( 5 ) * e( 5 ) ) );
      }

      double normStress( const Vector6d& stress )
      {
        return ContinuumMechanics::VoigtNotation::voigtToStress( stress ).norm();
      }

      double StrainVolumetricNegative( const Vector6d& strain )
      {
        Vector3d dEpPrincipal = principalStrains( strain );

        return Math::macauly( -dEpPrincipal( 0 ) ) + Math::macauly( -dEpPrincipal( 1 ) ) +
               Math::macauly( -dEpPrincipal( 2 ) );
      }

      double I1Strain( const Vector6d& strain )
      {
        return strain.head( 3 ).sum();
      }

      double I2Strain( const Vector6d& strain ) // you could also use normal I2, but
                                                // with epsilon12 instead of 2*epsilon12
      {
        const Vector6d& e = strain;

        return e( 0 ) * e( 1 ) + e( 1 ) * e( 2 ) + e( 2 ) * e( 0 ) - e( 3 ) / 2. * e( 3 ) / 2. -
               e( 4 ) / 2. * e( 4 ) / 2. - e( 5 ) / 2. * e( 5 ) / 2.;
      }

      double I3Strain( const Vector6d& strain ) // you could also use normal I3, but
                                                // with epsilon12 instead of 2*epsilon12
      {
        return voigtToStrain( strain ).determinant();
      }

      double J2Strain( const Vector6d& strain )
      {
        const double res = 1. / 3. * std::pow( I1( strain ), 2. ) - I2Strain( strain );
        return res > 0 ? res : 0;
      }

      double J3Strain( const Vector6d& strain ) // determinant of the deviatoric strain tensor
      {
        return voigtToStrain< double >( IDev * strain ).determinant();
      }

    } // namespace Invariants

    namespace Derivatives {
      using namespace Invariants;

      Vector6d dStressMean_dStress()
      {
        return 1. / 3 * I;
      }

      Vector6d dTheta_dStress( double theta, const Vector6d& stress )
      {
        if ( theta <= 1e-15 || theta >= Pi / 3 - 1e-15 )
          return Vector6d::Zero();

        // const double J2_ = J2(stress);
        // const double J3_ = J3(stress);

        const double dThetadJ2 = dTheta_dJ2( stress );
        const double dThetadJ3 = dTheta_dJ3( stress );

        if ( Math::isNaN( dThetadJ2 ) || Math::isNaN( dThetadJ3 ) )
          return Vector6d::Zero();

        return dThetadJ2 * dJ2_dStress( stress ) + dThetadJ3 * dJ3_dStress( stress );
      }

      double dTheta_dJ2( const Vector6d& stress )
      {
        const HaighWestergaardCoordinates hw    = haighWestergaard( stress );
        const double&                     theta = hw.theta;

        if ( theta <= 1e-14 || theta >= Pi / 3 - 1e-14 )
          return 1e16;

        const double J2_ = J2( stress );
        const double J3_ = J3( stress );

        const double cos2_3theta = std::cos( 3 * theta ) * std::cos( 3 * theta );
        const double dThetadJ2   = 3 * sqrt3 / 4 * J3_ / ( std::pow( J2_, 2.5 ) * std::sqrt( 1.0 - cos2_3theta ) );
        return dThetadJ2;
      }

      double dTheta_dJ3( const Vector6d& stress )
      {
        const HaighWestergaardCoordinates hw    = haighWestergaard( stress );
        const double&                     theta = hw.theta;

        if ( theta <= 1e-14 || theta >= Pi / 3 - 1e-14 )
          return -1e16;

        const double J2_ = J2( stress );
        // const double J3_ = J3(stress);

        const double cos2_3theta = std::cos( 3 * theta ) * std::cos( 3 * theta );
        const double dThetadJ3   = -sqrt3 / 2. * 1. / ( std::pow( J2_, 1.5 ) * std::sqrt( 1.0 - cos2_3theta ) );
        return dThetadJ3;
      }

      double dThetaStrain_dJ2Strain( const Vector6d& strain )
      {
        const HaighWestergaardCoordinates hw    = haighWestergaardFromStrain( strain );
        const double                      theta = hw.theta;

        if ( theta <= 1e-15 || theta >= Pi / 3 - 1e-15 )
          return 1e16;
        else
          return 3. * std::sqrt( 3. ) / 4. * J3Strain( strain ) /
                 ( std::pow( J2Strain( strain ), 5. / 2 ) * std::sqrt( 1. - std::pow( std::cos( 3. * theta ), 2. ) ) );
      }

      double dThetaStrain_dJ3Strain( const Vector6d& strain )
      {
        const HaighWestergaardCoordinates hw    = haighWestergaardFromStrain( strain );
        const double&                     theta = hw.theta;

        if ( theta <= 1e-15 || theta >= Pi / 3 - 1e-15 )
          return -1e16;
        else
          return -std::sqrt( 3. ) / 2. * 1. /
                 ( std::pow( J2Strain( strain ), 3. / 2 ) * std::sqrt( 1. - std::pow( std::cos( 3. * theta ), 2. ) ) );
      }

      Vector6d dJ2_dStress( const Vector6d& stress )
      {
        return P.array() * ( IDev * stress ).array();
      }

      Vector6d dJ3_dStress( const Vector6d& stress )
      {
        Vector6d s = IDev * stress;
        return ( P.array() * stressToVoigt< double >( voigtToStress( s ) * voigtToStress( s ) ).array() ).matrix() -
               2. / 3. * J2( stress ) * I;
      }

      Vector6d dJ2Strain_dStrain( const Vector6d& strain )
      {
        return PInv.array() * ( IDev * strain ).array();
      }

      Vector6d dJ3Strain_dStrain( const Vector6d& strain )
      {
        Vector6d e = IDev * strain;
        return ( PInv.array() * strainToVoigt( voigtToStrain( e ) * voigtToStrain( e ) ).array() ).matrix() -
               2. / 3. * J2Strain( strain ) * I;
      }

      Vector6d dThetaStrain_dStrain( const Vector6d& strain )

      {
        return dThetaStrain_dJ2Strain( strain ) * dJ2Strain_dStrain( strain ) +
               dThetaStrain_dJ3Strain( strain ) * dJ3Strain_dStrain( strain );
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

      Vector3d dStrainVolumetricNegative_dStrainPrincipal( const Vector6d& strain )
      {
        Vector3d       dEvdEpPrinc  = Vector3d::Zero();
        const Vector3d deltaEpPrinc = Invariants::sortedPrincipalStrains( strain );

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

      Matrix36 dSortedStrainPrincipal_dStrain(
        const Vector6d& dEp ) // equations from page 218-219 PhD Thesis David Unteregger
      {
        Vector3d dEpPrinc_dEpvol   = Vector3d::Zero();
        Vector3d dEpPrinc_dEprho   = Vector3d::Zero();
        Vector3d dEPprinc_dEptheta = Vector3d::Zero();

        const double                      sqrt2_3 = std::sqrt( 2. / 3. );
        const HaighWestergaardCoordinates hw      = haighWestergaardFromStrain( dEp );
        // const double& epsM =		hw(0);
        const double& rhoE = hw.rho;
        // const double& thetaE =		hw(2);

        dEpPrinc_dEpvol = 1. / 3. * Vector3d::Ones();
        dEpPrinc_dEprho << sqrt2_3 * std::cos( hw.theta ), sqrt2_3 * std::cos( hw.theta - 2. * Constants::Pi / 3. ),
          sqrt2_3 * std::cos( hw.theta + 2. * Constants::Pi / 3. );

        dEPprinc_dEptheta << -sqrt2_3 * hw.rho * std::sin( hw.theta ),
          -sqrt2_3 * hw.rho * std::sin( hw.theta - 2. * Constants::Pi / 3. ),
          -sqrt2_3 * hw.rho * std::sin( hw.theta + 2. * Constants::Pi / 3. );

        RowVector6d dEpvol_dEp   = RowVector6d::Zero();
        dEpvol_dEp               = I;
        RowVector6d dEprho_dEp   = RowVector6d::Zero();
        RowVector6d dEptheta_dEp = RowVector6d::Zero();

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

        Matrix36 dEpPrincdEpAna;
        dEpPrincdEpAna = ( dEpPrinc_dEpvol * dEpvol_dEp ) + ( dEpPrinc_dEprho * dEprho_dEp ) +
                         ( dEPprinc_dEptheta * dEptheta_dEp );

        return dEpPrincdEpAna;
      }

      RowVector6d dDeltaEpvneg_dE( const Vector6d& dEp, const Matrix6d& CelInv, const Matrix6d& Cep )
      {
        return dStrainVolumetricNegative_dStrainPrincipal( dEp ).transpose() * dSortedStrainPrincipal_dStrain( dEp ) *
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
                    pow(N(0,0),2), pow(N(0,1),2), pow(N(0,2),2), 2*N(0,0)*N(0,1), 2*N(0,0)*N(0,2), 2*N(0,2)*N(0,1),
                    pow(N(1,0),2), pow(N(1,1),2), pow(N(1,2),2), 2*N(1,0)*N(1,1), 2*N(1,0)*N(1,2), 2*N(1,2)*N(1,1),
                    pow(N(2,0),2), pow(N(2,1),2), pow(N(2,2),2), 2*N(2,0)*N(2,1), 2*N(2,0)*N(2,2), 2*N(2,2)*N(2,1),
                    N(0,0)*N(1,0), N(0,1)*N(1,1), N(0,2)*N(1,2), N(0,0)*N(1,1)+N(0,1)*N(1,0), N(0,0)*N(1,2)+N(0,2)*N(1,0), N(0,1)*N(1,2)+N(0,2)*N(1,1),
                    N(2,0)*N(0,0), N(2,1)*N(0,1), N(2,2)*N(0,2), N(2,0)*N(0,1)+N(2,1)*N(0,0), N(2,0)*N(0,2)+N(2,2)*N(0,0), N(2,1)*N(0,2)+N(2,2)*N(0,1),
                    N(1,0)*N(2,0), N(1,1)*N(2,1), N(1,2)*N(2,2), N(1,0)*N(2,1)+N(1,1)*N(2,0), N(1,0)*N(2,2)+N(1,2)*N(2,0), N(1,1)*N(2,2)+N(1,2)*N(2,1);
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

      Marmot::Vector6d rotateVoigtStress( const Eigen::Matrix3d& Q, const Marmot::Vector6d& voigtStress )
      {
        const Matrix3d& T  = voigtToStress( voigtStress );
        const Matrix3d& TR = Q * T * Q.transpose();
        return stressToVoigt( TR );
      }

      Marmot::Vector6d transformStressToLocalSystem( const Marmot::Vector6d& stress,
                                                     const Matrix3d&         transformedCoordinateSystem )
      {
        const Matrix3d s_prime = Math::transformToLocalSystem( voigtToStress( stress ), transformedCoordinateSystem );

        return stressToVoigt( s_prime );
      }

      Marmot::Vector6d transformStrainToLocalSystem( const Marmot::Vector6d& strain,
                                                     const Matrix3d&         transformedCoordinateSystem )
      {
        const Matrix3d e_prime = Math::transformToLocalSystem( voigtToStrain( strain ), transformedCoordinateSystem );

        return strainToVoigt( e_prime );
      }

      Marmot::Vector6d transformStressToGlobalSystem( const Marmot::Vector6d& stress,
                                                      const Matrix3d&         transformedCoordinateSystem )
      {
        const Matrix3d s_prime = Math::transformToGlobalSystem( voigtToStress( stress ), transformedCoordinateSystem );

        return stressToVoigt( s_prime );
      }
      Marmot::Vector6d transformStrainToGlobalSystem( const Marmot::Vector6d& strain,
                                                      const Matrix3d&         transformedCoordinateSystem )
      {
        const Matrix3d e_prime = Math::transformToGlobalSystem( voigtToStrain( strain ), transformedCoordinateSystem );

        return strainToVoigt( e_prime );
      }

      Matrix6d transformStiffnessToGlobalSystem( const Marmot::Matrix6d& stiffness,
                                                 const Matrix3d&         transformedCoordinateSystem )
      {
        const EigenTensors::Tensor3333d stiffnessTensorLocal = voigtToStiffness( stiffness );
        EigenTensors::Tensor3333d       stiffnessTensorGlobal;
        stiffnessTensorGlobal.setZero();
        Matrix3d N = transformedCoordinateSystem.transpose();

        for ( size_t i = 0; i < 3; i++ )
          for ( size_t j = 0; j < 3; j++ )
            for ( size_t k = 0; k < 3; k++ )
              for ( size_t l = 0; l < 3; l++ )
                for ( size_t m = 0; m < 3; m++ )
                  for ( size_t n = 0; n < 3; n++ )
                    for ( size_t o = 0; o < 3; o++ )
                      for ( size_t p = 0; p < 3; p++ )
                        stiffnessTensorGlobal( i, j, k, l ) += N( i, m ) * N( j, n ) * N( k, o ) * N( l, p ) *
                                                               stiffnessTensorLocal( m, n, o, p );

        Matrix6d res = stiffnessToVoigt( stiffnessTensorGlobal );
        return res;
      }

    } // namespace Transformations
  }   // namespace ContinuumMechanics::VoigtNotation
} // namespace Marmot
