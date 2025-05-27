#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include <iostream>
#include <ostream>

namespace Marmot::Testing {

  std::string getString( const double a )
  {
    return std::to_string( a );
  }

  std::string getString( const autodiff::dual a )
  {
    return "(" + std::to_string( a.val ) + ",  " + std::to_string( a.grad ) + ")";
  }

  bool checkIfEqual( const double a, const double b, const double tol )
  {
    // check for NaN values and return false if any of the values is NaN
    if ( std::isnan( a ) || std::isnan( b ) ) {
      std::cout << " Hint: a = " << a << " or "
                << " b = " << b << " is NaN" << std::endl;
      return false;
    }

    // check for inf values and return false if any of the values is inf
    if ( std::isinf( a ) || std::isinf( b ) ) {
      std::cout << " Hint: a = " << a << " or "
                << " b = " << b << " is inf" << std::endl;
      return false;
    }

    // check if values are equal with respect to the tolerance
    // relative error is used for non-zero values
    if ( std::abs( a - b ) / ( std::abs( a ) > 1e-10 ? std::abs( a ) : 1.0 ) > tol ) {
      std::cout << " Hint: a = " << a << " != "
                << " b = " << b << " ( tol = " << tol << " )" << std::endl;
      return false;
    }
    else
      return true;
  }

  bool checkIfEqual( const autodiff::dual a, const autodiff::dual b, const double tol )
  {
    if ( checkIfEqual( a.val, b.val, tol ) && checkIfEqual( a.grad, b.grad, tol ) ) {
      return true;
    }
    return false;
  }

  bool checkIfEqual( const std::complex< double > a, const std::complex< double > b, const double tol )
  {
    // Check equality for the real and imaginary parts independently
    if ( checkIfEqual( a.real(), b.real(), tol ) && checkIfEqual( a.imag(), b.imag(), tol ) ) {
      return true;
    }
    return false;
  }

  void throwExceptionOnFailure( const bool condition, const std::string& message )
  {
    if ( !condition ) {
      throw std::runtime_error( message );
    }
  }

  void executeTestsAndCollectExceptions( const std::vector< std::function< void() > >& testFunctions )
  {

    const auto length = testFunctions.size();

    auto exceptions = std::vector< std::string >( length );

    bool allPassed = true;

    for ( const auto& testFunction : testFunctions ) {
      try {
        testFunction();
      }
      catch ( const std::exception& e ) {
        allPassed = false;
        exceptions.push_back( e.what() );
      }
    }

    for ( const auto& exception : exceptions ) {
      if ( !exception.empty() ) {
        std::cout << "Exception: " << exception << std::endl;
      }
    }

    if ( !allPassed ) {
      throw std::runtime_error( "some tests failed" );
    }
  }

  bool spinTurbokreisel( const std::unique_ptr< MarmotMaterialHypoElastic >& material,
                         double*                                             stress,
                         double*                                             dStress_dStrain,
                         const double*                                       dStrain,
                         const double*                                       timeOld,
                         const double                                        dT,
                         double&                                             pNewDT,
                         const double                                        stressTol,
                         const double                                        stiffnessTol )
  {
    using namespace Eigen;
    using namespace ContinuumMechanics::VoigtNotation;

    const Vector6d initialStress( stress );
    const Vector6d dE( dStrain );

    // get material state
    int                         nStateVars = material->getNumberOfAssignedStateVars();
    double*                     stateVars  = material->getAssignedStateVars();
    const std::vector< double > initialStateVarVector( stateVars, stateVars + nStateVars );

    // evaluate material in reference coordinates to obtain reference solution
    material->computeStress( stress, dStress_dStrain, dStrain, timeOld, dT, pNewDT );
    const Vector6d refStress( stress );
    const Matrix6d refStiffness( dStress_dStrain );

    std::cout << "reference stress:\n" << refStress.transpose() << "\n" << std::endl;
    std::cout << "reference tangent:\n" << refStiffness.transpose() << "\n" << std::endl;

    const int                     N   = 100;
    Eigen::Matrix< double, N, 2 > pts = fibonacciLatticeHemisphere< N >();

    Eigen::Vector2d pt;
    double          phi, theta;

    for ( int i1 = 0; i1 < N; i1++ ) {
      pt    = pts.row( i1 );
      phi   = pt[0];
      theta = pt[1];

      // transformation matrix
      AngleAxisd phiX( theta, Vector3d::UnitX() );
      AngleAxisd phiY( phi, Vector3d::UnitY() );
      AngleAxisd phiZ( 0, Vector3d::UnitZ() );

      Matrix3d e  = ( phiX * phiY * phiZ ).matrix().transpose();
      Matrix3d e_ = e.transpose();

      // reset material state
      std::vector< double > stateVarVector( initialStateVarVector );
      material->assignStateVars( &stateVarVector[0], nStateVars );

      // transform strain increment to new coordinate system
      Matrix3d dStrain = voigtToStrain( dE );

      Matrix3d dStrain_ = Matrix3d::Zero();
      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          for ( int k = 0; k < 3; k++ )
            for ( int l = 0; l < 3; l++ ) {
              dStrain_( i, j ) += e( i, k ) * e( j, l ) * dStrain( k, l );
            };

      Vector6d dE_ = strainToVoigt( dStrain_ );

      // transform initial stress to new coordinate system
      Vector6d S( initialStress );
      Matrix3d Stress = voigtToStress( S );

      Matrix3d Stress_ = Matrix3d::Zero();
      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          for ( int k = 0; k < 3; k++ )
            for ( int l = 0; l < 3; l++ ) {
              Stress_( i, j ) += e( i, k ) * e( j, l ) * Stress( k, l );
            };

      Vector6d S_ = stressToVoigt( Stress_ );

      // evaluate material in new coordinate system
      Matrix6d dS_dE_ = Matrix6d::Zero();
      material->computeStress( S_.data(), dS_dE_.data(), dE_.data(), timeOld, dT, pNewDT );

      // transform stress to reference coordinate system
      Stress_ = voigtToStress( S_ );

      Stress.setZero();
      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          for ( int k = 0; k < 3; k++ )
            for ( int l = 0; l < 3; l++ ) {
              Stress( i, j ) += e_( i, k ) * e_( j, l ) * Stress_( k, l );
            };

      S = stressToVoigt( Stress );

      // check stress
      if ( !checkIfEqual( MatrixXd( S ), MatrixXd( refStress ), stressTol ) ) {
        std::cout << "\nTURBOKREISEL STRESS CHECK FAILED DURING ITERATION " << i1 << "\n" << std::endl;
        std::cout << "stress in new coordinate system:\n" << S_.transpose() << "\n" << std::endl;
        std::cout << "stress in reference coordinate system:\n" << S.transpose() << "\n" << std::endl;

        return false;
      };

      // transform tangent to reference coordinate system
      const Eigen::Tensor< double, 4 > Stiffness_ = voigtToStiffness( dS_dE_ );
      EigenTensors::Tensor3333d        Stiffness;

      Stiffness.setZero();
      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          for ( int k = 0; k < 3; k++ )
            for ( int l = 0; l < 3; l++ )
              for ( int m = 0; m < 3; m++ )
                for ( int n = 0; n < 3; n++ )
                  for ( int o = 0; o < 3; o++ )
                    for ( int p = 0; p < 3; p++ ) {
                      Stiffness( i, j, k, l ) += e_( i, m ) * e_( j, n ) * e_( k, o ) * e_( l, p ) *
                                                 Stiffness_( m, n, o, p );
                    };

      Matrix6d dS_dE = stiffnessToVoigt( Stiffness );

      // check tangent
      if ( !checkIfEqual( MatrixXd( dS_dE ), MatrixXd( refStiffness ), stiffnessTol ) ) {
        std::cout << "\nTURBOKREISEL TANGENT CHECK FAILED DURING ITERATION " << i1 << "\n" << std::endl;
        std::cout << "tangent in new coordinate system:\n" << dS_dE_.transpose() << "\n" << std::endl;
        std::cout << "tangent in reference coordinate system:\n" << dS_dE.transpose() << "\n" << std::endl;

        return false;
      };
    };
    return true;
  };

} // namespace Marmot::Testing
