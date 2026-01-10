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
    if ( std::abs( a - b ) / ( std::abs( a ) > 1 ? std::abs( a ) : 1.0 ) > tol ) {
      std::cout << "absolute error: " << std::abs( a - b ) << " > " << tol << std::endl;
      std::cout << "relative error: " << std::abs( a - b ) / std::abs( a ) << " > " << tol << std::endl;
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

  bool spinTurbokreisel( Marmot::Solvers::MarmotMaterialPointSolverHypoElastic& solver,
                         const double                                           stressTol,
                         const double                                           stiffnessTol )
  {
    using namespace Eigen;
    using namespace ContinuumMechanics::VoigtNotation;

    solver.solve();
    auto           history = solver.getHistory();
    const Vector6d refStress( history.back().stress );
    const Matrix6d refStiffness( history.back().dStressdStrain );

    const int                     N   = 100;
    Eigen::Matrix< double, N, 2 > pts = fibonacciLatticeHemisphere< N >();

    Eigen::Vector2d pt;
    double          phi, theta;

    // modify steps to account for rotation
    const auto steps = solver.getSteps();
    // check if at least one step exists
    if ( steps.size() == 0 ) {
      std::cout << "No steps found in solver." << std::endl;
      return false;
    }

    for ( auto& step : steps ) {
      // must be pure strain control, i.e. all strain increment components are controlled
      if ( step.isStrainComponentControlled.any() == false ) {
        std::cout << "TURBOKREISEL TEST REQUIRES PURE STRAIN CONTROLLED STEPS." << std::endl;
        return false;
      };
    }

    for ( int i1 = 0; i1 < N; i1++ ) {
      pt    = pts.row( i1 );
      phi   = pt[0];
      theta = pt[1];

      // transformation matrix
      AngleAxisd phiX( theta, Vector3d::UnitX() );
      AngleAxisd phiY( phi, Vector3d::UnitY() );
      AngleAxisd phiZ( 0, Vector3d::UnitZ() );

      Matrix3d e = ( phiX * phiY * phiZ ).matrix().transpose();
      /* Matrix3d e_ = e.transpose(); */

      std::cout << "coordinate system: iteration " << i1 << "\n" << e << "\n" << std::endl;

      // reset material state
      solver.resetToInitialState();

      // clear previous steps
      solver.clearSteps();

      // modify steps to account for rotation
      for ( auto& step : steps ) {
        auto modifiedStep = step;
        modifiedStep
          .strainIncrementTarget = Transformations::transformStrainToLocalSystem( Vector6d(
                                                                                    step.strainIncrementTarget ),
                                                                                  e );
        std::cout << "modified strain increment target:\n"
                  << Vector6d( modifiedStep.strainIncrementTarget ).transpose() << "\n"
                  << std::endl;
        solver.addStep( modifiedStep );
      };

      // solve in new coordinate system
      solver.solve();

      auto history_ = solver.getHistory();

      const Vector6d S  = Vector6d( history_.back().stress );
      const Vector6d S_ = Transformations::transformStressToGlobalSystem( S, e );

      // check stress
      if ( !checkIfEqual( MatrixXd( S_ ), MatrixXd( refStress ), stressTol ) ) {
        std::cout << "\nTURBOKREISEL STRESS CHECK FAILED DURING ITERATION " << i1 << "\n" << std::endl;
        std::cout << "stress in new coordinate system:\n" << S_.transpose() << "\n" << std::endl;
        std::cout << "stress in reference coordinate system:\n" << S.transpose() << "\n" << std::endl;

        return false;
      };

      // get tangent in new coordinate system
      Matrix6d dS_dE  = Matrix6d( history_.back().dStressdStrain );
      Matrix6d dS_dE_ = Transformations::transformStiffnessToGlobalSystem( dS_dE, e );

      // check tangent
      if ( !checkIfEqual( MatrixXd( dS_dE_ ), MatrixXd( refStiffness ), stiffnessTol ) ) {
        std::cout << "\nTURBOKREISEL TANGENT CHECK FAILED DURING ITERATION " << i1 << "\n" << std::endl;
        std::cout << "tangent in new coordinate system:\n" << dS_dE.transpose() << "\n" << std::endl;
        std::cout << "tangent in reference coordinate system:\n" << dS_dE_.transpose() << "\n" << std::endl;

        return false;
      };
    };
    return true;
  };

} // namespace Marmot::Testing
