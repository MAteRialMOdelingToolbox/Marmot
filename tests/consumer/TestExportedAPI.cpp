/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck,
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
 *
 * This file is part of the MAteRialMOdellingToolbox (marmot).
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of marmot.
 * ---------------------------------------------------------------------
 */

/** @file TestExportedAPI.cpp
 * @brief Uses Marmot the way a consumer (EdelweissFE, the Abaqus interface) does: linked against the shared library
 * only, through the exported API only.
 *
 * @details Every other test links Marmot's object files and so reaches internals that the shared library does not
 * export on Windows (or with MARMOT_EXPORT_API_ONLY), see MARMOT_API in MarmotPortability.h. This test is the one
 * that fails when the exported API is incomplete: a missing MARMOT_API, or a header whose inline definitions depend
 * on a symbol the library does not export. It therefore uses no test helper from the library (MarmotTesting.h is
 * not exported) and checks against closed-form results only.
 */

#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotKinematics.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotVoigt.h"
#include <Eigen/Dense>
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

  /** @brief Throws if @p condition is false; the consumer-side replacement for MarmotTesting's helpers.
   * @param condition The condition that must hold.
   * @param message Reported when it does not.
   */
  void require( bool condition, const std::string& message )
  {
    if ( !condition )
      throw std::runtime_error( message );
  }

  /** @brief Checks that the header-defined constants are usable without any symbol from the library. */
  void testHeaderDefinedConstants()
  {
    using namespace Marmot::ContinuumMechanics;

    // dOmega/dl (i,j,k,l) = 1/2 ( d_ik d_jl - d_kj d_il )
    require( Kinematics::VelocityGradient::dOmega_dVelocityGradient( 0, 1, 0, 1 ) == 0.5,
             "dOmega_dVelocityGradient( 0, 1, 0, 1 ) must be 1/2" );
    require( Kinematics::VelocityGradient::dOmega_dVelocityGradient( 0, 1, 1, 0 ) == -0.5,
             "dOmega_dVelocityGradient( 0, 1, 1, 0 ) must be -1/2" );
    // dD/dl in engineering Voigt notation: the shear row (0,1) -> 3 carries the factor 2
    require( Kinematics::VelocityGradient::dStretchingRate_dVelocityGradient( 3, 0, 1 ) == 1.0,
             "dStretchingRate_dVelocityGradient( 3, 0, 1 ) must be 1" );
    require( Kinematics::VelocityGradient::dStretchingRate_dVelocityGradient( 0, 0, 0 ) == 1.0,
             "dStretchingRate_dVelocityGradient( 0, 0, 0 ) must be 1" );

    require( VoigtNotation::P( 3 ) == 2.0 && VoigtNotation::PInv( 3 ) == 0.5, "VoigtNotation::P / PInv" );
    require( VoigtNotation::I( 0 ) == 1.0 && VoigtNotation::I( 3 ) == 0.0, "VoigtNotation::I" );
    require( VoigtNotation::IDev.trace() == 5.0, "VoigtNotation::IDev must have trace 5 (= 6 - 1)" );
  }

  /** @brief Creates a material through the factory and checks its stress response against the elastic tangent. */
  void testMaterialThroughFactory()
  {
    const double                E = 210e3, nu = 0.3;
    const std::vector< double > properties{ E, nu };

    std::unique_ptr< MarmotMaterialHypoElastic > material(
      MarmotLibrary::MarmotMaterialHypoElasticFactory::createMaterial( "LINEARELASTIC",
                                                                       properties.data(),
                                                                       static_cast< int >( properties.size() ),
                                                                       1 ) );
    require( material != nullptr, "the factory must create LINEARELASTIC" );

    const int             nStateVars = material->getNumberOfRequiredStateVars();
    std::vector< double > stateVars( nStateVars > 0 ? nStateVars : 1, 42.0 );
    material->initializeYourself( stateVars.data(), nStateVars );
    material->setCharacteristicElementLength( 1.0 ); // the one exported member of the material interface

    MarmotMaterialHypoElastic::state3D state;
    state.stateVars                  = stateVars.data();
    Marmot::Matrix6d dStress_dStrain = Marmot::Matrix6d::Zero();
    Marmot::Vector6d dStrain         = Marmot::Vector6d::Zero();
    dStrain( 0 )                     = 1e-3;
    material->computeStress( state, dStress_dStrain, dStrain, { 0.0, 1.0 } );

    const double lambda = E * nu / ( ( 1 + nu ) * ( 1 - 2 * nu ) );
    const double mu     = E / ( 2 * ( 1 + nu ) );
    const double tol    = 1e-12 * E;
    require( std::abs( state.stress( 0 ) - ( lambda + 2 * mu ) * dStrain( 0 ) ) < tol,
             "uniaxial strain: sigma_11 must be (lambda + 2 mu) eps_11" );
    require( std::abs( state.stress( 1 ) - lambda * dStrain( 0 ) ) < tol,
             "uniaxial strain: sigma_22 must be lambda eps_11" );
    require( std::abs( state.stress( 3 ) ) < tol, "uniaxial strain: no shear stress" );
    require( std::abs( dStress_dStrain( 0, 0 ) - ( lambda + 2 * mu ) ) < tol, "tangent C_1111" );
  }

  /** @brief Creates an element through the factory, as EdelweissFE does, and checks the assembled kernels. */
  void testElementThroughFactory()
  {
    std::unique_ptr< MarmotElement > element( MarmotLibrary::MarmotElementFactory::createElement( "CPS4", 1 ) );
    require( element != nullptr, "the factory must create CPS4" );
    require( element->getNNodes() == 4 && element->getNSpatialDimensions() == 2, "CPS4 is a 4-node 2D element" );

    const std::vector< double > coordinates{ 0, 0, 1, 0, 1, 1, 0, 1 }; // unit square
    const std::vector< double > elementProperties{ 1.0 };              // thickness
    const std::vector< double > materialProperties{ 210e3, 0.3 };
    element->assignNodeCoordinates( coordinates.data() );
    element->assignProperty(
      ElementProperties( elementProperties.data(), static_cast< int >( elementProperties.size() ) ) );
    element->assignProperty( MarmotMaterialSection( "LINEARELASTIC",
                                                    materialProperties.data(),
                                                    static_cast< int >( materialProperties.size() ) ) );
    std::vector< double > stateVars( element->getNumberOfRequiredStateVars(), 0.0 );
    element->assignStateVars( stateVars.data(), static_cast< int >( stateVars.size() ) );
    element->initializeYourself();

    const int nDof = element->getNDofPerElement();
    require( nDof == 8, "CPS4 has 8 degrees of freedom" );

    // Stiffness at the undeformed state
    Eigen::VectorXd Q = Eigen::VectorXd::Zero( nDof ), P = Eigen::VectorXd::Zero( nDof );
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero( nDof, nDof );
    element->computeKernels( Q.data(), Q.data(), P.data(), K.data(), 0.0, 1.0 );

    const double tol = 1e-12 * K.cwiseAbs().maxCoeff();
    require( K.cwiseAbs().maxCoeff() > 0.0, "the stiffness must not vanish" );
    require( P.cwiseAbs().maxCoeff() == 0.0, "no internal forces without deformation" );
    require( ( K - K.transpose() ).cwiseAbs().maxCoeff() < tol, "the stiffness must be symmetric" );
    require( ( K * Eigen::VectorXd::Ones( nDof ) ).cwiseAbs().maxCoeff() < tol,
             "a rigid translation must not produce internal forces" );

    // Linear elasticity: the internal forces of a fresh element are P = K Q
    std::fill( stateVars.begin(), stateVars.end(), 0.0 );
    element->initializeYourself();
    for ( int i = 0; i < nDof; i++ )
      Q( i ) = 1e-3 * std::sin( 1.0 + i );
    Eigen::VectorXd PQ = Eigen::VectorXd::Zero( nDof );
    Eigen::MatrixXd KQ = Eigen::MatrixXd::Zero( nDof, nDof );
    element->computeKernels( Q.data(), Q.data(), PQ.data(), KQ.data(), 0.0, 1.0 );
    require( ( PQ - K * Q ).cwiseAbs().maxCoeff() < tol * 1e-3,
             "internal forces must equal K Q for linear elasticity" );
    require( ( KQ - K ).cwiseAbs().maxCoeff() < tol, "the stiffness must not depend on the displacement" );
  }

} // namespace

/** @brief Runs the consumer-side checks; all are run, every failure is reported.
 * @return 0 if all checks pass, 1 otherwise.
 */
int main()
{
  MarmotJournal::setMSGOutputDirection( std::cout );

  const std::vector< std::pair< std::string, std::function< void() > > > tests{
    { "header-defined constants", testHeaderDefinedConstants },
    { "material through the factory", testMaterialThroughFactory },
    { "element through the factory", testElementThroughFactory },
  };

  int nFailed = 0;
  for ( const auto& [name, test] : tests ) {
    try {
      test();
      std::cout << "passed: " << name << "\n";
    }
    catch ( const std::exception& e ) {
      std::cout << "FAILED: " << name << ": " << e.what() << "\n";
      nFailed++;
    }
  }
  return nFailed == 0 ? 0 : 1;
}
