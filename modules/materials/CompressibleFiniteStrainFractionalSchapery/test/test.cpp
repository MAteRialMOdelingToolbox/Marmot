/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck
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

#include "Marmot/CompressibleFiniteStrainFractionalSchapery.h"
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
#include "Marmot/MarmotMaterialPointSolverFiniteStrain.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTesting.h"
#include <memory>
#include <string>
#include <vector>

using namespace Marmot::Testing;
using namespace Marmot::Materials;
using namespace Marmot::Math::Tensors::FastorTensors;
using namespace Marmot::Math::Tensors::FastorIndices;
using namespace Marmot::ContinuumMechanics;

// Helper property builders
// Layout: [hyperelasticBase, onlyShearCreep, <elastic props>, alpha, tau0, g_inf, a1, a2, b_S]

static std::vector< double > getNeoHookeProps()
{
  // Base 0 (NeoHooke): mu=100.0, kappa=300.0, alpha=0.3, tau0=5.0, g_inf=0.2, a1=0.1, a2=0.05, b_S=0.02
  return { 0.0, 1.0, 100.0, 300.0, 0.3, 5.0, 0.2, 0.1, 0.05, 0.02 };
}

static std::vector< double > getMooneyRivlinProps()
{
  // Base 2 (MooneyRivlin): C10=40.0, C01=10.0, kappa=300.0, alpha=0.25, tau0=10.0, g_inf=0.15, a1=0.0, a2=0.0, b_S=0.0
  return { 2.0, 1.0, 40.0, 10.0, 300.0, 0.25, 10.0, 0.15, 0.0, 0.0, 0.0 };
}

static std::vector< double > getYeohProps()
{
  // Base 1 (Yeoh): C10=50.0, C20=-2.0, C30=0.5, kappa=300.0, alpha=0.35, tau0=2.0, g_inf=0.1, a1=0.0, a2=0.0, b_S=0.0
  return { 1.0, 1.0, 50.0, -2.0, 0.5, 300.0, 0.35, 2.0, 0.1, 0.0, 0.0, 0.0 };
}

static Marmot::Materials::MaterialPointSolverFiniteStrain makeSolver( const std::string&     matName,
                                                                      std::vector< double >& matProps )
{
  auto solveropts = Marmot::Materials::MaterialPointSolverFiniteStrain::SolverOptions();
  return Marmot::Materials::MaterialPointSolverFiniteStrain( matName,
                                                             matProps.data(),
                                                             static_cast< int >( matProps.size() ),
                                                             solveropts );
}

static Marmot::Materials::MaterialPointSolverFiniteStrain::Step makeStep( const Tensor33d& gradUIncrement,
                                                                          double           timeStart,
                                                                          double           timeEnd,
                                                                          double           dT )
{
  Marmot::Materials::MaterialPointSolverFiniteStrain::Step step;
  step.gradUIncrementTarget        = gradUIncrement;
  step.stressIncrementTarget       = Tensor33d( 0.0 );
  step.isGradUComponentControlled  = Tensor33t< bool >( true );
  step.isStressComponentControlled = Tensor33t< bool >( false );
  step.timeStart                   = timeStart;
  step.timeEnd                     = timeEnd;
  step.dTStart                     = dT;
  step.dTMax                       = dT;
  return step;
}

// Test 1: Factory registration
void testFactoryRegistration()
{
  auto                                                       props = getNeoHookeProps();
  const std::string                                          name  = "COMPRESSIBLEFINITESTRAINFRACTIONALSCHAPERY";
  std::unique_ptr< Marmot::Materials::MaterialFiniteStrain > mat(
    Marmot::Registry::MaterialFiniteStrainFactory::createMaterial( name,
                                                                   props.data(),
                                                                   static_cast< int >( props.size() ),
                                                                   1 ) );
  throwExceptionOnFailure( mat != nullptr, "Material " + name + " could not be created from factory" );

  const std::string nameSub = "COMPRESSIBLEFINITESTRAINFRACTIONALSCHAPERY_SUBSTEPPED";
  std::unique_ptr< Marmot::Materials::MaterialFiniteStrain > matSub(
    Marmot::Registry::MaterialFiniteStrainFactory::createMaterial( nameSub,
                                                                   props.data(),
                                                                   static_cast< int >( props.size() ),
                                                                   2 ) );
  throwExceptionOnFailure( matSub != nullptr, "Substepped material " + nameSub + " could not be created from factory" );
}

// Test 2: Stress symmetry under general deformation
void testStressSymmetry()
{
  const std::string matName  = "COMPRESSIBLEFINITESTRAINFRACTIONALSCHAPERY";
  auto              matProps = getNeoHookeProps();
  auto              solver   = makeSolver( matName, matProps );

  Tensor33d gradU( 0.0 );
  gradU( 0, 0 ) = 0.05;
  gradU( 0, 1 ) = 0.02;
  gradU( 1, 0 ) = 0.02;
  gradU( 1, 1 ) = -0.01;
  gradU( 2, 2 ) = -0.02;

  solver.addStep( makeStep( gradU, 0.0, 1.0, 0.1 ) );
  solver.solve();

  auto tau = solver.getHistory().back().stress;
  throwExceptionOnFailure( checkIfEqual( tau, Fastor::transpose( tau ), 1e-10 ), "Kirchhoff stress must be symmetric" );
}

// Test 3: Pure rotation produces zero stress
void testPureRotationZeroStress()
{
  const std::string matName  = "COMPRESSIBLEFINITESTRAINFRACTIONALSCHAPERY";
  auto              matProps = getNeoHookeProps();

  for ( int phi_deg = 0; phi_deg <= 180; phi_deg += 45 ) {
    const double phi = Marmot::Math::degToRad( phi_deg );

    Tensor33d F( 0.0 );
    F( 0, 0 ) = std::cos( phi );
    F( 0, 1 ) = -std::sin( phi );
    F( 1, 0 ) = std::sin( phi );
    F( 1, 1 ) = std::cos( phi );
    F( 2, 2 ) = 1.0;

    auto solver = makeSolver( matName, matProps );
    solver.addStep( makeStep( F - Spatial3D::I, 0.0, 1.0, 1.0 ) );
    solver.solve();

    throwExceptionOnFailure( checkIfEqual( solver.getHistory().back().stress, Tensor33d( 0.0 ), 1e-10 ),
                             "Pure rotation must produce zero stress" );
  }
}

// Test 4: Objectivity: tau(Q * F) = Q * tau(F) * Q^T
void testObjectivity()
{
  const std::string matName  = "COMPRESSIBLEFINITESTRAINFRACTIONALSCHAPERY";
  auto              matProps = getNeoHookeProps();

  Tensor33d F( 0.0 );
  F( 0, 0 ) = 1.05;
  F( 0, 1 ) = 0.02;
  F( 1, 0 ) = 0.02;
  F( 1, 1 ) = 0.98;
  F( 2, 2 ) = 0.96;

  auto solverRef = makeSolver( matName, matProps );
  solverRef.addStep( makeStep( F - Spatial3D::I, 0.0, 1.0, 1.0 ) );
  solverRef.solve();
  Tensor33d tauRef = solverRef.getHistory().back().stress;

  for ( int phi_deg = 30; phi_deg <= 90; phi_deg += 30 ) {
    const double phi = Marmot::Math::degToRad( phi_deg );

    Tensor33d Q( 0.0 );
    Q( 0, 0 ) = std::cos( phi );
    Q( 0, 1 ) = -std::sin( phi );
    Q( 1, 0 ) = std::sin( phi );
    Q( 1, 1 ) = std::cos( phi );
    Q( 2, 2 ) = 1.0;

    Tensor33d F_rot  = einsum< ik, kj, to_ij >( Q, F );
    auto      solver = makeSolver( matName, matProps );
    solver.addStep( makeStep( F_rot - Spatial3D::I, 0.0, 1.0, 1.0 ) );
    solver.solve();
    Tensor33d tauRot = solver.getHistory().back().stress;

    Tensor33d tauExpected = einsum< iI, IJ, jJ, to_ij >( Q, tauRef, Q );

    throwExceptionOnFailure( checkIfEqual( tauRot, tauExpected, 1e-8 ),
                             "Objectivity failed for phi=" + std::to_string( phi_deg ) );
  }
}

// Test 5: Viscoelastic stress relaxation during hold
void testViscoelasticRelaxation()
{
  const std::string matName  = "COMPRESSIBLEFINITESTRAINFRACTIONALSCHAPERY";
  auto              matProps = getNeoHookeProps();
  auto              solver   = makeSolver( matName, matProps );

  // Ramp phase: 0 to 1s
  Tensor33d gradURamp( 0.0 );
  gradURamp( 0, 0 ) = 0.10;
  solver.addStep( makeStep( gradURamp, 0.0, 1.0, 0.1 ) );

  // Hold phase: 1 to 50s
  Tensor33d gradUHold( 0.0 );
  solver.addStep( makeStep( gradUHold, 1.0, 50.0, 1.0 ) );

  solver.solve();

  const auto& history = solver.getHistory();

  // Find peak stress at ramp end (t=1.0)
  double peakStress = 0.0;
  for ( const auto& h : history ) {
    if ( std::abs( h.time - 1.0 ) < 1e-6 ) {
      peakStress = h.stress( 0, 0 );
      break;
    }
  }

  double endStress = history.back().stress( 0, 0 );

  throwExceptionOnFailure( peakStress > 0.0, "Peak stress must be positive" );
  throwExceptionOnFailure( endStress < peakStress, "Stress must relax during hold" );
  throwExceptionOnFailure( endStress > 0.0, "Relaxed stress must remain positive" );
}

// Test 6: Multiple hyperelastic bases (NeoHooke, MooneyRivlin, Yeoh)
void testBases()
{
  const std::string matName = "COMPRESSIBLEFINITESTRAINFRACTIONALSCHAPERY";

  std::vector< std::vector< double > > allProps = {
    getNeoHookeProps(),
    getMooneyRivlinProps(),
    getYeohProps(),
  };

  Tensor33d gradU( 0.0 );
  gradU( 0, 0 ) = 0.05;

  for ( auto& props : allProps ) {
    auto solver = makeSolver( matName, props );
    solver.addStep( makeStep( gradU, 0.0, 1.0, 0.2 ) );
    solver.solve();

    double s11 = solver.getHistory().back().stress( 0, 0 );
    throwExceptionOnFailure( s11 > 0.0, "Base must produce positive axial stress" );
  }
}

int main()
{
  auto tests = std::vector< std::function< void() > >{
    testFactoryRegistration,
    testStressSymmetry,
    testPureRotationZeroStress,
    testObjectivity,
    testViscoelasticRelaxation,
    testBases,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
