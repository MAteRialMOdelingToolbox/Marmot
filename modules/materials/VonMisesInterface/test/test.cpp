// clang-format off
/**
 * Unit tests for VonMisesInterface.
 *
 * Shared parameters for all tests:
 *   E = 1e5,  nu = 0.3,  h = 1e-3
 *   G = E / (2*(1+nu)) = 38461.538...
 *   normal = {0, 0, 1}
 *
 * Voigt ordering: [S11, S22, S33, S12, S13, S23]
 * averageStress: flat row-major 3x3 buffer (9 doubles).
 *   voigtToStress fills [[v0,v3,v4],[v3,v1,v5],[v4,v5,v2]]
 *   stored row-major:    {v0,v3,v4, v3,v1,v5, v4,v5,v2}
 *
 * ── Elastic tests  (yieldStress = 1e10) ─────────────────────────────────────
 *
 *  A.  Displacement-jump only
 *      dU_top={0,1e-3,0}, dU_bot=0, normal={0,0,1}
 *      dJumpU = {0,1e-3,0}
 *      (1/h)*outer(dJumpU,n) -> dU_kl[1,2]=1  (h=1e-3, dU[1]=1e-3)
 *      Symmetrize: eps[1,2]=eps[2,1]=0.5
 *      Voigt gamma_23 (index 5) = 2*eps[1,2] = 1  ->  S23 = G*1 = G
 *      force = {0,G,0},  avgStress flat = {0,0,0, 0,0,G, 0,G,0}
 *
 *  B.  Surface-strain only
 *      dSurface_strain_top (row-major 3x3): [0,1]=[1,0]=1e-2, rest 0; bottom=0
 *      avg: [0,1]=[1,0]=5e-3  ->  Voigt gamma_12=1e-2  ->  S12=G*1e-2
 *      force = {0,0,0},  avgStress flat = {0,S12,0, S12,0,0, 0,0,0}
 *
 * ── Plastic tests  (yieldStress=100, HLin=10, deltaYieldStress=0) ───────────
 *
 *  C.  Displacement-jump only (same kinematics as A but smaller strain)
 *      dU_top={0,1e-4,0} -> gamma_23 = 1e-4/1e-3 = 0.1
 *      trial S23 = G*0.1 = 3846.15 >> tau_yield = 100/sqrt(3) = 57.74  ->  PLASTIC
 *      Newton return mapping:  S23 = 58.0633280
 *
 *  D.  Surface-strain only (same kinematics as B)
 *      dSurface_strain_top: [0,1]=[1,0]=1e-2 -> gamma_12 = 1e-2
 *      trial S12 = G*1e-2 = 384.615 >> tau_yield = 57.74  ->  PLASTIC
 *      Newton return mapping:  S12 = 57.7633540
 */
// clang-format on

#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotInterfaceMaterialHelperFunctions.h"
#include "Marmot/MarmotMaterialHypoElasticInterface.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotWiechertInterface.h"
#include "Marmot/VonMisesInterface.h"
#include <Eigen/Dense>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <iterator>

using namespace Marmot::Testing;
using namespace Marmot::ContinuumMechanics::Elasticity::Isotropic;
using namespace Marmot::ContinuumMechanics::Elasticity::TransverseIsotropic;
using namespace Marmot::Materials::InterfaceMaterialHelperFunctions;

// ---------------------------------------------------------------------------
// Helper: create a VonMisesInterface material object
// ---------------------------------------------------------------------------
std::unique_ptr< MarmotMaterialHypoElasticInterface > createMaterial( const double* props, int nProps )
{
  const int elLabel = 1;
  auto mat = std::unique_ptr< MarmotMaterialHypoElasticInterface >( dynamic_cast< MarmotMaterialHypoElasticInterface* >(
    MarmotLibrary::MarmotMaterialHypoElasticInterfaceFactory::createMaterial( "VONMISESINTERFACE",
                                                                              props,
                                                                              nProps,
                                                                              elLabel ) ) );
  if ( !mat )
    throw std::runtime_error( "VonMisesInterface: material creation failed" );
  return mat;
}

// ---------------------------------------------------------------------------
// Helper: initialise a fresh material and apply a single loading increment.
//   props[7]           = { E, nu, h, yieldStress, HLin, deltaYieldStress, delta }
//   dU[6]              = displacement increments (top[0:3], bottom[3:6])
//   dSurfaceStrain[18] = surface-strain increments (top[0:9], bottom[9:18])
//   normal[3]          = interface normal
//   Out: force[3], averageStress[9], pNewDT
// ---------------------------------------------------------------------------
void runSingleIncrement( const double* props,
                         const double* dU,
                         const double* dSurfaceStrain,
                         const double* normal,
                         double*       force,
                         double*       averageStress,
                         double&       pNewDT )
{
  auto mat = createMaterial( props, 7 );

  const int       nStateVars = mat->getNumberOfRequiredStateVars();
  Eigen::VectorXd stateVar( nStateVars );
  stateVar.setZero();
  mat->assignStateVars( stateVar.data(), nStateVars );

  double H_inv_ij[9]          = { 0 };
  double Z_ijkl[81]           = { 0 };
  double H_inv_nF_ijk[27]     = { 0 };
  double Yn_H_inv_Fn_ijkl[81] = { 0 };

  const double timeOld = 0.0;
  const double dT      = 1.0;

  mat->computeStress( force,
                      averageStress,
                      H_inv_ij,
                      Z_ijkl,
                      H_inv_nF_ijk,
                      Yn_H_inv_Fn_ijkl,
                      dU,
                      dSurfaceStrain,
                      normal,
                      &timeOld,
                      dT,
                      pNewDT );
}

// ===========================================================================
// TEST A – elastic response to a displacement jump
// ===========================================================================
void testElasticDisplacementJump()
{
  // props: [E, nu, h, yieldStress, HLin, deltaYieldStress, delta]
  // yieldStress = 1e10 ensures a purely elastic response
  const double props[7]           = { 1e5, 0.3, 1e-3, 1e10, 0., 0., 1. };
  const double dU[6]              = { 0., 1e-3, 0., 0., 0., 0. };
  const double dSurfaceStrain[18] = { 0. };
  const double normal[3]          = { 0., 0., 1. };

  double force[3]         = { 0. };
  double averageStress[9] = { 0. };
  double pNewDT           = 1.0;

  runSingleIncrement( props, dU, dSurfaceStrain, normal, force, averageStress, pNewDT );

  // G = E/(2*(1+nu)) = 1e5/2.6
  // gamma_23 = dJumpU[1]/h = 1e-3/1e-3 = 1  ->  S23 = G
  const double G   = 1e5 / ( 2. * 1.3 );
  const double S23 = G;

  // force[i] = S[i,2]:  {0, S23, 0}
  const double forceTarget[3] = { 0., S23, 0. };

  // voigtToStress with v[5]=S23:  [[0,0,0],[0,0,S23],[0,S23,0]]
  // flat row-major 3x3:            {0,0,0, 0,0,S23, 0,S23,0}
  const double avgTarget[9] = { 0., 0., 0., 0., 0., S23, 0., S23, 0. };

  Eigen::Map< const Eigen::Vector3d > forceVec( force );
  Eigen::Map< const Eigen::Vector3d > forceTgt( forceTarget );
  Eigen::Map< const Eigen::VectorXd > avgVec( averageStress, 9 );
  Eigen::Map< const Eigen::VectorXd > avgTgt( avgTarget, 9 );

  throwExceptionOnFailure( pNewDT == 1.0, "Newton diverged in " + std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( checkIfEqual< double >( forceVec, forceTgt, 1e-6 ),
                           "force mismatch in " + std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( checkIfEqual< double >( avgVec, avgTgt, 1e-6 ),
                           "averageStress mismatch in " + std::string( __PRETTY_FUNCTION__ ) );
}

// ===========================================================================
// TEST B – elastic response to a surface-strain increment
// ===========================================================================
void testElasticSurfaceStrain()
{
  const double props[7] = { 1e5, 0.3, 1e-3, 1e10, 0., 0., 1. };
  const double dU[6]    = { 0. };
  // top 3x3 row-major: [0,1]=1e-2, [1,0]=1e-2, rest 0;  bottom = 0
  const double dSurfaceStrain[18] = { 0., 1e-2, 0., 1e-2, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
  const double normal[3]          = { 0., 0., 1. };

  double force[3]         = { 0. };
  double averageStress[9] = { 0. };
  double pNewDT           = 1.0;

  runSingleIncrement( props, dU, dSurfaceStrain, normal, force, averageStress, pNewDT );

  // avg[0,1]=avg[1,0]=0.5*1e-2=5e-3
  // Voigt gamma_12=2*5e-3=1e-2  ->  S12 = G*1e-2
  const double G   = 1e5 / ( 2. * 1.3 );
  const double S12 = G * 1e-2;

  // force = S . n = S[:,2] = {0,0,0}
  const double forceTarget[3] = { 0., 0., 0. };

  // voigtToStress with v[3]=S12:  [[0,S12,0],[S12,0,0],[0,0,0]]
  // flat row-major 3x3:            {0,S12,0, S12,0,0, 0,0,0}
  const double avgTarget[9] = { 0., S12, 0., S12, 0., 0., 0., 0., 0. };

  Eigen::Map< const Eigen::Vector3d > forceVec( force );
  Eigen::Map< const Eigen::Vector3d > forceTgt( forceTarget );
  Eigen::Map< const Eigen::VectorXd > avgVec( averageStress, 9 );
  Eigen::Map< const Eigen::VectorXd > avgTgt( avgTarget, 9 );

  throwExceptionOnFailure( pNewDT == 1.0, "Newton diverged in " + std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( checkIfEqual< double >( forceVec, forceTgt, 1e-6 ),
                           "force mismatch in " + std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( checkIfEqual< double >( avgVec, avgTgt, 1e-6 ),
                           "averageStress mismatch in " + std::string( __PRETTY_FUNCTION__ ) );
}

// ===========================================================================
// TEST C – plastic response to a displacement jump
//
// dU_top={0,1e-4,0} -> gamma_23 = 1e-4/1e-3 = 0.1
// trial S23 = G*0.1 = 3846.15 >> tau_yield = 100/sqrt(3) = 57.74  ->  PLASTIC
// Newton return mapping (HLin=10):  S23 = 58.0633280
// ===========================================================================
void testPlasticDisplacementJump()
{
  // yieldStress=100, HLin=10, deltaYieldStress=0  (linear isotropic hardening)
  const double props[7]           = { 1e5, 0.3, 1e-3, 100., 10., 0., 1. };
  const double dU[6]              = { 0., 1e-4, 0., 0., 0., 0. };
  const double dSurfaceStrain[18] = { 0. };
  const double normal[3]          = { 0., 0., 1. };

  double force[3]         = { 0. };
  double averageStress[9] = { 0. };
  double pNewDT           = 1.0;

  runSingleIncrement( props, dU, dSurfaceStrain, normal, force, averageStress, pNewDT );

  // Analytically computed (see file header):
  const double S23 = 58.0633280;

  const double forceTarget[3] = { 0., S23, 0. };
  const double avgTarget[9]   = { 0., 0., 0., 0., 0., S23, 0., S23, 0. };

  Eigen::Map< const Eigen::Vector3d > forceVec( force );
  Eigen::Map< const Eigen::Vector3d > forceTgt( forceTarget );
  Eigen::Map< const Eigen::VectorXd > avgVec( averageStress, 9 );
  Eigen::Map< const Eigen::VectorXd > avgTgt( avgTarget, 9 );

  throwExceptionOnFailure( pNewDT == 1.0, "Newton diverged in " + std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( checkIfEqual< double >( forceVec, forceTgt, 1e-6 ),
                           "force mismatch in " + std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( checkIfEqual< double >( avgVec, avgTgt, 1e-6 ),
                           "averageStress mismatch in " + std::string( __PRETTY_FUNCTION__ ) );
}

// ===========================================================================
// TEST D – plastic response to a surface-strain increment
//
// dSurface_strain_top: [0,1]=[1,0]=1e-2 -> avg[0,1]=5e-3 -> gamma_12=1e-2
// trial S12 = G*1e-2 = 384.615 >> tau_yield = 57.74  ->  PLASTIC
// Newton return mapping (HLin=10):  S12 = 57.7633540
// ===========================================================================
void testPlasticSurfaceStrain()
{
  const double props[7]           = { 1e5, 0.3, 1e-3, 100., 10., 0., 1. };
  const double dU[6]              = { 0. };
  const double dSurfaceStrain[18] = { 0., 1e-2, 0., 1e-2, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
  const double normal[3]          = { 0., 0., 1. };

  double force[3]         = { 0. };
  double averageStress[9] = { 0. };
  double pNewDT           = 1.0;

  runSingleIncrement( props, dU, dSurfaceStrain, normal, force, averageStress, pNewDT );

  // Analytically computed (see file header):
  const double S12 = 57.7633540;

  const double forceTarget[3] = { 0., 0., 0. };
  const double avgTarget[9]   = { 0., S12, 0., S12, 0., 0., 0., 0., 0. };

  Eigen::Map< const Eigen::Vector3d > forceVec( force );
  Eigen::Map< const Eigen::Vector3d > forceTgt( forceTarget );
  Eigen::Map< const Eigen::VectorXd > avgVec( averageStress, 9 );
  Eigen::Map< const Eigen::VectorXd > avgTgt( avgTarget, 9 );

  throwExceptionOnFailure( pNewDT == 1.0, "Newton diverged in " + std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( checkIfEqual< double >( forceVec, forceTgt, 1e-6 ),
                           "force mismatch in " + std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( checkIfEqual< double >( avgVec, avgTgt, 1e-6 ),
                           "averageStress mismatch in " + std::string( __PRETTY_FUNCTION__ ) );
}

// ===========================================================================
int main()
{
  auto tests = std::vector< std::function< void() > >{
    testElasticDisplacementJump, // A: elastic, displacement jump
    testElasticSurfaceStrain,    // B: elastic, surface strain
    testPlasticDisplacementJump, // C: plastic, displacement jump
    testPlasticSurfaceStrain,    // D: plastic, surface strain
  };

  executeTestsAndCollectExceptions( tests );
  return 0;
}
