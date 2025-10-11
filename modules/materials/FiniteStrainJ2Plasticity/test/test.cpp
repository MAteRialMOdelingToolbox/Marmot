#include "Marmot/FiniteStrainJ2Plasticity.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotTesting.h"
#include <string>

using namespace Marmot::Testing;
using namespace Marmot::Materials;
using namespace Marmot::FastorStandardTensors;
using namespace Marmot::FastorIndices;

// NOTE: Scalar return mapping not implemented yet

std::string getAlgorithmName( int k )
{
  switch ( k ) {
  case 1: return "FullReturnMapping";
  case 2: return "FDAF";
  case 3: return "FDAC";
  case 4: return "CSDA";
  default: return "Unknown";
  }
}

void testSetup( const std::string& testName,
                const Tensor33d&   inputF,
                const Tensor33d&   targetStress,
                const Tensor33d&   targetFp,
                const double       targetAlphaP,
                bool               checkTangent      = false,
                const Tensor3333d& targetTangent     = Tensor3333d( 0.0 ),
                const Tensor3333d& targetTangentFDAF = Tensor3333d( 0.0 ),
                const Tensor3333d& targetTangentFDAC = Tensor3333d( 0.0 ),
                bool               ObjectivityCheck  = false,
                bool               IsotropyCheck     = false )
{

  for ( int k = 1; k < 5; k++ ) {
    std::string algorithmName = getAlgorithmName( k );

    // Material properties: K, G, fy, fyInf, eta, H, implementation type, (density)
    // Implementation type: 0 - scalar return mapping (not implemented yet), 1 - full return mapping, 2 - FDAF, 3 -
    // FDAC, 4 - CSDA; Not relevant for the current tests, as the computeStress routine directly called.
    std::array< double, 7 > materialProperties_ = { 175000, 80800, 260, 580, 9, 70, static_cast< double >( k ) };
    const double            nMaterialProperties = 7;
    const int               elLabel             = 1;

    // Create material instance
    FiniteStrainJ2Plasticity mat = FiniteStrainJ2Plasticity( &materialProperties_[0], nMaterialProperties, elLabel );

    if ( mat.getNumberOfRequiredStateVars() > 10 ) {
      throw std::runtime_error( "Number of required state vars changed!" );
    }

    // Set initial state variables
    double                   alphaP = 0.0;                                     // initial equivalent plastic strain
    Tensor33d                Fp = Marmot::FastorStandardTensors::Spatial3D::I; // initial plastic deformation gradient
    std::array< double, 10 > stateVars_;
    std::memcpy( stateVars_.data(), Fp.data(), 9 * sizeof( double ) );
    stateVars_[9] = alphaP;

    // Snapshot initial state for objectivity test
    std::array< double, 10 > stateVarsInitial = stateVars_;

    mat.assignStateVars( stateVars_.data(), 10 );

    // Create deformation, time increment, response and tangent objects required for stress computation
    FiniteStrainJ2Plasticity::Deformation< 3 > def;
    FiniteStrainJ2Plasticity::TimeIncrement    timeInc = { 0, 0.1 };

    FiniteStrainJ2Plasticity::ConstitutiveResponse< 3 > response;
    FiniteStrainJ2Plasticity::AlgorithmicModuli< 3 >    tangent;

    // Prescribe a deformation gradient tensor F for the considered load case
    def.F = inputF;

    // Compute stress response
    mat.computeStress( response, tangent, def, timeInc );

    // Check objectivity and isotropy if requested
    if ( ObjectivityCheck == false && IsotropyCheck == false ) {

      // Compare computed stress to target stress values
      throwExceptionOnFailure( checkIfEqual( response.tau, targetStress, 1e-10 ),
                               testName + " - Kirchhoff stress tensor (tau) computation failed for " + algorithmName +
                                 " for FiniteStrainJ2Plasticity material in " + std::string( __PRETTY_FUNCTION__ ) );

      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )

          throwExceptionOnFailure( checkIfEqual( response.tau( i, j ), response.tau( j, i ), 1e-10 ),
                                   testName + " - Kirchhoff stress tensor symmetry check failed for " + algorithmName +
                                     " for FiniteStrainJ2Plasticity material in " +
                                     std::string( __PRETTY_FUNCTION__ ) );

      // Extract updated state variables
      Tensor33d FpCurrent( 0.0 );
      std::memcpy( FpCurrent.data(), stateVars_.data(), 9 * sizeof( double ) );

      // Compare updated state variables to target values
      throwExceptionOnFailure( checkIfEqual( FpCurrent, targetFp, 1e-10 ),
                               testName +
                                 " - Plastic deformation gradient tensor (Fp) computation "
                                 "failed " +
                                 algorithmName + " for FiniteStrainJ2Plasticity material in " +
                                 std::string( __PRETTY_FUNCTION__ ) );

      throwExceptionOnFailure( checkIfEqual( stateVars_[9], targetAlphaP, 1e-10 ),
                               testName + " - Strain-like hardening variable (alphaP) computation failed for " +
                                 algorithmName + " for FiniteStrainJ2Plasticity material in " +
                                 std::string( __PRETTY_FUNCTION__ ) );

      if ( checkTangent && algorithmName == "FDAF" ) {
        // Compare algorithmic tangent for FDAF case to target tangent values
        throwExceptionOnFailure( checkIfEqual( tangent.dTau_dF, targetTangentFDAF, 1e-8 ),
                                 testName + " - Algorithmic tangent tensor computation failed for " + algorithmName +
                                   " for FiniteStrainJ2Plasticity material in " + std::string( __PRETTY_FUNCTION__ ) );
      }
      else if ( checkTangent && algorithmName == "FDAC" ) {
        // Compare algorithmic tangent for FDAC case to target tangent values
        throwExceptionOnFailure( checkIfEqual( tangent.dTau_dF, targetTangentFDAC, 1e-8 ),
                                 testName + " - Algorithmic tangent tensor computation failed for " + algorithmName +
                                   " for FiniteStrainJ2Plasticity material in " + std::string( __PRETTY_FUNCTION__ ) );
      }
      else if ( checkTangent ) {
        // Compare algorithmic tangent for FullReturnMapping and CSDA to target tangent values
        throwExceptionOnFailure( checkIfEqual( tangent.dTau_dF, targetTangent, 1e-10 ),
                                 testName + " - Algorithmic tangent tensor computation failed for " + algorithmName +
                                   " for FiniteStrainJ2Plasticity material in " + std::string( __PRETTY_FUNCTION__ ) );
      }
    }
    // Check objectivity and isotropy if requested
    if ( ObjectivityCheck ) {
      // Use already computed stress response and current F
      Tensor33d stressUnrotated = response.tau;
      Tensor33d F_unrotated     = def.F;

      for ( int phi_deg = 0; phi_deg <= 180; phi_deg += 30 ) {
        // Restore same initial state for each rotation
        stateVars_ = stateVarsInitial;
        mat.assignStateVars( stateVars_.data(), 10 );

        double phi = Marmot::Math::degToRad( phi_deg );

        Tensor33d Q( 0.0 );
        Q( 0, 0 ) = cos( phi );
        Q( 0, 1 ) = -sin( phi );
        Q( 1, 0 ) = sin( phi );
        Q( 1, 1 ) = cos( phi );
        Q( 2, 2 ) = 1;

        // Fr = Q * F -> Fr_ij = Q_ik F_kj
        Tensor33d F_rotated = einsum< ik, kj, to_ij >( Q, F_unrotated );
        def.F               = F_rotated;

        mat.computeStress( response, tangent, def, timeInc );

        Tensor33d stressNew = response.tau;

        // Tau (Q*F) = Q * Tau(F) * Q^T -> Tau(Q*F)_ij = Q_iI Tau(F)_IJ Q_jJ
        Tensor33d stressRotated = einsum< iI, IJ, jJ, to_ij >( Q, stressUnrotated, Q );

        throwExceptionOnFailure( checkIfEqual( stressNew, stressRotated, 1e-10 ),
                                 testName + " - Objectivity test failed for " + algorithmName +
                                   " (phi_deg=" + std::to_string( phi_deg ) +
                                   ") for FiniteStrainJ2Plasticity material in " + std::string( __PRETTY_FUNCTION__ ) );
      }
    }

    if ( IsotropyCheck ) {
      // Use already computed stress response and current deformed state
      Tensor33d stressUnrotated = response.tau;
      Tensor33d F_unrotated     = def.F;
      Tensor33d Fp_unrotated( 0.0 );
      std::memcpy( Fp_unrotated.data(), stateVars_.data(), 9 * sizeof( double ) );
      double alphaP_unrotated = stateVars_[9];

      for ( int phi_deg = 0; phi_deg <= 180; phi_deg += 30 ) {
        double phi = Marmot::Math::degToRad( phi_deg );

        Tensor33d Q( 0.0 );
        Q( 0, 0 ) = cos( phi );
        Q( 0, 1 ) = -sin( phi );
        Q( 1, 0 ) = sin( phi );
        Q( 1, 1 ) = cos( phi );
        Q( 2, 2 ) = 1;

        // Fr = F * Q -> Fr_ij = F_ik Q_kj
        Tensor33d F_rotated  = einsum< ik, kj, to_ij >( F_unrotated, Q );
        Tensor33d Fp_rotated = einsum< ik, kj, to_ij >( Fp_unrotated, Q );

        def.F = F_rotated;
        std::memcpy( stateVars_.data(), Fp_rotated.data(), 9 * sizeof( double ) );
        stateVars_[9] = alphaP_unrotated;
        mat.assignStateVars( stateVars_.data(), 10 );

        mat.computeStress( response, tangent, def, timeInc );

        Tensor33d stressNew = response.tau;

        throwExceptionOnFailure( checkIfEqual( stressNew, stressUnrotated, 1e-10 ),
                                 testName + " - Isotropy test failed for " + algorithmName +
                                   " (phi_deg=" + std::to_string( phi_deg ) +
                                   ") for FiniteStrainJ2Plasticity material in " + std::string( __PRETTY_FUNCTION__ ) );
      }
    }
  }
}

// Test I-1: Undeformed configuration
void testUndeformedResponse()
{
  Tensor33d inputF       = Marmot::FastorStandardTensors::Spatial3D::I;
  Tensor33d targetStress = Tensor33d( 0.0 );
  Tensor33d targetFp     = Marmot::FastorStandardTensors::Spatial3D::I;
  double    targetAlphaP = 0.0;
  testSetup( "I-1: Undeformed configuration", inputF, targetStress, targetFp, targetAlphaP );
}

// Test I-2: Deformation response
void testDeformationResponse()
{
  // Test I-2a: Simple shear deformation
  {
    Tensor33d inputF = Marmot::FastorStandardTensors::Spatial3D::I;
    inputF( 1, 0 ) += 0.02;

    Tensor33d stressTarget( 0.0 );
    stressTarget( 0, 0 ) = -1.66365138295391;
    stressTarget( 0, 1 ) = 166.959293480769;
    stressTarget( 1, 0 ) = 166.959293480769;
    stressTarget( 1, 1 ) = 1.67553448665712;
    stressTarget( 2, 2 ) = -0.0118831037141365;

    Tensor33d FpTarget( 0.0 );
    FpTarget( 0, 0 ) = 1.00013018240255;
    FpTarget( 0, 1 ) = 0.00896629252627631;
    FpTarget( 1, 0 ) = 0.00896629252627631;
    FpTarget( 1, 1 ) = 0.999950856552022;
    FpTarget( 2, 2 ) = 0.999999361845117;

    double alphaPTarget = 0.0103537584382;

    testSetup( "I-2a: Simple shear", inputF, stressTarget, FpTarget, alphaPTarget );
  }

  // Test I-2b: Hydrostatic deformation
  {
    Tensor33d inputF = Marmot::FastorStandardTensors::Spatial3D::I;
    inputF( 0, 0 ) += 0.002;
    inputF( 1, 1 ) += 0.002;
    inputF( 2, 2 ) += 0.002;

    Tensor33d stressTarget( 0.0 );
    stressTarget( 0, 0 ) = 1048.97652265991;
    stressTarget( 1, 1 ) = 1048.97652265991;
    stressTarget( 2, 2 ) = 1048.97652265991;

    Tensor33d FpTarget     = Marmot::FastorStandardTensors::Spatial3D::I;
    double    alphaPTarget = 0.0;

    testSetup( "I-2b: Hydrostatic", inputF, stressTarget, FpTarget, alphaPTarget );
  }

  // Test I-2c: Arbitrary deformation
  {
    Tensor33d inputF = Marmot::FastorStandardTensors::Spatial3D::I;
    inputF( 0, 0 )   = 1.01;
    inputF( 0, 1 )   = 0.06;
    inputF( 0, 2 )   = -0.03;
    inputF( 1, 0 )   = 0.06;
    inputF( 1, 1 )   = 1.02;
    inputF( 1, 2 )   = 0.04;
    inputF( 2, 0 )   = -0.03;
    inputF( 2, 1 )   = 0.04;
    inputF( 2, 2 )   = 0.95;

    Tensor33d stressTarget( 0.0 );
    stressTarget( 0, 0 ) = -4815.15772271678;
    stressTarget( 0, 1 ) = 179.062345389071;
    stressTarget( 0, 2 ) = -95.4629813550829;
    stressTarget( 1, 0 ) = 179.062345389071;
    stressTarget( 1, 1 ) = -4786.68287740693;
    stressTarget( 1, 2 ) = 124.546217296781;
    stressTarget( 2, 0 ) = -95.4629813550829;
    stressTarget( 2, 1 ) = 124.546217296781;
    stressTarget( 2, 2 ) = -4995.74104468547;

    Tensor33d FpTarget( 0.0 );
    FpTarget( 0, 0 ) = 1.01899554888122;
    FpTarget( 0, 1 ) = 0.0594232892187332;
    FpTarget( 0, 2 ) = -0.0297463602782848;
    FpTarget( 1, 0 ) = 0.0594232892187332;
    FpTarget( 1, 1 ) = 1.02889141910947;
    FpTarget( 1, 2 ) = 0.0396457910879559;
    FpTarget( 2, 0 ) = -0.0297463602782848;
    FpTarget( 2, 1 ) = 0.0396457910879559;
    FpTarget( 2, 2 ) = 0.959563358208919;

    double alphaPTarget = 0.0998876084740522;

    testSetup( "I-2c: Arbitrary deformation", inputF, stressTarget, FpTarget, alphaPTarget );
  }
}
// Test I-3: Algorithmic tangent
void testAlgorithmicTangent()
{
  Tensor33d inputF = Marmot::FastorStandardTensors::Spatial3D::I;
  inputF( 0, 0 ) += 0.001;
  inputF( 1, 1 ) += 0.002;
  inputF( 2, 2 ) += 0.003;

  Tensor33d stressTarget( 0.0 );
  stressTarget( 0, 0 ) = 898.575312428168;
  stressTarget( 1, 1 ) = 1048.76543449885;
  stressTarget( 2, 2 ) = 1199.06587694072;

  Tensor33d FpTarget( 0.0 );
  FpTarget( 0, 0 ) = 0.99993174180691;
  FpTarget( 1, 1 ) = 0.999999983290772;
  FpTarget( 2, 2 ) = 1.00006827956296;

  double alphaPTarget = 7.88301112597152e-05;

  Tensor3333d tangentTarget( 0.0 );
  tangentTarget( 0, 0, 0, 0 ) = 200890.525305335;
  tangentTarget( 0, 0, 1, 1 ) = 124564.513352951;
  tangentTarget( 0, 0, 2, 2 ) = 198537.043113999;
  tangentTarget( 0, 1, 0, 1 ) = 75132.5523288948;
  tangentTarget( 0, 1, 1, 0 ) = 75057.5697417401;
  tangentTarget( 0, 2, 0, 2 ) = 75197.6138238636;
  tangentTarget( 0, 2, 2, 0 ) = 75047.6684323903;
  tangentTarget( 1, 0, 0, 1 ) = 75132.5523288948;
  tangentTarget( 1, 0, 1, 0 ) = 75057.5697417401;
  tangentTarget( 1, 1, 0, 0 ) = 124688.953426225;
  tangentTarget( 1, 1, 1, 1 ) = 274826.645041077;
  tangentTarget( 1, 1, 2, 2 ) = 124474.348696745;
  tangentTarget( 1, 2, 1, 2 ) = 75187.7026280292;
  tangentTarget( 1, 2, 2, 1 ) = 75112.7398138437;
  tangentTarget( 2, 0, 0, 2 ) = 75197.6138238636;
  tangentTarget( 2, 0, 2, 0 ) = 75047.6684323903;
  tangentTarget( 2, 1, 1, 2 ) = 75187.7026280292;
  tangentTarget( 2, 1, 2, 1 ) = 75112.7398138437;
  tangentTarget( 2, 2, 0, 0 ) = 198933.720522818;
  tangentTarget( 2, 2, 1, 1 ) = 124598.574593654;
  tangentTarget( 2, 2, 2, 2 ) = 200455.918711321;

  Tensor3333d tangentTargetFDAF( 0.0 );
  tangentTargetFDAF( 0, 0, 0, 0 ) = 200891.207248906;
  tangentTargetFDAF( 0, 0, 0, 1 ) = 0.293344510107489;
  tangentTargetFDAF( 0, 0, 0, 2 ) = 0.29385296412801;
  tangentTargetFDAF( 0, 0, 1, 0 ) = 0.292965005011075;
  tangentTargetFDAF( 0, 0, 1, 1 ) = 124564.921473659;
  tangentTargetFDAF( 0, 0, 1, 2 ) = 0.294008477858679;
  tangentTargetFDAF( 0, 0, 2, 0 ) = 0.293582659210213;
  tangentTargetFDAF( 0, 0, 2, 1 ) = 0.294117661344026;
  tangentTargetFDAF( 0, 0, 2, 2 ) = 198536.580316643;
  tangentTargetFDAF( 0, 1, 0, 1 ) = 75132.5523290578;
  tangentTargetFDAF( 0, 1, 1, 0 ) = 75057.5697419029;
  tangentTargetFDAF( 0, 2, 0, 2 ) = 75197.6138240272;
  tangentTargetFDAF( 0, 2, 2, 0 ) = 75047.6684325537;
  tangentTargetFDAF( 1, 0, 0, 1 ) = 75132.5523290578;
  tangentTargetFDAF( 1, 0, 1, 0 ) = 75057.5697419029;
  tangentTargetFDAF( 1, 1, 0, 0 ) = 124689.184484852;
  tangentTargetFDAF( 1, 1, 0, 1 ) = -0.001343538153446;
  tangentTargetFDAF( 1, 1, 0, 2 ) = -0.00134265592944414;
  tangentTargetFDAF( 1, 1, 1, 0 ) = -0.00134406872000709;
  tangentTargetFDAF( 1, 1, 1, 1 ) = 274826.634674348;
  tangentTargetFDAF( 1, 1, 1, 2 ) = -0.00134226715396347;
  tangentTargetFDAF( 1, 1, 2, 0 ) = -0.00134303384086018;
  tangentTargetFDAF( 1, 1, 2, 1 ) = -0.00134211451069891;
  tangentTargetFDAF( 1, 1, 2, 2 ) = 124474.117290285;
  tangentTargetFDAF( 1, 2, 1, 2 ) = 75187.7026281932;
  tangentTargetFDAF( 1, 2, 2, 1 ) = 75112.7398140078;
  tangentTargetFDAF( 2, 0, 0, 2 ) = 75197.6138240272;
  tangentTargetFDAF( 2, 0, 2, 0 ) = 75047.6684325537;
  tangentTargetFDAF( 2, 1, 1, 2 ) = 75187.7026281932;
  tangentTargetFDAF( 2, 1, 2, 1 ) = 75112.7398140078;
  tangentTargetFDAF( 2, 2, 0, 0 ) = 198934.179430311;
  tangentTargetFDAF( 2, 2, 0, 1 ) = -0.299419075895853;
  tangentTargetFDAF( 2, 2, 0, 2 ) = -0.29992769650696;
  tangentTargetFDAF( 2, 2, 1, 0 ) = -0.299039033663026;
  tangentTargetFDAF( 2, 2, 1, 1 ) = 124598.157567648;
  tangentTargetFDAF( 2, 2, 1, 2 ) = -0.30008287730396;
  tangentTargetFDAF( 2, 2, 2, 0 ) = -0.299657009019649;
  tangentTargetFDAF( 2, 2, 2, 1 ) = -0.300192215324636;
  tangentTargetFDAF( 2, 2, 2, 2 ) = 200455.224503959;

  Tensor3333d tangentTargetFDAC( 0.0 );
  tangentTargetFDAC( 0, 0, 0, 0 ) = 200890.509333366;
  tangentTargetFDAC( 0, 0, 1, 1 ) = 124564.488040833;
  tangentTargetFDAC( 0, 0, 2, 2 ) = 198537.025177391;
  tangentTargetFDAC( 0, 1, 0, 1 ) = 75132.5793713203;
  tangentTargetFDAC( 0, 1, 1, 0 ) = 75057.5967841801;
  tangentTargetFDAC( 1, 0, 0, 1 ) = 75132.5793713203;
  tangentTargetFDAC( 1, 0, 1, 0 ) = 75057.5967841801;
  tangentTargetFDAC( 0, 2, 0, 2 ) = 75197.6409636925;
  tangentTargetFDAC( 0, 2, 2, 0 ) = 75047.69557222;
  tangentTargetFDAC( 2, 0, 0, 2 ) = 75197.6409636925;
  tangentTargetFDAC( 2, 0, 2, 0 ) = 75047.69557222;
  tangentTargetFDAC( 1, 1, 0, 0 ) = 124688.985901659;
  tangentTargetFDAC( 1, 1, 1, 1 ) = 274826.696170781;
  tangentTargetFDAC( 1, 1, 2, 2 ) = 124474.382115553;
  tangentTargetFDAC( 1, 2, 1, 2 ) = 75187.7298654265;
  tangentTargetFDAC( 1, 2, 2, 1 ) = 75112.7670512278;
  tangentTargetFDAC( 2, 1, 1, 2 ) = 75187.7298654265;
  tangentTargetFDAC( 2, 1, 2, 1 ) = 75112.7670512278;
  tangentTargetFDAC( 2, 2, 0, 0 ) = 198933.702151485;
  tangentTargetFDAC( 2, 2, 1, 1 ) = 124598.548812737;
  tangentTargetFDAC( 2, 2, 2, 2 ) = 200455.905150422;

  bool checkTangent = true;

  testSetup( "I-3: Algorithmic tangent",
             inputF,
             stressTarget,
             FpTarget,
             alphaPTarget,
             checkTangent,
             tangentTarget,
             tangentTargetFDAF,
             tangentTargetFDAC );
}

// Test I-4: Rotation
void testRotation()
{
  // Test I-4a: Pure rotation about the z-axis
  {
    for ( int phi_deg = 0; phi_deg <= 180; phi_deg++ ) {
      double phi = Marmot::Math::degToRad( phi_deg );

      Tensor33d inputF( 0.0 );
      inputF( 0, 0 ) = cos( phi );
      inputF( 0, 1 ) = -sin( phi );
      inputF( 0, 2 ) = 0;
      inputF( 1, 0 ) = sin( phi );
      inputF( 1, 1 ) = cos( phi );
      inputF( 1, 2 ) = 0;
      inputF( 2, 0 ) = 0;
      inputF( 2, 1 ) = 0;
      inputF( 2, 2 ) = 1;

      Tensor33d stressTarget( 0.0 );
      Tensor33d FpTarget     = Marmot::FastorStandardTensors::Spatial3D::I;
      double    alphaPTarget = 0.0;

      testSetup( "I-4a: Pure rotation about z-axis", inputF, stressTarget, FpTarget, alphaPTarget );
    }
  }

  // Test I-4b: Objectivity test for arbitrary deformation
  // Test I-4c: Isotropy test for arbitrary deformation
  {
    Tensor33d inputF( 0.0 );
    inputF( 0, 0 ) = 1.01;
    inputF( 0, 1 ) = 0.06;
    inputF( 0, 2 ) = -0.03;
    inputF( 1, 0 ) = 0.06;
    inputF( 1, 1 ) = 1.02;
    inputF( 1, 2 ) = 0.04;
    inputF( 2, 0 ) = -0.03;
    inputF( 2, 1 ) = 0.04;
    inputF( 2, 2 ) = 0.95;

    testSetup( "I-4b: ",
               inputF,
               Tensor33d( 0.0 ),
               Tensor33d( 0.0 ),
               0.0,
               false,
               Tensor3333d( 0.0 ),
               Tensor3333d( 0.0 ),
               Tensor3333d( 0.0 ),
               true,
               false );
    testSetup( "I-4c: ",
               inputF,
               Tensor33d( 0.0 ),
               Tensor33d( 0.0 ),
               0.0,
               false,
               Tensor3333d( 0.0 ),
               Tensor3333d( 0.0 ),
               Tensor3333d( 0.0 ),
               false,
               true );
  }
}

int main()
{
  auto tests = std::vector< std::function< void() > >{ testUndeformedResponse,
                                                       testDeformationResponse,
                                                       testAlgorithmicTangent,
                                                       testRotation };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
