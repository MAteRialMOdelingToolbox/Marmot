#include "Marmot/FiniteStrainJ2Plasticity.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotTesting.h"

using namespace Marmot::Testing;
using namespace Marmot::Materials;
using namespace Marmot::FastorStandardTensors;
using namespace Marmot::FastorIndices;

// Scalar return mapping not implemented yet

////////////////////////////////////////
// I test group - full return mapping //
////////////////////////////////////////

// Test I-1: Undeformed configuration
void testUndeformedResponseFullReturnMapping()
{
  // Material properties: K, G, fy, fyInf, eta, H, implementation type, (density)
  // Implementation type: 0 - scalar return mapping (not implemented yet), 1 - full return mapping, 2 - FDAF, 3 - FDAC,
  // 4 - CSDA; Not relevant for the current tests, as the computeStress routine directly called. 
  std::array< double, 7 > materialProperties_ = { 175000, 80800, 260, 580, 9, 70, 1 };
  const double            nMaterialProperties = 7;
  const int               elLabel             = 1;

  // Create material instance
  FiniteStrainJ2Plasticity mat = FiniteStrainJ2Plasticity( &materialProperties_[0], nMaterialProperties, elLabel );

  if ( mat.getNumberOfRequiredStateVars() > 10 ) {
    throw std::runtime_error( "Number of required state vars changed!" );
  }

  // Set initial state variables
  double                   alphaP = 0.0;                                         // initial equivalent plastic strain
  Tensor33d                Fp     = Marmot::FastorStandardTensors::Spatial3D::I; // initial plastic deformation gradient
  std::array< double, 10 > stateVars_;
  std::memcpy( stateVars_.data(), Fp.data(), 9 * sizeof( double ) );
  stateVars_[9] = alphaP;

  mat.assignStateVars( stateVars_.data(), 10 );

  // Create deformation, time increment, response and tangent objects required for stress computation
  FiniteStrainJ2Plasticity::Deformation< 3 > def;
  FiniteStrainJ2Plasticity::TimeIncrement    timeInc = { 0, 0.1 };

  FiniteStrainJ2Plasticity::ConstitutiveResponse< 3 > response;
  FiniteStrainJ2Plasticity::AlgorithmicModuli< 3 >    tangent;

  // Prescribe a deformation gradient tensor F for the considered load case
  def.F = Marmot::FastorStandardTensors::Spatial3D::I;

  // Compute stress response
  mat.computeStressWithFullReturnMapping( response, tangent, def, timeInc );

  // Target Kirchhoff stress values for the given deformation gradient and material parameters
  Tensor33d stressTarget( 0.0 );

  // Compare computed stress to target stress values
  throwExceptionOnFailure( checkIfEqual( response.tau, stressTarget, 1e-10 ),
                           "I-1: Undeformed configuration - Kirchhoff stress tensor (tau) computation failed for "
                           "FiniteStrainJ2Plasticity material in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // Extract updated state variables
  Tensor33d FpCurrent( 0.0 );
  std::memcpy( FpCurrent.data(), stateVars_.data(), 9 * sizeof( double ) );

  Tensor33d FpTarget = Marmot::FastorStandardTensors::Spatial3D::I;

  // Compare updated state variables to target values
  throwExceptionOnFailure( checkIfEqual( FpCurrent, FpTarget, 1e-10 ),
                           "I-1: Undeformed configuration - Plastic deformation gradient tensor (Fp) computation "
                           "failed for "
                           "FiniteStrainJ2Plasticity material in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  double alphaPTarget = 0.0;

  throwExceptionOnFailure( checkIfEqual( stateVars_[9], alphaPTarget, 1e-10 ),
                           "I-1: Undeformed configuration - Strain-like hardening variable (alphaP) computation failed "
                           "for "
                           "FiniteStrainJ2Plasticity material in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Test I-2: Deformation response

void testDeformationResponseFullReturnMapping()
{

  // Test I-2a: Simple shear deformation
  {
    std::array< double, 7 > materialProperties_ = { 175000, 80800, 260, 580, 9, 70, 1 };
    const double            nMaterialProperties = 7;
    const int               elLabel             = 1;

    FiniteStrainJ2Plasticity mat = FiniteStrainJ2Plasticity( &materialProperties_[0], nMaterialProperties, elLabel );

    if ( mat.getNumberOfRequiredStateVars() > 10 ) {
      throw std::runtime_error( "Number of required state vars changed!" );
    }

    double                   alphaP = 0.0;
    Tensor33d                Fp     = Marmot::FastorStandardTensors::Spatial3D::I;
    std::array< double, 10 > stateVars_;
    std::memcpy( stateVars_.data(), Fp.data(), 9 * sizeof( double ) );
    stateVars_[9] = alphaP;

    mat.assignStateVars( stateVars_.data(), 10 );

    FiniteStrainJ2Plasticity::Deformation< 3 > def;
    FiniteStrainJ2Plasticity::TimeIncrement    timeInc = { 0, 0.1 };

    FiniteStrainJ2Plasticity::ConstitutiveResponse< 3 > response;
    FiniteStrainJ2Plasticity::AlgorithmicModuli< 3 >    tangent;

    def.F = Marmot::FastorStandardTensors::Spatial3D::I;
    def.F( 1, 0 ) += 0.02;

    mat.computeStressWithFullReturnMapping( response, tangent, def, timeInc );

    Tensor33d stressTarget( 0.0 );
    stressTarget( 0, 0 ) = -1.66365138295391;
    stressTarget( 0, 1 ) = 166.959293480769;
    stressTarget( 1, 0 ) = 166.959293480769;
    stressTarget( 1, 1 ) = 1.67553448665712;
    stressTarget( 2, 2 ) = -0.0118831037141365;

    throwExceptionOnFailure( checkIfEqual( response.tau, stressTarget, 1e-10 ),
                             "I-2a: Simple shear - Kirchhoff stress tensor (tau) computation failed for "
                             "FiniteStrainJ2Plasticity material in " +
                               std::string( __PRETTY_FUNCTION__ ) );

    Tensor33d FpCurrent( 0.0 );
    std::memcpy( FpCurrent.data(), stateVars_.data(), 9 * sizeof( double ) );

    Tensor33d FpTarget( 0.0 );
    FpTarget( 0, 0 ) = 1.00013018240255;
    FpTarget( 0, 1 ) = 0.00896629252627631;
    FpTarget( 1, 0 ) = 0.00896629252627631;
    FpTarget( 1, 1 ) = 0.999950856552022;
    FpTarget( 2, 2 ) = 0.999999361845117;

    throwExceptionOnFailure( checkIfEqual( FpCurrent, FpTarget, 1e-10 ),
                             "I-2a: Simple shear - Plastic deformation gradient tensor (Fp) computation failed for "
                             "FiniteStrainJ2Plasticity material in " +
                               std::string( __PRETTY_FUNCTION__ ) );

    double alphaPTarget = 0.0103537584382;

    throwExceptionOnFailure( checkIfEqual( stateVars_[9], alphaPTarget, 1e-10 ),
                             "I-2a: Simple shear - Strain-like hardening variable (alphaP) computation failed for "
                             "FiniteStrainJ2Plasticity material in " +
                               std::string( __PRETTY_FUNCTION__ ) );
  }

  // Test I-2b: Hydrostatic deformation
  {
    std::array< double, 7 > materialProperties_ = { 175000, 80800, 260, 580, 9, 70, 1 };
    const double            nMaterialProperties = 7;
    const int               elLabel             = 1;

    FiniteStrainJ2Plasticity mat = FiniteStrainJ2Plasticity( &materialProperties_[0], nMaterialProperties, elLabel );

    if ( mat.getNumberOfRequiredStateVars() > 10 ) {
      throw std::runtime_error( "Number of required state vars changed!" );
    }

    double                   alphaP = 0.0;
    Tensor33d                Fp     = Marmot::FastorStandardTensors::Spatial3D::I;
    std::array< double, 10 > stateVars_;
    std::memcpy( stateVars_.data(), Fp.data(), 9 * sizeof( double ) );
    stateVars_[9] = alphaP;

    mat.assignStateVars( stateVars_.data(), 10 );

    FiniteStrainJ2Plasticity::Deformation< 3 > def;
    FiniteStrainJ2Plasticity::TimeIncrement    timeInc = { 0, 0.1 };

    FiniteStrainJ2Plasticity::ConstitutiveResponse< 3 > response;
    FiniteStrainJ2Plasticity::AlgorithmicModuli< 3 >    tangent;

    def.F = Marmot::FastorStandardTensors::Spatial3D::I;
    def.F( 0, 0 ) += 0.002;
    def.F( 1, 1 ) += 0.002;
    def.F( 2, 2 ) += 0.002;

    mat.computeStressWithFullReturnMapping( response, tangent, def, timeInc );

    Tensor33d stressTarget( 0.0 );
    stressTarget( 0, 0 ) = 1048.97652265991;
    stressTarget( 1, 1 ) = 1048.97652265991;
    stressTarget( 2, 2 ) = 1048.97652265991;

    throwExceptionOnFailure( checkIfEqual( response.tau, stressTarget, 1e-10 ),
                             "I-2b: Hydrostatic - Kirchhoff stress tensor (tau) computation failed for "
                             "FiniteStrainJ2Plasticity material in " +
                               std::string( __PRETTY_FUNCTION__ ) );

    Tensor33d FpCurrent( 0.0 );
    std::memcpy( FpCurrent.data(), stateVars_.data(), 9 * sizeof( double ) );

    Tensor33d FpTarget = Marmot::FastorStandardTensors::Spatial3D::I;

    throwExceptionOnFailure( checkIfEqual( FpCurrent, FpTarget, 1e-10 ),
                             "I-2b: Hydrostatic - Plastic deformation gradient tensor (Fp) computation failed for "
                             "FiniteStrainJ2Plasticity material in " +
                               std::string( __PRETTY_FUNCTION__ ) );

    double alphaPTarget = 0.0;

    throwExceptionOnFailure( checkIfEqual( stateVars_[9], alphaPTarget, 1e-10 ),
                             "I-2b: Hydrostatic - Strain-like hardening variable (alphaP) computation failed for "
                             "FiniteStrainJ2Plasticity material in " +
                               std::string( __PRETTY_FUNCTION__ ) );
  }

  // Test I-2c: Arbitrary deformation
  {
    std::array< double, 7 > materialProperties_ = { 175000, 80800, 260, 580, 9, 70, 1 };
    const double            nMaterialProperties = 7;
    const int               elLabel             = 1;

    FiniteStrainJ2Plasticity mat = FiniteStrainJ2Plasticity( &materialProperties_[0], nMaterialProperties, elLabel );

    if ( mat.getNumberOfRequiredStateVars() > 10 ) {
      throw std::runtime_error( "Number of required state vars changed!" );
    }

    double                   alphaP = 0.0;
    Tensor33d                Fp     = Marmot::FastorStandardTensors::Spatial3D::I;
    std::array< double, 10 > stateVars_;
    std::memcpy( stateVars_.data(), Fp.data(), 9 * sizeof( double ) );
    stateVars_[9] = alphaP;

    mat.assignStateVars( stateVars_.data(), 10 );

    FiniteStrainJ2Plasticity::Deformation< 3 > def;
    FiniteStrainJ2Plasticity::TimeIncrement    timeInc = { 0, 0.1 };

    FiniteStrainJ2Plasticity::ConstitutiveResponse< 3 > response;
    FiniteStrainJ2Plasticity::AlgorithmicModuli< 3 >    tangent;

    def.F         = Marmot::FastorStandardTensors::Spatial3D::I;
    def.F( 0, 0 ) = 1.01;
    def.F( 0, 1 ) = 0.06;
    def.F( 0, 2 ) = -0.03;
    def.F( 1, 0 ) = 0.06;
    def.F( 1, 1 ) = 1.02;
    def.F( 1, 2 ) = 0.04;
    def.F( 2, 0 ) = -0.03;
    def.F( 2, 1 ) = 0.04;
    def.F( 2, 2 ) = 0.95;

    mat.computeStressWithFullReturnMapping( response, tangent, def, timeInc );

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

    throwExceptionOnFailure( checkIfEqual( response.tau, stressTarget, 1e-10 ),
                             "I-2c: Arbitrary deformation - Kirchhoff stress tensor (tau) computation failed for "
                             "FiniteStrainJ2Plasticity material in " +
                               std::string( __PRETTY_FUNCTION__ ) );

    for ( int i = 0; i < 3; i++ )
      for ( int j = 0; j < 3; j++ )

        throwExceptionOnFailure( checkIfEqual( response.tau( i, j ), response.tau( j, i ), 1e-10 ),
                                 "I-2c: Kirchhoff stress tensor symmetry check for the arbitrary deformation load "
                                 "case failed for FiniteStrainJ2Plasticity "
                                 "material in " +
                                   std::string( __PRETTY_FUNCTION__ ) );

    Tensor33d FpCurrent( 0.0 );
    std::memcpy( FpCurrent.data(), stateVars_.data(), 9 * sizeof( double ) );

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

    throwExceptionOnFailure( checkIfEqual( FpCurrent, FpTarget, 1e-10 ),
                             "I-2c: Arbitrary deformation - Plastic deformation gradient tensor (Fp) computation "
                             "failed for "
                             "FiniteStrainJ2Plasticity material in " +
                               std::string( __PRETTY_FUNCTION__ ) );

    double alphaPTarget = 0.0998876084740522;

    throwExceptionOnFailure( checkIfEqual( stateVars_[9], alphaPTarget, 1e-10 ),
                             "I-2c: Arbitrary deformation - Strain-like hardening variable (alphaP) computation failed "
                             "for "
                             "FiniteStrainJ2Plasticity material in " +
                               std::string( __PRETTY_FUNCTION__ ) );
  }
}

// Test I-3: Algorithmic tangent
void testAlgorithmicTangentFullReturnMapping()

{
  std::array< double, 7 > materialProperties_ = { 175000, 80800, 260, 580, 9, 70, 1 };
  const double            nMaterialProperties = 7;
  const int               elLabel             = 1;

  FiniteStrainJ2Plasticity mat = FiniteStrainJ2Plasticity( &materialProperties_[0], nMaterialProperties, elLabel );

  if ( mat.getNumberOfRequiredStateVars() > 10 ) {
    throw std::runtime_error( "Number of required state vars changed!" );
  }

  double                   alphaP = 0.0;
  Tensor33d                Fp     = Marmot::FastorStandardTensors::Spatial3D::I;
  std::array< double, 10 > stateVars_;
  std::memcpy( stateVars_.data(), Fp.data(), 9 * sizeof( double ) );
  stateVars_[9] = alphaP;

  mat.assignStateVars( stateVars_.data(), 10 );

  FiniteStrainJ2Plasticity::Deformation< 3 > def;
  FiniteStrainJ2Plasticity::TimeIncrement    timeInc = { 0, 0.1 };

  FiniteStrainJ2Plasticity::ConstitutiveResponse< 3 > response;
  FiniteStrainJ2Plasticity::AlgorithmicModuli< 3 >    tangent;

  def.F = Marmot::FastorStandardTensors::Spatial3D::I;
  def.F( 0, 0 ) += 0.001;
  def.F( 1, 1 ) += 0.002;
  def.F( 2, 2 ) += 0.003;

  mat.computeStressWithFullReturnMapping( response, tangent, def, timeInc );

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

  throwExceptionOnFailure( checkIfEqual( tangent.dTau_dF, tangentTarget, 1e-10 ),
                           "I-3: Algorithmic tangent computation failed for FiniteStrainJ2Plasticity material in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Test I-4: Rotation
void testRotationFullReturnMapping()
{
  // Test I-4a: Pure rotation about the z-axis
  {
    std::array< double, 7 > materialProperties_ = { 175000, 80800, 260, 580, 9, 70, 1 };
    const double            nMaterialProperties = 7;
    const int               elLabel             = 1;

    FiniteStrainJ2Plasticity mat = FiniteStrainJ2Plasticity( &materialProperties_[0], nMaterialProperties, elLabel );

    if ( mat.getNumberOfRequiredStateVars() > 10 ) {
      throw std::runtime_error( "Number of required state vars changed!" );
    }

    double                   alphaP = 0.0;
    Tensor33d                Fp     = Marmot::FastorStandardTensors::Spatial3D::I;
    std::array< double, 10 > stateVars_;
    std::memcpy( stateVars_.data(), Fp.data(), 9 * sizeof( double ) );
    stateVars_[9] = alphaP;

    mat.assignStateVars( stateVars_.data(), 10 );

    FiniteStrainJ2Plasticity::Deformation< 3 > def;
    FiniteStrainJ2Plasticity::TimeIncrement    timeInc = { 0, 0.1 };

    FiniteStrainJ2Plasticity::ConstitutiveResponse< 3 > response;
    FiniteStrainJ2Plasticity::AlgorithmicModuli< 3 >    tangent;

    for ( int phi_deg = 0; phi_deg <= 180; phi_deg++ ) {
      double phi = Marmot::Math::degToRad( phi_deg );

      def.F( 0, 0 ) = cos( phi );
      def.F( 0, 1 ) = -sin( phi );
      def.F( 0, 2 ) = 0;
      def.F( 1, 0 ) = sin( phi );
      def.F( 1, 1 ) = cos( phi );
      def.F( 1, 2 ) = 0;
      def.F( 2, 0 ) = 0;
      def.F( 2, 1 ) = 0;
      def.F( 2, 2 ) = 1;

      mat.computeStressWithFullReturnMapping( response, tangent, def, timeInc );

      Tensor33d stressTarget( 0.0 );

      throwExceptionOnFailure( checkIfEqual( response.tau, stressTarget, 1e-10 ),
                               "I-4a: Rotation about the z-axis - Kirchhoff stress tensor (tau) computation failed for "
                               "FiniteStrainJ2Plasticity material in " +
                                 std::string( __PRETTY_FUNCTION__ ) );

      Tensor33d FpCurrent( 0.0 );
      std::memcpy( FpCurrent.data(), stateVars_.data(), 9 * sizeof( double ) );

      Tensor33d FpTarget = Marmot::FastorStandardTensors::Spatial3D::I;

      throwExceptionOnFailure( checkIfEqual( FpCurrent, FpTarget, 1e-10 ),
                               "I-4a: Pure rotation about the z-axis - Plastic deformation gradient tensor (Fp) "
                               "computation "
                               "failed "
                               "for "
                               "FiniteStrainJ2Plasticity material in " +
                                 std::string( __PRETTY_FUNCTION__ ) );

      double alphaPTarget = 0.0;

      throwExceptionOnFailure( checkIfEqual( stateVars_[9], alphaPTarget, 1e-10 ),
                               "I-4a: Rotation about the z-axis - Strain-like hardening variable (alphaP) computation "
                               "failed "
                               "for "
                               "FiniteStrainJ2Plasticity material in " +
                                 std::string( __PRETTY_FUNCTION__ ) );
    }
  }

  // Test I-4b: Objectivity test for arbitrary deformation and rotation about the z-axis
  {
    std::array< double, 7 > materialProperties_ = { 175000, 80800, 260, 580, 9, 70, 1 };
    const double            nMaterialProperties = 7;
    const int               elLabel             = 1;

    FiniteStrainJ2Plasticity mat = FiniteStrainJ2Plasticity( &materialProperties_[0], nMaterialProperties, elLabel );

    if ( mat.getNumberOfRequiredStateVars() > 10 ) {
      throw std::runtime_error( "Number of required state vars changed!" );
    }

    double                   alphaP = 0.0;
    Tensor33d                Fp     = Marmot::FastorStandardTensors::Spatial3D::I;
    std::array< double, 10 > stateVars_;
    std::memcpy( stateVars_.data(), Fp.data(), 9 * sizeof( double ) );
    stateVars_[9] = alphaP;

    mat.assignStateVars( stateVars_.data(), 10 );

    FiniteStrainJ2Plasticity::Deformation< 3 > def;
    FiniteStrainJ2Plasticity::TimeIncrement    timeInc = { 0, 0.1 };

    FiniteStrainJ2Plasticity::ConstitutiveResponse< 3 > response;
    FiniteStrainJ2Plasticity::AlgorithmicModuli< 3 >    tangent;

    def.F( 0, 0 ) = 1.01;
    def.F( 0, 1 ) = 0.06;
    def.F( 0, 2 ) = -0.03;
    def.F( 1, 0 ) = 0.06;
    def.F( 1, 1 ) = 1.02;
    def.F( 1, 2 ) = 0.04;
    def.F( 2, 0 ) = -0.03;
    def.F( 2, 1 ) = 0.04;
    def.F( 2, 2 ) = 0.95;

    mat.computeStressWithFullReturnMapping( response, tangent, def, timeInc );

    Tensor33d stressUnrotated = response.tau;
    Tensor33d F_unrotated     = def.F;

    for ( int phi_deg = 0; phi_deg <= 180; phi_deg++ ) {
      double phi = Marmot::Math::degToRad( phi_deg );

      Tensor33d Q( 0.0 );
      Q( 0, 0 ) = cos( phi );
      Q( 0, 1 ) = -sin( phi );
      Q( 1, 0 ) = sin( phi );
      Q( 1, 1 ) = cos( phi );
      Q( 2, 2 ) = 1;

      // Fr = Q * F -> Fr_ij = Q_ik F_kj

      Tensor33d F_rotated = einsum< ik, kj, to_ij >( Q, F_unrotated );

      def.F = F_rotated;

      mat.computeStressWithFullReturnMapping( response, tangent, def, timeInc );

      Tensor33d stressNew = response.tau;

      // Tau (Q*F) = Q * Tau(F) * Q^T -> Tau(Q*F)_ij = Q_iI Tau(F)_IJ Q_jJ
      Tensor33d stressRotated = einsum< iI, IJ, jJ, to_ij >( Q, stressUnrotated, Q );

      throwExceptionOnFailure( checkIfEqual( stressNew, stressRotated, 1e-10 ),
                               "I-4b: Objectivity test for arbitrary deformation and rotation around z-axis failed "
                               "for "
                               "FiniteStrainJ2Plasticity "
                               "material in " +
                                 std::string( __PRETTY_FUNCTION__ ) );
    }
  }
  // Test I-4c: Isotropy test for arbitrary deformation
  {
    std::array< double, 7 > materialProperties_ = { 175000, 80800, 260, 580, 9, 70, 1 };
    const double            nMaterialProperties = 7;
    const int               elLabel             = 1;

    FiniteStrainJ2Plasticity mat = FiniteStrainJ2Plasticity( &materialProperties_[0], nMaterialProperties, elLabel );

    if ( mat.getNumberOfRequiredStateVars() > 10 ) {
      throw std::runtime_error( "Number of required state vars changed!" );
    }

    double    alphaP = 0.0;
    Tensor33d Fp     = Marmot::FastorStandardTensors::Spatial3D::I;

    std::array< double, 10 > stateVars_;
    std::memcpy( stateVars_.data(), Fp.data(), 9 * sizeof( double ) );
    stateVars_[9] = alphaP;

    mat.assignStateVars( stateVars_.data(), 10 );

    FiniteStrainJ2Plasticity::Deformation< 3 > def;
    FiniteStrainJ2Plasticity::TimeIncrement    timeInc = { 0, 0.1 };

    FiniteStrainJ2Plasticity::ConstitutiveResponse< 3 > response;
    FiniteStrainJ2Plasticity::AlgorithmicModuli< 3 >    tangent;

    def.F( 0, 0 ) = 1.01;
    def.F( 0, 1 ) = 0.06;
    def.F( 0, 2 ) = -0.03;
    def.F( 1, 0 ) = 0.06;
    def.F( 1, 1 ) = 1.02;
    def.F( 1, 2 ) = 0.04;
    def.F( 2, 0 ) = -0.03;
    def.F( 2, 1 ) = 0.04;
    def.F( 2, 2 ) = 0.95;

    mat.computeStressWithFullReturnMapping( response, tangent, def, timeInc );

    Tensor33d stressUnrotated = response.tau;
    Tensor33d F_unrotated     = def.F;
    Tensor33d Fp_unrotated( 0.0 );
    std::memcpy( Fp_unrotated.data(), stateVars_.data(), 9 * sizeof( double ) );

    for ( int phi_deg = 0; phi_deg <= 180; phi_deg++ ) {
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

      mat.computeStressWithFullReturnMapping( response, tangent, def, timeInc );

      Tensor33d stressNew = response.tau;

      throwExceptionOnFailure( checkIfEqual( stressNew, stressUnrotated, 1e-10 ),
                               "I-4c: Isotropy test for arbitrary deformation failed "
                               "for "
                               "FiniteStrainJ2Plasticity "
                               "material in " +
                                 std::string( __PRETTY_FUNCTION__ ) );
    }
  }
}

////////////////////////////////////////
//        II test group - FDAF        //
////////////////////////////////////////

// Test II-1: Deformation response
void testDeformationResponseFDAF()
{
  // Test II-1a: Hydrostatic deformation
  {
    std::array< double, 7 > materialProperties_ = { 175000, 80800, 260, 580, 9, 70, 2 };
    const double            nMaterialProperties = 7;
    const int               elLabel             = 1;

    FiniteStrainJ2Plasticity mat = FiniteStrainJ2Plasticity( &materialProperties_[0], nMaterialProperties, elLabel );

    if ( mat.getNumberOfRequiredStateVars() > 10 ) {
      throw std::runtime_error( "Number of required state vars changed!" );
    }

    double                   alphaP = 0.0;
    Tensor33d                Fp     = Marmot::FastorStandardTensors::Spatial3D::I;
    std::array< double, 10 > stateVars_;
    std::memcpy( stateVars_.data(), Fp.data(), 9 * sizeof( double ) );
    stateVars_[9] = alphaP;

    mat.assignStateVars( stateVars_.data(), 10 );

    FiniteStrainJ2Plasticity::Deformation< 3 > def;
    FiniteStrainJ2Plasticity::TimeIncrement    timeInc = { 0, 0.1 };

    FiniteStrainJ2Plasticity::ConstitutiveResponse< 3 > response;
    FiniteStrainJ2Plasticity::AlgorithmicModuli< 3 >    tangent;

    def.F = Marmot::FastorStandardTensors::Spatial3D::I;
    def.F( 0, 0 ) += 0.002;
    def.F( 1, 1 ) += 0.002;
    def.F( 2, 2 ) += 0.002;

    mat.computeStressFDAF( response, tangent, def, timeInc );

    Tensor33d stressTarget( 0.0 );
    stressTarget( 0, 0 ) = 1048.97652265991;
    stressTarget( 1, 1 ) = 1048.97652265991;
    stressTarget( 2, 2 ) = 1048.97652265991;

    throwExceptionOnFailure( checkIfEqual( response.tau, stressTarget, 1e-10 ),
                             "II-1a: Hydrostatic FDAF - Kirchhoff stress tensor (tau) computation failed for "
                             "FiniteStrainJ2Plasticity material in " +
                               std::string( __PRETTY_FUNCTION__ ) );

    Tensor33d FpCurrent( 0.0 );
    std::memcpy( FpCurrent.data(), stateVars_.data(), 9 * sizeof( double ) );

    Tensor33d FpTarget = Marmot::FastorStandardTensors::Spatial3D::I;

    throwExceptionOnFailure( checkIfEqual( FpCurrent, FpTarget, 1e-10 ),
                             "II-1a: Hydrostatic FDAF - Plastic deformation gradient tensor (Fp) computation failed "
                             "for "
                             "FiniteStrainJ2Plasticity material in " +
                               std::string( __PRETTY_FUNCTION__ ) );

    double alphaPTarget = 0.0;

    throwExceptionOnFailure( checkIfEqual( stateVars_[9], alphaPTarget, 1e-10 ),
                             "II-1a: Hydrostatic FDAF - Strain-like hardening variable (alphaP) computation failed for "
                             "FiniteStrainJ2Plasticity material in " +
                               std::string( __PRETTY_FUNCTION__ ) );
  }

  // Test II-1b: Arbitrary deformation
  {
    std::array< double, 7 > materialProperties_ = { 175000, 80800, 260, 580, 9, 70, 2 };
    const double            nMaterialProperties = 7;
    const int               elLabel             = 1;

    FiniteStrainJ2Plasticity mat = FiniteStrainJ2Plasticity( &materialProperties_[0], nMaterialProperties, elLabel );

    if ( mat.getNumberOfRequiredStateVars() > 10 ) {
      throw std::runtime_error( "Number of required state vars changed!" );
    }

    double                   alphaP = 0.0;
    Tensor33d                Fp     = Marmot::FastorStandardTensors::Spatial3D::I;
    std::array< double, 10 > stateVars_;
    std::memcpy( stateVars_.data(), Fp.data(), 9 * sizeof( double ) );
    stateVars_[9] = alphaP;

    mat.assignStateVars( stateVars_.data(), 10 );

    FiniteStrainJ2Plasticity::Deformation< 3 > def;
    FiniteStrainJ2Plasticity::TimeIncrement    timeInc = { 0, 0.1 };

    FiniteStrainJ2Plasticity::ConstitutiveResponse< 3 > response;
    FiniteStrainJ2Plasticity::AlgorithmicModuli< 3 >    tangent;

    def.F         = Marmot::FastorStandardTensors::Spatial3D::I;
    def.F( 0, 0 ) = 1.01;
    def.F( 0, 1 ) = 0.06;
    def.F( 0, 2 ) = -0.03;
    def.F( 1, 0 ) = 0.06;
    def.F( 1, 1 ) = 1.02;
    def.F( 1, 2 ) = 0.04;
    def.F( 2, 0 ) = -0.03;
    def.F( 2, 1 ) = 0.04;
    def.F( 2, 2 ) = 0.95;

    mat.computeStressFDAF( response, tangent, def, timeInc );

    Tensor33d stressTarget( 0.0 );
    stressTarget( 0, 0 ) = -4815.15772271671;
    stressTarget( 0, 1 ) = 179.062345389071;
    stressTarget( 0, 2 ) = -95.4629813550829;
    stressTarget( 1, 0 ) = 179.062345389071;
    stressTarget( 1, 1 ) = -4786.68287740687;
    stressTarget( 1, 2 ) = 124.546217296781;
    stressTarget( 2, 0 ) = -95.4629813550829;
    stressTarget( 2, 1 ) = 124.546217296781;
    stressTarget( 2, 2 ) = -4995.74104468542;

    throwExceptionOnFailure( checkIfEqual( response.tau, stressTarget, 1e-10 ),
                             "II-2b: Arbitrary deformation FDAF - Kirchhoff stress tensor (tau) computation failed for "
                             "FiniteStrainJ2Plasticity material in " +
                               std::string( __PRETTY_FUNCTION__ ) );

    for ( int i = 0; i < 3; i++ )
      for ( int j = 0; j < 3; j++ )

        throwExceptionOnFailure( checkIfEqual( response.tau( i, j ), response.tau( j, i ), 1e-10 ),
                                 "II-2b: Kirchhoff stress tensor symmetry check for the arbitrary deformation load "
                                 "case failed for FiniteStrainJ2Plasticity "
                                 "material in " +
                                   std::string( __PRETTY_FUNCTION__ ) );

    Tensor33d FpCurrent( 0.0 );
    std::memcpy( FpCurrent.data(), stateVars_.data(), 9 * sizeof( double ) );

    Tensor33d FpTarget( 0.0 );
    FpTarget( 0, 0 ) = 1.01899554888122;
    FpTarget( 0, 1 ) = 0.0594232892187332;
    FpTarget( 0, 2 ) = -0.0297463602782848;
    FpTarget( 1, 0 ) = 0.0594232892187332;
    FpTarget( 1, 1 ) = 1.02889141910946;
    FpTarget( 1, 2 ) = 0.0396457910879559;
    FpTarget( 2, 0 ) = -0.0297463602782848;
    FpTarget( 2, 1 ) = 0.0396457910879559;
    FpTarget( 2, 2 ) = 0.959563358208919;

    throwExceptionOnFailure( checkIfEqual( FpCurrent, FpTarget, 1e-10 ),
                             "II-2b: Arbitrary deformation FDAF - Plastic deformation gradient tensor (Fp) computation "
                             "failed for "
                             "FiniteStrainJ2Plasticity material in " +
                               std::string( __PRETTY_FUNCTION__ ) );

    double alphaPTarget = 0.0998876084740521;

    throwExceptionOnFailure( checkIfEqual( stateVars_[9], alphaPTarget, 1e-10 ),
                             "II-2b: Arbitrary deformation FDAF - Strain-like hardening variable (alphaP) computation "
                             "failed "
                             "for "
                             "FiniteStrainJ2Plasticity material in " +
                               std::string( __PRETTY_FUNCTION__ ) );
  }
}

// Test II-2: Algorithmic tangent
void testAlgorithmicTangentFDAF()

{
  std::array< double, 7 > materialProperties_ = { 175000, 80800, 260, 580, 9, 70, 2 };
  const double            nMaterialProperties = 7;
  const int               elLabel             = 1;

  FiniteStrainJ2Plasticity mat = FiniteStrainJ2Plasticity( &materialProperties_[0], nMaterialProperties, elLabel );

  if ( mat.getNumberOfRequiredStateVars() > 10 ) {
    throw std::runtime_error( "Number of required state vars changed!" );
  }

  double                   alphaP = 0.0;
  Tensor33d                Fp     = Marmot::FastorStandardTensors::Spatial3D::I;
  std::array< double, 10 > stateVars_;
  std::memcpy( stateVars_.data(), Fp.data(), 9 * sizeof( double ) );
  stateVars_[9] = alphaP;

  mat.assignStateVars( stateVars_.data(), 10 );

  FiniteStrainJ2Plasticity::Deformation< 3 > def;
  FiniteStrainJ2Plasticity::TimeIncrement    timeInc = { 0, 0.1 };

  FiniteStrainJ2Plasticity::ConstitutiveResponse< 3 > response;
  FiniteStrainJ2Plasticity::AlgorithmicModuli< 3 >    tangent;

  def.F = Marmot::FastorStandardTensors::Spatial3D::I;
  def.F( 0, 0 ) += 0.001;
  def.F( 1, 1 ) += 0.002;
  def.F( 2, 2 ) += 0.003;

  mat.computeStressFDAF( response, tangent, def, timeInc );

  Tensor3333d tangentTarget( 0.0 );

  tangentTarget( 0, 0, 0, 0 ) = 200891.207248906;
  tangentTarget( 0, 0, 0, 1 ) = 0.293344510107489;
  tangentTarget( 0, 0, 0, 2 ) = 0.29385296412801;
  tangentTarget( 0, 0, 1, 0 ) = 0.292965005011075;
  tangentTarget( 0, 0, 1, 1 ) = 124564.921473659;
  tangentTarget( 0, 0, 1, 2 ) = 0.294008477858679;
  tangentTarget( 0, 0, 2, 0 ) = 0.293582659210213;
  tangentTarget( 0, 0, 2, 1 ) = 0.294117661344026;
  tangentTarget( 0, 0, 2, 2 ) = 198536.580316643;
  tangentTarget( 0, 1, 0, 1 ) = 75132.5523290578;
  tangentTarget( 0, 1, 1, 0 ) = 75057.5697419029;
  tangentTarget( 0, 2, 0, 2 ) = 75197.6138240272;
  tangentTarget( 0, 2, 2, 0 ) = 75047.6684325537;
  tangentTarget( 1, 0, 0, 1 ) = 75132.5523290578;
  tangentTarget( 1, 0, 1, 0 ) = 75057.5697419029;
  tangentTarget( 1, 1, 0, 0 ) = 124689.184484852;
  tangentTarget( 1, 1, 0, 1 ) = -0.001343538153446;
  tangentTarget( 1, 1, 0, 2 ) = -0.00134265592944414;
  tangentTarget( 1, 1, 1, 0 ) = -0.00134406872000709;
  tangentTarget( 1, 1, 1, 1 ) = 274826.634674348;
  tangentTarget( 1, 1, 1, 2 ) = -0.00134226715396347;
  tangentTarget( 1, 1, 2, 0 ) = -0.00134303384086018;
  tangentTarget( 1, 1, 2, 1 ) = -0.00134211451069891;
  tangentTarget( 1, 1, 2, 2 ) = 124474.117290285;
  tangentTarget( 1, 2, 1, 2 ) = 75187.7026281932;
  tangentTarget( 1, 2, 2, 1 ) = 75112.7398140078;
  tangentTarget( 2, 0, 0, 2 ) = 75197.6138240272;
  tangentTarget( 2, 0, 2, 0 ) = 75047.6684325537;
  tangentTarget( 2, 1, 1, 2 ) = 75187.7026281932;
  tangentTarget( 2, 1, 2, 1 ) = 75112.7398140078;
  tangentTarget( 2, 2, 0, 0 ) = 198934.179430311;
  tangentTarget( 2, 2, 0, 1 ) = -0.299419075895853;
  tangentTarget( 2, 2, 0, 2 ) = -0.29992769650696;
  tangentTarget( 2, 2, 1, 0 ) = -0.299039033663026;
  tangentTarget( 2, 2, 1, 1 ) = 124598.157567648;
  tangentTarget( 2, 2, 1, 2 ) = -0.30008287730396;
  tangentTarget( 2, 2, 2, 0 ) = -0.299657009019649;
  tangentTarget( 2, 2, 2, 1 ) = -0.300192215324636;
  tangentTarget( 2, 2, 2, 2 ) = 200455.224503959;

  throwExceptionOnFailure( checkIfEqual( tangent.dTau_dF, tangentTarget, 1e-8 ),
                           "II-2: FDAF Algorithmic tangent computation failed for FiniteStrainJ2Plasticity material "
                           "in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

////////////////////////////////////////
//       III test group - FDAC        //
////////////////////////////////////////

// Test III-1: Deformation response
void testDeformationResponseFDAC()
{
  // Test III-1a: Hydrostatic deformation
  {
    std::array< double, 7 > materialProperties_ = { 175000, 80800, 260, 580, 9, 70, 3 };
    const double            nMaterialProperties = 7;
    const int               elLabel             = 1;

    FiniteStrainJ2Plasticity mat = FiniteStrainJ2Plasticity( &materialProperties_[0], nMaterialProperties, elLabel );

    if ( mat.getNumberOfRequiredStateVars() > 10 ) {
      throw std::runtime_error( "Number of required state vars changed!" );
    }

    double                   alphaP = 0.0;
    Tensor33d                Fp     = Marmot::FastorStandardTensors::Spatial3D::I;
    std::array< double, 10 > stateVars_;
    std::memcpy( stateVars_.data(), Fp.data(), 9 * sizeof( double ) );
    stateVars_[9] = alphaP;

    mat.assignStateVars( stateVars_.data(), 10 );

    FiniteStrainJ2Plasticity::Deformation< 3 > def;
    FiniteStrainJ2Plasticity::TimeIncrement    timeInc = { 0, 0.1 };

    FiniteStrainJ2Plasticity::ConstitutiveResponse< 3 > response;
    FiniteStrainJ2Plasticity::AlgorithmicModuli< 3 >    tangent;

    def.F = Marmot::FastorStandardTensors::Spatial3D::I;
    def.F( 0, 0 ) += 0.002;
    def.F( 1, 1 ) += 0.002;
    def.F( 2, 2 ) += 0.002;

    mat.computeStressFDAC( response, tangent, def, timeInc );

    Tensor33d stressTarget( 0.0 );
    stressTarget( 0, 0 ) = 1048.97652265991;
    stressTarget( 1, 1 ) = 1048.97652265991;
    stressTarget( 2, 2 ) = 1048.97652265991;

    throwExceptionOnFailure( checkIfEqual( response.tau, stressTarget, 1e-10 ),
                             "III-1a: Hydrostatic FDAC - Kirchhoff stress tensor (tau) computation failed for "
                             "FiniteStrainJ2Plasticity material in " +
                               std::string( __PRETTY_FUNCTION__ ) );

    Tensor33d FpCurrent( 0.0 );
    std::memcpy( FpCurrent.data(), stateVars_.data(), 9 * sizeof( double ) );

    Tensor33d FpTarget = Marmot::FastorStandardTensors::Spatial3D::I;

    throwExceptionOnFailure( checkIfEqual( FpCurrent, FpTarget, 1e-10 ),
                             "III-1a: Hydrostatic FDAC - Plastic deformation gradient tensor (Fp) computation failed "
                             "for "
                             "FiniteStrainJ2Plasticity material in " +
                               std::string( __PRETTY_FUNCTION__ ) );

    double alphaPTarget = 0.0;

    throwExceptionOnFailure( checkIfEqual( stateVars_[9], alphaPTarget, 1e-10 ),
                             "III-1a: Hydrostatic FDAC - Strain-like hardening variable (alphaP) computation failed "
                             "for "
                             "FiniteStrainJ2Plasticity material in " +
                               std::string( __PRETTY_FUNCTION__ ) );
  }

  // Test III-1b: Arbitrary deformation
  {
    std::array< double, 7 > materialProperties_ = { 175000, 80800, 260, 580, 9, 70, 3 };
    const double            nMaterialProperties = 7;
    const int               elLabel             = 1;

    FiniteStrainJ2Plasticity mat = FiniteStrainJ2Plasticity( &materialProperties_[0], nMaterialProperties, elLabel );

    if ( mat.getNumberOfRequiredStateVars() > 10 ) {
      throw std::runtime_error( "Number of required state vars changed!" );
    }

    double                   alphaP = 0.0;
    Tensor33d                Fp     = Marmot::FastorStandardTensors::Spatial3D::I;
    std::array< double, 10 > stateVars_;
    std::memcpy( stateVars_.data(), Fp.data(), 9 * sizeof( double ) );
    stateVars_[9] = alphaP;

    mat.assignStateVars( stateVars_.data(), 10 );

    FiniteStrainJ2Plasticity::Deformation< 3 > def;
    FiniteStrainJ2Plasticity::TimeIncrement    timeInc = { 0, 0.1 };

    FiniteStrainJ2Plasticity::ConstitutiveResponse< 3 > response;
    FiniteStrainJ2Plasticity::AlgorithmicModuli< 3 >    tangent;

    def.F         = Marmot::FastorStandardTensors::Spatial3D::I;
    def.F( 0, 0 ) = 1.01;
    def.F( 0, 1 ) = 0.06;
    def.F( 0, 2 ) = -0.03;
    def.F( 1, 0 ) = 0.06;
    def.F( 1, 1 ) = 1.02;
    def.F( 1, 2 ) = 0.04;
    def.F( 2, 0 ) = -0.03;
    def.F( 2, 1 ) = 0.04;
    def.F( 2, 2 ) = 0.95;

    mat.computeStressFDAC( response, tangent, def, timeInc );

    Tensor33d stressTarget( 0.0 );
    stressTarget( 0, 0 ) = -4815.15772271685;
    stressTarget( 0, 1 ) = 179.062345389071;
    stressTarget( 0, 2 ) = -95.4629813550829;
    stressTarget( 1, 0 ) = 179.062345389071;
    stressTarget( 1, 1 ) = -4786.682877407;
    stressTarget( 1, 2 ) = 124.546217296781;
    stressTarget( 2, 0 ) = -95.4629813550829;
    stressTarget( 2, 1 ) = 124.546217296781;
    stressTarget( 2, 2 ) = -4995.74104468542;

    throwExceptionOnFailure( checkIfEqual( response.tau, stressTarget, 1e-10 ),
                             "III-2b: Arbitrary deformation FDAC - Kirchhoff stress tensor (tau) computation failed "
                             "for "
                             "FiniteStrainJ2Plasticity material in " +
                               std::string( __PRETTY_FUNCTION__ ) );

    for ( int i = 0; i < 3; i++ )
      for ( int j = 0; j < 3; j++ )

        throwExceptionOnFailure( checkIfEqual( response.tau( i, j ), response.tau( j, i ), 1e-10 ),
                                 "III-2b: Kirchhoff stress tensor symmetry check for the arbitrary deformation load "
                                 "case failed for FiniteStrainJ2Plasticity "
                                 "material in " +
                                   std::string( __PRETTY_FUNCTION__ ) );

    Tensor33d FpCurrent( 0.0 );
    std::memcpy( FpCurrent.data(), stateVars_.data(), 9 * sizeof( double ) );

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

    throwExceptionOnFailure( checkIfEqual( FpCurrent, FpTarget, 1e-10 ),
                             "III-2b: Arbitrary deformation FDAC - Plastic deformation gradient tensor (Fp) "
                             "computation "
                             "failed for "
                             "FiniteStrainJ2Plasticity material in " +
                               std::string( __PRETTY_FUNCTION__ ) );

    double alphaPTarget = 0.0998876084740522;

    throwExceptionOnFailure( checkIfEqual( stateVars_[9], alphaPTarget, 1e-10 ),
                             "III-2b: Arbitrary deformation FDAC - Strain-like hardening variable (alphaP) computation "
                             "failed "
                             "for "
                             "FiniteStrainJ2Plasticity material in " +
                               std::string( __PRETTY_FUNCTION__ ) );
  }
}

// Test III-2: Algorithmic tangent
void testAlgorithmicTangentFDAC()

{
  std::array< double, 7 > materialProperties_ = { 175000, 80800, 260, 580, 9, 70, 3 };
  const double            nMaterialProperties = 7;
  const int               elLabel             = 1;

  FiniteStrainJ2Plasticity mat = FiniteStrainJ2Plasticity( &materialProperties_[0], nMaterialProperties, elLabel );

  if ( mat.getNumberOfRequiredStateVars() > 10 ) {
    throw std::runtime_error( "Number of required state vars changed!" );
  }

  double                   alphaP = 0.0;
  Tensor33d                Fp     = Marmot::FastorStandardTensors::Spatial3D::I;
  std::array< double, 10 > stateVars_;
  std::memcpy( stateVars_.data(), Fp.data(), 9 * sizeof( double ) );
  stateVars_[9] = alphaP;

  mat.assignStateVars( stateVars_.data(), 10 );

  FiniteStrainJ2Plasticity::Deformation< 3 > def;
  FiniteStrainJ2Plasticity::TimeIncrement    timeInc = { 0, 0.1 };

  FiniteStrainJ2Plasticity::ConstitutiveResponse< 3 > response;
  FiniteStrainJ2Plasticity::AlgorithmicModuli< 3 >    tangent;

  def.F = Marmot::FastorStandardTensors::Spatial3D::I;
  def.F( 0, 0 ) += 0.001;
  def.F( 1, 1 ) += 0.002;
  def.F( 2, 2 ) += 0.003;

  mat.computeStressFDAC( response, tangent, def, timeInc );

  Tensor3333d tangentTarget( 0.0 );
  tangentTarget( 0, 0, 0, 0 ) = 200890.509333366;
  tangentTarget( 0, 0, 1, 1 ) = 124564.488040833;
  tangentTarget( 0, 0, 2, 2 ) = 198537.025177391;
  tangentTarget( 0, 1, 0, 1 ) = 75132.5793713203;
  tangentTarget( 0, 1, 1, 0 ) = 75057.5967841801;
  tangentTarget( 1, 0, 0, 1 ) = 75132.5793713203;
  tangentTarget( 1, 0, 1, 0 ) = 75057.5967841801;
  tangentTarget( 0, 2, 0, 2 ) = 75197.6409636925;
  tangentTarget( 0, 2, 2, 0 ) = 75047.69557222;
  tangentTarget( 2, 0, 0, 2 ) = 75197.6409636925;
  tangentTarget( 2, 0, 2, 0 ) = 75047.69557222;
  tangentTarget( 1, 1, 0, 0 ) = 124688.985901659;
  tangentTarget( 1, 1, 1, 1 ) = 274826.696170781;
  tangentTarget( 1, 1, 2, 2 ) = 124474.382115553;
  tangentTarget( 1, 2, 1, 2 ) = 75187.7298654265;
  tangentTarget( 1, 2, 2, 1 ) = 75112.7670512278;
  tangentTarget( 2, 1, 1, 2 ) = 75187.7298654265;
  tangentTarget( 2, 1, 2, 1 ) = 75112.7670512278;
  tangentTarget( 2, 2, 0, 0 ) = 198933.702151485;
  tangentTarget( 2, 2, 1, 1 ) = 124598.548812737;
  tangentTarget( 2, 2, 2, 2 ) = 200455.905150422;

  throwExceptionOnFailure( checkIfEqual( tangent.dTau_dF, tangentTarget, 1e-8 ),
                           "III-2: FDAC Algorithmic tangent computation failed for FiniteStrainJ2Plasticity material "
                           "in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

////////////////////////////////////////
//        IV test group - CSDA        //
////////////////////////////////////////

// Test IV-1: Deformation response
void testDeformationResponseCSDA()
{
  // Test IV-1a: Hydrostatic deformation
  {
    std::array< double, 7 > materialProperties_ = { 175000, 80800, 260, 580, 9, 70, 4 };
    const double            nMaterialProperties = 7;
    const int               elLabel             = 1;

    FiniteStrainJ2Plasticity mat = FiniteStrainJ2Plasticity( &materialProperties_[0], nMaterialProperties, elLabel );

    if ( mat.getNumberOfRequiredStateVars() > 10 ) {
      throw std::runtime_error( "Number of required state vars changed!" );
    }

    double                   alphaP = 0.0;
    Tensor33d                Fp     = Marmot::FastorStandardTensors::Spatial3D::I;
    std::array< double, 10 > stateVars_;
    std::memcpy( stateVars_.data(), Fp.data(), 9 * sizeof( double ) );
    stateVars_[9] = alphaP;

    mat.assignStateVars( stateVars_.data(), 10 );

    FiniteStrainJ2Plasticity::Deformation< 3 > def;
    FiniteStrainJ2Plasticity::TimeIncrement    timeInc = { 0, 0.1 };

    FiniteStrainJ2Plasticity::ConstitutiveResponse< 3 > response;
    FiniteStrainJ2Plasticity::AlgorithmicModuli< 3 >    tangent;

    def.F = Marmot::FastorStandardTensors::Spatial3D::I;
    def.F( 0, 0 ) += 0.002;
    def.F( 1, 1 ) += 0.002;
    def.F( 2, 2 ) += 0.002;

    mat.computeStressCSDA( response, tangent, def, timeInc );

    Tensor33d stressTarget( 0.0 );
    stressTarget( 0, 0 ) = 1048.97652265991;
    stressTarget( 1, 1 ) = 1048.97652265991;
    stressTarget( 2, 2 ) = 1048.97652265991;

    throwExceptionOnFailure( checkIfEqual( response.tau, stressTarget, 1e-10 ),
                             "IV-1a: Hydrostatic CSDA - Kirchhoff stress tensor (tau) computation failed for "
                             "FiniteStrainJ2Plasticity material in " +
                               std::string( __PRETTY_FUNCTION__ ) );

    Tensor33d FpCurrent( 0.0 );
    std::memcpy( FpCurrent.data(), stateVars_.data(), 9 * sizeof( double ) );

    Tensor33d FpTarget = Marmot::FastorStandardTensors::Spatial3D::I;

    throwExceptionOnFailure( checkIfEqual( FpCurrent, FpTarget, 1e-10 ),
                             "IV-1a: Hydrostatic CSDA - Plastic deformation gradient tensor (Fp) computation failed "
                             "for "
                             "FiniteStrainJ2Plasticity material in " +
                               std::string( __PRETTY_FUNCTION__ ) );

    double alphaPTarget = 0.0;

    throwExceptionOnFailure( checkIfEqual( stateVars_[9], alphaPTarget, 1e-10 ),
                             "IV-1a: Hydrostatic CSDA - Strain-like hardening variable (alphaP) computation failed for "
                             "FiniteStrainJ2Plasticity material in " +
                               std::string( __PRETTY_FUNCTION__ ) );
  }

  // Test IV-1b: Arbitrary deformation
  {
    std::array< double, 7 > materialProperties_ = { 175000, 80800, 260, 580, 9, 70, 4 };
    const double            nMaterialProperties = 7;
    const int               elLabel             = 1;

    FiniteStrainJ2Plasticity mat = FiniteStrainJ2Plasticity( &materialProperties_[0], nMaterialProperties, elLabel );

    if ( mat.getNumberOfRequiredStateVars() > 10 ) {
      throw std::runtime_error( "Number of required state vars changed!" );
    }

    double                   alphaP = 0.0;
    Tensor33d                Fp     = Marmot::FastorStandardTensors::Spatial3D::I;
    std::array< double, 10 > stateVars_;
    std::memcpy( stateVars_.data(), Fp.data(), 9 * sizeof( double ) );
    stateVars_[9] = alphaP;

    mat.assignStateVars( stateVars_.data(), 10 );

    FiniteStrainJ2Plasticity::Deformation< 3 > def;
    FiniteStrainJ2Plasticity::TimeIncrement    timeInc = { 0, 0.1 };

    FiniteStrainJ2Plasticity::ConstitutiveResponse< 3 > response;
    FiniteStrainJ2Plasticity::AlgorithmicModuli< 3 >    tangent;

    def.F         = Marmot::FastorStandardTensors::Spatial3D::I;
    def.F( 0, 0 ) = 1.01;
    def.F( 0, 1 ) = 0.06;
    def.F( 0, 2 ) = -0.03;
    def.F( 1, 0 ) = 0.06;
    def.F( 1, 1 ) = 1.02;
    def.F( 1, 2 ) = 0.04;
    def.F( 2, 0 ) = -0.03;
    def.F( 2, 1 ) = 0.04;
    def.F( 2, 2 ) = 0.95;

    mat.computeStressCSDA( response, tangent, def, timeInc );

    Tensor33d stressTarget( 0.0 );
    stressTarget( 0, 0 ) = -4815.15772271683;
    stressTarget( 0, 1 ) = 179.062345389071;
    stressTarget( 0, 2 ) = -95.4629813550829;
    stressTarget( 1, 0 ) = 179.062345389071;
    stressTarget( 1, 1 ) = -4786.68287740697;
    stressTarget( 1, 2 ) = 124.546217296781;
    stressTarget( 2, 0 ) = -95.4629813550829;
    stressTarget( 2, 1 ) = 124.546217296781;
    stressTarget( 2, 2 ) = -4995.74104468553;

    throwExceptionOnFailure( checkIfEqual( response.tau, stressTarget, 1e-10 ),
                             "IV-2b: Arbitrary deformation CSDA - Kirchhoff stress tensor (tau) computation failed "
                             "for "
                             "FiniteStrainJ2Plasticity material in " +
                               std::string( __PRETTY_FUNCTION__ ) );

    for ( int i = 0; i < 3; i++ )
      for ( int j = 0; j < 3; j++ )

        throwExceptionOnFailure( checkIfEqual( response.tau( i, j ), response.tau( j, i ), 1e-10 ),
                                 "IV-2b: Kirchhoff stress tensor symmetry check for the arbitrary deformation load "
                                 "case failed for FiniteStrainJ2Plasticity "
                                 "material in " +
                                   std::string( __PRETTY_FUNCTION__ ) );

    Tensor33d FpCurrent( 0.0 );
    std::memcpy( FpCurrent.data(), stateVars_.data(), 9 * sizeof( double ) );

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

    throwExceptionOnFailure( checkIfEqual( FpCurrent, FpTarget, 1e-10 ),
                             "IV-2b: Arbitrary deformation CSDA - Plastic deformation gradient tensor (Fp) "
                             "computation "
                             "failed for "
                             "FiniteStrainJ2Plasticity material in " +
                               std::string( __PRETTY_FUNCTION__ ) );

    double alphaPTarget = 0.0998876084740522;

    throwExceptionOnFailure( checkIfEqual( stateVars_[9], alphaPTarget, 1e-10 ),
                             "IV-2b: Arbitrary deformation CSDA - Strain-like hardening variable (alphaP) computation "
                             "failed "
                             "for "
                             "FiniteStrainJ2Plasticity material in " +
                               std::string( __PRETTY_FUNCTION__ ) );
  }
}

// Test IV-2: Algorithmic tangent
void testAlgorithmicTangentCSDA()

{
  std::array< double, 7 > materialProperties_ = { 175000, 80800, 260, 580, 9, 70, 4 };
  const double            nMaterialProperties = 7;
  const int               elLabel             = 1;

  FiniteStrainJ2Plasticity mat = FiniteStrainJ2Plasticity( &materialProperties_[0], nMaterialProperties, elLabel );

  if ( mat.getNumberOfRequiredStateVars() > 10 ) {
    throw std::runtime_error( "Number of required state vars changed!" );
  }

  double                   alphaP = 0.0;
  Tensor33d                Fp     = Marmot::FastorStandardTensors::Spatial3D::I;
  std::array< double, 10 > stateVars_;
  std::memcpy( stateVars_.data(), Fp.data(), 9 * sizeof( double ) );
  stateVars_[9] = alphaP;

  mat.assignStateVars( stateVars_.data(), 10 );

  FiniteStrainJ2Plasticity::Deformation< 3 > def;
  FiniteStrainJ2Plasticity::TimeIncrement    timeInc = { 0, 0.1 };

  FiniteStrainJ2Plasticity::ConstitutiveResponse< 3 > response;
  FiniteStrainJ2Plasticity::AlgorithmicModuli< 3 >    tangent;

  def.F = Marmot::FastorStandardTensors::Spatial3D::I;
  def.F( 0, 0 ) += 0.001;
  def.F( 1, 1 ) += 0.002;
  def.F( 2, 2 ) += 0.003;

  mat.computeStressCSDA( response, tangent, def, timeInc );

  Tensor3333d tangentTarget( 0.0 );
  tangentTarget( 0, 0, 0, 0 ) = 200890.525305334;
  tangentTarget( 0, 0, 1, 1 ) = 124564.513352951;
  tangentTarget( 0, 0, 2, 2 ) = 198537.043113999;
  tangentTarget( 0, 1, 0, 1 ) = 75132.5523288947;
  tangentTarget( 0, 1, 1, 0 ) = 75057.5697417401;
  tangentTarget( 1, 0, 0, 1 ) = 75132.5523288947;
  tangentTarget( 1, 0, 1, 0 ) = 75057.5697417401;
  tangentTarget( 0, 2, 0, 2 ) = 75197.6138238636;
  tangentTarget( 0, 2, 2, 0 ) = 75047.6684323903;
  tangentTarget( 2, 0, 0, 2 ) = 75197.6138238636;
  tangentTarget( 2, 0, 2, 0 ) = 75047.6684323903;
  tangentTarget( 1, 1, 0, 0 ) = 124688.953426225;
  tangentTarget( 1, 1, 1, 1 ) = 274826.645041077;
  tangentTarget( 1, 1, 2, 2 ) = 124474.348696745;
  tangentTarget( 1, 2, 1, 2 ) = 75187.7026280292;
  tangentTarget( 1, 2, 2, 1 ) = 75112.7398138437;
  tangentTarget( 2, 1, 1, 2 ) = 75187.7026280292;
  tangentTarget( 2, 1, 2, 1 ) = 75112.7398138437;
  tangentTarget( 2, 2, 0, 0 ) = 198933.720522818;
  tangentTarget( 2, 2, 1, 1 ) = 124598.574593654;
  tangentTarget( 2, 2, 2, 2 ) = 200455.918711321;

  throwExceptionOnFailure( checkIfEqual( tangent.dTau_dF, tangentTarget, 1e-8 ),
                           "IV-2: CSDA Algorithmic tangent computation failed for FiniteStrainJ2Plasticity material "
                           "in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

int main()
{

  auto tests = std::vector< std::function< void() > >{ testUndeformedResponseFullReturnMapping,
                                                       testDeformationResponseFullReturnMapping,
                                                       testAlgorithmicTangentFullReturnMapping,
                                                       testRotationFullReturnMapping,
                                                       testDeformationResponseFDAF,
                                                       testAlgorithmicTangentFDAF,
                                                       testDeformationResponseFDAC,
                                                       testAlgorithmicTangentFDAC,
                                                       testDeformationResponseCSDA,
                                                       testAlgorithmicTangentCSDA };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
