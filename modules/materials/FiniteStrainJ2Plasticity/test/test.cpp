#include "Marmot/FiniteStrainJ2Plasticity.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotTesting.h"

using namespace Marmot::Testing;
using namespace Marmot::Materials;
using namespace Marmot::FastorStandardTensors;
using namespace Marmot::FastorIndices;

////////////////////////////////////////
// I test group - full return mapping //
////////////////////////////////////////

// Test I-1: Undeformed configuration
void testUndeformedResponseFullReturnMapping()
{
  // Material properties: K, G, fy, fyInf, eta, H, implementation type, (density)
  // Implementation type: 0 - scalar return mapping (not implemented yet), 1 - full return mapping, 2 - FDAF, 3 - FDAC,
  // 4 - CSDA
  std::array< double, 7 > materialProperties_ = { 175000, 80800, 260, 580, 9, 70, 1 };
  const double            nMaterialProperties = 7;
  const int               elLabel             = 1;

  // Create material instance
  FiniteStrainJ2Plasticity mat = FiniteStrainJ2Plasticity( &materialProperties_[0], nMaterialProperties, elLabel );

  // Set state vars
  if ( mat.getNumberOfRequiredStateVars() > 10 ) {
    throw std::runtime_error( "Number of required state vars changed!" );
  }

  double                   alphaP = 0.0;                                         // initial equivalent plastic strain
  Tensor33d                Fp     = Marmot::FastorStandardTensors::Spatial3D::I; // initial plastic deformation gradient
  std::array< double, 10 > stateVars_ = { Fp.data()[0],
                                          Fp.data()[1],
                                          Fp.data()[2],
                                          Fp.data()[3],
                                          Fp.data()[4],
                                          Fp.data()[5],
                                          Fp.data()[6],
                                          Fp.data()[7],
                                          Fp.data()[8],
                                          alphaP };

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

  // Target Kirchoff stress values for the given deformation gradient and material parameters
  Tensor33d stressTarget( 0.0 );

  // Compare computed stress to target stress values
  throwExceptionOnFailure( checkIfEqual( response.tau, stressTarget, 1e-10 ),
                           "I-1 Test: Undeformed configuration failed for FiniteStrainJ2Plasticity material in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Test I-2: Deformation response

void testDeformationResponseFullReturnMapping()
{

  // Test I-2a: Uniaxial stretch deformation
  {
    std::array< double, 7 > materialProperties_ = { 175000, 80800, 260, 580, 9, 70, 1 };
    const double            nMaterialProperties = 7;
    const int               elLabel             = 1;

    FiniteStrainJ2Plasticity mat = FiniteStrainJ2Plasticity( &materialProperties_[0], nMaterialProperties, elLabel );

    if ( mat.getNumberOfRequiredStateVars() > 10 ) {
      throw std::runtime_error( "Number of required state vars changed!" );
    }

    double                   alphaP = 0.0;          
    Tensor33d                Fp = Marmot::FastorStandardTensors::Spatial3D::I;
    std::array< double, 10 > stateVars_ = { Fp.data()[0],
                                            Fp.data()[1],
                                            Fp.data()[2],
                                            Fp.data()[3],
                                            Fp.data()[4],
                                            Fp.data()[5],
                                            Fp.data()[6],
                                            Fp.data()[7],
                                            Fp.data()[8],
                                            alphaP };

    mat.assignStateVars( stateVars_.data(), 10 );

    FiniteStrainJ2Plasticity::Deformation< 3 > def;
    FiniteStrainJ2Plasticity::TimeIncrement    timeInc = { 0, 0.1 };

    FiniteStrainJ2Plasticity::ConstitutiveResponse< 3 > response;
    FiniteStrainJ2Plasticity::AlgorithmicModuli< 3 >    tangent;

    def.F = Marmot::FastorStandardTensors::Spatial3D::I;
    def.F( 0, 0 ) += 0.008;

    mat.computeStressWithFullReturnMapping( response, tangent, def, timeInc );

  Tensor33d stressTarget( 0.0 );
  stressTarget( 0, 0 ) = 109.197470795418;
  stressTarget( 0, 1 ) = 0.;
  stressTarget( 0, 2 ) = 0.;
  stressTarget( 1, 0 ) = 0.;
  stressTarget( 1, 1 ) = 49.3922392274353;
  stressTarget( 1, 2 ) = 0.;
  stressTarget( 2, 0 ) = 0.;
  stressTarget( 2, 1 ) = 0.;
  stressTarget( 2, 2 ) = 49.3922392274353;

   /* throwExceptionOnFailure( checkIfEqual( response.tau, stressTarget, 1e-10 ),
                             "I-2a: Uniaxial stretch failed for FiniteStrainJ2Plasticity material in " +
                               std::string( __PRETTY_FUNCTION__ ) );

    */                           

    std::cout << "1. Test deformation=\n" << def.F << std::endl;
    std::cout << "1. Test stress = \n" << response.tau << std::endl;
    // std::cout << "1. Test tangent = \n" << tangent.dTau_dF << std::endl;
    std::cout << "StateVars" << stateVars_[0] << " " << stateVars_[1] << " " << stateVars_[2] << " " << stateVars_[3]
              << " " << stateVars_[4] << " " << stateVars_[5] << " " << stateVars_[6] << " " << stateVars_[7] << " "
              << stateVars_[8] << " " << stateVars_[9] << std::endl;
  }


 // Test I-2b: Simple shear deformation
  {
    std::array< double, 7 > materialProperties_ = { 175000, 80800, 260, 580, 9, 70, 1 };
    const double            nMaterialProperties = 7;
    const int               elLabel             = 1;

    FiniteStrainJ2Plasticity mat = FiniteStrainJ2Plasticity( &materialProperties_[0], nMaterialProperties, elLabel );

    if ( mat.getNumberOfRequiredStateVars() > 10 ) {
      throw std::runtime_error( "Number of required state vars changed!" );
    }

    double                   alphaP = 0.0;          
    Tensor33d                Fp = Marmot::FastorStandardTensors::Spatial3D::I;
    std::array< double, 10 > stateVars_ = { Fp.data()[0],
                                            Fp.data()[1],
                                            Fp.data()[2],
                                            Fp.data()[3],
                                            Fp.data()[4],
                                            Fp.data()[5],
                                            Fp.data()[6],
                                            Fp.data()[7],
                                            Fp.data()[8],
                                            alphaP };

    mat.assignStateVars( stateVars_.data(), 10 );

    FiniteStrainJ2Plasticity::Deformation< 3 > def;
    FiniteStrainJ2Plasticity::TimeIncrement    timeInc = { 0, 0.1 };

    FiniteStrainJ2Plasticity::ConstitutiveResponse< 3 > response;
    FiniteStrainJ2Plasticity::AlgorithmicModuli< 3 >    tangent;

    def.F = Marmot::FastorStandardTensors::Spatial3D::I;
    def.F( 1, 0 ) += 0.02;

    mat.computeStressWithFullReturnMapping( response, tangent, def, timeInc );

  Tensor33d stressTarget( 0.0 );
  stressTarget( 0, 0 ) = 109.197470795418;
  stressTarget( 0, 1 ) = 0.;
  stressTarget( 0, 2 ) = 0.;
  stressTarget( 1, 0 ) = 0.;
  stressTarget( 1, 1 ) = 49.3922392274353;
  stressTarget( 1, 2 ) = 0.;
  stressTarget( 2, 0 ) = 0.;
  stressTarget( 2, 1 ) = 0.;
  stressTarget( 2, 2 ) = 49.3922392274353;

   /* throwExceptionOnFailure( checkIfEqual( response.tau, stressTarget, 1e-10 ),
                             "I-2a: Uniaxial stretch failed for FiniteStrainJ2Plasticity material in " +
                               std::string( __PRETTY_FUNCTION__ ) );

    */                           

    std::cout << "1. Test deformation=\n" << def.F << std::endl;
    std::cout << "1. Test stress = \n" << response.tau << std::endl;
    // std::cout << "1. Test tangent = \n" << tangent.dTau_dF << std::endl;
    std::cout << "StateVars" << stateVars_[0] << " " << stateVars_[1] << " " << stateVars_[2] << " " << stateVars_[3]
              << " " << stateVars_[4] << " " << stateVars_[5] << " " << stateVars_[6] << " " << stateVars_[7] << " "
              << stateVars_[8] << " " << stateVars_[9] << std::endl;
  }




 // Test I-2c: Hydrostatic deformation
  {
    std::array< double, 7 > materialProperties_ = { 175000, 80800, 260, 580, 9, 70, 1 };
    const double            nMaterialProperties = 7;
    const int               elLabel             = 1;

    FiniteStrainJ2Plasticity mat = FiniteStrainJ2Plasticity( &materialProperties_[0], nMaterialProperties, elLabel );

    if ( mat.getNumberOfRequiredStateVars() > 10 ) {
      throw std::runtime_error( "Number of required state vars changed!" );
    }

    double                   alphaP = 0.0;          
    Tensor33d                Fp = Marmot::FastorStandardTensors::Spatial3D::I;
    std::array< double, 10 > stateVars_ = { Fp.data()[0],
                                            Fp.data()[1],
                                            Fp.data()[2],
                                            Fp.data()[3],
                                            Fp.data()[4],
                                            Fp.data()[5],
                                            Fp.data()[6],
                                            Fp.data()[7],
                                            Fp.data()[8],
                                            alphaP };

    mat.assignStateVars( stateVars_.data(), 10 );

    FiniteStrainJ2Plasticity::Deformation< 3 > def;
    FiniteStrainJ2Plasticity::TimeIncrement    timeInc = { 0, 0.1 };

    FiniteStrainJ2Plasticity::ConstitutiveResponse< 3 > response;
    FiniteStrainJ2Plasticity::AlgorithmicModuli< 3 >    tangent;

    def.F = Marmot::FastorStandardTensors::Spatial3D::I;
    def.F( 0, 0 ) += 0.002;
    def.F( 1, 1 ) += 0.002;
    def.F( 1, 1 ) += 0.002;

    mat.computeStressWithFullReturnMapping( response, tangent, def, timeInc );

  Tensor33d stressTarget( 0.0 );
  stressTarget( 0, 0 ) = 109.197470795418;
  stressTarget( 0, 1 ) = 0.;
  stressTarget( 0, 2 ) = 0.;
  stressTarget( 1, 0 ) = 0.;
  stressTarget( 1, 1 ) = 49.3922392274353;
  stressTarget( 1, 2 ) = 0.;
  stressTarget( 2, 0 ) = 0.;
  stressTarget( 2, 1 ) = 0.;
  stressTarget( 2, 2 ) = 49.3922392274353;

   /* throwExceptionOnFailure( checkIfEqual( response.tau, stressTarget, 1e-10 ),
                             "I-2a: Uniaxial stretch failed for FiniteStrainJ2Plasticity material in " +
                               std::string( __PRETTY_FUNCTION__ ) );

    */                           

    std::cout << "1. Test deformation=\n" << def.F << std::endl;
    std::cout << "1. Test stress = \n" << response.tau << std::endl;
    // std::cout << "1. Test tangent = \n" << tangent.dTau_dF << std::endl;
    std::cout << "StateVars" << stateVars_[0] << " " << stateVars_[1] << " " << stateVars_[2] << " " << stateVars_[3]
              << " " << stateVars_[4] << " " << stateVars_[5] << " " << stateVars_[6] << " " << stateVars_[7] << " "
              << stateVars_[8] << " " << stateVars_[9] << std::endl;
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
    Tensor33d                Fp = Marmot::FastorStandardTensors::Spatial3D::I;
    std::array< double, 10 > stateVars_ = { Fp.data()[0],
                                            Fp.data()[1],
                                            Fp.data()[2],
                                            Fp.data()[3],
                                            Fp.data()[4],
                                            Fp.data()[5],
                                            Fp.data()[6],
                                            Fp.data()[7],
                                            Fp.data()[8],
                                            alphaP };

    mat.assignStateVars( stateVars_.data(), 10 );

    FiniteStrainJ2Plasticity::Deformation< 3 > def;
    FiniteStrainJ2Plasticity::TimeIncrement    timeInc = { 0, 0.1 };

    FiniteStrainJ2Plasticity::ConstitutiveResponse< 3 > response;
    FiniteStrainJ2Plasticity::AlgorithmicModuli< 3 >    tangent;

    def.F = Marmot::FastorStandardTensors::Spatial3D::I;
    def.F( 0, 0 ) += 0.008;

    mat.computeStressWithFullReturnMapping( response, tangent, def, timeInc );

  Tensor33d stressTarget( 0.0 );
  stressTarget( 0, 0 ) = 109.197470795418;
  stressTarget( 0, 1 ) = 0.;
  stressTarget( 0, 2 ) = 0.;
  stressTarget( 1, 0 ) = 0.;
  stressTarget( 1, 1 ) = 49.3922392274353;
  stressTarget( 1, 2 ) = 0.;
  stressTarget( 2, 0 ) = 0.;
  stressTarget( 2, 1 ) = 0.;
  stressTarget( 2, 2 ) = 49.3922392274353;

   /* throwExceptionOnFailure( checkIfEqual( response.tau, stressTarget, 1e-10 ),
                             "I-2a: Uniaxial stretch failed for FiniteStrainJ2Plasticity material in " +
                               std::string( __PRETTY_FUNCTION__ ) );

    */                           

    std::cout << "1. Test deformation=\n" << def.F << std::endl;
    std::cout << "1. Test stress = \n" << response.tau << std::endl;
    // std::cout << "1. Test tangent = \n" << tangent.dTau_dF << std::endl;
    std::cout << "StateVars" << stateVars_[0] << " " << stateVars_[1] << " " << stateVars_[2] << " " << stateVars_[3]
              << " " << stateVars_[4] << " " << stateVars_[5] << " " << stateVars_[6] << " " << stateVars_[7] << " "
              << stateVars_[8] << " " << stateVars_[9] << std::endl;
  }


// Test I-4: Rotation



int main()
{

  auto tests = std::vector< std::function< void() > >{ testUndeformedResponseFullReturnMapping,
                                                       testDeformationResponseFullReturnMapping, testAlgorithmicTangentFullReturnMapping };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
