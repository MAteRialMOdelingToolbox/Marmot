#include "Marmot/CompressibleNeoHooke.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotTesting.h"

using namespace Marmot::Testing;
using namespace Marmot::Materials;
using namespace Marmot::FastorStandardTensors;


// 1. Test: Undeformed configuration
void testUndeformedResponse()
{
  // idx 0 - Bulk modulus K, idx 1 - Shear modulus G in MPa
  std::array< double, 2 > materialProperties_ = { 3500, 1500 };
  const double            nMaterialProperties = 2;
  const int               elLabel             = 1;

  // Create material instance
  CompressibleNeoHooke mat = CompressibleNeoHooke( &materialProperties_[0], nMaterialProperties, elLabel );

  // Create deformation, time increment, response and tangent objects required for stress computation
  CompressibleNeoHooke::Deformation< 3 > def;
  CompressibleNeoHooke::TimeIncrement    timeInc = { 0, 0.1 };

  CompressibleNeoHooke::ConstitutiveResponse< 3 > response;
  CompressibleNeoHooke::AlgorithmicModuli< 3 >    tangent;


  // Prescribe a deformation gradient tensor F for the considered load case
  def.F = Marmot::FastorStandardTensors::Spatial3D::I;

  // Compute stress response
  mat.computeStress( response, tangent, def, timeInc );

  // Target Kirchoff stress values for the given deformation gradient and material parameters
  Tensor33d stressTarget;
  stressTarget( 0, 0 ) = 0.;
  stressTarget( 0, 1 ) = 0.;
  stressTarget( 0, 2 ) = 0.;
  stressTarget( 1, 0 ) = 0.;
  stressTarget( 1, 1 ) = 0.;
  stressTarget( 1, 2 ) = 0.;
  stressTarget( 2, 0 ) = 0.;
  stressTarget( 2, 1 ) = 0.;
  stressTarget( 2, 2 ) = 0.;

  

  // Compare computed stress to target stress values
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      throwExceptionOnFailure(
        checkIfEqual( response.tau( i, j ), stressTarget( i, j ), 1e-10 ),
        "1. Test: Undeformed configuration failed for CompressibleNeoHooke material in " + std::string( __PRETTY_FUNCTION__ ) );
        
  
   /*std::cout << "1. Test deformation=\n" << def.F << std::endl; 
   std::cout << "1. Test stress = \n" << response.tau << std::endl; 
   std::cout << "1. Test tangent = \n" << tangent.dTau_dF << std::endl; */
}





// Test the uniaxial response of the CompressibleNeoHooke material model
void testFiniteStrainUniaxialResponse()
{
  // idx 0 - Bulk modulus K, idx 1 - Shear modulus G in MPa
  std::array< double, 2 > materialProperties_ = { 3500, 1500 };
  const double            nMaterialProperties = 2;
  const int               elLabel             = 1;

  // Create material instance
  CompressibleNeoHooke mat = CompressibleNeoHooke( &materialProperties_[0], nMaterialProperties, elLabel );

  // Create deformation, time increment, response and tangent objects required for stress computation
  CompressibleNeoHooke::Deformation< 3 > def;
  CompressibleNeoHooke::TimeIncrement    timeInc = { 0, 0.1 };

  CompressibleNeoHooke::ConstitutiveResponse< 3 > response;
  CompressibleNeoHooke::AlgorithmicModuli< 3 >    tangent;


  // Prescribe a deformation gradient tensor F for the considered load case
  def.F = Marmot::FastorStandardTensors::Spatial3D::I;
  def.F( 0, 0 ) += 0.2;

  // Compute stress response
  mat.computeStress( response, tangent, def, timeInc );

  // Target Kirchoff stress values for the given deformation gradient and material parameters
  Tensor33d stressTarget;
  stressTarget( 0, 0 ) = 1042.00258647807;
  stressTarget( 0, 1 ) = 0.;
  stressTarget( 0, 2 ) = 0.;
  stressTarget( 1, 0 ) = 0.;
  stressTarget( 1, 1 ) = 457.540373427632;
  stressTarget( 1, 2 ) = 0.;
  stressTarget( 2, 0 ) = 0.;
  stressTarget( 2, 1 ) = 0.;
  stressTarget( 2, 2 ) = 457.540373427632;

  // Compare computed stress to target stress values
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      throwExceptionOnFailure(
        checkIfEqual( response.tau( i, j ), stressTarget( i, j ), 1e-10 ),
        "Test 2: Finite Strain Uniaxial deformation computation failed for CompressibleNeoHooke material in " + std::string( __PRETTY_FUNCTION__ ) );
  
  // std::cout << "2. Test = \n" << def.F << std::endl; 
  // std::cout << "2. Test = \n" << response.tau << std::endl; 
  // std::cout << "2. Test Uniaxial = \n" << tangent.dTau_dF << std::endl; 
}


void testShearResponse()
{

  std::array< double, 2 > materialProperties_ = { 3500, 1500 };
  const double            nMaterialProperties = 2;
  const int               elLabel             = 1;

  CompressibleNeoHooke mat = CompressibleNeoHooke( &materialProperties_[0], nMaterialProperties, elLabel );

  CompressibleNeoHooke::Deformation< 3 > def;
  CompressibleNeoHooke::TimeIncrement    timeInc = { 0, 0.1 };

  CompressibleNeoHooke::ConstitutiveResponse< 3 > response;
  CompressibleNeoHooke::AlgorithmicModuli< 3 >    tangent;

  def.F = Marmot::FastorStandardTensors::Spatial3D::I;
  def.F( 1, 0 ) += 0.2;

  mat.computeStress( response, tangent, def, timeInc );

  Tensor33d stressTarget;
  stressTarget( 0, 0 ) = -120;
  stressTarget( 0, 1 ) = 300;
  stressTarget( 0, 2 ) = 0.;
  stressTarget( 1, 0 ) = 300;
  stressTarget( 1, 1 ) = 40;
  stressTarget( 1, 2 ) = 0.;
  stressTarget( 2, 0 ) = 0.;
  stressTarget( 2, 1 ) = 0.;
  stressTarget( 2, 2 ) = -20;

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      throwExceptionOnFailure(
        checkIfEqual( response.tau( i, j ), stressTarget( i, j ), 1e-10 ),
        "Uniaxial stress computation failed for CompressibleNeoHooke material in " + std::string( __PRETTY_FUNCTION__ ) );
 
   // std::cout << "3. Test = \n" << def.F << std::endl; 
   // std::cout << "3. Test = \n" << response.tau << std::endl; 
   // std::cout << "3. Test Uniaxial = \n" << tangent.dTau_dF << std::endl; 
}


void testHydrostaticResponse()
{

  std::array< double, 2 > materialProperties_ = { 3500, 1500 };
  const double            nMaterialProperties = 2;
  const int               elLabel             = 1;

  CompressibleNeoHooke mat = CompressibleNeoHooke( &materialProperties_[0], nMaterialProperties, elLabel );

  CompressibleNeoHooke::Deformation< 3 > def;
  CompressibleNeoHooke::TimeIncrement    timeInc = { 0, 0.1 };

  CompressibleNeoHooke::ConstitutiveResponse< 3 > response;
  CompressibleNeoHooke::AlgorithmicModuli< 3 >    tangent;

  def.F = Marmot::FastorStandardTensors::Spatial3D::I;
  def.F( 0, 0 ) += 0.2;
  def.F( 1, 1 ) += 0.2;
  def.F( 2, 2 ) += 0.2;

  mat.computeStress( response, tangent, def, timeInc );

/*  Tensor33d stressTarget;
  stressTarget( 0, 0 ) = ;
  stressTarget( 0, 1 ) = 0.;
  stressTarget( 0, 2 ) = 0.;
  stressTarget( 1, 0 ) = 0.;
  stressTarget( 1, 1 ) = ;
  stressTarget( 1, 2 ) = 0.;
  stressTarget( 2, 0 ) = 0.;
  stressTarget( 2, 1 ) = 0.;
  stressTarget( 2, 2 ) = ;

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      throwExceptionOnFailure(
        checkIfEqual( response.tau( i, j ), stressTarget( i, j ), 1e-10 ),
        "Uniaxial stress computation failed for CompressibleNeoHooke material in " + std::string( __PRETTY_FUNCTION__ ) );
        */
   std::cout << "4. Test = \n" << def.F << std::endl; 
   std::cout << "4. Test = \n" << response.tau << std::endl; 
  // std::cout << "4. Test = \n" << tangent.dTau_dF << std::endl;  
}


void testGeneralDeformationResponse()
{

  std::array< double, 2 > materialProperties_ = { 3500, 1500 };
  const double            nMaterialProperties = 2;
  const int               elLabel             = 1;

  CompressibleNeoHooke mat = CompressibleNeoHooke( &materialProperties_[0], nMaterialProperties, elLabel );

  CompressibleNeoHooke::Deformation< 3 > def;
  CompressibleNeoHooke::TimeIncrement    timeInc = { 0, 0.1 };

  CompressibleNeoHooke::ConstitutiveResponse< 3 > response;
  CompressibleNeoHooke::AlgorithmicModuli< 3 >    tangent;

  def.F = Marmot::FastorStandardTensors::Spatial3D::I;
  def.F( 0, 0 ) += 0.1;
  def.F( 1, 1 ) += 0.2;
  def.F( 2, 2 ) += 0.3;

  mat.computeStress( response, tangent, def, timeInc );

/*  Tensor33d stressTarget;
  stressTarget( 0, 0 ) = ;
  stressTarget( 0, 1 ) = ;
  stressTarget( 0, 2 ) = ;
  stressTarget( 1, 0 ) = ;
  stressTarget( 1, 1 ) = ;
  stressTarget( 1, 2 ) = ;
  stressTarget( 2, 0 ) = ;
  stressTarget( 2, 1 ) = ;
  stressTarget( 2, 2 ) = ;

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      throwExceptionOnFailure(
        checkIfEqual( response.tau( i, j ), stressTarget( i, j ), 1e-10 ),
        "Uniaxial stress computation failed for CompressibleNeoHooke material in " + std::string( __PRETTY_FUNCTION__ ) );
        */
   std::cout << "5. Test = \n" << def.F << std::endl; 
   std::cout << "5. Test = \n" << response.tau << std::endl; 
  // std::cout << "5. Test Uniaxial = \n" << tangent.dTau_dF << std::endl;   
}

void testAlgorithmicTangent()
{

  std::array< double, 2 > materialProperties_ = { 3500, 1500 };
  const double            nMaterialProperties = 2;
  const int               elLabel             = 1;

  CompressibleNeoHooke mat = CompressibleNeoHooke( &materialProperties_[0], nMaterialProperties, elLabel );

  CompressibleNeoHooke::Deformation< 3 > def;
  CompressibleNeoHooke::TimeIncrement    timeInc = { 0, 0.1 };

  CompressibleNeoHooke::ConstitutiveResponse< 3 > response;
  CompressibleNeoHooke::AlgorithmicModuli< 3 >    tangent;

  def.F = Marmot::FastorStandardTensors::Spatial3D::I;
  def.F( 0, 0 ) += 0.1;
  def.F( 1, 1 ) += 0.2;
  def.F( 2, 2 ) += 0.3;

  mat.computeStress( response, tangent, def, timeInc );

/*  Tensor33d stressTarget;
  stressTarget( 0, 0 ) = ;
  stressTarget( 0, 1 ) = ;
  stressTarget( 0, 2 ) = ;
  stressTarget( 1, 0 ) = ;
  stressTarget( 1, 1 ) = ;
  stressTarget( 1, 2 ) = ;
  stressTarget( 2, 0 ) = ;
  stressTarget( 2, 1 ) = ;
  stressTarget( 2, 2 ) = ;

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      throwExceptionOnFailure(
        checkIfEqual( response.tau( i, j ), stressTarget( i, j ), 1e-10 ),
        "Uniaxial stress computation failed for CompressibleNeoHooke material in " + std::string( __PRETTY_FUNCTION__ ) );
        */
  
}


int main()
{

  auto tests = std::vector< std::function< void() > >{ testUndeformedResponse,
                                                       testFiniteStrainUniaxialResponse,
                                                       testShearResponse,
                                                       testHydrostaticResponse,
                                                       testGeneralDeformationResponse,
                                                       testAlgorithmicTangent };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
