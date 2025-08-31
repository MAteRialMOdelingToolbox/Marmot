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





// Test 2.
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
        "Test 2: Finite Strain Uniaxial deformation load case computation failed for CompressibleNeoHooke material in " + std::string( __PRETTY_FUNCTION__ ) );
  
  // std::cout << "2. Test = \n" << def.F << std::endl; 
  // std::cout << "2. Test = \n" << response.tau << std::endl; 
  // std::cout << "2. Test Uniaxial = \n" << tangent.dTau_dF << std::endl; 
}


// Test 3.
void testFiniteStrainShearResponse()
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
  stressTarget( 0, 0 ) = -20;
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
        "Test 3. Finite strain simple shear load case failed for CompressibleNeoHooke material in " + std::string( __PRETTY_FUNCTION__ ) );
 
    //std::cout << "3. Test = \n" << def.F << std::endl; 
    //std::cout << "3. Test = \n" << response.tau << std::endl; 
    //std::cout << "3. Test = \n" << tangent.dTau_dF << std::endl; 
}

// Test 4.
void testSmallStrainShearResponse()
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
  def.F( 1, 0 ) += 1e-06;

  mat.computeStress( response, tangent, def, timeInc );

  Tensor33d stressTarget;
  stressTarget( 0, 0 ) = -4.99994712299667e-10;
  stressTarget( 0, 1 ) = 0.0015;
  stressTarget( 0, 2 ) = 0.;
  stressTarget( 1, 0 ) = 0.0015;
  stressTarget( 1, 1 ) = 9.99793777126387e-10;
  stressTarget( 1, 2 ) = 0.;
  stressTarget( 2, 0 ) = 0.;
  stressTarget( 2, 1 ) = 0.;
  stressTarget( 2, 2 ) = -4.99994712299667e-10;

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      throwExceptionOnFailure(
        checkIfEqual( response.tau( i, j ), stressTarget( i, j ), 1e-10 ),
        "Test 4. Small strain simple shear failed for CompressibleNeoHooke material in " + std::string( __PRETTY_FUNCTION__ ) );
 
    // std::cout << "4. Test = \n" << def.F << std::endl; 
    // std::cout << "4. Test = \n" << response.tau << std::endl; 
    // std::cout << "4. Test = \n" << tangent.dTau_dF << std::endl; 
}

// Test 5.
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
  def.F( 0, 0 ) += 0.02;
  def.F( 1, 1 ) += 0.02;
  def.F( 2, 2 ) += 0.02;

  mat.computeStress( response, tangent, def, timeInc );

  Tensor33d stressTarget;
  stressTarget( 0, 0 ) = 208.417157443082;
  stressTarget( 0, 1 ) = 0.;
  stressTarget( 0, 2 ) = 0.;
  stressTarget( 1, 0 ) = 0.;
  stressTarget( 1, 1 ) = 208.417157443082;
  stressTarget( 1, 2 ) = 0.;
  stressTarget( 2, 0 ) = 0.;
  stressTarget( 2, 1 ) = 0.;
  stressTarget( 2, 2 ) = 208.417157443082;

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      throwExceptionOnFailure(
        checkIfEqual( response.tau( i, j ), stressTarget( i, j ), 1e-10 ),
        "Test 5. Hydrostatic load case computation failed for CompressibleNeoHooke material in " + std::string( __PRETTY_FUNCTION__ ) );
    
  // std::cout << "5. Test = \n" << def.F << std::endl; 
  // std::cout << "5. Test = \n" << response.tau << std::endl; 
  // std::cout << "5. Test = \n" << tangent.dTau_dF << std::endl;  
}

// Test 6.
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
  def.F( 0, 0 ) += 0.01;
  def.F( 1, 1 ) += 0.02;
  def.F( 2, 2 ) += 0.03;

  mat.computeStress( response, tangent, def, timeInc );

  Tensor33d stressTarget;
  stressTarget( 0, 0 ) = 178.712770583994;
  stressTarget( 0, 1 ) = 0.;
  stressTarget( 0, 2 ) = 0.;
  stressTarget( 1, 0 ) = 0.;
  stressTarget( 1, 1 ) = 207.982235529133;
  stressTarget( 1, 2 ) = 0.;
  stressTarget( 2, 0 ) = 0.;
  stressTarget( 2, 1 ) = 0.;
  stressTarget( 2, 2 ) = 237.540069587033;

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      throwExceptionOnFailure(
        checkIfEqual( response.tau( i, j ), stressTarget( i, j ), 1e-10 ),
        "Test 6: General deformation response computation failed for CompressibleNeoHooke material in " + std::string( __PRETTY_FUNCTION__ ) );
        
  //std::cout << "6. Test = \n" << def.F << std::endl; 
  //std::cout << "6. Test = \n" << response.tau << std::endl; 
  // std::cout << "6. Test Uniaxial = \n" << tangent.dTau_dF << std::endl;   
}

// Test 7.
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
  def.F( 0, 0 ) += 0.01;
  def.F( 1, 1 ) += 0.02;
  def.F( 2, 2 ) += 0.03;

  mat.computeStress( response, tangent, def, timeInc );

Tensor3333d tangentTarget(0.0);

tangentTarget(0,0,0,0) = 5450.8251046444;
tangentTarget(0,0,1,1) = 2494.28143117259;
tangentTarget(0,0,2,2) = 2450.93382241823;

tangentTarget(0,1,0,1) = 1470.68247507597;
tangentTarget(0,1,1,0) = 1456.26401943797;

tangentTarget(0,2,0,2) = 1485.10093071397;
tangentTarget(0,2,2,0) = 1456.26401943797;

tangentTarget(1,0,0,1) = 1470.68247507597;
tangentTarget(1,0,1,0) = 1456.26401943797;

tangentTarget(1,1,0,0) = 2518.97728692678;
tangentTarget(1,1,1,1) = 5416.51601207936;
tangentTarget(1,1,2,2) = 2431.98918491329;

tangentTarget(1,2,1,2) = 1485.10093071397;
tangentTarget(1,2,2,1) = 1470.68247507597;

tangentTarget(2,0,0,2) = 1485.10093071397;
tangentTarget(2,0,2,0) = 1456.26401943797;

tangentTarget(2,1,1,2) = 1485.10093071397;
tangentTarget(2,1,2,1) = 1470.68247507597;

tangentTarget(2,2,0,0) = 2499.46716543642;
tangentTarget(2,2,1,1) = 2455.83221613793;
tangentTarget(2,2,2,2) = 5383.05976216137;


  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
        for (int k = 0; k < 3; k++)
            for (int l = 0; l < 3; l++)
      throwExceptionOnFailure(
        checkIfEqual( tangent.dTau_dF( i, j, k, l ), tangentTarget( i, j, k, l ), 1e-10 ),
        "Test 7: Algorithmic tangent computation failed for CompressibleNeoHooke material in " + std::string( __PRETTY_FUNCTION__ ) );
        

 // std::cout << "7. Test = \n" << tangent.dTau_dF << std::endl;
  
}

// Test 8.
void testPureRotation()
{

  std::array< double, 2 > materialProperties_ = { 3500, 1500 };
  const double            nMaterialProperties = 2;
  const int               elLabel             = 1;

  CompressibleNeoHooke mat = CompressibleNeoHooke( &materialProperties_[0], nMaterialProperties, elLabel );

  CompressibleNeoHooke::Deformation< 3 > def;
  CompressibleNeoHooke::TimeIncrement    timeInc = { 0, 0.1 };

  CompressibleNeoHooke::ConstitutiveResponse< 3 > response;
  CompressibleNeoHooke::AlgorithmicModuli< 3 >    tangent;
 
  double phi = std::numbers::pi / 6; // 30 degrees

  def.F( 0, 0 ) = cos( phi );
  def.F( 0, 1 ) = -sin( phi );
  def.F( 0, 2 ) = 0;
  def.F( 1, 0 ) = sin( phi );
  def.F( 1, 1 ) = cos( phi );
  def.F( 1, 2 ) = 0;
  def.F( 2, 0 ) = 0;
  def.F( 2, 1 ) = 0;
  def.F( 2, 2 ) = 1;

  mat.computeStress( response, tangent, def, timeInc );

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

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      throwExceptionOnFailure(
        checkIfEqual( response.tau( i, j ), stressTarget( i, j ), 1e-10 ),
        "Test 7: Algorithmic tangent computation failed for CompressibleNeoHooke material in " + std::string( __PRETTY_FUNCTION__ ) );
        
  //std::cout << "8. Test Pure Rotation deformation = \n" << def.F << std::endl;
  //std::cout << "8. Test Pure Rotation stress = \n" << response.tau << std::endl;
  
}


int main()
{

  auto tests = std::vector< std::function< void() > >{ testUndeformedResponse,
                                                       testFiniteStrainUniaxialResponse,
                                                       testFiniteStrainShearResponse,
                                                       testSmallStrainShearResponse,
                                                       testHydrostaticResponse,
                                                       testGeneralDeformationResponse,
                                                       testAlgorithmicTangent,
                                                        testPureRotation };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
