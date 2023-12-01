#include "Marmot/CompressibleNeoHooke.h"
#include "Marmot/Marmot.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotMicromorphicTensorBasics.h"

int main()
{

  using namespace Marmot::Materials;
  std::array< double, 2 > materialProperties_ = { 3500, 1500 };
  const double            nMaterialProperties = 2;
  const int               elLabel             = 1;

  CompressibleNeoHooke mat = CompressibleNeoHooke( &materialProperties_[0], nMaterialProperties, elLabel );

  CompressibleNeoHooke::Deformation< 3 > def;
  CompressibleNeoHooke::TimeIncrement    timeInc = { 0, 0.1 };

  CompressibleNeoHooke::ConstitutiveResponse< 3 > response;
  CompressibleNeoHooke::AlgorithmicModuli< 3 >    tangent;

  def.F = Marmot::FastorStandardTensors::Spatial3D::I;
  def.F( 0, 0 ) += 1e-4;
  def.F( 1, 1 ) += 1e-4;
  def.F( 2, 2 ) += 1e-4;

  mat.computeStress( response, tangent, def, timeInc );

  std::cout << "deformation=\n" << def.F << std::endl;
  std::cout << "stress = \n" << response.tau << std::endl;
  std::cout << "tangent = \n" << tangent.dTau_dF << std::endl;

  return 0;
}
