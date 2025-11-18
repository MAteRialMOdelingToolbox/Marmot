#include "Marmot/MarmotElementFactory.h"

using namespace MarmotLibrary;

MarmotElementFactory::ElementFactoryMap& MarmotElementFactory::elementFactoryFunctionByName()
{
  static ElementFactoryMap map;
  return map;
}
