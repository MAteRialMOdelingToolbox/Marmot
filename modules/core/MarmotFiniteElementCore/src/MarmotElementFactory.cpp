#include "Marmot/MarmotElementFactory.h"

using namespace Marmot::Factory;

MarmotElementFactory::ElementFactoryMap& MarmotElementFactory::elementFactoryFunctionByName()
{
  static ElementFactoryMap map;
  return map;
}
