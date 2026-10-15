#include "Marmot/MarmotElementFactory.h"

namespace Marmot::Registry {

ElementFactory::ElementFactoryMap& ElementFactory::elementFactoryFunctionByName()
{
  static ElementFactoryMap map;
  return map;
}

} // namespace Marmot::Registry
