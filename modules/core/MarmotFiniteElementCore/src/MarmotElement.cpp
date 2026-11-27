#include "Marmot/MarmotElement.h"

Marmot::Discretizations::Elements::Element::~Element() {}

Marmot::Utilities::StateView Marmot::Discretizations::Elements::Element::getStateView( const std::string& stateName,
                                                                                       int quadraturePoint )
{
  const int nQuadraturePoints = getNumberOfQuadraturePoints();
  if ( quadraturePoint < 0 || quadraturePoint >= nQuadraturePoints )
    throw std::invalid_argument( Marmot::Utilities::MakeString()
                                 << "Element::getStateView: state '" << stateName << "' asked for quadrature point "
                                 << quadraturePoint << " of a " << getElementShape() << " element, which has "
                                 << nQuadraturePoints << " (0 to " << nQuadraturePoints - 1 << ")" );

  return getStateViewAtQuadraturePoint( stateName, quadraturePoint );
}

void Marmot::Discretizations::Elements::Element::assignProperty(
  const Marmot::Discretizations::ElementProperties& property )
{
}

void Marmot::Discretizations::Elements::Element::assignProperty(
  const Marmot::Discretizations::MaterialSection& property )
{
}
