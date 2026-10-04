#include "Marmot/MarmotGeometryInterfaceElement.h"

using namespace Eigen;

/* -------------------------------------------------------------------------- */
/* Shape functions                                                            */
/* -------------------------------------------------------------------------- */

// ILine2:
// 2-node line interface embedded in 2D,
// with 4 total interface element nodes.
template <>
MarmotGeometryInterfaceElement< 2, 4 >::NSized MarmotGeometryInterfaceElement< 2, 4 >::N( const XiSized& xi ) const
{
  return Marmot::FiniteElement::Spatial1D::Bar2::N( xi( 0 ) );
}

// IQuad4:
// 4-node quadrilateral interface surface embedded in 3D,
// with 8 total interface element nodes.
template <>
MarmotGeometryInterfaceElement< 3, 8 >::NSized MarmotGeometryInterfaceElement< 3, 8 >::N( const XiSized& xi ) const
{
  return Marmot::FiniteElement::Spatial2D::Quad4::N( xi );
}

/* -------------------------------------------------------------------------- */
/* Shape-function derivatives wrt interface coordinates                       */
/* -------------------------------------------------------------------------- */

// ILine2:
// derivative wrt line coordinate xi.
template <>
MarmotGeometryInterfaceElement< 2, 4 >::dNdXiSized MarmotGeometryInterfaceElement< 2, 4 >::dNdXi(
  const XiSized& xi ) const
{
  return Marmot::FiniteElement::Spatial1D::Bar2::dNdXi( xi( 0 ) );
}

// IQuad4:
// derivatives wrt surface coordinates xi, eta.
template <>
MarmotGeometryInterfaceElement< 3, 8 >::dNdXiSized MarmotGeometryInterfaceElement< 3, 8 >::dNdXi(
  const XiSized& xi ) const
{
  return Marmot::FiniteElement::Spatial2D::Quad4::dNdXi( xi );
}

/* -------------------------------------------------------------------------- */
/* Explicit template instantiation                                            */
/* -------------------------------------------------------------------------- */

template class MarmotGeometryInterfaceElement< 2, 4 >;
template class MarmotGeometryInterfaceElement< 3, 8 >;