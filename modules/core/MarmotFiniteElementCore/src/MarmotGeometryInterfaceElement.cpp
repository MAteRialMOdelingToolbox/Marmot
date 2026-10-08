#include "Marmot/MarmotGeometryInterfaceElement.h"

using namespace Eigen;

/* -------------------------------------------------------------------------- */
/* Shape functions                                                            */
/* -------------------------------------------------------------------------- */

// ILine2:
// 2-node line interface embedded in 2D,
// with 4 total interface element nodes.
/**
 * @brief Shape functions of the two-node line interface at a parametric coordinate.
 *
 * @param[in] xi parametric coordinate on the line
 * @return the two shape functions of the line
 */
template <>
MarmotGeometryInterfaceElement< 2, 4 >::NSized MarmotGeometryInterfaceElement< 2, 4 >::N( const XiSized& xi ) const
{
  return Marmot::FiniteElement::Spatial1D::Bar2::N( xi( 0 ) );
}

// IQuad4:
// 4-node quadrilateral interface surface embedded in 3D,
// with 8 total interface element nodes.
/**
 * @brief Shape functions of the four-node quadrilateral interface surface at a parametric point.
 *
 * @param[in] xi parametric coordinates on the surface
 * @return the four shape functions of the quadrilateral
 */
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
/**
 * @brief Derivatives of the shape functions of the two-node line interface.
 *
 * @param[in] xi parametric coordinate on the line
 * @return the derivatives of the two shape functions with respect to the line coordinate
 */
template <>
MarmotGeometryInterfaceElement< 2, 4 >::dNdXiSized MarmotGeometryInterfaceElement< 2, 4 >::dNdXi(
  const XiSized& xi ) const
{
  return Marmot::FiniteElement::Spatial1D::Bar2::dNdXi( xi( 0 ) );
}

// IQuad4:
// derivatives wrt surface coordinates xi, eta.
/**
 * @brief Derivatives of the shape functions of the four-node quadrilateral interface surface.
 *
 * @param[in] xi parametric coordinates on the surface
 * @return the derivatives of the four shape functions with respect to the surface coordinates
 */
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