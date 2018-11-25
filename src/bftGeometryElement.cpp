#include "bftGeometryElement.h"

/* Template Specialization for Shape functions.
 * These (full) specializations are precompiled into the current static bftMechanics.
 * Thus, they act as static links to the corresponding FiniteElement functions.
 *
 * For a new element, please specialize ALL functions.
 *
 * */

// Truss2
template <>
BftGeometryElement<1, 2>::NSized BftGeometryElement<1, 2>::N( const XiSized& xi )
{
    return bft::FiniteElement::Spatial1D::Truss2::N( xi( 0 ) );
}
// Quad4
template <>
BftGeometryElement<2, 4>::NSized BftGeometryElement<2, 4>::N( const XiSized& xi )
{
    return bft::FiniteElement::Spatial2D::Quad4::N( xi );
}
// Quad8
template <>
BftGeometryElement<2, 8>::NSized BftGeometryElement<2, 8>::N( const XiSized& xi )
{
    return bft::FiniteElement::Spatial2D::Quad8::N( xi );
}
// Hexa8
template <>
BftGeometryElement<3, 8>::NSized BftGeometryElement<3, 8>::N( const XiSized& xi )
{
    return bft::FiniteElement::Spatial3D::Hexa8::N( xi );
}
// Hexa20
template <>
BftGeometryElement<3, 20>::NSized BftGeometryElement<3, 20>::N( const XiSized& xi )
{
    return bft::FiniteElement::Spatial3D::Hexa20::N( xi );
}

/* Template Specialization for Shape functions derivatives.
 * */
// Truss2
template <> BftGeometryElement<1, 2>::dNdXiSized BftGeometryElement<1, 2>::dNdXi( const XiSized& xi )
{
    return bft::FiniteElement::Spatial1D::Truss2::dNdXi( xi( 0 ) );
}
// Quad4
template <>
BftGeometryElement<2, 4>::dNdXiSized BftGeometryElement<2, 4>::dNdXi( const XiSized& xi )
{
    return bft::FiniteElement::Spatial2D::Quad4::dNdXi( xi );
}
// Quad8
template <>
BftGeometryElement<2, 8>::dNdXiSized BftGeometryElement<2, 8>::dNdXi( const XiSized& xi )
{
    return bft::FiniteElement::Spatial2D::Quad8::dNdXi( xi );
}
// Hexa8
template <>
BftGeometryElement<3, 8>::dNdXiSized BftGeometryElement<3, 8>::dNdXi( const XiSized& xi )
{
    return bft::FiniteElement::Spatial3D::Hexa8::dNdXi( xi );
}
// Hexa20
template <>
BftGeometryElement<3, 20>::dNdXiSized BftGeometryElement<3, 20>::dNdXi( const XiSized& xi )
{
    return bft::FiniteElement::Spatial3D::Hexa20::dNdXi( xi );
}

/* Template Specialization for B Operator
 * Note that partial template specialization (i.e., for all nDim=2 and all nDim=3) is not valid for
 * methods (only for classes) in C++, and, hence, the B Operator must be specialized for each
 * nDim/nNodes combination individually.
 * */

/* TODO: Replace by constexpr if expression for each dimension -> get rid of partial specializations; full specialization should be possible
 *
 * */

// Truss2
template <>
typename BftGeometryElement<1, 2>::BSized BftGeometryElement<1, 2>::B( const dNdXiSized& dNdX )
{
    return dNdX;
}
// Quad4
template <>
typename BftGeometryElement<2, 4>::BSized BftGeometryElement<2, 4>::B( const dNdXiSized& dNdX )
{
    return bft::FiniteElement::Spatial2D::B<4>( dNdX );
}
// Quad8
template <>
typename BftGeometryElement<2, 8>::BSized BftGeometryElement<2, 8>::B( const dNdXiSized& dNdX )
{
    return bft::FiniteElement::Spatial2D::B<8>( dNdX );
}
// Hexa8
template <>
typename BftGeometryElement<3, 8>::BSized BftGeometryElement<3, 8>::B( const dNdXiSized& dNdX )
{
    return bft::FiniteElement::Spatial3D::B<8>( dNdX );
}
// Hexa20
template <>
typename BftGeometryElement<3, 20>::BSized BftGeometryElement<3, 20>::B( const dNdXiSized& dNdX )
{
    return bft::FiniteElement::Spatial3D::B<20>( dNdX );
}

// Truss2
template <>
typename BftGeometryElement<1, 2>::BSized BftGeometryElement<1, 2>::BGreen( const dNdXiSized& dNdX, const JacobianSized& F)
{
    return dNdX * F(0,0);
}
// Quad4
template <>
typename BftGeometryElement<2, 4>::BSized BftGeometryElement<2, 4>::BGreen( const dNdXiSized& dNdX,  const JacobianSized& F)
{
    return bft::FiniteElement::Spatial2D::BGreen<4>( dNdX, F);
}
// Quad8
template <>
typename BftGeometryElement<2, 8>::BSized BftGeometryElement<2, 8>::BGreen( const dNdXiSized& dNdX,  const JacobianSized& F)
{
    return bft::FiniteElement::Spatial2D::BGreen<8>( dNdX, F);
}
// Hexa8
template <>
typename BftGeometryElement<3, 8>::BSized BftGeometryElement<3, 8>::BGreen( const dNdXiSized& dNdX,  const JacobianSized& F)
{
    return bft::FiniteElement::Spatial3D::BGreen<8>( dNdX, F);
}
// Hexa20
template <>
typename BftGeometryElement<3, 20>::BSized BftGeometryElement<3, 20>::BGreen( const dNdXiSized& dNdX, const JacobianSized& F)
{
    return bft::FiniteElement::Spatial3D::BGreen<20>( dNdX, F);
}
