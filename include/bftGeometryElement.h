#pragma once 
#include "bftTypedefs.h"
#include "bftFiniteElement.h"

template< int nDim, int  nNodes>
class BftGeometryElement
{
    /* This is the Geometry Base element, which serves as a bases for all BftUels.
     * It corresponds to the GeometryElement in mpFEM,
     * although this as a static templated version.
     *
     * BftUels (corresponding do DofElements in mpFEM) can inherit from this element,
     * and access shape functions, derivatives and B Operator
     *
     * The element automatically determines its shape by the given  nDimension and number of nodes
     * */

    public:
        
        /*Typedefs*/
        static constexpr int    VoigtSize = (((nDim*nDim) +  nDim) / 2);

        typedef Matrix<double,  nDim, 1 >                                   XiSized;
        typedef Matrix<double,  nDim*nNodes, 1 >                            CoordinateVector;
        typedef Matrix<double,  nDim,  nDim>                                JacobianSized;

        typedef Matrix<double, 1, nNodes>                                   NSized;
        typedef Matrix<double, nDim, nNodes*nDim>                           NBSized;
        typedef Matrix<double, nDim, nNodes>                                dNdXiSized;
        typedef Matrix<double, VoigtSize, nNodes*nDim >                     BSized;


        /*Properties*/
        const Map<const CoordinateVector>                                   coordinates;
        const bft::FiniteElement::ElementShapes                             shape;

        /*Methods*/
        BftGeometryElement(const double* coords):
            coordinates(coords),
            shape ( bft::FiniteElement::getElementShapeByMetric(nDim, nNodes) ) 
            {};

        /*Please specialize these functions for each element individially
         *.cpp file. 
         *Fully specialized templates are precompiled in bftMechanics (rather than the unspecialized and partially specialized templates)
         * */
        NSized   N( const Ref< const XiSized>&   xi);
        dNdXiSized dNdXi( const Ref< const XiSized>&   xi);
        BSized   B( const Ref< const dNdXiSized>&  dNdX);

        /*These functions are equal for each element and independent of node number and  nDimension*/
        NBSized NB(const Ref<const NSized>&  N) {
            return bft::FiniteElement::NB<nDim, nNodes>(N);}
        JacobianSized Jacobian( const Ref< const dNdXiSized>& dNdXi) { 
            return bft::FiniteElement::Jacobian<nDim, nNodes> (dNdXi, coordinates); }
        dNdXiSized dNdX(const Ref<const dNdXiSized>& dNdXi, const Ref<const JacobianSized>& JacobianInverse) {
            return (dNdXi.transpose() * JacobianInverse).transpose();}

};
