#pragma once
#include "bftFiniteElement.h"
#include "bftTypedefs.h"
#include <map>
#include <iostream>

template <int nDim, int nNodes>
class BftGeometryElement {
    /* This is the Geometry Base element, which serves as a base for all BftElements.
     * It corresponds to the GeometryElement in mpFEM,
     * although this as a static templated version.
     *
     * BftElements (corresponding do DofElements in mpFEM) can inherit from this element,
     * and access shape functions, derivatives and B Operator
     *
     * The element automatically determines its shape by the given  nDimension and number of nodes
     * */

  public:
    /*Typedefs*/
    static constexpr int VoigtSize = ( ( ( nDim * nDim ) + nDim ) / 2 );

    typedef Eigen::Matrix<double, nDim, 1>                  XiSized;
    typedef Eigen::Matrix<double, nDim * nNodes, 1>         CoordinateVector;
    typedef Eigen::Matrix<double, nDim, nDim>               JacobianSized;
    typedef Eigen::Matrix<double, 1, nNodes>                NSized;
    typedef Eigen::Matrix<double, nDim, nNodes * nDim>      NBSized;
    typedef Eigen::Matrix<double, nDim, nNodes>             dNdXiSized;
    typedef Eigen::Matrix<double, VoigtSize, nNodes * nDim> BSized;

    /*Properties*/
    Eigen::Map<const CoordinateVector>      coordinates;
    const bft::FiniteElement::ElementShapes shape;

    /*Methods*/
    BftGeometryElement()
        : coordinates( nullptr ), shape( bft::FiniteElement::getElementShapeByMetric( nDim, nNodes ) ){};

    std::string getElementShape()
    {
        using namespace bft::FiniteElement;
        static std::map<ElementShapes, std::string> shapes = {{Bar2, "bar2"},
                                                              {Quad4, "quad4"},
                                                              {Quad8, "quad8"},
                                                              {Tetra4, "tetra4"},
                                                              {Hexa8, "hexa8"},
                                                              {Hexa20, "hexa20"}};

        return shapes[this->shape];
    }

    void initializeYourself( const double* coords )
    {
        new ( &coordinates ) Eigen::Map<const CoordinateVector>( coords );
    }

    /*Please specialize these functions for each element individially
     *.cpp file.
     *Fully specialized templates are precompiled in bftMechanics (rather than the unspecialized and
     *partially specialized templates)
     * */
    NSized     N( const XiSized& xi );
    dNdXiSized dNdXi( const XiSized& xi );
    BSized     B( const dNdXiSized& dNdX );
    BSized     BGreen( const dNdXiSized& dNdX, const JacobianSized& F );

    /*These functions are equal for each element and independent of node number and  nDimension*/
    NBSized NB( const NSized& N ) { return bft::FiniteElement::NB<nDim, nNodes>( N ); }

    JacobianSized Jacobian( const dNdXiSized& dNdXi )
    {
        return bft::FiniteElement::Jacobian<nDim, nNodes>( dNdXi, coordinates );
    }

    dNdXiSized dNdX( const dNdXiSized& dNdXi, const JacobianSized& JacobianInverse )
    {
        return ( dNdXi.transpose() * JacobianInverse ).transpose();
    }

    JacobianSized F( const dNdXiSized& dNdX, const CoordinateVector& Q )
    {
        return bft::FiniteElement::Jacobian<nDim, nNodes>( dNdX, Q ) + JacobianSized::Identity();
    }
};
