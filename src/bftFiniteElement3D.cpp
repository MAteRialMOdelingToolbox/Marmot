#include "bftFiniteElement.h"
#include "bftTypedefs.h"
#include <iostream>

namespace bft{

    //****************************************************
    namespace FiniteElement
    {
        namespace Spatial3D
        {
            namespace Hexa8{
                NSized N(const Ref<const Vector3d>& xi)
                { /*
                   *        (8)________(7)   
                   *          /|      /|   
                   *      (5)/_|_____(6)    
                   *        |  |     | |
                   *        |(4)_____|_|(3)
                   *        | /      | /
                   *        |/_______|/
                   *     (1)        (2) 
                   *
                   *
                   *     x3  x2
                   *     | /
                   *     |/___x1 
                   *
                   *  */

                    NSized N_;

                    N_ << (1-xi(0)) * (1-xi(1)) * (1-xi(2)),
                       (1+xi(0)) * (1-xi(1)) * (1-xi(2)),
                       (1+xi(0)) * (1+xi(1)) * (1-xi(2)),
                       (1-xi(0)) * (1+xi(1)) * (1-xi(2)),
                       (1-xi(0)) * (1-xi(1)) * (1+xi(2)),
                       (1+xi(0)) * (1-xi(1)) * (1+xi(2)),
                       (1+xi(0)) * (1+xi(1)) * (1+xi(2)),
                       (1-xi(0)) * (1+xi(1)) * (1+xi(2));

                    N_ *= 1./8;
                    return N_;

                }
                dNdXiSized dNdXi(const Ref<const Vector3d>& xi)
                {
                    dNdXiSized dNdXi_;

                    dNdXi_ <<
                        -1 * (1-xi(1)) * (1-xi(2)),
                        +1 * (1-xi(1)) * (1-xi(2)),
                        +1 * (1+xi(1)) * (1-xi(2)),
                        -1 * (1+xi(1)) * (1-xi(2)),
                        -1 * (1-xi(1)) * (1+xi(2)),
                        +1 * (1-xi(1)) * (1+xi(2)),
                        +1 * (1+xi(1)) * (1+xi(2)),
                        -1 * (1+xi(1)) * (1+xi(2)),

                        (1-xi(0)) * -1 * (1-xi(2)),
                        (1+xi(0)) * -1 * (1-xi(2)),
                        (1+xi(0)) * +1 * (1-xi(2)),
                        (1-xi(0)) * +1 * (1-xi(2)),
                        (1-xi(0)) * -1 * (1+xi(2)),
                        (1+xi(0)) * -1 * (1+xi(2)),
                        (1+xi(0)) * +1 * (1+xi(2)),
                        (1-xi(0)) * +1 * (1+xi(2)),

                        (1-xi(0)) * (1-xi(1)) * -1,
                        (1+xi(0)) * (1-xi(1)) * -1,
                        (1+xi(0)) * (1+xi(1)) * -1,
                        (1-xi(0)) * (1+xi(1)) * -1,
                        (1-xi(0)) * (1-xi(1)) * +1,
                        (1+xi(0)) * (1-xi(1)) * +1,
                        (1+xi(0)) * (1+xi(1)) * +1,
                        (1-xi(0)) * (1+xi(1)) * +1;

                    dNdXi_ *= 1./8;

                    return dNdXi_;
                }

                Matrix3d Jacobian(const Ref<const dNdXiSized >& dNdXi, const Ref<const Matrix<double, nNodes*nDim, 1>>& coordinates)
                {
                    // convenience wrapper to the templated version of the Jacobian
                    return FiniteElement::Jacobian<nDim, nNodes> ( dNdXi, coordinates );
                }

                BSized B(const Ref<const  dNdXiSized>& dNdX)
                {
                    // convenience wrapper to the templated version of the Spatial2D B Operator
                    return FiniteElement::Spatial3D::B<nNodes> ( dNdX ) ;
                }

            } // end of namespace Hexa8 

            namespace Hexa20{
                NSized N(const Ref<const Vector3d>& xi)
                { /*
                   *           (8)_____(15)_____(7)   
                   *            /|             /|   
                   *       (16)/ |            /(14)
                   *          /  |           /  |
                   *      (5)/___|__(13)___(6)  |
                   *        |  (18)         | (17)
                   *        |    |          |   |
                   *        |  (4)___(11)___|___|(3)
                   *      (15)  /         (16)  /
                   *        |  /(12)        |  /(10)
                   *        | /             | /
                   *        |/______________|/
                   *     (1)       (9)      (2) 
                   *
                   *
                   *     x3  x2
                   *     | /
                   *     |/___x1 
                   *
                   *  */

                    NSized N_;
                    // Taken from mpFEM - Peter Gamnitzer
                    N_ << 
                    -0.125*(1-xi(0))*(1-xi(1))*(1-xi(2))*(2+xi(0)+xi(1)+xi(2)),
                    -0.125*(1+xi(0))*(1-xi(1))*(1-xi(2))*(2-xi(0)+xi(1)+xi(2)),
                    -0.125*(1+xi(0))*(1+xi(1))*(1-xi(2))*(2-xi(0)-xi(1)+xi(2)),
                    -0.125*(1-xi(0))*(1+xi(1))*(1-xi(2))*(2+xi(0)-xi(1)+xi(2)),
                    -0.125*(1-xi(0))*(1-xi(1))*(1+xi(2))*(2+xi(0)+xi(1)-xi(2)),
                    -0.125*(1+xi(0))*(1-xi(1))*(1+xi(2))*(2-xi(0)+xi(1)-xi(2)),
                    -0.125*(1+xi(0))*(1+xi(1))*(1+xi(2))*(2-xi(0)-xi(1)-xi(2)),
                    -0.125*(1-xi(0))*(1+xi(1))*(1+xi(2))*(2+xi(0)-xi(1)-xi(2)),
                    0.25*(1-xi(0)*xi(0))*(1-xi(1))*(1-xi(2)),
                    0.25*(1-xi(1)*xi(1))*(1+xi(0))*(1-xi(2)),
                    0.25*(1-xi(0)*xi(0))*(1+xi(1))*(1-xi(2)),
                    0.25*(1-xi(1)*xi(1))*(1-xi(0))*(1-xi(2)),
                    0.25*(1-xi(0)*xi(0))*(1-xi(1))*(1+xi(2)),
                    0.25*(1-xi(1)*xi(1))*(1+xi(0))*(1+xi(2)),
                    0.25*(1-xi(0)*xi(0))*(1+xi(1))*(1+xi(2)),
                    0.25*(1-xi(1)*xi(1))*(1-xi(0))*(1+xi(2)),
                    0.25*(1-xi(0))*(1-xi(1))*(1-xi(2)*xi(2)),
                    0.25*(1+xi(0))*(1-xi(1))*(1-xi(2)*xi(2)),
                    0.25*(1+xi(0))*(1+xi(1))*(1-xi(2)*xi(2)),
                    0.25*(1-xi(0))*(1+xi(1))*(1-xi(2)*xi(2));
                    return N_;

                }
                dNdXiSized dNdXi(const Ref<const Vector3d>& xi)
                {
                    dNdXiSized dNdXi_;

                    dNdXi_ <<
                   0.125*(1-xi[1])*(1-xi[2])*(1+2*xi[0]+xi[1]+xi[2]),
                  -0.125*(1-xi[1])*(1-xi[2])*(1-2*xi[0]+xi[1]+xi[2]),
                  -0.125*(1+xi[1])*(1-xi[2])*(1-2*xi[0]-xi[1]+xi[2]),
                   0.125*(1+xi[1])*(1-xi[2])*(1+2*xi[0]-xi[1]+xi[2]),
                   0.125*(1-xi[1])*(1+xi[2])*(1+2*xi[0]+xi[1]-xi[2]),
                  -0.125*(1-xi[1])*(1+xi[2])*(1-2*xi[0]+xi[1]-xi[2]),
                  -0.125*(1+xi[1])*(1+xi[2])*(1-2*xi[0]-xi[1]-xi[2]),
                   0.125*(1+xi[1])*(1+xi[2])*(1+2*xi[0]-xi[1]-xi[2]),
                  -0.5*xi[0]*(1-xi[1])*(1-xi[2])                    ,
                   0.25*(1-xi[1]*xi[1])*(1-xi[2])                   ,
                  -0.5*xi[0]*(1+xi[1])*(1-xi[2])                    ,
                  -0.25*(1-xi[1]*xi[1])*(1-xi[2])                   ,
                  -0.5*xi[0]*(1-xi[1])*(1+xi[2])                    ,
                   0.25*(1-xi[1]*xi[1])*(1+xi[2])                   ,
                  -0.5*xi[0]*(1+xi[1])*(1+xi[2])                    ,
                  -0.25*(1-xi[1]*xi[1])*(1+xi[2])                   ,
                  -0.25*(1-xi[1])*(1-xi[2]*xi[2])                   ,
                   0.25*(1-xi[1])*(1-xi[2]*xi[2])                   ,
                   0.25*(1+xi[1])*(1-xi[2]*xi[2])                   ,
                  -0.25*(1+xi[1])*(1-xi[2]*xi[2])                   ,
                   0.125*(1-xi[0])*(1-xi[2])*(1+xi[0]+2*xi[1]+xi[2]),
                   0.125*(1+xi[0])*(1-xi[2])*(1-xi[0]+2*xi[1]+xi[2]),
                  -0.125*(1+xi[0])*(1-xi[2])*(1-xi[0]-2*xi[1]+xi[2]),
                  -0.125*(1-xi[0])*(1-xi[2])*(1+xi[0]-2*xi[1]+xi[2]),
                   0.125*(1-xi[0])*(1+xi[2])*(1+xi[0]+2*xi[1]-xi[2]),
                   0.125*(1+xi[0])*(1+xi[2])*(1-xi[0]+2*xi[1]-xi[2]),
                  -0.125*(1+xi[0])*(1+xi[2])*(1-xi[0]-2*xi[1]-xi[2]),
                  -0.125*(1-xi[0])*(1+xi[2])*(1+xi[0]-2*xi[1]-xi[2]),
                  -0.25*(1-xi[0]*xi[0])*(1-xi[2])                   ,
                  -0.5*xi[1]*(1+xi[0])*(1-xi[2])                    ,
                   0.25*(1-xi[0]*xi[0])*(1-xi[2])                   ,
                  -0.5*xi[1]*(1-xi[0])*(1-xi[2])                    ,
                  -0.25*(1-xi[0]*xi[0])*(1+xi[2])                   ,
                  -0.5*xi[1]*(1+xi[0])*(1+xi[2])                    ,
                   0.25*(1-xi[0]*xi[0])*(1+xi[2])                   ,
                  -0.5*xi[1]*(1-xi[0])*(1+xi[2])                    ,
                  -0.25*(1-xi[0])*(1-xi[2]*xi[2])                   ,
                  -0.25*(1+xi[0])*(1-xi[2]*xi[2])                   ,
                   0.25*(1+xi[0])*(1-xi[2]*xi[2])                   ,
                   0.25*(1-xi[0])*(1-xi[2]*xi[2])                   ,
                   0.125*(1-xi[0])*(1-xi[1])*(1+xi[0]+xi[1]+2*xi[2]),
                   0.125*(1+xi[0])*(1-xi[1])*(1-xi[0]+xi[1]+2*xi[2]),
                   0.125*(1+xi[0])*(1+xi[1])*(1-xi[0]-xi[1]+2*xi[2]),
                   0.125*(1-xi[0])*(1+xi[1])*(1+xi[0]-xi[1]+2*xi[2]),
                  -0.125*(1-xi[0])*(1-xi[1])*(1+xi[0]+xi[1]-2*xi[2]),
                  -0.125*(1+xi[0])*(1-xi[1])*(1-xi[0]+xi[1]-2*xi[2]),
                  -0.125*(1+xi[0])*(1+xi[1])*(1-xi[0]-xi[1]-2*xi[2]),
                  -0.125*(1-xi[0])*(1+xi[1])*(1+xi[0]-xi[1]-2*xi[2]),
                  -0.25*(1-xi[0]*xi[0])*(1-xi[1])                   ,
                  -0.25*(1-xi[1]*xi[1])*(1+xi[0])                   ,
                  -0.25*(1-xi[0]*xi[0])*(1+xi[1])                   ,
                  -0.25*(1-xi[1]*xi[1])*(1-xi[0])                   ,
                   0.25*(1-xi[0]*xi[0])*(1-xi[1])                   ,
                   0.25*(1-xi[1]*xi[1])*(1+xi[0])                   ,
                   0.25*(1-xi[0]*xi[0])*(1+xi[1])                   ,
                   0.25*(1-xi[1]*xi[1])*(1-xi[0])                   ,
                  -0.5*(1-xi[0])*(1-xi[1])*xi[2]                    ,
                  -0.5*(1+xi[0])*(1-xi[1])*xi[2]                    ,
                  -0.5*(1+xi[0])*(1+xi[1])*xi[2]                    ,
                  -0.5*(1-xi[0])*(1+xi[1])*xi[2];

                    return dNdXi_;
                }

                Matrix3d Jacobian(const Ref<const dNdXiSized >& dNdXi, const Ref<const Matrix<double, nNodes*nDim, 1>>& coordinates)
                {
                    // convenience wrapper to the templated version of the Jacobian
                    return FiniteElement::Jacobian<nDim, nNodes> ( dNdXi, coordinates );
                }

                BSized B(const Ref<const  dNdXiSized>& dNdX)
                {
                    // convenience wrapper to the templated version of the Spatial2D B Operator
                    return FiniteElement::Spatial3D::B<nNodes> ( dNdX ) ;
                }
            } // end of namespace Hexa20
        } // end of namespace Spatial3D
    }
} // end of namespace bft
