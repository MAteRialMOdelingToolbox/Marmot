#include "bftFiniteElement.h"
#include "bftTypedefs.h"
#include <iostream>

namespace bft{

    //****************************************************
    namespace FiniteElement
    {
        namespace Spatial2D
        {
            namespace Quad4
            {

                NSized N(const Ref<const Vector2d>& xi){

                    /* 
                     *   Shape functions
                     *    (1) _______(0)
                     *       |       |  
                     *       |       |  
                     *       |_______|  
                     *    (2)        (3) 
                     */

                    NSized N_;
                    N_ <<   1./4 * (1.+xi(0))*(1.+xi(1)),
                       1./4 * (1.-xi(0))*(1.+xi(1)),
                       1./4 * (1.-xi(0))*(1.-xi(1)),
                       1./4 * (1.+xi(0))*(1.-xi(1));
                    return N_;
                }

                dNdXiSized dNdXi(const Ref<const Vector2d>& xi){

                    dNdXiSized result;

                    result <<   /* 1                2                 3                 4*/
                        /* ,xi1 */  +1./4*(1+xi(1)),    -1./4*(1+xi(1)), -1./4*(1-xi(1)), +1./4*(1-xi(1)), 
                        /*, xi2 */  +1./4*(1+xi(0)),    +1./4*(1-xi(0)), -1./4*(1-xi(0)), -1./4*(1+xi(0));

                    return result;}

                NSized get2DCoordinateIndicesOfBoundaryTruss(int elementFace){
                    NSized truss2IndicesInQuad4;
                    switch(elementFace){
                        case 1: truss2IndicesInQuad4 <<     0,1,    2,3 /*  4,5,    6,7,*/  ; break;
                        case 2: truss2IndicesInQuad4 << /*  0,1,*/  2,3,    4,5/*,  6,7,*/  ; break;
                        case 3: truss2IndicesInQuad4 << /*  0,1,    2,3,  */4,5,    6,7     ; break;
                        case 4: truss2IndicesInQuad4 <<     6,7,/*  2,3,    4,5, */ 0,1     ; break;
                    }
                    return truss2IndicesInQuad4;
                }

                Matrix2d Jacobian(const Ref<const dNdXiSized >& dNdXi, const Ref<const Vector8d>& coordinates)
                {
                    // convenience wrapper to the templated version of the Jacobian
                    return FiniteElement::Jacobian<nDim, nNodes> ( dNdXi, coordinates );
                }

                BSized B(const Ref<const  dNdXiSized>& dNdX)
                {
                    // convenience wrapper to the templated version of the Spatial2D B Operator
                    return FiniteElement::Spatial2D::B<nNodes> ( dNdX ) ;
                }


                //****************************************************
                //namespace Boundary2
                //{

                    /*    Create Edge of 2d element
                     *
                     *               face 1
                     *           (1) _______(0)
                     *              |       |  
                     *       face 2 |       | face 4  
                     *              |_______|  
                     *           (2)        (3) 
                     *               face 3    
                     * */

                    //const Matrix2d gaussPts1d_2(int elementFace){
                        //// create gausspoints in order to use shapefunctions of quad4 element
                        //Matrix2d gp;
                        //double xi = 0.577350269189625764509;
                        //switch(elementFace){
                            //case 1: { gp << -xi, 1,
                                        //xi, 1;  break;}
                            //case 2: { gp << -1, -xi,
                                        //-1,  +xi; break;}
                            //case 3: { gp <<-xi, -1,
                                        //+xi, -1; break;}
                            //case 4: { gp << 1, -xi,
                                        //1,  +xi;  break;}
                        //}
                        //return gp;
                    //}

                    //NSized dNdXi(int elementFace, const Ref<const Vector2d>& xi){
                        //// derivative of shapeFunction of element face and direction of gradient 
                        //// following counterclockwise element numbering
                        //switch(elementFace){
                            //case 1: { return -Quad4::dNdXi(xi).row(0); break;}
                            //case 2: { return -Quad4::dNdXi(xi).row(1); break;}
                            //case 3: { return Quad4::dNdXi(xi).row(0); break;}
                            //case 4: { return Quad4::dNdXi(xi).row(1); break;}
                            //default: {throw std::invalid_argument("Boundary 2: invalid face ID specifed");}
                        //}
                    //}
                //} 


                //****************************************************
            } // end of namespace Quad4
            namespace Quad8
            {
                NSized N(const Ref<const Vector2d>& xi){
                    /* Shape functions
                     *
                     *         (6)
                     *   (3) _______(2)
                     *      |       |  
                     *  (7) |       | (5)  
                     *      |_______|  
                     *   (0)        (1) 
                     *         (4)
                     *
                     * */

                    const double xi0 = xi(0);
                    const double xi1 = xi(1);
                    Matrix<double, nNodes, 1> N_;
                    N_ <<                (1.-xi0)*(1.-xi1) * (-1-xi0-xi1)/4,
                       (1.+xi0)*(1.-xi1) * (-1+xi0-xi1)/4,
                       (1.+xi0)*(1.+xi1) * (-1+xi0+xi1)/4,
                       (1.-xi0)*(1.+xi1) * (-1-xi0+xi1)/4,   
                       (1-xi0*xi0)       * (1-xi1      )/2,
                       (1+xi0      )     * (1-xi1*xi1)/2,
                       (1-xi0*xi0)       * (1+xi1     )/2,
                       (1-xi0      )     * (1-xi1*xi1)/2;

                    return N_; 
                }


                dNdXiSized dNdXi(const Ref<const Vector2d>& xi){

                    const double xi0 = xi(0);
                    const double xi1 = xi(1);
                    Matrix<double, nDim, nNodes> result;
                    result <<       
                        /*                  1                           2                           3                           4
                         *                  5                           6                           7                           8 */
                        /* ,Xi1 */     (1-xi1)*(2*xi0+xi1)/4,      (1-xi1)*(2*xi0-xi1)/4,      (1+xi1)*(2*xi0+xi1)/4,      (1+xi1)*(2*xi0-xi1)/4,
                        /* ,Xi1 */     -xi0*(1-xi1),               (1-xi1*xi1)/2,              -xi0*(1+xi1),               -(1-xi1*xi1)/2,
                        /*                  1                           2                           3                           4
                         *                  5                           6                           7                           8*/
                        /* ,Xi2 */     (1-xi0)*(xi0+2*xi1)/4,      -(1+xi0)*(xi0-2*xi1)/4,     (1+xi0)*(xi0+2*xi1)/4,      (1-xi0)*(-xi0+2*xi1)/4,
                        /* ,Xi2 */     -(1-xi0*xi0)/2,             -xi1*(1+xi0),               +(1-xi0*xi0)/2,             -xi1*(1-xi0);


                    return result;
                }

                Matrix2d Jacobian(const Ref<const dNdXiSized >& dNdXi, const Ref<const Matrix<double, nNodes*nDim, 1>>& coordinates)
                {
                    // convenience wrapper to the templated version of the Jacobian
                    return FiniteElement::Jacobian<nDim, nNodes> ( dNdXi, coordinates );
                }

                BSized B(const Ref<const  dNdXiSized>& dNdX)
                {
                    // convenience wrapper to the templated version of the Spatial2D B Operator
                    return FiniteElement::Spatial2D::B<nNodes> ( dNdX ) ;
                }

                std::array<int,3> getNodesOfFace(int elementFace){
                    // create gausspoints in order to use shapefunctions of quad4 element
                    switch(elementFace){
                        case 1: { return {0, 1, 4};  break;}
                        case 2: { return {1, 2, 5};  break;}
                        case 3: { return {2, 3, 6};  break;}
                        case 4: { return {3, 0, 7};  break;} 
                        default: {throw std::invalid_argument("Quad8: invalid face ID specifed");}
                    }
                }

                Vector6 get2DCoordinateIndicesOfBoundaryTruss(int elementFace){
                    Vector6 truss3IndicesInQuad8;
                    switch(elementFace){
                        case 1: truss3IndicesInQuad8 <<     0,1,    2,3, /* 4,5,    6,7, */ 8,9/*   10,11,  12,13,   14,15 */; break;
                        case 2: truss3IndicesInQuad8 << /*  0,1,*/  2,3,    4,5/*,  6,7,    8,9*/,  10,11/*,12,13,   14,15 */; break;
                        case 3: truss3IndicesInQuad8 << /*  0,1,    2,3,  */4,5,    6,7,/*, 8,9     10,11,*/12,13/*, 14,15 */; break;
                        case 4: truss3IndicesInQuad8 <<     6,7,/*  2,3,    4,5, */ 0,1,/*  8,9     10,11,  12,13,*/ 14,15 ; break;
                        default: {throw std::invalid_argument("Quad8: invalid face ID specifed");}

                    }
                    return truss3IndicesInQuad8;
                }
            } // end of namespace Quad8
            namespace Truss2
            {
                /*  (0)            (1)
                 *    o------------o
                 *
                 *    <------|-----> xi
                 *     -1    0     1
                 */

                Vector2d N(double  xi)
                {
                    Vector2d N;
                    N <<    (1-xi) * 0.5,
                      (1+xi) * 0.5;
                    return N;
                }

                Vector2d dNdXi(double xi)
                {
                    Vector2d dNdXi;
                    dNdXi << -0.5, 0.5;  
                    return dNdXi;
                }

                Vector2d Jacobian(const Vector2d& dNdXi, const Vector4d& coordinates)
                {
                    Vector2d J = Vector2d::Zero();

                    for(int i = 0; i < nDim; i++)		
                        for(int k=0; k<nNodes; k++) 
                            J(i) += dNdXi(k) * coordinates(i + k*nDim);
                    return J;
                }

                Vector2d TangentialVector(const Vector2d& Jacobian)
                {
                    return Jacobian / Jacobian.norm();
                }

                Vector2d NormalVector(const Vector2d& Jacobian)
                {
                    Vector2d n;
                    n <<    Jacobian(1), -Jacobian(0); 
                    return n / n.norm();
                }
            } // end of namespace Truss2


            //****************************************************
            namespace Truss3
            {
                /* 
                 * (0)    (2)     (1)
                 *  o------o-------o
                 *
                 *   <-----.-----> xi
                 *   -1    0    1 
                 * 
                 * */

                Vector3d N(double  xi)
                {
                    Vector3d N;
                    N << (xi - 1 ) * xi / 2,
                      (xi + 1 ) * xi / 2,
                      (1 - xi*xi );
                    return N;
                }

                Vector3d dNdXi(double xi)
                {
                    Vector3d dNdXi;
                    dNdXi <<    xi - 0.5, 
                          xi + 0.5,
                          -2*xi;

                    return dNdXi;
                }

                Vector2d Jacobian(const Vector3d& dNdXi, const Vector6& coordinates)
                {
                    Vector2d J = Vector2d::Zero();

                    for(int i = 0; i < nDim; i++)		
                        for(int k=0; k<nNodes; k++) 
                            J(i) += dNdXi(k) * coordinates(i + k*nDim);
                    return J;
                }

                //Vector2d TangentialVector(const Vector2d& Jacobian)
                //{
                    //return Jacobian / Jacobian.norm();
                //}

                //Vector2d NormalVector(const Vector2d& Jacobian)
                //{
                    //Vector2d n;
                    //n <<    Jacobian(1),
                      //-Jacobian(0); 

                    //return n / n.norm();
                //}
            } // end of namespace Truss3

            void modifyCharElemLengthAbaqusLike(double& charElemLength, int intPoint){
                switch(intPoint){
                    case 0:         // central node
                        charElemLength *= 2./3.; break;
                    case 1:         // corner nodes 
                    case 2: 
                    case 3: 
                    case 4:
                        charElemLength *= 5./12; break;
                    case 5:         // middle nodes
                    case 6:
                    case 7:
                    case 8:
                        charElemLength *= std::sqrt(5./18.); break; }
            }
        } // end of namespace Spatial2D 
    } // end of namespace FiniteElement
    //****************************************************
} // end of namespace bft
