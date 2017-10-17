#include "bftFiniteElement.h"
#include "bftTypedefs.h"
#include <iostream>

namespace bft{

    //****************************************************
    namespace FiniteElement
    {
        ElementShapes getElementShapeByMetric(int nDim, int nNodes){

            switch(nDim){

                case (2):{ switch(nNodes){
                             case 4:    { return ElementShapes::Quad4;}
                             case 8:    { return ElementShapes::Quad8;}
                             default:   {std::cout << "Undefined element shape in nDim" << nDim << ", nNodes" << nNodes << std::endl; exit(-1); }
                         }}
                case (3):{ switch(nNodes){
                             case 8:    { return ElementShapes::Hexa8;}
                             case 20:   { return ElementShapes::Hexa20;}
                             default:   {std::cout << "Undefined element shape in nDim" << nDim << ", nNodes" << nNodes << std::endl; exit(-1); }
                         }}

                default:{std::cout << "Undefined element shape in nDim" << nDim << ", nNodes" << nNodes << std::endl; exit(-1); }
            }
        }
        // create a large N matrix dependent on the degrees of freedom per node
        MatrixXd NB(const Ref<const VectorXd>& N, const int nDoFPerNode){

            MatrixXd N_(nDoFPerNode, N.size()*nDoFPerNode);
            N_ = MatrixXd::Zero(nDoFPerNode, N.size()*nDoFPerNode);

            for (int i=0; i<N.size(); i++){
                for (int j=0; j<nDoFPerNode; j++){
                    N_(j,nDoFPerNode*i+j) = N(i);
                }
            }

            return N_;} 

        MatrixXd Jacobian(const MatrixXd& dNdXi, const VectorXd& coordinates)
        {
            int nDim = dNdXi.rows();
            int nNodes = dNdXi.cols();

            MatrixXd J_ = MatrixXd::Zero(nDim, nDim);

            for(int i = 0; i < nDim; i++)		// loop over global dimensions
                for(int j=0; j < nDim; j++)		// loop over local dimensions
                    for(int k=0; k<nNodes; k++) // Loop over nodes
                        J_(i, j) += dNdXi(j, k) * coordinates(i + k*nDim);
            return J_;
        }


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
                namespace Boundary2
                {

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

                    const Matrix2d gaussPts1d_2(int elementFace){
                        // create gausspoints in order to use shapefunctions of quad4 element
                        Matrix2d gp;
                        double xi = 0.577350269189625764509;
                        switch(elementFace){
                            case 1: { gp << -xi, 1,
                                        xi, 1;  break;}
                            case 2: { gp << -1, -xi,
                                        -1,  +xi; break;}
                            case 3: { gp <<-xi, -1,
                                        +xi, -1; break;}
                            case 4: { gp << 1, -xi,
                                        1,  +xi;  break;}
                        }
                        return gp;
                    }

                    NSized dNdXi(int elementFace, const Ref<const Vector2d>& xi){
                        // derivative of shapeFunction of element face and direction of gradient 
                        // following counterclockwise element numbering
                        switch(elementFace){
                            case 1: { return -Quad4::dNdXi(xi).row(0); break;}
                            case 2: { return -Quad4::dNdXi(xi).row(1); break;}
                            case 3: { return Quad4::dNdXi(xi).row(0); break;}
                            case 4: { return Quad4::dNdXi(xi).row(1); break;}
                            default: {std::cout << "Boundary 2: invalid face ID specifed" << std::endl; exit(-1);}
                        }
                    }
                } 


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
                        default: {std::cout << "Quad8: invalid face ID specifed" << std::endl; exit(-1);}
                    }
                }

                Vector6 get2DCoordinateIndicesOfBoundaryTruss(int elementFace){
                    Vector6 truss3IndicesInQuad8;
                    switch(elementFace){
                        case 1: truss3IndicesInQuad8 <<     0,1,    2,3, /* 4,5,    6,7, */ 8,9/*   10,11,  12,13,   14,15 */; break;
                        case 2: truss3IndicesInQuad8 << /*  0,1,*/  2,3,    4,5/*,  6,7,    8,9*/,  10,11/*,12,13,   14,15 */; break;
                        case 3: truss3IndicesInQuad8 << /*  0,1,    2,3,  */4,5,    6,7,/*, 8,9     10,11,*/12,13/*, 14,15 */; break;
                        case 4: truss3IndicesInQuad8 <<     6,7,/*  2,3,    4,5, */ 0,1,/*  8,9     10,11,  12,13,*/ 14,15 ; break;
                        default: {std::cout << "Quad8: invalid face ID specifed" << std::endl; exit(-1);}

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
                    n <<    Jacobian(1),
                      -Jacobian(0); 
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

                Vector2d TangentialVector(const Vector2d& Jacobian)
                {
                    return Jacobian / Jacobian.norm();
                }

                Vector2d NormalVector(const Vector2d& Jacobian)
                {
                    Vector2d n;
                    n <<    Jacobian(1),
                      -Jacobian(0); 

                    return n / n.norm();
                }
            } // end of namespace Truss3
        } // end of namespace Spatial2D

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

    namespace BoundaryElementFactory{
        
        using bft::FiniteElement::ElementShapes;
       
        VectorXd getBoundaryNodeList(bft::FiniteElement::ElementShapes shape, const int& elementFace){
            switch(shape)
            {
                case(ElementShapes::Quad4):  { using namespace FiniteElement::Spatial2D::Quad4;
                                               return get2DCoordinateIndicesOfBoundaryTruss(elementFace); }
                case(ElementShapes::Quad8):  { using namespace FiniteElement::Spatial2D::Quad8;
                                               return get2DCoordinateIndicesOfBoundaryTruss(elementFace); } 
                default: {std::cout << "NodeIdxList: Invalid shape combination for boundary element" << std::endl; exit(-1);}
            }
        }
        
        MatrixXd getNormalVector(bft::FiniteElement::ElementShapes shape, const Ref<const VectorXd>& coords,  const Ref<const VectorXd>& gp){
            switch(shape)
            {
                case(ElementShapes::Quad4): { using namespace FiniteElement::Spatial2D::Truss2;
                                              return NormalVector(Jacobian(dNdXi(gp(0)), coords)); }
                case(ElementShapes::Quad8): { using namespace FiniteElement::Spatial2D::Truss3; 
                                              return NormalVector(Jacobian(dNdXi(gp(0)), coords)); }
                default: {std::cout << "Jacobian: Invalid shape combination for boundary element" << std::endl; exit(-1);}
            }
        }

        MatrixXd getGaussPointList(bft::FiniteElement::ElementShapes shape)
        {
            switch(shape)
            {
                case(ElementShapes::Quad4): {  return gaussPtList2; }
                case(ElementShapes::Quad8): {  return gaussPtList3; } 
                default: {std::cout << "Gausspoints: Invalid shape/integrationType combination boundary eement" << std::endl; exit(-1);}
            }
        }
        
        VectorXd getGaussWeights(bft::FiniteElement::ElementShapes shape)
        {
            switch(shape)
            {
                case(ElementShapes::Quad4): {  return gaussPtList2Weights; }
                case(ElementShapes::Quad8): {  return gaussPtList3Weights; } 
                default: {std::cout << "Boundary element: invalid gauss weights" << std::endl; exit(-1);}
            }
        }
        
        MatrixXd getNB(bft::FiniteElement::ElementShapes shape, const Ref<const VectorXd>& gp)
        {
            using namespace bft::FiniteElement;
            switch(shape)
            {
                case(ElementShapes::Quad4): {  return NB( Spatial2D::Truss2::N(gp(0)), 2); }
                case(ElementShapes::Quad8): {  return NB( Spatial2D::Truss3::N(gp(0)), 2); } 
                default: {std::cout << "Boundary element: invalid NB shape" << std::endl; exit(-1);}
            }
        }
    
        double getIntVol(bft::FiniteElement::ElementShapes shape, const Ref<const VectorXd>& coords, const Ref<const VectorXd>& gp)
        {
            switch(shape)
            {
                case(ElementShapes::Quad4): {  using namespace Spatial2D::Truss2; 
                                               return Jacobian(dNdXi(gp(0)), coords).norm(); }
                case(ElementShapes::Quad8): {  using namespace Spatial2D::Truss3; 
                                               return Jacobian(dNdXi(gp(0)), coords).norm(); } 
                default: {std::cout << "Boundary element: invalid integration volume" << std::endl; exit(-1);}
            }
        }
        
        }// end of namespace BoundaryElementFactory


    } // end of namespace FiniteElement




    //****************************************************
    namespace NumIntegration
    {
        MatrixXd getGaussPointList(bft::FiniteElement::ElementShapes shape, IntegrationTypes integrationType)
        {
            using bft::FiniteElement::ElementShapes;
            switch(shape)
            {
                case(ElementShapes::Quad4): {   if(integrationType == IntegrationTypes::FullIntegration)
                                                return Spatial2D::gaussPtList2x2;
                                                else
                                                    return Spatial2D::gaussPtList1x1;}

                case(ElementShapes::Quad8): {   if(integrationType == IntegrationTypes::FullIntegration)
                                                return Spatial2D::gaussPtList3x3;
                                                else
                                                    return Spatial2D::gaussPtList2x2;}

                case(ElementShapes::Hexa8): {   if(integrationType == IntegrationTypes::FullIntegration)
                                                return Spatial3D::gaussPtList2x2x2;
                                                else
                                                    return Spatial3D::gaussPtList1x1x1;}

                case(ElementShapes::Hexa20): {   if(integrationType == IntegrationTypes::FullIntegration)
                                                    return Spatial3D::gaussPtList3x3x3;
                                                 else
                                                    return Spatial3D::gaussPtList2x2x2;}

                default: {std::cout << "Invalid shape/integrationType combination" << std::endl; exit(-1);}
            }
        }

        VectorXd getGaussWeights (bft::FiniteElement::ElementShapes shape, IntegrationTypes integrationType){
            using bft::FiniteElement::ElementShapes;
            switch(shape)
            {
                case(ElementShapes::Quad4): {   if(integrationType == IntegrationTypes::FullIntegration)
                                                    return Spatial2D::gaussPtList2x2Weights;
                                                else
                                                    return Spatial2D::gaussPtList1x1Weights;}

                case(ElementShapes::Quad8): {   if(integrationType == IntegrationTypes::FullIntegration)
                                                    return Spatial2D::gaussPtList3x3Weights;
                                                else
                                                    return Spatial2D::gaussPtList2x2Weights;}

                case(ElementShapes::Hexa8): {   if(integrationType == IntegrationTypes::FullIntegration)
                                                    return Spatial3D::gaussPtList2x2x2Weights;
                                                else
                                                    return Spatial3D::gaussPtList1x1x1Weights;}

                case(ElementShapes::Hexa20): {   if(integrationType == IntegrationTypes::FullIntegration)
                                                    return Spatial3D::gaussPtList3x3x3Weights;
                                                else
                                                    return Spatial3D::gaussPtList2x2x2Weights;}

                default: {std::cout << "Invalid shape/integrationType combination" << std::endl; exit(-1);}
            }
        }
        namespace Spatial2D{
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
    } // end of namespace NumIntegration
} // end of namespace bft
