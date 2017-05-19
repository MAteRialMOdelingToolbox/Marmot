#include "bftFiniteElement.h"
#include "bftTypedefs.h"

namespace bft{
    
    //****************************************************
	namespace FiniteElement
    {
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
        
        
        namespace Spatial2D
        {
        
            //****************************************************
            namespace Quad4
            {
                
                Vector4d N(const Ref<const Vector2d>& xi){
                
              /* Shape functions
                 (1) _______(0)
                    |       |  
                    |       |  
                    |_______|  
                 (2)        (3) */
                
                    Vector4d N_;
                    N_ <<   1./4 * (1.+xi(0))*(1.+xi(1)),
                            1./4 * (1.-xi(0))*(1.+xi(1)),
                            1./4 * (1.-xi(0))*(1.-xi(1)),
                            1./4 * (1.+xi(0))*(1.-xi(1));
                    return N_;
                }
                
                Matrix<double, nNodes,nDim> dNdXi(const Ref<const Vector2d>& xi){

                        Matrix<double, nNodes, nDim> result;
                        // 					dN(i)/dxi1		dN(i)/dxi2
                        result <<
                                        +1./4*(1+xi(1)), +1./4*(1+xi(0)),   
                                        -1./4*(1+xi(1)), +1./4*(1-xi(0)),
                                        -1./4*(1-xi(1)), -1./4*(1-xi(0)),
                                        +1./4*(1-xi(1)), -1./4*(1+xi(0));
                        return result;}

                 Vector4d get2DCoordinateIndicesOfBoundaryTruss(int elementFace){
                        Vector4d truss2IndicesInQuad4;
                        switch(elementFace){
                            case 1: truss2IndicesInQuad4 <<     0,1,    2,3 /*  4,5,    6,7,*/  ; break;
                            case 2: truss2IndicesInQuad4 << /*  0,1,*/  2,3,    4,5/*,  6,7,*/  ; break;
                            case 3: truss2IndicesInQuad4 << /*  0,1,    2,3,  */4,5,    6,7     ; break;
                            case 4: truss2IndicesInQuad4 <<     6,7,/*  2,3,    4,5, */ 0,1     ; break;
                        }
                        return truss2IndicesInQuad4;
                    }
            } // end of namespace Quad4
            
            
            //****************************************************
            namespace Quad8
            {
                Matrix<double, nNodes, 1> N(const Ref<const Vector2d>& xi){
                
              /* Shape functions
                       (6)
                 (3) _______(2)
                    |       |  
                (7) |       | (5)  
                    |_______|  
                 (0)   (4)   (1) */
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
            
                
                Matrix<double, nNodes, nDim> dNdXi(const Ref<const Vector2d>& xi){

                        const double xi0 = xi(0);
                        const double xi1 = xi(1);
                        Matrix<double, nNodes, nDim> result;
                        // 					dN(i)/dxi1		                    dN(i)/dxi2
                        result <<           (1-xi1) * (2*xi0 + xi1)/4,       (1-xi0)*(xi0+ 2*xi1)/4,
                                            (1-xi1) * (2*xi0 - xi1)/4,      -(1+xi0)*(xi0- 2*xi1)/4,
                                            (1+xi1) * (2*xi0 + xi1)/4,       (1+xi0)*(xi0+ 2*xi1)/4,
                                            (1+xi1) * (2*xi0 - xi1)/4,      (1-xi0)*(-xi0+ 2*xi1)/4,
                                                  -xi0 * (1-xi1),                -(1-xi0*xi0)/2,
                                                   (1-xi1*xi1)/2,              -xi1 * (1 + xi0),
                                                  -xi0 * (1+xi1),                +(1-xi0*xi0)/2,
                                                 -(1-xi1*xi1)/2,               -xi1 * (1 - xi0);
                        return result;
                        }
                        
                std::array<int,3> getNodesOfFace(int elementFace){
                        // create gausspoints in order to use shapefunctions of quad4 element
                        //std::array<int,3> nodes;
                        switch(elementFace){
                            case 1: { return {0, 1, 4};  break;}
                            case 2: { return {1, 2, 5};  break;}
                            case 3: { return {2, 3, 6};  break;}
                            case 4: { return {3, 0, 7};  break;} }
                        //return nodes;
                        }

                    Vector6 get2DCoordinateIndicesOfBoundaryTruss(int elementFace){
                        Vector6 truss3IndicesInQuad8;
                        switch(elementFace){
                            case 1: truss3IndicesInQuad8 <<     0,1,    2,3, /* 4,5,    6,7, */ 8,9/*   10,11,  12,13,   14,15 */; break;
                            case 2: truss3IndicesInQuad8 << /*  0,1,*/  2,3,    4,5/*,  6,7,    8,9*/,  10,11/*,12,13,   14,15 */; break;
                            case 3: truss3IndicesInQuad8 << /*  0,1,    2,3,  */4,5,    6,7,/*, 8,9     10,11,*/12,13/*, 14,15 */; break;
                            case 4: truss3IndicesInQuad8 <<     6,7,/*  2,3,    4,5, */ 0,1,/*  8,9     10,11,  12,13,*/ 14,15 ; break;
                        }
                        return truss3IndicesInQuad8;
                    }
            } // end of namespace Quad8
            
            
            //****************************************************
            namespace Truss2
            {
                /*  (0)            (1)
                //   o-------------o
                     
                     * ------|------> xi
                     * -1    0      1
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
                /* (0)      (2)     (1)
                    o------o-------o

                     <-----.-----> xi
                     -1    0    1 */

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
                
            
            //****************************************************
            namespace Boundary2
            {
                
              /* Create Edge of 2d element
               
               *     face 1
                 (1) _______(0)
                    |       |  
             face 2 |       | face 4  
                    |_______|  
                 (2)        (3) */
                  //  face 3    
                
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
            
                Vector4d dNdXi(int elementFace, const Ref<const Vector2d>& xi){
                        // derivative of shapeFunction of element face and direction of gradient 
                        // following counterclockwise element numbering
                        switch(elementFace){
                            case 1: { return -Quad4::dNdXi(xi).col(0); break;}
                            case 2: { return -Quad4::dNdXi(xi).col(1); break;}
                            case 3: { return Quad4::dNdXi(xi).col(0); break;}
                            case 4: { return Quad4::dNdXi(xi).col(1); break;}
                        }
                    }
            
            } // end of namespace Boundary2
        } // end of namespace Spatial2D
   } // end of namespace FiniteElement
 
  
   //****************************************************
   namespace NumIntegration
   {
            const Vector2d gaussPts1d_2()
            {
    		    double gp = 0.577350269189625764509;
    		    Eigen::Matrix<double, 2,  1> gaussPts;
                gaussPts << -gp, +gp; 
                return gaussPts; 
            }
            
            const Matrix<double, 4, 2> gaussPts2d_2x2()
            {
    		    double gp = 0.577350269189625764509;
    		    Eigen::Matrix<double, 4,  2> gaussPts;
                gaussPts << +gp,  +gp, 
                            -gp,  +gp,
    						-gp,  -gp, 
                            +gp,  -gp;
                return gaussPts; 
            }
    } // end of namespace NumIntegration
} // end of namespace bft
