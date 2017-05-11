#include "bftVoigt.h"
#include "bftFunctions.h"
#include "bftConstants.h"
#include <iostream>
#include <map>
#include "bftFiniteElement.h"

namespace bft{
	//****************************************************
	namespace FiniteElement
    {
        //****************************************************
        namespace Quad4
        {
            
            Matrix<double, nNodes, 1> shapeFunctions(const Ref<const Vector2d>& xi){
            
          /* Shape functions
             (1) _______(0)
                |       |  
                |       |  
                |_______|  
             (2)        (3) */
            
                Matrix<double, nNodes, 1> N_;
                N_ <<   1./4 * (1.+xi(0))*(1.+xi(1)),
                        1./4 * (1.-xi(0))*(1.+xi(1)),
                        1./4 * (1.-xi(0))*(1.-xi(1)),
                        1./4 * (1.+xi(0))*(1.-xi(1));
                return N_;
            }
            
            Matrix<double, nNodes, nDim> dNdXi(const Ref<const Vector2d>& xi){

                    Matrix<double, nNodes, nDim> result;
                	// 					dN(i)/dxi1		dN(i)/dxi2
                    result <<
                                    +1./4*(1+xi(1)), +1./4*(1+xi(0)),   
                                    -1./4*(1+xi(1)), +1./4*(1-xi(0)),
                                    -1./4*(1-xi(1)), -1./4*(1-xi(0)),
                                    +1./4*(1-xi(1)), -1./4*(1+xi(0));
                    return result;}

        } // end of namespace Quad4
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
        //****************************************************
        namespace NumIntegration
        {

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
    } // end of namespace FiniteElement
} // end of namespace bft
