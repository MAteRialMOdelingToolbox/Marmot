#include "bftVoigt.h"
#include "bftFunctions.h"
#include "bftConstants.h"
#include <iostream>
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
        namespace NumIntegration
        {

            const Eigen::Matrix<double, 4, 2> gaussPts2d_2x2()
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

