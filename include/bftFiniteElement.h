#include "bftVoigt.h"
#include "bftFunctions.h"
#include "bftConstants.h"
#include <iostream>

namespace bft{

	//****************************************************
	namespace FiniteElement
    {
        namespace Quad4
        {
            static constexpr int nNodes = 4;
            static constexpr int nDim = 2;
            
            Matrix<double, nNodes, 1> shapeFunctions(const Ref<const Vector2d>& xi);           
            Matrix<double, nNodes, nDim> dNdXi(const Ref <const Vector2d>& xi);
        } 

        namespace NumIntegration
        {
            const Eigen::Matrix<double, 4, 2> gaussPts2d_2x2();
        } 
    } 
} 

