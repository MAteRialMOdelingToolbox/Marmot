#include <iostream>
#include "Eigen/Core"
#include "../include/bftVoigt.h"
#include "../include/bftTypedefs.h"
#include "../include/bftFiniteElement.h"

using namespace bft;

Matrix3d test_PlaneStressTangentAnalytically(double E, double nu)
{

    Matrix3d CAnalytic = Matrix3d::Zero();
    CAnalytic(0,0) = E/(1-nu*nu);
    CAnalytic(1,1) = E/(1-nu*nu);
    CAnalytic(2,2) = E/(1-nu*nu) * (1-nu)/2.;
    CAnalytic(0,1) = E/(1-nu*nu) * nu;
    CAnalytic(1,0) = E/(1-nu*nu) * nu;

    Matrix6 C = mechanics::Cel(E,nu); 
    // Test Plane Stress Tangent
    Matrix3d CPlaneStress = mechanics::getPlaneStressTangent(C);

    std::cout << "PLANE STRESS TANGENT WITH ANALYTIC FORMULAS " << std::endl;
    std::cout << CAnalytic << std::endl;
    std::cout << "PLANE STRESS TANGENT OF PLANE STRESS WRAPPER " << std::endl;
    std::cout << CPlaneStress << std::endl;
}

void test_NBold(const int nDoF)
{
    const Matrix<double, 4, 2> gp = bft::NumIntegration::gaussPts2d_2x2();
    //Vector3d NLinear = bft::FiniteElement::Truss3::shapeFunctions(gp(0,0));
    Vector3d NLinear = bft::FiniteElement::Truss3::shapeFunctions(gp(0,0));
    std::cout << "N " << std::endl;
    std::cout << bft::FiniteElement::createNBold(NLinear, nDoF) << std::endl;
}

int main(void){

    std::cout << "this is the test area for bftMechanics" << std::endl; 
    
    test_PlaneStressTangentAnalytically(20000, 0.3); 

    test_NBold(3);

    return 0;

}
