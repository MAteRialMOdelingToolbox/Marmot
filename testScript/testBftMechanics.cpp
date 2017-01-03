#include <iostream>
#include "Eigen/Core"
#include "../include/bftVoigt.h"
#include "../include/bftTypedefs.h"

using namespace bft;

Matrix3d testPlaneStressTangentAnalytically(double E, double nu)
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


int main(void){

    std::cout << "this is the test area for bftMechanics" << std::endl; 
    testPlaneStressTangentAnalytically(20000, 0.3); 

    return 0;

}
