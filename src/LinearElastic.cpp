#include "LinearElastic.h"
#include "bftTypedefs.h"
#include "bftFunctions.h"
#include "bftAbaqusUtility.h"
#include "bftVoigt.h"

using namespace bft;

void LinearElastic::computeStress(
        double *stress, 
        double* dStressDDStrain,  
        const double *strainOld,
        const double *dStrain,
        const double* timeOld,
        const double  dT,
        double& pNewDT)
{		
    //material properties
    const double& E  =                      materialProperties[0];
    const double& nu =                      materialProperties[1];

    const int nTensor = 6;

    Matrix6 Cel = 							mechanics::Cel(E,nu);

    Map<VectorXd>                           abqNomStress(stress, nTensor);
    Map<const VectorXd>                     abqDE(dStrain, nTensor);
    Map<MatrixXd>                           abqC(dStressDDStrain, nTensor, nTensor);   

    bft::Vector6 oldStress =                Vector6::Zero(); oldStress.head(nTensor)= abqNomStress;
    bft::Vector6 dE =                       Vector6::Zero(); dE.head(nTensor)       = abqDE;

    // Zero strain  increment check
    if((dE.array() == 0).all())
        return backToAbaqus(Cel, abqC, oldStress, abqNomStress, nTensor);               

    Vector6 newStress = oldStress + Cel * dE;                                                                         

    return backToAbaqus(Cel, abqC,newStress, abqNomStress, nTensor);
}

