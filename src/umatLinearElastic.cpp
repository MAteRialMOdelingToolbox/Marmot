#include "umatLinearElastic.h"
#include <iostream>
#include <string>
#include <cstring>
#include "bftTypedefs.h"
#include "bftFunctions.h"
#include "bftAbaqusUtility.h"
#include "bftVoigt.h"

using namespace bft;

extern "C" void umatLinearElastic(
        /*to be def.*/  double stress[],                // stress vector in order: S11, S22, (S33), S12, (S13), (S23) 
        /*to be def.*/  double stateVariables[],        // solution dependent state variables; passed in values @ beginning of increment -> set to values @ end of increment
        /*to be def.*/  double jacobianSigmaEpsilon[],  // material Jacobian matrix ddSigma/ddEpsilon
        /*to be def.*/  double &sSE,                    // specific elastic strain energy  |-
        /*to be def.*/  double &sPD,                    // specific plastic dissipation    |---> Should be defined in Abaqus/Standard
        /*to be def.*/  double &sCD,                    // specific creep dissipation      |-
        //FOLLOWING ARGUMENTS: only in a fully coupled thermal-stress or thermal-electrical-structural analysis
        /*to be def.*/  double &rpl,                    // volumetric heat generation per unit time @ end of inc. caused by mech. working of material            
        /*to be def.*/  double ddSigma_ddTemp[],        // variation of stress with respect to temperature
        /*to be def.*/  double dRpl_dEpsilon[],         // variation of rpl with respect to strain increments
        /*to be def.*/  double &dRpl_dTemp,             // variation of rpl with respect to temperature
        const   double strain[],                        // array containing total strains @ beginning of increment     (only mechanical, no thermal); Shear Strain in engineering: e.g. gamma12 = 2*epsilon12
        const   double dStrain[],                       // array containing strain increment                           (only mechanical, no thermal)
        const   double time[2],                         // time[1]: value of step time @ beginning of current inc. or frequency
        //                                              // time[2]: value of total time @ beginning of current inc.
        const   double &dtime,                          // time increment
        const   double &temp,                           // temp @ start of inc.
        const   double &dTemp,                          // temp increment
        const   double preDef[],                        // array of interpolated values of predefined field variables @ this point @ start of inc., based on values read in nodes
        const   double dPreDef[],                       // array of inc. of pre. def. field variables
        const   char matName[80],                       // user defined material name, 
        const   int &nDirect,                           // number of direct stress components @ this point
        const   int &nShear,                            // number of engineering shear stress components @ this point
        const   int &nTensor,                           // size of stress and strain component array (nDirect + nShear)
        const   int &nStatV,                            // number of solution dependent state variables associated with this mat. type
        const   double props[],                         // user def. array of mat. constants associated with this material
        const   int &nProps,                            // number of user def. variables
        const   double coords[3],                       // coordinates of this point
        const   double dRot[9],                         // rotation increment matrix 3x3
        /*may be def.*/ double &pNewdT,                 // propagation for new time increment
        const   double &charElemLength,                 // characteristic element Length
        const   double dfGrd0[9],                       // deformation gradient @ beginning of increment      3x3 |
        const   double dfGrd1[9],                       // deformation gradient @ end of increment            3x3 |--> always stored as 3D-matrix
        const   int &noEl,                              // element number
        const   int &nPt,                               // integration Point number
        const   int &layer,                             // layer number (composite shells @ layered solids)
        const   int &kSectPt,                           // section point number within current layer
        const   int jStep[4],                           // step number **NOTE:** documentation and course material differ here, kStep[1] <-> jStep[4]
        //                                              // additional fields *may* be: procedure key, large deformation flag, perturbation step flag
        const   int &kInc,                              // increment Number
        const   int matNameLength                       // length of Material Name := 80, passed in when FORTRAN calls c/c++: Microsoft C compiler
        )                        
{		

        //material properties
        const double& E  =                      props[0];
        const double& nu =                      props[1];

	    // state variables
        bool& umatIsInitialized =               reinterpret_cast<bool&>(stateVariables[0]);

        Matrix6 Cel = 							mechanics::Cel(E,nu);
       // Matrix6 CelInv = 						mechanics::CelInverse(E,nu);
		Matrix6 C = 							Cel;

        Map<VectorXd>                           abqNomStress(stress, nTensor);
        Map<const VectorXd>                     abqDE(dStrain, nTensor);
        Map<MatrixXd>                           abqC(jacobianSigmaEpsilon, nTensor, nTensor);   
		
        bft::Vector6 nomStress =                     Vector6::Zero(); nomStress.head(nTensor)= abqNomStress;
        bft::Vector6 dE =                            Vector6::Zero(); dE.head(nTensor)       = abqDE;

        if(umatIsInitialized == false){
				umatIsInitialized =             true;}

        // Zero strain  increment check
        if((dE.array() == 0).all())
                return backToAbaqus(C, abqC, nomStress, abqNomStress, nTensor);               

        Vector6                                 effStress = nomStress;
        Vector6                                 trialStress;

        Matrix6                                 Cep;

		trialStress = effStress + Cel * dE;                                                                         
        nomStress = trialStress;                                                                                                        
        C =  Cel;                          

        return backToAbaqus(C, abqC, nomStress, abqNomStress, nTensor);
}

