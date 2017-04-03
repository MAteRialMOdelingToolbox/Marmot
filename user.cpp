#include <aba_for_c.h>
#include "userLibrary.h"
#include <iostream>
#include <map>
#include <tuple>
#include <string>

//These functions are provided to the 'sub' UMATs for easy printing Messages and Warnings. 
#ifndef FOR_NAME
    #define FOR_NAME(a) a##_
#endif 
namespace MainConstants
{
	bool printWarnings = true;    
	bool printMessages = false;
}

extern "C" 
{
        void FOR_NAME(stdb_abqerr)(const int *lop, const char* stringZT, const int *intArray, const double *realArray, const char *appendix, const int lengthString, const int lengthAppendix);
        void FOR_NAME(xit)();
}

extern "C" bool warningToMSG(const std::string& message)
{
    // return always false(!)
    const int lop = -1;
    if(MainConstants::printWarnings)
            FOR_NAME(stdb_abqerr)(&lop, message.c_str(), nullptr , nullptr , nullptr, message.length(), 0);
    return false;
}

extern "C" bool notificationToMSG(const std::string& message)
{
    // return always true(!)
    const int lop = 1;
    if(MainConstants::printMessages)
			FOR_NAME(stdb_abqerr)(&lop, message.c_str(), nullptr , nullptr , nullptr, message.length(), 0);
    return true;
}


extern "C" void FOR_NAME(uel)(
        double rightHandSide[/*lVarx , nRightHandSide*/],           // right hand side load vector(s) 1: common, 2: additional for RIKS (see documentation)
        double KMatrix[/*nDof * nDof*/],                            // stiffness matrix 
        double stateVars[],                                         // solution dependent state variables; passed in values @ beginning of increment -> set to values @ end of increment
        double energies[8],                                         // may be updated: energies
        const int &nDegreesOfFreedom ,
        const int &nRightHandSide,                
        const int &nStateVars, 
        const double properties[/*nProperties*/],
        const int &nProperties,
        const double coordinates[/*mcrd, nNodes*/],                    // undeformed coordinates of the node respective DOFs
        const int &maxNCoords,                                        // max number of coordinates (see documentation)
        const int &nNodes,                  
        const double U_[/*nDof*/],                                   // current solution (end of increment)
        const double dU_[/*mlvarx=nDof(?), nRightHandSide*/],        // increment of solutions
        const double UDot_[/*nDof*/],                                // first derivative (velocity..)
        const double UDotDot_[/*nDof*/],                             // second derivative (acceleration..)
        const int &elementType,                                     // user defined element type id
        const double time[2],                                       // 1: time of step, 2: total time
        const double &dTime,                                        // time increment
        const int& stepNumber,
        const int& incrementNumber,
        const int& elementNumber,                  
        const double solutionParams[],                              // solution procedure dependent parameters 
        const int& nDLoadActive,                                    // id of load / flux currently active on this element
        const int distributedLoadTypes[/*mDLoads, * */],            // An array containing the integers used to define distributed load types for the element.  
        const double distributedLoadMags[/*mDloads, * */],          // magnitudes @ beginnning of increment
        const double Predef[/*nPredef*/], 
        const int &nPredef,
        const int lFlags[],
        const int mlvarx[],
        const double dDistributedLoadMags[/*mDloads, * */],         // increment of maginitudes
        const int &mDload,                                          // total number of distributed loads and fluxes defined on this element
        double &pNewdT,         
        const int integerProperties[],
        const int &nIntegerProperties, 
        const double &period)
{    
        // get umatPointer and number of stateVars for umat
        const int& materialID =           integerProperties[0];
        const int nStateVarsUmat =        integerProperties[1];

        bft::pUmatType umatPointer =      userLibrary::getUmatById(materialID);
        bft::pSimpleUelWithUmatType simpleUel = userLibrary::getSimpleUelWithUmatById(elementType);

        simpleUel(  rightHandSide, KMatrix, stateVars, nStateVars, 
                    properties, nProperties, coordinates, U_, dU_, 
                    time, dTime, elementNumber, pNewdT, integerProperties, 
                    nIntegerProperties, umatPointer, nStateVarsUmat);  
}

extern "C" void FOR_NAME(umat)(
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
        const   char matName[80],                       // user defined material name, attention: If Intel Compiler is used, matNameLength *may*
                                                        // be passed directly after this argument
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
        const   int matNameLength                       // length of Material Name := 80, passed in when FORTRAN calls c/c++: Microsoft C compiler AND GCC (it *may* differ for IntelC++)
        ){       
           
        const std::string materialName(matName);

        bft::pUmatType umat = userLibrary::getUmatByName(materialName.substr(0, materialName.find_first_of(' ')));

        if(nDirect == 3) 
            umat(stress, stateVariables, jacobianSigmaEpsilon, sSE, sPD, sCD, rpl, ddSigma_ddTemp, 
                    dRpl_dEpsilon, dRpl_dTemp, strain, dStrain, time, dtime, temp, dTemp,
                    preDef, dPreDef, matName, nDirect, nShear, nTensor, nStatV, props, 
                    nProps, coords, dRot, pNewdT, charElemLength, dfGrd0, dfGrd1, noEl, nPt,
                    layer, kSectPt, jStep, kInc, matNameLength);	

        else if(nDirect == 2)
            userLibrary::umatPlaneStressWrapped(umat, stress, stateVariables, jacobianSigmaEpsilon, sSE, sPD, sCD, rpl, ddSigma_ddTemp, 
                    dRpl_dEpsilon, dRpl_dTemp, strain, dStrain, time, dtime, temp, dTemp,
                    preDef, dPreDef, matName, nDirect, nShear, nTensor, nStatV, props, 
                    nProps, coords, dRot, pNewdT, charElemLength, dfGrd0, dfGrd1, noEl, nPt,
                    layer, kSectPt, jStep, kInc, matNameLength);	
}

