#include <aba_for_c.h>
#include "userLibrary.h"
#include <iostream>
#include <string>
#include "bftUel.h"
#include "bftMaterialHypoElastic.h"

//These functions are provided to the 'sub' UMATs for easy printing Messages and Warnings. 
//
#ifndef FOR_NAME
    #define FOR_NAME(a) a##_
#endif 
namespace MainConstants
{
	bool printWarnings = false;    
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
        const double U[/*nDof*/],                                   // current solution (end of increment)
        const double dU[/*mlvarx=nDof(?), nRightHandSide*/],        // increment of solutions
        const double UDot[/*nDof*/],                                // first derivative (velocity..)
        const double UDotDot[/*nDof*/],                             // second derivative (acceleration..)
        const int &elementType,                                     // user defined element type id
        const double time[2],                                       // 1: time of step, 2: total time
        const double &dTime,                                        // time increment
        const int& stepNumber,
        const int& incrementNumber,
        const int& elementNumber,                  
        const double solutionParams[],                              // solution procedure dependent parameters 
        const int& nDLoadActive,                                    // id of load / flux currently active on this element
        const int distributedLoadTypes[/*mDLoads, * */],            // An array containing the integers used to define distributed load types for the element.  
        const double distributedLoadMags[/*mDloads, * */],          // magnitudes @ end of increment
        const double Predef[/*nPredef*/], 
        const int &nPredef,
        const int lFlags[],
        const int mlvarx[],
        const double dDistributedLoadMags[/*mDloads, * */],         // increment of magnitudes 
        const int &mDload,                                          // total number of distributed loads and fluxes defined on this element
        double &pNewDT,         
        const int integerProperties[],
        const int &nIntegerProperties, 
        const double &period)
{    
        // get umatPointer and number of stateVars for umat
        userLibrary::MaterialCode materialID =  static_cast<userLibrary::MaterialCode>( integerProperties[0] );
        const int nStateVarsUmat =          integerProperties[1];
        const int nPropertiesUmat =         integerProperties[2];
        const int nPropertiesElement =      integerProperties[3];

        // get additional definitions: [active Geostatic stress], 
        const bool activeGeostatic = nIntegerProperties > 4 && integerProperties[4] > 0; 
        
        int sumProps = 0;
        for (int i = 2; i< nIntegerProperties; i++)
            sumProps += integerProperties[i];
        sumProps = activeGeostatic ? sumProps+5 : sumProps;

        if( sumProps < nProperties){
            std::cout << "insufficient properties defined" << sumProps  << " / " << nProperties << std::endl;
            pNewDT = 1e-36;
            return;}

        //bft::pUmatType umatPointer = userLibrary::getUmatById(materialID);

        const double* propertiesUmat =    &properties[0];
        const double* propertiesElement = &properties[nPropertiesUmat];

        BftUel* myUel = userLibrary::UelFactory(elementType, 
                                                coordinates,
                                                stateVars,
                                                nStateVars,
                                                propertiesElement,
                                                nPropertiesElement, 
                                                elementNumber,
                                                materialID,
                                                nStateVarsUmat,
                                                propertiesUmat, 
                                                nPropertiesUmat);

        // apply geostatic stress by setting values to statevars corresponding to stress 
        switch(lFlags[0]) {
            case Abaqus::UelFlags1::GeostaticStress: {
                myUel->setInitialConditions(BftUel::GeostaticStress, &properties[nPropertiesUmat+nPropertiesElement] ); 
                break;}
            default: break;}       

        // compute K and P 
        myUel->computeYourself(U , dU, rightHandSide, KMatrix, time, dTime, pNewDT); 
                     if (pNewDT < 0.25)
                         pNewDT = 0.25;

        // recompute distributed loads in nodal forces and add it to P 
        for (int i =0; i<mDload; i++){
            if (distributedLoadMags[i]<1.e-16)
                continue;
            myUel->computeDistributedLoad(BftUel::Pressure, rightHandSide, distributedLoadTypes[i], &distributedLoadMags[i], time, dTime);}

        delete myUel;
}

extern "C" void FOR_NAME(umat)(
        /*to be def.*/  double stress[],                // stress vector in order: S11, S22, (S33), S12, (S13), (S23) 
        /*to be def.*/  double stateVars[],        // solution dependent state variables; passed in values @ beginning of increment -> set to values @ end of increment
        /*to be def.*/  double dStressDDStrain[],  // material Jacobian matrix ddSigma/ddEpsilon
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
        const   int &nStateVars,                            // number of solution dependent state variables associated with this mat. type
        const   double materialProperties[],                         // user def. array of mat. constants associated with this material
        const   int &nMaterialProperties,                            // number of user def. variables
        const   double coords[3],                       // coordinates of this point
        const   double dRot[9],                         // rotation increment matrix 3x3
        /*may be def.*/ double &pNewDT,                 // propagation for new time increment
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
          
        userLibrary::MaterialCode materialCode = static_cast<userLibrary::MaterialCode> ( stateVars[nStateVars-1] );
        if ( materialCode <= 0){
            const std::string materialName(matName);
            materialCode = userLibrary::getMaterialCodeFromName ( materialName.substr(0, materialName.find_first_of(' ')). substr(0, materialName.find_first_of('-'))  ); 
            stateVars[nStateVars-1] = static_cast<double> (materialCode);}

        BftMaterialHypoElastic* material = dynamic_cast<BftMaterialHypoElastic*> (bftMaterialFactory( 
                                    materialCode, stateVars, nStateVars-1, 
                                    materialProperties, nMaterialProperties, noEl, nPt));

        material->setCharacteristicElementLength(charElemLength);
        
        double stress6[6], strain6[6], dStrain6[6], dStressDDStrain66[36];
        userLibrary::extendAbaqusToVoigt(stress6, stress, strain6, strain, dStrain6, dStrain, nDirect, nShear);
        if(nDirect == 3) 
            material->computeStress(stress6, dStressDDStrain66, strain6, dStrain6, time, dtime, pNewDT);
        else if(nDirect == 2)
            material->computePlaneStress(stress6, dStressDDStrain66, strain6, dStrain6, time, dtime, pNewDT);

        if(pNewDT < 0.25){
            pNewDT = 0.25;
            return;}

        userLibrary::backToAbaqus(stress, stress6, dStressDDStrain, dStressDDStrain66, nDirect, nShear);
}

