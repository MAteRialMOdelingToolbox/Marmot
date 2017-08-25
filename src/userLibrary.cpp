#include <map>
#include <tuple>
#include "bftUel.h"
#include "bftTypedefs.h"
#include <iostream>
#include <string>

#ifdef LinearElastic 
    #include "umatLinearElastic.h"
#endif
#ifdef ModLeon 
    #include "umatModLeon.h"
#endif
#ifdef ShotLeon 
    #include "umatShotLeon.h"
#endif
#ifdef ShotLeonV2
    #include "umatShotLeonV2.h"
#endif
#ifdef ShotLeonNonLocal
    #include "umatShotLeonNonLocal.h"
#endif
#ifdef ShotLeonV2NonLocal
    #include "umatShotLeonV2NonLocal.h"
#endif
#ifdef ModLeonNonLocal 
    #include "umatModLeonNonLocal.h"
#endif
#ifdef Meschke 
    #include "umatMeschke.h"
#endif
#ifdef SchaedlichSchweiger 
    #include "umatSchaedlichSchweiger.h"
#endif
#ifdef HoekBrown 
    #include "umatHoekBrown.h"
#endif
#ifdef UntereggerRockMass 
    #include "umatUntereggerRockMass.h"
#endif
#ifdef MohrCoulomb 
    #include "umatMohrCoulomb.h"
#endif
#ifdef UntereggerRockMassNonLocal 
    #include "umatUntereggerRockMassNonLocal.h"
#endif

namespace userLibrary{
    bft::pUmatType getUmatById(int id){
        static std::map <int, bft::pUmatType> userMaterials= { 
            #ifdef LinearElastic 
            {0,   umatLinearElastic},  
            #endif
            #ifdef ModLeon 
            {1,   umatModLeon},  
            #endif
            #ifdef ShotLeon 
            {2,   umatShotLeon},  
            #endif
            #ifdef Meschke 
            {3,   umatMeschke},  
            #endif
            #ifdef SchaedlichSchweiger 
            {4,   umatSchaedlichSchweiger},  
            #endif
            #ifdef ModLeonNonLocal 
            {5,   umatModLeonNonLocal},  
            #endif 
            #ifdef HoekBrown 
            {6,   umatHoekBrown},  
            #endif 
            #ifdef UntereggerRockMass 
            {7,   umatUntereggerRockMass},  
            #endif 
            #ifdef MohrCoulomb 
            {8,   umatMohrCoulomb},  
            #endif 
            #ifdef UntereggerRockMassNonLocal 
            {9,   umatUntereggerRockMassNonLocal},  
            #endif 
            #ifdef ShotLeonNonLocal
            {10,   umatShotLeonNonLocal},  
            #endif
            #ifdef ShotLeonV2
            {11,   umatShotLeonV2},  
            #endif
            #ifdef ShotLeonV2NonLocal
            {12,   umatShotLeonV2NonLocal},  
            #endif
            };

    try{
        return userMaterials.at(id);}
    catch (const std::exception& ) {
        std::cout << "Material with ID "<<id<<" not found!" << std::endl;
        return nullptr; }
    }

    bft::pUmatType getUmatByName(const std::string& nameUpperCase)
    {
        static std::map<std::string, bft::pUmatType> userMaterials= { 
            #ifdef LinearElastic 
            {"LINEARELASTIC",   umatLinearElastic},  
            #endif
            #ifdef ModLeon 
            {"MODLEON",   umatModLeon},  
            #endif
            #ifdef ShotLeon 
            {"SHOTLEON",   umatShotLeon},  
            #endif
            #ifdef ShotLeonV2
            {"SHOTLEONV2",   umatShotLeonV2},  
            #endif
            #ifdef ShotLeonNonLocal 
            {"SHOTLEONNONLOCAL",   umatShotLeonNonLocal},  
            #endif
            #ifdef ShotLeonV2NonLocal 
            {"SHOTLEONV2NONLOCAL",   umatShotLeonV2NonLocal},  
            #endif
            #ifdef Meschke 
            {"MESCHKE",   umatMeschke},  
            #endif
            #ifdef SchaedlichSchweiger 
            {"SCHAEDLICHSCHWEIGER",   umatSchaedlichSchweiger},  
            #endif
            #ifdef ModLeonNonLocal 
            {"MODLEONNONLOCAL",   umatModLeonNonLocal},  
            #endif
            #ifdef HoekBrown 
            {"HOEKBROWN",   umatHoekBrown},  
            #endif
            #ifdef UntereggerRockMass 
            {"UNTEREGGERROCKMASS",   umatUntereggerRockMass},  
            #endif
            #ifdef UntereggerRockMassNonLocal 
            {"UNTEREGGERROCKMASSNONLOCAL",   umatUntereggerRockMassNonLocal},  
            #endif
            #ifdef MohrCoulomb 
            {"MOHRCOULOMB",   umatMohrCoulomb},  
            #endif
            };
    try{
        return userMaterials.at(nameUpperCase);}
    catch (const std::exception& ) {
        std::cout << "Material "<<nameUpperCase<<" not found!" << std::endl;
        return nullptr; }
    }
} 

#ifdef uelDisplacement
    #include "uelDisplacementFactory.h"
#endif
#ifdef uelNonLocal
    #include "uelNonLocalFactory.h"
#endif

/* UEL ID System
 *
 * XXXX
 * ||||_ 4: type of element 
 * |||__ 3: active fields 
 * ||___ 2: number of nodes
 * |____ 1: (number of nodes)
 *
 *
 * active fields:   0: mechanical (=displacement),
 *                  1: mechanical + nonlocal damage,
 *
 * type of element: 1: 1D full integration, 
 *                  2: 2D full integration, plane stress 
 *                  3: 3D full integration, 
 *                  4: 1D red. integration, 
 *                  5: 2D red. integration, plane stress
 *                  6: 3D red. integration
 *                  7: 2D full integration, plane strain 
 *                  8: 2D red. integration, plane strain 
 *
 * examples:
 *
 * C3D20:           2003
 * C3D8:            803
 * C3D8R:           806
 * CPS4NonLocal:    412
 * CPE4NonLocal:    472
 * */

namespace userLibrary{

    BftUel* UelFactory(int id, const double* elementCoordinates, double* stateVars, int nStateVars, 
            const double* propertiesElement, int nPropertiesElement, int elementNumber, 
            bft::pUmatType umat, int nStateVarsUmat, const double* propertiesUmat, int nPropertiesUmat)
    {
        switch(id){
            #ifdef uelDisplacement
            case 402: {return UelDisplacementFactory:: generateUelCPS4 (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, umat, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 407: {return UelDisplacementFactory:: generateUelCPE4 (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, umat, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 802: {return UelDisplacementFactory:: generateUelCPS8 (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, umat, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 805: {return UelDisplacementFactory:: generateUelCPS8R (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, umat, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 803: {return UelDisplacementFactory:: generateUelC3D8 (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, umat, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 806: {return UelDisplacementFactory:: generateUelC3D8R (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, umat, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 808: {return UelDisplacementFactory:: generateUelCPE8R (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, umat, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 807: {return UelDisplacementFactory:: generateUelCPE8 (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, umat, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 2003: {return UelDisplacementFactory:: generateUelC3D20 (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, umat, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 2006: {return UelDisplacementFactory:: generateUelC3D20R (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, umat, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            #endif 
            #ifdef uelNonLocal
            case 412: {return UelNonLocalFactory:: generateUelCPS4NonLocal (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, umat, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 417: {return UelNonLocalFactory:: generateUelCPE4NonLocal (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, umat, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 418: {return UelNonLocalFactory:: generateUelCPE4RNonLocal (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, umat, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 812: {return UelNonLocalFactory:: generateUelCPS8NonLocal (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, umat, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 815: {return UelNonLocalFactory:: generateUelCPS8RNonLocal (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, umat, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 813: {return UelNonLocalFactory:: generateUelC3D8NonLocal (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, umat, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 816: {return UelNonLocalFactory:: generateUelC3D8RNonLocal (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, umat, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 817: {return UelNonLocalFactory:: generateUelCPE8NonLocal (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, umat, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 818: {return UelNonLocalFactory:: generateUelCPE8RNonLocal (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, umat, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 2013: {return UelNonLocalFactory:: generateUelC3D20NonLocal (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, umat, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 2016: {return UelNonLocalFactory:: generateUelC3D20RNonLocal (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, umat, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            #endif 

            default:{ std::cout << "bftUserLibrary: Element with ID " << id << " not found" << std::endl; exit(-1);}
        }
    }

} 
