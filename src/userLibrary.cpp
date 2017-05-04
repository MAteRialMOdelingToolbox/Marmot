#include <map>
#include <tuple>
//#include "userLibrary.h"
#include "bftUel.h"
#include "bftTypedefs.h"
#include <iostream>
#include <string>

#ifdef ModLeon 
    #include "umatModLeon.h"
#endif
#ifdef ShotLeon 
    #include "umatShotLeon.h"
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
#ifdef linearElastic 
    #include "umatLinearElastic.h"
#endif

namespace userLibrary{
    bft::pUmatType getUmatById(int id){
        static std::map <int, bft::pUmatType> userMaterials= { 
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
            #ifdef linearElastic 
            {9,   umatLinearElastic},  
            #endif 
            };

    return userMaterials.at(id);
    }

    bft::pUmatType getUmatByName(const std::string& nameUpperCase)
    {
        static std::map<std::string, bft::pUmatType> userMaterials= { 
            #ifdef ModLeon 
            {"MODLEON",   umatModLeon},  
            #endif
            #ifdef ShotLeon 
            {"SHOTLEON",   umatShotLeon},  
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
            #ifdef MohrCoulomb 
            {"MOHRCOULOMB",   umatMohrCoulomb},  
            #endif
            #ifdef linearElastic 
            {"LINEARELASTIC",   umatLinearElastic},  
            #endif
            };

    return userMaterials.at(nameUpperCase);
    }
} 

#ifdef uelCPE4 
    #include "cpe4.h"
#endif
#ifdef uelCPS4
    #include "cps4.h"
#endif
#ifdef uelCPS4NonLocal
    #include "cps4NonLocal.h"
#endif
#ifdef uelCPS8R
    #include "cps8r.h"
#endif
#ifdef uelCPS8RNonLocal 
    #include "cps8rNonLocal.h"
#endif


/* UEL ID System
 *
 * XXXX
 * ||||_ 4: active fields: 
 * |||__ 3: type of element 
 * ||___ 2:  number of nodes
 * |____ 1:  (number of nodes)
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
 *
 * */

namespace userLibrary{

    BftUel* UelFactory(int id, const double* elementCoordinates, double* stateVars, int nStateVars, 
            const double* propertiesElement, int nPropertiesElement, int elementNumber, 
            bft::pUmatType umat, int nStateVarsUmat, const double* propertiesUmat, int nPropertiesUmat)
    {
        switch(id){
            #ifdef uelCPS4
            case 402:{return new CPS4(elementCoordinates, stateVars, nStateVars, propertiesElement, nPropertiesElement, elementNumber, umat , nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            #endif
            #ifdef uelCPE4
            case 417:{return new CPE4(elementCoordinates, stateVars, nStateVars, propertiesElement, nPropertiesElement, elementNumber, umat , nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            #endif
            #ifdef uelCPS4NonLocal 
            case 412:{return new CPS4NonLocal(elementCoordinates, stateVars, nStateVars, propertiesElement, nPropertiesElement, elementNumber, umat , nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            #endif
            #ifdef uelCPS8R
            case 805:{return new CPS8R(elementCoordinates, stateVars, nStateVars, propertiesElement, nPropertiesElement, elementNumber, umat , nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            #endif
            #ifdef uelCPS8RNonLocal 
            case 815:{return new CPS8RNonLocal(elementCoordinates, stateVars, nStateVars, propertiesElement, nPropertiesElement, elementNumber, umat , nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            #endif
            default: return nullptr;
        }
    }

} 
