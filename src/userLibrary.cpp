#include <map>
#include <tuple>
#include <iostream>
#include "userLibrary.h"

#ifdef LinearElastic 
    #include "umatLinearElastic.h"
#endif
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
            };

    return userMaterials.at(id);
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

    return userMaterials.at(nameUpperCase);
    }
} 

#ifdef uelCPE4 
    #include "uelCPE4SimpleUmat.h"
#endif
#ifdef uelCPS4
    #include "uelCPS4SimpleUmat.h"
#endif
#ifdef uelCPS4NonLocal
    #include "uelCPS4NonLocalSimpleUmat.h"
#endif
#ifdef uelCPE4NonLocal
    #include "uelCPE4NonLocalSimpleUmat.h"
#endif
#ifdef uelCPS8R
    #include "uelCPS8RSimpleUmat.h"
#endif
#ifdef uelCPS8NonLocal 
    #include "uelCPS8NonLocalSimpleUmat.h"
#endif
#ifdef uelCPS8RNonLocal 
    #include "uelCPS8RNonLocalSimpleUmat.h"
#endif
#ifdef uelCPE8NonLocal 
    #include "uelCPE8NonLocalSimpleUmat.h"
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
    bft::pSimpleUelWithUmatType getSimpleUelWithUmatById(int id){
        static std::map<int, bft::pSimpleUelWithUmatType> userElements= { 
            #ifdef uelCPE4
            {407, uelCPE4SimpleUmat} ,
            #endif
            #ifdef uelCPS4 
            {402, uelCPS4SimpleUmat} ,
            #endif
            #ifdef uelCPS4NonLocal 
            {412, uelCPS4NonLocalSimpleUmat} ,
            #endif
            #ifdef uelCPE4NonLocal 
            {417, uelCPE4NonLocalSimpleUmat} ,
            #endif
            #ifdef uelCPS8R
            {805, uelCPS8RSimpleUmat} ,
            #endif
            #ifdef uelCPS8NonLocal 
            {812, uelCPS8NonLocalSimpleUmat} ,
            #endif
            #ifdef uelCPS8RNonLocal 
            {815, uelCPS8RNonLocalSimpleUmat} ,
            #endif
            #ifdef uelCPE8NonLocal 
            {817, uelCPE8NonLocalSimpleUmat} ,
            #endif
        };
    return userElements.at(id);
    }

} 
