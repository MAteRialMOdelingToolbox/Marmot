#include <map>
#include <tuple>
#include "userLibrary.h"

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

namespace userLibrary{
    bft::pUmatType getUmatById(int id){
        static std::map < int, bft::pUmatType> userMaterials= { 
            #ifdef ModLeon 
            {1,   umatModLeon},  
            #endif
            #ifdef ShotLeon 
            {2,   umatShotLeon},  
            #endif
            #ifdef Meschke 
            {3,   umatMeschke },  
            #endif
            #ifdef SchaedlichSchweiger 
            {4,   umatSchaedlichSchweiger},  
            #endif
            #ifdef ModLeonNonLocal 
            {5,   umatModLeonNonLocal},  
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

#ifdef uelCPS8R
    #include "uelCPS8RSimpleUmat.h"
#endif
#ifdef uelCPS8RNonLocal 
    #include "uelCPS8RNonLocalSimpleUmat.h"
#endif

namespace userLibrary{
    bft::pSimpleUelWithUmatType getSimpleUelWithUmatById(int id){

        static std::map<int, bft::pSimpleUelWithUmatType> userElements= { 
            #ifdef uelCPE4
            {412, uelCPE4SimpleUmat} ,
            #endif

            #ifdef uelCPS4 
            {402, uelCPS4SimpleUmat} ,
            #endif
            #ifdef uelCPS4NonLocal 
            {403, uelCPS4NonLocalSimpleUmat} ,
            #endif
            #ifdef uelCPS8R
            {802, uelCPS8RSimpleUmat} ,
            #endif
            #ifdef uelCPS8RNonLocal 
            {803, uelCPS8RNonLocalSimpleUmat} ,
            #endif
        };

    return userElements.at(id);
    }

} 
