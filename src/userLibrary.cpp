#include <map>
#include <tuple>
#include "userLibrary.h"
#include "bftUel.h"
#include "bftTypedefs.h"
#include "bftMaterial.h"
#include "bftVoigt.h"
#include <iostream>
#include <string>

#ifdef LINEARELASTIC
    #include "LinearElastic.h"
#endif
#ifdef LINEARELASTICSOLIDIFICATIONCREEP
    #include "LESolidificationCreep.h"
#endif
#ifdef MODLEON
    #include "ModLeon.h"
#endif
#ifdef MODLEONSEMIEXPLICIT
    #include "ModLeonSemiExplicit.h"
#endif
#ifdef MODLEONADAPTIVE
    #include "ModLeonAdaptive.h"
#endif
#ifdef MODLEONSEMIEXPLICITADAPTIVE
    #include "ModLeonSemiExplicitAdaptive.h"
#endif
#ifdef MODLEONPLANESTRESS
    #include "ModLeonPS.h"
#endif
#ifdef SHOTLEON
    #include "ShotLeon.h"
#endif
//#ifdef ShotLeonV2
    //#include "ShotLeonV2.h"
//#endif
//#ifdef SHOTLEONNONLOCAL
    //#include "ShotLeonNonLocal.h"
//#endif
//#ifdef ShotLeonV2NonLocal
    //#include "ShotLeonV2NonLocal.h"
//#endif
#ifdef MODLEONNONLOCAL
    #include "ModLeonNonLocal.h"
#endif
#ifdef MESCHKE
    #include "Meschke.h"
#endif
#ifdef SCHAEDLICHSCHWEIGER
    #include "SchaedlichSchweiger.h"
#endif
//#ifdef HoekBrown 
    //#include "materialCodeHoekBrown.h"
//#endif
//#ifdef UntereggerRockMass 
    //#include "materialCodeUntereggerRockMass.h"
//#endif
//#ifdef MohrCoulomb 
    //#include "materialCodeMohrCoulomb.h"
//#endif
//#ifdef UntereggerRockMassNonLocal 
    //#include "materialCodeUntereggerRockMassNonLocal.h"
//#endif
namespace userLibrary{

    MaterialCode getMaterialCodeFromName(const std::string& materialCode)
    {
        if(     materialCode == "LINEARELASTIC" )                   return MaterialCode::LinearElastic;
        else if(materialCode == "LESOLIDIFICATIONCREEP" )           return MaterialCode::LinearElasticSolidificationCreep;
        else if(materialCode == "MODLEON" )                         return MaterialCode::ModLeon;
        else if(materialCode == "SHOTLEON" )                        return MaterialCode::ShotLeon;
        else if(materialCode == "MODLEONPLANESTRESS" )              return MaterialCode::ModLeonPlaneStress;
        else if(materialCode == "MODLEONSEMIEXPLICIT" )             return MaterialCode::ModLeonSemiExplicit;
        else if(materialCode == "MODLEONADAPTIVE" )                 return MaterialCode::ModLeonAdaptive;
        else if(materialCode == "MODLEONSEMIEXPLICITADAPTIVE" )     return MaterialCode::ModLeonSemiExplicitAdaptive;
        else if(materialCode == "MODLEONNONLOCAL" )                 return MaterialCode::ModLeonNonLocal;
        else if(materialCode == "MESCHKE" )                         return MaterialCode::Meschke;
        else if(materialCode == "SCHAEDLICHSCHWEIGER" )             return MaterialCode::SchaedlichSchweiger;

        else{ throw std::invalid_argument("bftUserLibrary: Material Not Found: "+materialCode);}
    }

    BftMaterial* bftMaterialFactory( MaterialCode materialCode,
                                    double *stateVars,
                                    int nStateVars,
                                    const double* materialProperties, 
                                    int nMaterialProperties,
                                    int element, 
                                    int gaussPt)
    {
        switch(materialCode)
        {
            #ifdef LINEARELASTIC
            case LinearElastic: { return new class LinearElastic(stateVars, nStateVars, materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef LINEARELASTICSOLIDIFICATIONCREEP
            case LinearElasticSolidificationCreep: { return new class LinearElasticSolidificationCreep(stateVars, nStateVars, materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef MODLEON
            case ModLeon: { return new class ModLeon(stateVars, nStateVars, materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef SHOTLEON
            case ShotLeon: { return new class ShotLeon(stateVars, nStateVars, materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef MODLEONSEMIEXPLICIT
            case ModLeonSemiExplicit: { return new class LinearElastic(stateVars, nStateVars, materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef MODLEONADAPTIVE
            case ModLeonAdaptive: { return new class ModLeonAdaptive(stateVars, nStateVars, materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef MODLEONSEMIEXPLICITADAPTIVE
            case ModLeonSemiExplicitAdaptive: { return new class ModLeonSemiExplicitAdaptive(stateVars, nStateVars, materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef MODLEONNONLOCAL
            case ModLeonNonLocal: { return new class ModLeonNonLocal(stateVars, nStateVars, materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef MESCHKE
            case Meschke: { return new class Meschke(stateVars, nStateVars, materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef SCHAEDLICHSCHWEIGER
            case SchaedlichSchweiger: { return new class SchaedlichSchweiger(stateVars, nStateVars, materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            default: throw std::invalid_argument("bftUserLibrary: Invalid Material Code Requested!");
        }
    }

} 

#ifdef UELDISPLACEMENT
    #include "uelDisplacementFactory.h"
#endif
#ifdef UELNONLOCAL
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

    BftUel* UelFactory(int elementCode, const double* elementCoordinates, double* stateVars, int nStateVars, 
            const double* propertiesElement, int nPropertiesElement, int elementNumber, 
            MaterialCode materialCode,
            int nStateVarsUmat, const double* propertiesUmat, int nPropertiesUmat)
    {
        switch(elementCode){
            #ifdef UELDISPLACEMENT
            case 402: {return UelDisplacementFactory:: generateUelCPS4 (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, materialCode , nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 407: {return UelDisplacementFactory:: generateUelCPE4 (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 802: {return UelDisplacementFactory:: generateUelCPS8 (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 805: {return UelDisplacementFactory:: generateUelCPS8R (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 803: {return UelDisplacementFactory:: generateUelC3D8 (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 806: {return UelDisplacementFactory:: generateUelC3D8R (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 808: {return UelDisplacementFactory:: generateUelCPE8R (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 807: {return UelDisplacementFactory:: generateUelCPE8 (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 2003: {return UelDisplacementFactory:: generateUelC3D20 (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 2006: {return UelDisplacementFactory:: generateUelC3D20R (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            #endif 
            #ifdef UELNONLOCAL
            case 412: {return UelNonLocalFactory:: generateUelCPS4NonLocal (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 417: {return UelNonLocalFactory:: generateUelCPE4NonLocal (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 418: {return UelNonLocalFactory:: generateUelCPE4RNonLocal (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 812: {return UelNonLocalFactory:: generateUelCPS8NonLocal (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 815: {return UelNonLocalFactory:: generateUelCPS8RNonLocal (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 813: {return UelNonLocalFactory:: generateUelC3D8NonLocal (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 816: {return UelNonLocalFactory:: generateUelC3D8RNonLocal (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 817: {return UelNonLocalFactory:: generateUelCPE8NonLocal (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 818: {return UelNonLocalFactory:: generateUelCPE8RNonLocal (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 2013: {return UelNonLocalFactory:: generateUelC3D20NonLocal (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case 2016: {return UelNonLocalFactory:: generateUelC3D20RNonLocal (elementCoordinates, stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            #endif 

            default:{throw std::invalid_argument("bftUserLibrary: Element Not Found");}
        }
    }

    void extendAbaqusToVoigt(double* stress6, double* stress, double* strain6, const double* strain, double* dStrain6, const double* dStrain, int nDirect, int nShear)
    {
        bft::mVector6 s(stress6); 
        Map< VectorXd> abqStress ( stress, nDirect + nShear );
        s.setZero();
        s.head(nDirect) =           abqStress.head(nDirect);
        s.segment(3, nShear) =      abqStress.tail(nShear);

        bft::mVector6 e(strain6); 
        Map< const VectorXd> abqStrain( strain, nDirect + nShear );
        e.setZero();
        e.head(nDirect) =           abqStrain.head(nDirect);
        e.segment(3, nShear) =      abqStrain.tail(nShear);

        bft::mVector6 de(dStrain6); 
        Map< const VectorXd> abqDStrain( dStrain, nDirect + nShear );
        de.setZero();
        de.head(nDirect) =      abqDStrain.head(nDirect);
        de.segment(3, nShear) = abqDStrain.tail(nShear);

    }

    void backToAbaqus(double* stress, double* stress6, double* dStressDDStrain, double* dStressDDStrain66, int nDirect, int nShear)
    {
        bft::mVector6 s(stress6); 
        Map< VectorXd> abqStress ( stress, nDirect + nShear );

        abqStress.head(nDirect) = s.head(nDirect);
        abqStress.tail(nShear) = s.segment(3, nShear);

        bft::mMatrix6 C(dStressDDStrain66);
        Map<MatrixXd> abqC ( dStressDDStrain, nDirect + nShear, nDirect + nShear);

        if(nDirect == 2 ){
            abqC = bft::mechanics::getPlaneStressTangent(C); }
        else{
            abqC.topLeftCorner(nDirect, nDirect) =      C.topLeftCorner(nDirect, nDirect);
            abqC.bottomRightCorner(nShear, nShear) =    C.block(3,3, nShear, nShear);
            abqC.topRightCorner(nDirect, nShear) =      C.block(0,3, nDirect, nShear);
            abqC.bottomLeftCorner(nShear, nDirect) =    C.block(0,3, nShear, nDirect);}
    }

} 
