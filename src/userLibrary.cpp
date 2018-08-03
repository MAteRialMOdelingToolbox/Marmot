#include <map>
#include <tuple>
#include "userLibrary.h"
#include "bftUel.h"
#include "bftTypedefs.h"
#include "bftMaterial.h"
#include "bftVoigt.h"
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
#ifdef SHOTLEONV2 
    #include "ShotLeonV2.h"
#endif
#ifdef SHOTLEONNONLOCAL
    #include "ShotLeonNonLocal.h"
#endif
#ifdef SHOTLEONV2NONLOCAL 
    #include "ShotLeonV2NonLocal.h"
#endif
#ifdef MODLEONNONLOCAL
    #include "ModLeonNonLocal.h"
#endif
#ifdef MESCHKE
    #include "Meschke.h"
#endif
#ifdef SCHAEDLICHSCHWEIGER
    #include "SchaedlichSchweiger.h"
#endif
#ifdef HOEKBROWN 
    #include "HoekBrown.h"
#endif
#ifdef UNTEREGGERROCKMASS
    #include "UntereggerRockMass.h"
#endif
#ifdef UNTEREGGERROCKMASSPLAXIS
    #include "UntereggerRockMassPlaxis.h"
#endif
#ifdef MohrCoulomb 
    #include "materialCodeMohrCoulomb.h"
#endif
#ifdef UNTEREGGERROCKMASSNONLOCAL 
    #include "UntereggerRockMassNonLocal.h"
#endif
#ifdef LINEARELASTICNONLOCAL
    #include "LinearElasticNonLocal.h"
#endif
namespace userLibrary{

    MaterialCode getMaterialCodeFromName(const std::string& materialCode)
    {
        static std::map<std::string, MaterialCode> materialCodeMap =
        {
            {"MODLEON", ModLeon},
            {"SHOTLEON", ShotLeon},
            {"MESCHKE", Meschke},
            {"SCHAEDLICHSCHWEIGER", SchaedlichSchweiger},
            {"MODLEONNONLOCAL", ModLeonNonLocal},
            {"HOEKBROWN", HoekBrown},
            {"UNTEREGGERROCKMASS", UntereggerRockMass},
            {"MOHRCOULOMB", MohrCoulomb},
            {"UNTEREGGERROCKMASSNONLOCAL", UntereggerRockMassNonLocal},
            {"SHOTLEONNONLOCAL", ShotLeonNonLocal},
            {"SHOTLEONV2", ShotLeonV2},
            {"LINEARELASTIC", LinearElastic},
            {"LINEARELASTICSOLIDIFICATIONCREEP", LinearElasticSolidificationCreep},
            {"MODLEONADAPTIVE", ModLeonAdaptive},
            {"MODLEONSEMIEXPLICIT", ModLeonSemiExplicit},
            {"MODLEONSEMIEXPLICITADAPTIVE", ModLeonSemiExplicitAdaptive},
            {"MODLEONPLANESTRESS", ModLeonPlaneStress},
            {"SHOTLEONV2NONLOCAL", ShotLeonV2NonLocal},
            {"UNTEREGGERROCKMASSPLAXIS", UntereggerRockMassPlaxis},
            {"LINEARELASTICNONLOCAL", LinearElasticNonLocal},
        };

        return materialCodeMap [ materialCode ];
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
            #ifdef HOEKBROWN 
            case HoekBrown: { return new class HoekBrown(stateVars, nStateVars, materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef MODLEON
            case ModLeon: { return new class ModLeon(stateVars, nStateVars, materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef SHOTLEON
            case ShotLeon: { return new class ShotLeon(stateVars, nStateVars, materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef SHOTLEONV2
            case ShotLeonV2: { return new class ShotLeonV2(stateVars, nStateVars, materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef SHOTLEONNONLOCAL
            case ShotLeonNonLocal: { return new class ShotLeonNonLocal(stateVars, nStateVars, materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef SHOTLEONV2NONLOCAL
            case ShotLeonV2NonLocal: { return new class ShotLeonV2NonLocal(stateVars, nStateVars, materialProperties, nMaterialProperties, element, gaussPt);}
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
            #ifdef UNTEREGGERROCKMASS
            case UntereggerRockMass: { return new class UntereggerRockMass(stateVars, nStateVars, materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef UNTEREGGERROCKMASSNONLOCAL
            case UntereggerRockMassNonLocal: { return new class UntereggerRockMassNonLocal(stateVars, nStateVars, materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef LINEARELASTICNONLOCAL
            case LinearElasticNonLocal: { return new class LinearElasticNonLocal(stateVars, nStateVars, materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef UNTEREGGERROCKMASSPLAXIS
            case UntereggerRockMassPlaxis: { return new class UntereggerRockMassPlaxis(stateVars, nStateVars, materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            default: std::cout << " MaterialCode " << materialCode << std::endl; throw std::invalid_argument("bftUserLibrary: Invalid Material Code Requested!");
        }
    }

} 

#ifdef UELDISPLACEMENT
    #include "uelDisplacementFactory.h"
#endif
#ifdef UELNONLOCAL
    #include "uelNonLocalFactory.h"
#endif
#ifdef UELDISPLACEMENTEAS
    #include "uelDisplacementEASFactory.h"
#endif
#ifdef UELNONLOCALEAS
    #include "uelNonLocalEASFactory.h"
#endif

namespace userLibrary{

    ElementCode getElementCodeFromName(const std::string& elementName)
    {
        static std::map<std::string, ElementCode> elementCodeMap =
        {
            { "UelT2D2", UelT2D2 }, 
            { "UelCPS4", UelCPS4 }, 
            { "UelCPS8", UelCPS8 }, 
            { "UelCPS8R", UelCPS8R }, 
            { "UelCPE4", UelCPE4 }, 
            { "UelCPE8", UelCPE8 }, 
            { "UelCPE8R", UelCPE8R }, 
            { "UelCPE4EAS2", UelCPE4EAS2 }, 
            { "UelCPE4EAS4", UelCPE4EAS4 }, 
            { "UelCPE4EAS5", UelCPE4EAS5 }, 
            { "UelC3D8", UelC3D8 }, 
            { "UelC3D8R", UelC3D8R }, 
            { "UelC3D20", UelC3D20 }, 
            { "UelC3D20R", UelC3D20R }, 
            { "UelC3D8EAS3", UelC3D8EAS3 }, 
            { "UelC3D8EAS9", UelC3D8EAS9 }, 
            { "UelCPS4NonLocal", UelCPS4NonLocal }, 
            { "UelCPS8NonLocal", UelCPS8NonLocal }, 
            { "UelCPS8RNonLocal", UelCPS8RNonLocal }, 
            { "UelCPS4NonLocalEAS2", UelCPS4NonLocalEAS2 }, 
            { "UelCPS4NonLocalEAS4", UelCPS4NonLocalEAS4 }, 
            { "UelCPE4NonLocal", UelCPE4NonLocal }, 
            { "UelCPE4RNonLocal", UelCPE4RNonLocal }, 
            { "UelCPE8NonLocal", UelCPE8NonLocal }, 
            { "UelCPE8RNonLocal", UelCPE8RNonLocal }, 
            { "UelCPE4NonLocalEAS2", UelCPE4NonLocalEAS2 }, 
            { "UelCPE4NonLocalEAS4", UelCPE4NonLocalEAS4 }, 
            { "UelCPE4NonLocalEAS5", UelCPE4NonLocalEAS5 }, 
            { "UelC3D8NonLocalEAS3", UelC3D8NonLocalEAS3 }, 
            { "UelC3D8NonLocalEAS6b", UelC3D8NonLocalEAS6b }, 
            { "UelC3D8NonLocalEAS9", UelC3D8NonLocalEAS9 }, 
            { "UelC3D8NonLocal", UelC3D8NonLocal }, 
            { "UelC3D8RNonLocal", UelC3D8RNonLocal }, 
            { "UelC3D20NonLocal", UelC3D20NonLocal }, 
            { "UelC3D20RNonLocal", UelC3D20RNonLocal }, 

        };
        return elementCodeMap[ elementName ];
    }

    BftUel* UelFactory(ElementCode  elementCode, const double* elementCoordinates, double* stateVars, int nStateVars, 
            const double* propertiesElement, int nPropertiesElement, int elementNumber, 
            MaterialCode materialCode,
            int nStateVarsUmat, const double* propertiesUmat, int nPropertiesUmat)
    {
        switch(elementCode){
            #ifdef UELDISPLACEMENT
            case UelT2D2: {return UelDisplacementFactory:: generateUelT2D2 (elementCoordinates,  stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, materialCode , nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case UelCPS4: {return UelDisplacementFactory:: generateUelCPS4 (elementCoordinates,  stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, materialCode , nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case UelCPE4: {return UelDisplacementFactory:: generateUelCPE4 (elementCoordinates,  stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case UelCPS8: {return UelDisplacementFactory:: generateUelCPS8 (elementCoordinates,  stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case UelCPS8R: {return UelDisplacementFactory:: generateUelCPS8R (elementCoordinates,  stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case UelC3D8: {return UelDisplacementFactory:: generateUelC3D8 (elementCoordinates,  stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case UelC3D8R: {return UelDisplacementFactory:: generateUelC3D8R (elementCoordinates,  stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case UelCPE8R: {return UelDisplacementFactory:: generateUelCPE8R (elementCoordinates,  stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case UelCPE8: {return UelDisplacementFactory:: generateUelCPE8 (elementCoordinates,  stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case UelC3D20: {return UelDisplacementFactory:: generateUelC3D20 (elementCoordinates,  stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case UelC3D20R: {return UelDisplacementFactory:: generateUelC3D20R (elementCoordinates,  stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            #endif 
            #ifdef UELNONLOCAL
            case UelCPS4NonLocal: {return UelNonLocalFactory:: generateUelCPS4NonLocal (elementCoordinates,  stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case UelCPE4NonLocal: {return UelNonLocalFactory:: generateUelCPE4NonLocal (elementCoordinates,  stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case UelCPE4RNonLocal: {return UelNonLocalFactory:: generateUelCPE4RNonLocal (elementCoordinates,  stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case UelCPS8NonLocal: {return UelNonLocalFactory:: generateUelCPS8NonLocal (elementCoordinates,  stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case UelCPS8RNonLocal: {return UelNonLocalFactory:: generateUelCPS8RNonLocal (elementCoordinates,  stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case UelC3D8NonLocal: {return UelNonLocalFactory:: generateUelC3D8NonLocal (elementCoordinates,  stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case UelC3D8RNonLocal: {return UelNonLocalFactory:: generateUelC3D8RNonLocal (elementCoordinates,  stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case UelCPE8NonLocal: {return UelNonLocalFactory:: generateUelCPE8NonLocal (elementCoordinates,  stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case UelCPE8RNonLocal: {return UelNonLocalFactory:: generateUelCPE8RNonLocal (elementCoordinates,  stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case UelC3D20RNonLocal: {return UelNonLocalFactory:: generateUelC3D20RNonLocal (elementCoordinates,  stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            #endif 
            #ifdef UELDISPLACEMENTEAS
            case UelCPE4EAS2: {return UelDisplacementEASFactory:: generateUelCPE4EAS2 (elementCoordinates,  stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case UelCPE4EAS4: {return UelDisplacementEASFactory:: generateUelCPE4EAS4 (elementCoordinates,  stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case UelCPE4EAS5: {return UelDisplacementEASFactory:: generateUelCPE4EAS5 (elementCoordinates,  stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case UelC3D8EAS3: {return UelDisplacementEASFactory:: generateUelC3D8EAS3 (elementCoordinates,  stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case UelC3D8EAS9: {return UelDisplacementEASFactory:: generateUelC3D8EAS9 (elementCoordinates,  stateVars, nStateVars, propertiesElement, 
                             nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            #endif
            #ifdef UELNONLOCALEAS 
            case UelCPE4NonLocalEAS2: {return UelNonLocalEASFactory:: generateUelCPE4NonLocalEAS2 (elementCoordinates,  stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case UelCPE4NonLocalEAS4: {return UelNonLocalEASFactory:: generateUelCPE4NonLocalEAS4 (elementCoordinates,  stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case UelCPE4NonLocalEAS5: {return UelNonLocalEASFactory:: generateUelCPE4NonLocalEAS4 (elementCoordinates,  stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}

            case UelCPS4NonLocalEAS4: {return UelNonLocalEASFactory:: generateUelCPS4NonLocalEAS4 (elementCoordinates,  stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case UelC3D8NonLocalEAS3: {return UelNonLocalEASFactory:: generateUelC3D8NonLocalEAS3(elementCoordinates,  stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case UelC3D8NonLocalEAS6b: {return UelNonLocalEASFactory:: generateUelC3D8NonLocalEAS6b(elementCoordinates,  stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            case UelC3D8NonLocalEAS9: {return UelNonLocalEASFactory:: generateUelC3D8NonLocalEAS9(elementCoordinates,  stateVars, nStateVars, propertiesElement, 
                            nPropertiesElement, elementNumber, materialCode, nStateVarsUmat, propertiesUmat, nPropertiesUmat);}
            #endif

            default:{throw std::invalid_argument("bftUserLibrary: Element Not Found");}
        }
    }

    ElementCode getElementCodeFromInfo(int nNodes, 
                                       int activeFields, 
                                       int dim,
                                       bool fullIntegration, 
                                       const std::string& formulation,
                                       int EASParam){
           int elInfo[4] = {0};
           
           elInfo[0] = nNodes;
           elInfo[1] = activeFields - 1;
           
           switch(dim)
               {
                case 1:
                    elInfo[2] = fullIntegration ? 1 : 4; break; 
                case 2:
                    {
                    if (formulation == "PLANESTRESS")
                        elInfo[2] = (fullIntegration) ? 2 : 5; 
                    else if (formulation == "PLANESTRAIN")
                        elInfo[2] = (fullIntegration) ? 7 : 8; 
                    else
                        throw std::invalid_argument("bftUserLibrary: formulation type cannot be identified for ElementCode");
                    break; }
               case 3:
                    elInfo[2] = (fullIntegration) ? 3 : 6; break;
               
               default:{throw std::invalid_argument("bftUserLibrary: cannot generate ElementCode from Info");}
               }
           
           elInfo[3] = (EASParam>0) ? EASParam : 0;

           // concatenate to elementcode as in defined list elementCode
           int elCode;
           if (elInfo[3]>0)
               elCode = elInfo[3] + elInfo[2] * 100 + elInfo[1] * 1000 + elInfo[0] * 10000;
           else
               elCode = elInfo[2] + elInfo[1] * 10 + elInfo[0] * 100;
   
           return static_cast<ElementCode>(elCode); 
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
