#include "bftElement.h"
#include "bftMaterial.h"
#include "bftTypedefs.h"
#include "bftVoigt.h"
#include "userLibrary.h"
#include <map>
#include <string>
#include <tuple>

#ifdef BARODESY
#    include "Barodesy.h"
#endif
#ifdef BARODESYGRADIENTVOID
#    include "BarodesyGradientVoid.h"
#endif
#ifdef BARODESYGRADIENTDEFORMATIONMODULUS
#    include "BarodesyGradientDeformationModulus.h"
#endif
#ifdef COSSERATLINEARELASTIC
#    include "CosseratLinearElastic.h"
#endif
#ifdef COSSERATDRUCKERPRAGER
#    include "CosseratDruckerPrager.h"
#endif
#ifdef GMCDP 
#    include "GMCDP.h"
#endif
#ifdef LINEARELASTIC
#    include "LinearElastic.h"
#endif
#ifdef LINEARELASTICSOLIDIFICATIONCREEP
#    include "LESolidificationCreep.h"
#endif
#ifdef MCDP 
#    include "MCDP.h"
#endif
#ifdef MODLEON
#    include "ModLeon.h"
#endif
#ifdef MODLEONSEMIEXPLICIT
#    include "ModLeonSemiExplicit.h"
#endif
#ifdef MODLEONADAPTIVE
#    include "ModLeonAdaptive.h"
#endif
#ifdef MODLEONSEMIEXPLICITADAPTIVE
#    include "ModLeonSemiExplicitAdaptive.h"
#endif
#ifdef MODLEONPLANESTRESS
#    include "ModLeonPS.h"
#endif
#ifdef POROUSELASTIC
#    include "PorousElastic.h"
#endif
#ifdef SHOTLEON
#    include "ShotLeon.h"
#endif
#ifdef SHOTLEONV2
#    include "ShotLeonV2.h"
#endif
#ifdef SHOTLEONNONLOCAL
#    include "ShotLeonNonLocal.h"
#endif
#ifdef SHOTLEONV2NONLOCAL
#    include "ShotLeonV2NonLocal.h"
#endif
#ifdef MODIFIEDCAMCLAY
#    include "ModifiedCamClay.h"
#endif
#ifdef MODLEONNONLOCAL
#    include "ModLeonNonLocal.h"
#endif
#ifdef MESCHKE
#    include "Meschke.h"
#endif
#ifdef SCHAEDLICHSCHWEIGER
#    include "SchaedlichSchweiger.h"
#endif
#ifdef HOEKBROWN
#    include "HoekBrown.h"
#endif
#ifdef ROCKDAMAGEPLASTICITY
#    include "RockDamagePlasticity.h"
#endif
#ifdef UNTEREGGERROCKMASSPLAXIS
#    include "UntereggerRockMassPlaxis.h"
#endif
#ifdef MOHRCOULOMB
#    include "materialCodeMohrCoulomb.h"
#endif
#ifdef ROCKDAMAGEPLASTICITYNONLOCAL
#    include "RockDamagePlasticityNonLocal.h"
#endif
#ifdef LINEARELASTICNONLOCAL
#    include "LinearElasticNonLocal.h"
#endif
#ifdef STVENANTKIRCHHOFFISOTROPIC
#    include "StVenantKirchhoffIsotropic.h"
#endif

using namespace bft;
using namespace Eigen;

namespace userLibrary {

    MaterialCode getMaterialCodeFromName( const std::string& materialCode )
    {
        static std::map<std::string, MaterialCode> materialCodeMap = {
            {"BARODESY", Barodesy},
            {"BARODESYGRADIENTVOID", BarodesyGradientVoid},
            {"BARODESYGRADIENTDEFORMATIONMODULUS", BarodesyGradientDeformationModulus},
            {"COSSERATLINEARELASTIC", CosseratLinearElastic},
            {"COSSERATDRUCKERPRAGER", CosseratDruckerPrager},
            {"GMCDP", GMCDPModel},
            {"MODLEON", ModLeon},
            {"SHOTLEON", ShotLeon},
            {"MESCHKE", Meschke},
            {"SCHAEDLICHSCHWEIGER", SchaedlichSchweiger},
            {"MODLEONNONLOCAL", ModLeonNonLocal},
            {"MCDP", MCDPModel},
            {"MODIFIEDCAMCLAY", ModifiedCamClay},
            {"HOEKBROWN", HoekBrown},
            {"ROCKDAMAGEPLASTICITY", RockDamagePlasticity},
            {"MOHRCOULOMB", MohrCoulomb},
            {"ROCKDAMAGEPLASTICITYNONLOCAL", RockDamagePlasticityNonLocal},
            {"SHOTLEONNONLOCAL", ShotLeonNonLocal},
            {"SHOTLEONV2", ShotLeonV2},
            {"LINEARELASTIC", LinearElastic},
            {"LINEARELASTICSOLIDIFICATIONCREEP", LinearElasticSolidificationCreep},
            {"MODLEONADAPTIVE", ModLeonAdaptive},
            {"MODLEONSEMIEXPLICIT", ModLeonSemiExplicit},
            {"MODLEONSEMIEXPLICITADAPTIVE", ModLeonSemiExplicitAdaptive},
            {"MODLEONPLANESTRESS", ModLeonPlaneStress},
            {"POROUSELASTIC", PorousElastic},
            {"SHOTLEONV2NONLOCAL", ShotLeonV2NonLocal},
            {"UNTEREGGERROCKMASSPLAXIS", UntereggerRockMassPlaxis},
            {"LINEARELASTICNONLOCAL", LinearElasticNonLocal},
            {"STVENANTKIRCHHOFFISOTROPIC", StVenantKirchhoffIsotropic},
        };

        return materialCodeMap[materialCode];
    }

    BftMaterial* bftMaterialFactory( MaterialCode  materialCode,
                                     const double* materialProperties,
                                     int           nMaterialProperties,
                                     int           element,
                                     int           gaussPt )
    {
        switch ( materialCode ) {
// clang-format off
            #ifdef BARODESY 
            case Barodesy: { return new class Barodesy(materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef BARODESYGRADIENTVOID
            case BarodesyGradientVoid: { return new class BarodesyGradientVoid(materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef BARODESYGRADIENTDEFORMATIONMODULUS
            case BarodesyGradientDeformationModulus: { return new class BarodesyGradientDeformationModulus(materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef COSSERATLINEARELASTIC
            case CosseratLinearElastic: { return new class CosseratLinearElastic(materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef COSSERATDRUCKERPRAGER 
            case CosseratDruckerPrager: { return new class CosseratDruckerPrager(materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef LINEARELASTIC
            case LinearElastic: { return new class LinearElastic(materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef GMCDP
            case GMCDPModel: { return new class GMCDPModel( materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef HOEKBROWN 
            case HoekBrown: { return new class HoekBrown(materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef MCDP
            case MCDPModel: { return new class MCDPModel( materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef MODLEON
            case ModLeon: { return new class ModLeon( materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef POROUSELASTIC
            case PorousElastic: { return new class PorousElastic( materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef MODIFIEDCAMCLAY
            case ModifiedCamClay: { return new class ModifiedCamClay( materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef SHOTLEON
            case ShotLeon: { return new class ShotLeon(materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef SHOTLEONV2
            case ShotLeonV2: { return new class ShotLeonV2(materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef SHOTLEONNONLOCAL
            case ShotLeonNonLocal: { return new class ShotLeonNonLocal(materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef SHOTLEONV2NONLOCAL
            case ShotLeonV2NonLocal: { return new class ShotLeonV2NonLocal(materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef MODLEONSEMIEXPLICIT
            case ModLeonSemiExplicit: { return new class ModLeonSemiExplicit(materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef MODLEONADAPTIVE
            case ModLeonAdaptive: { return new class ModLeonAdaptive(materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef MODLEONSEMIEXPLICITADAPTIVE
            case ModLeonSemiExplicitAdaptive: { return new class ModLeonSemiExplicitAdaptive(materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef MODLEONNONLOCAL
            case ModLeonNonLocal: { return new class ModLeonNonLocal(materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef MESCHKE
            case Meschke: { return new class Meschke(materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef SCHAEDLICHSCHWEIGER
            case SchaedlichSchweiger: { return new class SchaedlichSchweiger(materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef ROCKDAMAGEPLASTICITY
            case RockDamagePlasticity: { return new class RockDamagePlasticity(materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef ROCKDAMAGEPLASTICITYNONLOCAL
            case RockDamagePlasticityNonLocal: { return new class RockDamagePlasticityNonLocal(materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef LINEARELASTICNONLOCAL
            case LinearElasticNonLocal: { return new class LinearElasticNonLocal(materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef UNTEREGGERROCKMASSPLAXIS
            case UntereggerRockMassPlaxis: { return new class UntereggerRockMassPlaxis(materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            #ifdef STVENANTKIRCHHOFFISOTROPIC
            case StVenantKirchhoffIsotropic: { return new class StVenantKirchhoffIsotropic(materialProperties, nMaterialProperties, element, gaussPt);}
            #endif
            default: {  std::ostringstream str; str<<"bftUserLibrary: Invalid material code "<< materialCode << " requested!" << std::endl; 
                         throw std::invalid_argument(str.str()); }
            // clang-format on
        }
    }

} // namespace userLibrary

#ifdef UELDISPLACEMENT
#    include "uelDisplacementFactory.h"
#endif
#ifdef UELDISPLACEMENTTL
#    include "uelDisplacementTLFactory.h"
#endif
#ifdef UELDISPLACEMENTUL
#    include "uelDisplacementULFactory.h"
#endif
#ifdef UELNONLOCAL
#    include "uelNonLocalFactory.h"
#endif
#ifdef UELDISPLACEMENTEAS
#    include "uelDisplacementEASFactory.h"
#endif
#ifdef UELNONLOCALEAS
#    include "uelNonLocalEASFactory.h"
#endif
#ifdef UELNONLOCALMIXED
#    include "uelNonLocalMixedFactory.h"
#endif
#ifdef UELNONLOCALUL
#    include "uelNonLocalULFactory.h"
#endif
#ifdef UELNONLOCALULFBAR
#    include "uelNonLocalULFBarFactory.h"
#endif
#ifdef UELCOSSERAT
#    include "uelCosseratFactory.h"
#endif
#ifdef UELNONLOCALCOSSERAT
#    include "uelNonLocalCosseratFactory.h"
#endif

namespace userLibrary {

    ElementCode getElementCodeFromName( const std::string& elementName )
    {
        static std::map<std::string, ElementCode> elementCodeMap = {
            {"UELT2D2", UelT2D2},
            {"UELCPS4", UelCPS4},
            {"UELCPS8", UelCPS8},
            {"UELCPS8R", UelCPS8R},
            {"UELCPE4", UelCPE4},
            {"UELCPE8", UelCPE8},
            {"UELCPE8R", UelCPE8R},
            {"UELC3D8", UelC3D8},
            {"UELC3D8R", UelC3D8R},
            {"UELC3D20", UelC3D20},
            {"UELC3D20R", UelC3D20R},
            {"UELT2D2TL", UelT2D2TL},
            {"UELCPS4TL", UelCPS4TL},
            {"UELCPS8TL", UelCPS8TL},
            {"UELCPS8RTL", UelCPS8RTL},
            {"UELCPE4TL", UelCPE4TL},
            {"UELCPE8TL", UelCPE8TL},
            {"UELCPE8RTL", UelCPE8RTL},
            {"UELC3D8TL", UelC3D8TL},
            {"UELC3D8RTL", UelC3D8RTL},
            {"UELC3D20TL", UelC3D20TL},
            {"UELC3D20RTL", UelC3D20RTL},
            {"UELCPE4UL", UelCPE4UL},
            {"UELCPE8RUL", UelCPE8RUL},
            {"UELC3D8UL", UelC3D8UL},
            {"UELCPE4EAS2", UelCPE4EAS2},
            {"UELCPE4EAS4", UelCPE4EAS4},
            {"UELCPE4EAS5", UelCPE4EAS5},
            {"UELC3D8EAS3", UelC3D8EAS3},
            {"UELC3D8EAS9", UelC3D8EAS9},
            {"UELCPS4NONLOCAL", UelCPS4NonLocal},
            {"UELCPS8NONLOCAL", UelCPS8NonLocal},
            {"UELCPS8RNONLOCAL", UelCPS8RNonLocal},
            {"UELCPS4NONLOCALEAS2", UelCPS4NonLocalEAS2},
            {"UELCPS4NONLOCALEAS4", UelCPS4NonLocalEAS4},
            {"UELCPE4NONLOCAL", UelCPE4NonLocal},
            {"UELCPE4RNONLOCAL", UelCPE4RNonLocal},
            {"UELCPE8NONLOCAL", UelCPE8NonLocal},
            {"UELCPE8RNONLOCAL", UelCPE8RNonLocal},
            {"UELCPE4NONLOCALEAS2", UelCPE4NonLocalEAS2},
            {"UELCPE4NONLOCALEAS4", UelCPE4NonLocalEAS4},
            {"UELCPE4NONLOCALEAS5", UelCPE4NonLocalEAS5},
            {"UELC3D8NONLOCALEAS3", UelC3D8NonLocalEAS3},
            {"UELC3D8NONLOCALEAS6B", UelC3D8NonLocalEAS6b},
            {"UELC3D8NONLOCALEAS9", UelC3D8NonLocalEAS9},
            {"UELC3D4NONLOCAL", UelC3D4NonLocal},
            {"UELC3D8NONLOCAL", UelC3D8NonLocal},
            {"UELC3D8RNONLOCAL", UelC3D8RNonLocal},
            {"UELC3D20NONLOCAL", UelC3D20NonLocal},
            {"UELC3D20RNONLOCAL", UelC3D20RNonLocal},
            //[>NonLocal, large strain<]
            {"UELC3D8NONLOCALUL", UelC3D8NonLocalUL},
            {"UELCPE4NONLOCALUL", UelCPE4NonLocalUL},
            {"UELCPE4RNONLOCALUL", UelCPE4RNonLocalUL},
            {"UELCPE8RNONLOCALUL", UelCPE8RNonLocalUL},
            //[>NonLocal, large strain, FBar<]
            {"UELC3D8NONLOCALULFBAR", UelC3D8NonLocalULFBar},
            {"UELCPE4NONLOCALULFBAR", UelCPE4NonLocalULFBar},
            // Cosserat
            {"UELCCPE4", UelCCPE4},
            {"UELCCPE8R", UelCCPE8R},
            {"UELCC3D20R", UelCC3D20R},
            //NonLocal Cosserat
            //{"UELNCCPE4", UelNCCPE4},
            {"UELNCCPE8R", UelNCCPE8R},
            {"UELNCCPS8R", UelNCCPS8R},
            {"UELNCC3D8R", UelNCC3D8R},
            {"UELCPS8RNONLOCALMIXED", UelCPS8RNonLocalMixed},
            {"UELCPE8RNONLOCALMIXED", UelCPE8RNonLocalMixed},
            {"UELC3D20RNONLOCALMIXED", UelC3D20RNonLocalMixed},

        };
        return elementCodeMap[elementName];
    }

    BftElement* bftElementFactory( ElementCode elementCode, int elementNumber )
    {
        switch ( elementCode ) {
// clang-format off
            #ifdef UELDISPLACEMENT
            case UelT2D2: {return UelDisplacementFactory:: generateUelT2D2(elementNumber );}
            case UelCPS4: {return UelDisplacementFactory:: generateUelCPS4(elementNumber);}
            case UelCPE4: {return UelDisplacementFactory:: generateUelCPE4(elementNumber);}
            case UelCPS8: {return UelDisplacementFactory:: generateUelCPS8(elementNumber);}
            case UelCPS8R: {return UelDisplacementFactory:: generateUelCPS8R(elementNumber);}
            case UelC3D8: {return UelDisplacementFactory:: generateUelC3D8(elementNumber);}
            case UelC3D8R: {return UelDisplacementFactory:: generateUelC3D8R(elementNumber);}
            case UelCPE8R: {return UelDisplacementFactory:: generateUelCPE8R(elementNumber);}
            case UelCPE8: {return UelDisplacementFactory:: generateUelCPE8(elementNumber);}
            case UelC3D20: {return UelDisplacementFactory:: generateUelC3D20(elementNumber);}
            case UelC3D20R: {return UelDisplacementFactory:: generateUelC3D20R(elementNumber);}
            #endif 
            #ifdef UELDISPLACEMENTTL
            //case UelT2D2: {return UelDisplacementFactory:: generateUelT2D2(elementNumber );}
            //case UelCPS4: {return UelDisplacementFactory:: generateUelCPS4(elementNumber);}
            case UelCPE4TL: {return UelDisplacementTLFactory:: generateUelCPE4TL(elementNumber);}
            //case UelCPS8: {return UelDisplacementFactory:: generateUelCPS8(elementNumber);}
            //case UelCPS8R: {return UelDisplacementFactory:: generateUelCPS8R(elementNumber);}
            case UelC3D8TL: {return UelDisplacementTLFactory:: generateUelC3D8TL(elementNumber);}
            //case UelC3D8R: {return UelDisplacementFactory:: generateUelC3D8R(elementNumber);}
            //case UelCPE8R: {return UelDisplacementFactory:: generateUelCPE8R(elementNumber);}
            //case UelCPE8: {return UelDisplacementFactory:: generateUelCPE8(elementNumber);}
            case UelC3D20TL: {return UelDisplacementFactory:: generateUelC3D20(elementNumber);}
            //case UelC3D20R: {return UelDisplacementFactory:: generateUelC3D20R(elementNumber);}
            #endif 
            #ifdef UELDISPLACEMENTUL
            case UelCPE4UL: {return UelDisplacementULFactory:: generateUelCPE4UL(elementNumber);}
            //case UelCPE8RUL: {return UelDisplacementULFactory:: generateUelCPE8RUL(elementNumber);}
            case UelC3D8UL: {return UelDisplacementULFactory:: generateUelC3D8UL(elementNumber);}
            #endif 
            #ifdef UELNONLOCAL
            case UelCPS4NonLocal: {return UelNonLocalFactory:: generateUelCPS4NonLocal(elementNumber);}
            case UelCPE4NonLocal: {return UelNonLocalFactory:: generateUelCPE4NonLocal(elementNumber);}
            case UelCPE4RNonLocal: {return UelNonLocalFactory:: generateUelCPE4RNonLocal(elementNumber);}
            case UelCPS8NonLocal: {return UelNonLocalFactory:: generateUelCPS8NonLocal(elementNumber);}
            case UelCPS8RNonLocal: {return UelNonLocalFactory:: generateUelCPS8RNonLocal(elementNumber);}
            case UelC3D4NonLocal: {return UelNonLocalFactory:: generateUelC3D4NonLocal(elementNumber);}
            case UelC3D8NonLocal: {return UelNonLocalFactory:: generateUelC3D8NonLocal(elementNumber);}
            case UelC3D8RNonLocal: {return UelNonLocalFactory:: generateUelC3D8RNonLocal(elementNumber);}
            case UelCPE8NonLocal: {return UelNonLocalFactory:: generateUelCPE8NonLocal(elementNumber);}
            case UelCPE8RNonLocal: {return UelNonLocalFactory:: generateUelCPE8RNonLocal(elementNumber);}
            case UelC3D20RNonLocal: {return UelNonLocalFactory:: generateUelC3D20RNonLocal(elementNumber);}
            case UelC3D20NonLocal: {return UelNonLocalFactory:: generateUelC3D20NonLocal(elementNumber);}
            #endif 
            #ifdef UELDISPLACEMENTEAS
            case UelCPE4EAS2: {return UelDisplacementEASFactory:: generateUelCPE4EAS2(elementNumber);}
            case UelCPE4EAS4: {return UelDisplacementEASFactory:: generateUelCPE4EAS4(elementNumber);}
            case UelCPE4EAS5: {return UelDisplacementEASFactory:: generateUelCPE4EAS5(elementNumber);}
            case UelC3D8EAS3: {return UelDisplacementEASFactory:: generateUelC3D8EAS3(elementNumber);}
            case UelC3D8EAS9: {return UelDisplacementEASFactory:: generateUelC3D8EAS9(elementNumber);}
            #endif
            #ifdef UELNONLOCALEAS 
            case UelCPE4NonLocalEAS2: {return UelNonLocalEASFactory:: generateUelCPE4NonLocalEAS2(elementNumber);}
            case UelCPE4NonLocalEAS4: {return UelNonLocalEASFactory:: generateUelCPE4NonLocalEAS4(elementNumber);}
            case UelCPE4NonLocalEAS5: {return UelNonLocalEASFactory:: generateUelCPE4NonLocalEAS4(elementNumber);}

            case UelCPS4NonLocalEAS4: {return UelNonLocalEASFactory:: generateUelCPS4NonLocalEAS4(elementNumber);}
            case UelC3D8NonLocalEAS3: {return UelNonLocalEASFactory:: generateUelC3D8NonLocalEAS3(elementNumber);}
            case UelC3D8NonLocalEAS6b: {return UelNonLocalEASFactory:: generateUelC3D8NonLocalEAS6b(elementNumber);}
            case UelC3D8NonLocalEAS9: {return UelNonLocalEASFactory:: generateUelC3D8NonLocalEAS9(elementNumber);}
            #endif
            #ifdef UELNONLOCALMIXED
            case UelCPS8RNonLocalMixed: {return UelNonLocalMixedFactory:: generateUelCPS8RNonLocalMixed(elementNumber);}
            case UelCPE8RNonLocalMixed: {return UelNonLocalMixedFactory:: generateUelCPE8RNonLocalMixed(elementNumber);}
            case UelC3D20RNonLocalMixed: {return UelNonLocalMixedFactory:: generateUelC3D20RNonLocalMixed(elementNumber);}
            #endif
            #ifdef UELNONLOCALUL
            case UelCPE4NonLocalUL: {return UelNonLocalULFactory:: generateUelCPE4NonLocalUL(elementNumber);}
            case UelCPE4RNonLocalUL: {return UelNonLocalULFactory:: generateUelCPE4RNonLocalUL(elementNumber);}
            case UelCPE8RNonLocalUL: {return UelNonLocalULFactory:: generateUelCPE8RNonLocalUL(elementNumber);}
            case UelC3D8NonLocalUL: {return UelNonLocalULFactory:: generateUelC3D8NonLocalUL(elementNumber);}
            #endif 
            #ifdef UELNONLOCALULFBAR
            case UelCPE4NonLocalULFBar: {return UelNonLocalULFBarFactory:: generateUelCPE4NonLocalULFBar(elementNumber);}
            case UelC3D8NonLocalULFBar: {return UelNonLocalULFBarFactory:: generateUelC3D8NonLocalULFBar(elementNumber);}
            #endif 
            #ifdef UELCOSSERAT
            case UelCCPE4 : {return UelCosseratFactory::generateUelCCPE4(elementNumber);}
            case UelCCPE8R : {return UelCosseratFactory::generateUelCCPE8R(elementNumber);}
            case UelCC3D20R : {return UelCosseratFactory::generateUelCC3D20R(elementNumber);}
            #endif
            #ifdef UELNONLOCALCOSSERAT
            //case UelNCCPE4 : {return UelNonLocalCosseratFactory::generateUelNCCPE4(elementNumber);}
            case UelNCCPE8R : {return UelNonLocalCosseratFactory::generateUelNCCPE8R(elementNumber);}
            case UelNCCPS8R : {return UelNonLocalCosseratFactory::generateUelNCCPS8R(elementNumber);}
            case UelNCC3D8R : {return UelNonLocalCosseratFactory::generateUelNCC3D8R(elementNumber);}
            #endif
        // clang-format on
        default: {
            std::ostringstream str;
            str << "bftUserLibrary: Invalid element " << elementCode << " requested!" << std::endl;
            throw std::invalid_argument( str.str() );
        }
        }
    }

    ElementCode getElementCodeFromInfo( int                nNodes,
                                        int                activeFields,
                                        int                dim,
                                        bool               fullIntegration,
                                        const std::string& formulation,
                                        int                EASParam )
    {
        int elInfo[4] = {0};

        elInfo[0] = nNodes;
        elInfo[1] = activeFields - 1;

        switch ( dim ) {
        case 1: elInfo[2] = fullIntegration ? 1 : 4; break;
        case 2: {
            if ( formulation == "PLANESTRESS" )
                elInfo[2] = ( fullIntegration ) ? 2 : 5;
            else if ( formulation == "PLANESTRAIN" )
                elInfo[2] = ( fullIntegration ) ? 7 : 8;
            else
                throw std::invalid_argument( "bftUserLibrary: formulation type cannot be identified for ElementCode" );
            break;
        }
        case 3: elInfo[2] = ( fullIntegration ) ? 3 : 6; break;

        default: {
            throw std::invalid_argument( "bftUserLibrary: cannot generate ElementCode from Info" );
        }
        }

        elInfo[3] = ( EASParam > 0 ) ? EASParam : 0;

        // concatenate to elementcode as in defined list elementCode
        int elCode;
        if ( elInfo[3] > 0 )
            elCode = elInfo[3] + elInfo[2] * 100 + elInfo[1] * 1000 + elInfo[0] * 10000;
        else
            elCode = elInfo[2] + elInfo[1] * 10 + elInfo[0] * 100;

        return static_cast<ElementCode>( elCode );
    }

} // namespace userLibrary
