#pragma once 
#include <string>
#include <map>
#include "bftUel.h"
#include "bftMaterial.h"

namespace userLibrary{

    enum MaterialCode{
        ModLeon=1,
        ShotLeon=2,
        Meschke=3,
        SchaedlichSchweiger=4,
        ModLeonNonLocal=5,
        HoekBrown=6,
        UntereggerRockMass=7,
        MohrCoulomb=8,
        UntereggerRockMassNonLocal=9,
        ShotLeonNonLocal=10,
        ShotLeonV2=11,
        LinearElastic=12,
        LinearElasticSolidificationCreep = 13,
        ModLeonAdaptive=14,
        ModLeonSemiExplicit=15,
        ModLeonSemiExplicitAdaptive=16,
        ModLeonPlaneStress=17,
        UntereggerRockMassPlaxis=18,
        ShotLeonV2NonLocal=19,
    };

    std::map<std::string, MaterialCode> materialCodeMap =
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
    };

    MaterialCode getMaterialCodeFromName(const std::string& materialName);
/*
     * XXXXXX
     * ||||||_    6: if EAS: number of EAS Parameters
     * |||||__    5: if EAS: number of EAS Parameters
     * ||||___    4: type of element 
     * |||____    3: active fields 
     * ||_____    2: number of nodes
     * |______    1: (number of nodes)
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
     * */
    enum ElementCode{
            
            //Displacement
            //Plane Stress
            UelCPS4 =       402,
            UelCPS8 =       802,
            UelCPS8R =      805,

            // Plane Strain
            UelCPE4 =       407,
            UelCPE8 =       807,
            UelCPE8R =       808,

            // Plane Strain - EAS
            UelCPE4EAS2 = 40702,
            UelCPE4EAS4 = 40704,
            UelCPE4EAS5 = 40705,
            
            // Solid
            UelC3D8 =       803,
            UelC3D8R =      806,
            UelC3D20 =      2003,
            UelC3D20R =     2006,

            // Solid EAS
            UelC3D8EAS9 =   80309,
                     
            // Nonlocal 
            // Plane Stress
            
            UelCPS4NonLocal = 412,
            UelCPS8NonLocal = 812,
            UelCPS8RNonLocal = 815,

            // Plane Stress - EAS
            UelCPS4NonLocalEAS2 =41202,
            UelCPS4NonLocalEAS4 =41204,

            // Plane Strain
            UelCPE4NonLocal =417,
            UelCPE4RNonLocal = 418,
            UelCPE8NonLocal =817,
            UelCPE8RNonLocal = 818,

            // Plane Strain - EAS
            UelCPE4NonLocalEAS2 =41702,
            UelCPE4NonLocalEAS4 =41704,

            // Solid
            UelC3D8NonLocal = 813,
            UelC3D8RNonLocal = 816,
            UelC3D20NonLocal = 2013,
            UelC3D20RNonLocal = 2016,
    };

    std::map<std::string, ElementCode> elementCodeMap =
    {
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
        { "UelC3D8NonLocal", UelC3D8NonLocal }, 
        { "UelC3D8RNonLocal", UelC3D8RNonLocal }, 
        { "UelC3D20NonLocal", UelC3D20NonLocal }, 
        { "UelC3D20RNonLocal", UelC3D20RNonLocal }, };

    ElementCode getElementCodeFromName(const std::string& elementName);

    BftUel* UelFactory( ElementCode elementCode,
                        const double* elementCoordinates, 
                        double* stateVars, int nStateVars, 
                        const double* propertiesElement, 
                        int nPropertiesElement, 
                        int elementNumber, 
                        MaterialCode material,
                        int nStateVarsUmat, 
                        const double* propertiesUmat, 
                        int nPropertiesUmat);

    BftMaterial* bftMaterialFactory(MaterialCode material,
                                    double* stateVars,
                                    int nStateVars,
                                    const double* materialProperties, 
                                    int nMaterialProperties,
                                    int element, 
                                    int gaussPt);

    void extendAbaqusToVoigt(double* stress6, double* stress, double* strain6, const double* strain, double* dStrain6, const double* dStrain, int nDirect, int nShear);
    void backToAbaqus(double* stress, double* stress6, double* dStressDDStrain, double* dStressDDStrain66, int nDirect, int nShear);

    static const int sizeGeostaticDefinition = 6;
}

namespace Abaqus{
    
    enum UelFlags1{
        GeostaticStress=61,   // Geostatic stress field according to Abaqus Analysis User's Guide Tab. 5.1.2-1 Keys to procedure types.
    };
}


