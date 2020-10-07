#pragma once
#include "bftElement.h"
#include "bftMaterial.h"
#include <functional>
#include <iostream>
#include <map>
#include <string>

namespace userLibrary {

    enum MaterialCode : int {
        ModLeon                          = 1,
        ShotLeon                         = 2,
        ViscoPlasticShotcreteModel       = 3,
        SchaedlichSchweiger              = 4,
        ModLeonNonLocal                  = 5,
        HoekBrown                        = 6,
        RockDamagePlasticity             = 7,
        MohrCoulomb                      = 8,
        RockDamagePlasticityNonLocal     = 9,
        ShotLeonNonLocal                 = 10,
        ShotLeonV2                       = 11,
        LinearElastic                    = 12,
        LinearElasticSolidificationCreep = 13,
        ModLeonAdaptive                  = 14,
        ModLeonSemiExplicit              = 15,
        ModLeonSemiExplicitAdaptive      = 16,
        ModLeonPlaneStress               = 17,
        UntereggerRockMassPlaxis         = 18,
        ShotLeonV2NonLocal               = 19,
        UntereggerRockMassAssocNonLocal  = 20,
        LinearElasticNonLocal            = 21,
        StVenantKirchhoffIsotropic       = 22,
        ModLeonAnalytical                = 23,
        Barodesy                         = 24,
        BarodesyGradientVoid             = 25,
        CosseratLinearElastic            = 26,
        CosseratDruckerPrager            = 27,
        MCDPModel                        = 28,
        GMCDPModel                       = 29,
        CDPModel                         = 30,
        CDPFibreReinforcedModel          = 31,
        PorousElastic                    = 32,
        BarodesyGradientDeformationModulus = 33,
        ModifiedCamClay                  = 34,
        JointedRock                      = 35,
        TransverseIsotropicLinearElastic = 36,
        MisesTI                          = 37,
        CDPM2                            = 38,
        DruckerPrager                    = 39,
        SandHypo                         = 40,
        SandHypoMicropolar               = 41,
        GMBiotElastic                    = 42,
        GMDruckerPrager                  = 43,
        DruckerPragerMD                  = 44,
        GMNeoHooke                       = 45,
    };

    enum ElementCode {

        /*
         * XXXXXX
         * ||||||_    6: if EAS: number of EAS Parameters
         * |||||__    5: if EAS: number of EAS Parameters (01 = FBar / Bbar)
         * ||||___    4: type of element
         * |||____    3: active fields
         * ||_____    2: number of nodes
         * |______    1: (number of nodes)
         *
         *
         * active fields:   0: mechanical (=displacement),
         *                  1: mechanical + nonlocal damage,
         *                  2: mechanical (=displacement) large strain TL,
         *                  3: mechanical (=displacement) large strain UL,
         *                  4: mechanical + nonlocal damage, large strain UL,
         *                  5: cosserat 
         *                  6: cosserat   + nonlocal damage, 
         *                  7: mechanical + nonlocal damage mixed,
         *                  9: Reserved for tests
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

        // Displacement
        // Truss
        UelT2D2 = 202,
        // Plane Stress
        UelCPS4 = 402,
        // UelCPS8  = 802,
        UelCPS8R = 805,

        // Plane Strain
        UelCPE4 = 407,
        // UelCPE8  = 807,
        UelCPE8R = 808,

        // Plane Strain - EAS
        UelCPE4EAS2 = 40702,
        UelCPE4EAS4 = 40704,
        UelCPE4EAS5 = 40705,

        // Solid
        UelC3D8   = 803,
        UelC3D8R  = 806,
        UelC3D20  = 2003,
        UelC3D20R = 2006,

        UelT2D2TL = 222,
        // Plane Stress
        UelCPS4TL  = 422,
        UelCPS8TL  = 822,
        UelCPS8RTL = 825,

        // Plane Strain
        UelCPE4TL  = 427,
        UelCPE8TL  = 827,
        UelCPE8RTL = 828,

        // Solid
        UelC3D8TL   = 823,
        UelC3D8RTL  = 826,
        UelC3D20TL  = 2023,
        UelC3D20RTL = 2026,

        UelCPE4UL  = 437,
        UelCPE8RUL = 838,
        UelC3D8UL  = 833,

        // Solid EAS
        UelC3D8EAS3  = 80303,
        UelC3D8EAS6b = 80361,
        UelC3D8EAS9  = 80309,

        // Nonlocal
        // Plane Stress

        UelCPS4NonLocal  = 412,
        UelCPS8NonLocal  = 812,
        UelCPS8RNonLocal = 815,
        
        // Nonlocal mixed
        UelCPS8RNonLocalMixed  = 875,
        UelCPE8RNonLocalMixed  = 877,
        UelC3D20RNonLocalMixed = 2076,

        // Plane Stress - EAS
        UelCPS4NonLocalEAS2 = 41202,
        UelCPS4NonLocalEAS4 = 41204,
        UelCPS4NonLocalEAS5 = 41205,

        // Plane Strain
        UelCPE4NonLocal  = 417,
        UelCPE4RNonLocal = 418,
        UelCPE8NonLocal  = 817,
        UelCPE8RNonLocal = 818,

        // Plane Strain - EAS
        UelCPE4NonLocalEAS2 = 41702,
        UelCPE4NonLocalEAS4 = 41704,
        UelCPE4NonLocalEAS5 = 41705,

        UelC3D8NonLocalEAS3  = 81303,
        UelC3D8NonLocalEAS6b = 81361,
        UelC3D8NonLocalEAS9  = 81309,

        // Solid
        UelC3D4NonLocal   = 413,
        UelC3D10NonLocal   = 1013,
        UelC3D8NonLocal   = 813,
        UelC3D8RNonLocal  = 816,
        UelC3D20NonLocal  = 2013,
        UelC3D20RNonLocal = 2016,

        // Nonlocal, Updated Lagrange
        UelC3D8NonLocalUL  = 843,
        UelCPE4NonLocalUL  = 447,
        UelCPE4RNonLocalUL = 448,
        UelCPE8RNonLocalUL = 848,
        // FBar versions
        UelC3D8NonLocalULFBar = 84301,
        UelCPE4NonLocalULFBar = 44701,

        // Cosserat
        UelCCPE4  = 458,
        UelCCPE8R  = 858,
        UelCC3D20R = 2056,
        UelCC3D8 = 853,

        // Nonlocal Cosserat
        UelNCCPS4  = 465,
        UelNCCPE8R  = 868,
        UelNCCPS8R  = 865,
        UelNCC3D20R = 2066,
        UelNCC3D8 = 866,

        UelGMCPE8R  = 898,
        UelGMC3D20R = 2096,
        UelGMC3D8  ,
        UelGMCPE8RUL ,
        UelGMC3D8UL ,
        UelGMC3D20RUL ,
    };

    // MaterialFactory
    //
    // - Allows materials to register themselve with their name and ID
    // - Allows the user to create instances of materials

    class BftMaterialFactory {
      public:
        using materialFactoryFunction = BftMaterial* (*)( const double* materialProperties,
                                                          int           nMaterialProperties,
                                                          int           element,
                                                          int           gaussPt );
        BftMaterialFactory()          = delete;

        static MaterialCode getMaterialCodeFromName( const std::string& materialName );

        static BftMaterial* createMaterial( MaterialCode  material,
                                            const double* materialProperties,
                                            int           nMaterialProperties,
                                            int           element,
                                            int           gaussPt );

        static bool registerMaterial( MaterialCode            materialCode,
                                      const std::string&      materialName,
                                      materialFactoryFunction factoryFunction );

      private:
        static std::map<std::string, MaterialCode>             materialNameToCodeAssociation;
        static std::map<MaterialCode, materialFactoryFunction> materialFactoryFunctionByCode;
    };

    // ElementFactory
    //
    // - Allows elements to register themselve with their name and ID
    // - Allows the user to create instances of elements

    class BftElementFactory {
      public:
        using elementFactoryFunction = BftElement* (*)( int elementNumber );
        BftElementFactory()          = delete;

        static ElementCode getElementCodeFromName( const std::string& elementName );

        static BftElement* createElement( ElementCode elementCode, int elementNumber );

        static bool registerElement( const std::string&     elementName,
                                     ElementCode            elementCode,
                                     elementFactoryFunction factoryFunction );

      private:
        static std::map<std::string, ElementCode>            elementNameToCodeAssociation;
        static std::map<ElementCode, elementFactoryFunction> elementFactoryFunctionByCode;
    };

} // namespace userLibrary
