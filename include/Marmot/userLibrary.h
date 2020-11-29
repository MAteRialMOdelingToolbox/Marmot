#pragma once
#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotMaterial.h"
#include <functional>
#include <iostream>
#include <map>
#include <string>

namespace userLibrary {

    enum MaterialCode : int {
        ModLeon                            = 1,
        ShotLeon                           = 2,
        ViscoPlasticShotcreteModel         = 3,
        SchaedlichSchweiger                = 4,
        ModLeonNonLocal                    = 5,
        HoekBrown                          = 6,
        RockDamagePlasticity               = 7,
        MohrCoulomb                        = 8,
        RockDamagePlasticityNonLocal       = 9,
        ShotLeonNonLocal                   = 10,
        ShotLeonV2                         = 11,
        LinearElastic                      = 12,
        LinearElasticSolidificationCreep   = 13,
        ModLeonAdaptive                    = 14,
        ModLeonSemiExplicit                = 15,
        ModLeonSemiExplicitAdaptive        = 16,
        ModLeonPlaneStress                 = 17,
        UntereggerRockMassPlaxis           = 18,
        ShotLeonV2NonLocal                 = 19,
        UntereggerRockMassAssocNonLocal    = 20,
        LinearElasticNonLocal              = 21,
        StVenantKirchhoffIsotropic         = 22,
        ModLeonAnalytical                  = 23,
        Barodesy                           = 24,
        BarodesyGradientVoid               = 25,
        CosseratLinearElastic              = 26,
        CosseratDruckerPrager              = 27,
        MCDPModel                          = 28,
        GMCDPModel                         = 29,
        CDPModel                           = 30,
        CDPFibreReinforcedModel            = 31,
        PorousElastic                      = 32,
        BarodesyGradientDeformationModulus = 33,
        ModifiedCamClay                    = 34,
        JointedRock                        = 35,
        TransverseIsotropicLinearElastic   = 36,
        MisesTI                            = 37,
        CDPM2                              = 38,
        DruckerPrager                      = 39,
        SandHypo                           = 40,
        SandHypoMicropolar                 = 41,
        GMBiotElastic                      = 42,
        GMDruckerPrager                    = 43,
        DruckerPragerMD                    = 44,
        GMNeoHooke                         = 45,
        GradientEnhancedDruckerPrager      = 46,
        GMCDPFiniteStrain                  = 47,
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
        T2D2 = 202,
        // Plane Stress
        CPS4 = 402,
        // CPS8  = 802,
        CPS8R = 805,

        // Plane Strain
        CPE4 = 407,
        // CPE8  = 807,
        CPE8R = 808,

        // Plane Strain - EAS
        CPE4EAS2 = 40702,
        CPE4EAS4 = 40704,
        CPE4EAS5 = 40705,

        // Solid
        C3D8   = 803,
        C3D8R  = 806,
        C3D20  = 2003,
        C3D20R = 2006,

        T2D2TL = 222,
        // Plane Stress
        CPS4TL  = 422,
        CPS8TL  = 822,
        CPS8RTL = 825,

        // Plane Strain
        CPE4TL  = 427,
        CPE8TL  = 827,
        CPE8RTL = 828,

        // Solid
        C3D8TL   = 823,
        C3D8RTL  = 826,
        C3D20TL  = 2023,
        C3D20RTL = 2026,

        CPE4UL  = 437,
        CPE8RUL = 838,
        C3D8UL  = 833,

        // Solid EAS
        C3D8EAS3  = 80303,
        C3D8EAS6b = 80361,
        C3D8EAS9  = 80309,

        // Nonlocal
        // Plane Stress

        CPS4NonLocal  = 412,
        CPS8NonLocal  = 812,
        CPS8RNonLocal = 815,

        // Nonlocal mixed
        CPS8RNonLocalMixed  = 875,
        CPE8RNonLocalMixed  = 877,
        C3D20RNonLocalMixed = 2076,

        // Plane Stress - EAS
        CPS4NonLocalEAS2 = 41202,
        CPS4NonLocalEAS4 = 41204,
        CPS4NonLocalEAS5 = 41205,

        // Plane Strain
        CPE4NonLocal  = 417,
        CPE4RNonLocal = 418,
        CPE8NonLocal  = 817,
        CPE8RNonLocal = 818,

        // Plane Strain - EAS
        CPE4NonLocalEAS2 = 41702,
        CPE4NonLocalEAS4 = 41704,
        CPE4NonLocalEAS5 = 41705,

        C3D8NonLocalEAS3  = 81303,
        C3D8NonLocalEAS6b = 81361,
        C3D8NonLocalEAS9  = 81309,

        // Solid
        C3D4NonLocal   = 413,
        C3D10NonLocal  = 1013,
        C3D8NonLocal   = 813,
        C3D8RNonLocal  = 816,
        C3D20NonLocal  = 2013,
        C3D20RNonLocal = 2016,

        // Nonlocal, Updated Lagrange
        C3D8NonLocalUL  = 843,
        CPE4NonLocalUL  = 447,
        CPE4RNonLocalUL = 448,
        CPE8RNonLocalUL = 848,
        // FBar versions
        C3D8NonLocalULFBar = 84301,
        CPE4NonLocalULFBar = 44701,

        // Cosserat
        CCPE4   = 458,
        CCPE8R  = 858,
        CC3D20R = 2056,
        CC3D8   = 853,

        // Nonlocal Cosserat
        NCCPS4   = 465,
        NCCPE8R  = 868,
        NCCPS8R  = 865,
        NCC3D20R = 2066,
        NCC3D8   = 866,

        GMCPE8R  = 898,
        GMC3D20R = 2096,
        GMC3D8,
        GMCPE8RUL,
        GMC3D8UL,
        GMC3D20RUL,
    };

    // MaterialFactory
    //
    // - Allows materials to register themselve with their name and ID
    // - Allows the user to create instances of materials

    class MarmotMaterialFactory {
      public:
        using materialFactoryFunction = MarmotMaterial* (*)( const double* materialProperties, int nMaterialProperties, int materialNumber );
        MarmotMaterialFactory()          = delete;

        static MaterialCode getMaterialCodeFromName( const std::string& materialName );

        static MarmotMaterial* createMaterial( MaterialCode  material, 
                const double* materialProperties,
                int nMaterialProperties, 
                int materialNumber);

        static bool registerMaterial( MaterialCode            materialCode,
                                      const std::string&      materialName,
                                      materialFactoryFunction factoryFunction );

      private:
        static std::map< std::string, MaterialCode >             materialNameToCodeAssociation;
        static std::map< MaterialCode, materialFactoryFunction > materialFactoryFunctionByCode;
    };

    // ElementFactory
    //
    // - Allows elements to register themselve with their name and ID
    // - Allows the user to create instances of elements

    class MarmotElementFactory {
      public:
        using elementFactoryFunction = MarmotElement* (*)( int elementNumber );
        MarmotElementFactory()          = delete;

        static ElementCode getElementCodeFromName( const std::string& elementName );

        static MarmotElement* createElement( ElementCode elementCode, int elementNumber );

        static bool registerElement( const std::string&     elementName,
                                     ElementCode            elementCode,
                                     elementFactoryFunction factoryFunction );

      private:
        static std::map< std::string, ElementCode >            elementNameToCodeAssociation;
        static std::map< ElementCode, elementFactoryFunction > elementFactoryFunctionByCode;
    };

} // namespace userLibrary
