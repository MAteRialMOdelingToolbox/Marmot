#pragma once
#include "bftMaterial.h"
#include "bftUel.h"
#include <iostream>
#include <map>
#include <string>

namespace userLibrary {

    enum MaterialCode : int {
        ModLeon                          = 1,
        ShotLeon                         = 2,
        Meschke                          = 3,
        SchaedlichSchweiger              = 4,
        ModLeonNonLocal                  = 5,
        HoekBrown                        = 6,
        UntereggerRockMass               = 7,
        MohrCoulomb                      = 8,
        UntereggerRockMassNonLocal       = 9,
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
    };

    MaterialCode getMaterialCodeFromName( const std::string& materialName );

    enum ElementCode {

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
         *                  2: mechanical (=displacement) large strain TL,
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
        UelCPS4  = 402,
        UelCPS8  = 802,
        UelCPS8R = 805,

        // Plane Strain
        UelCPE4  = 407,
        UelCPE8  = 807,
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

        // Solid EAS
        UelC3D8EAS3  = 80303,
        UelC3D8EAS6b = 80361,
        UelC3D8EAS9  = 80309,

        // Nonlocal
        // Plane Stress

        UelCPS4NonLocal  = 412,
        UelCPS8NonLocal  = 812,
        UelCPS8RNonLocal = 815,

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
        UelC3D8NonLocal   = 813,
        UelC3D8RNonLocal  = 816,
        UelC3D20NonLocal  = 2013,
        UelC3D20RNonLocal = 2016,
    };

    ElementCode getElementCodeFromInfo( int                nNodes,
                                        int                activeFields,
                                        int                dim,
                                        bool               fullIntegration,
                                        const std::string& formulation,
                                        int                EASParam );

    ElementCode getElementCodeFromName( const std::string& elementName );

    BftUel* UelFactory( ElementCode elementCode, int elementNumber );

    BftMaterial* bftMaterialFactory( MaterialCode  material,
                                     const double* materialProperties,
                                     int           nMaterialProperties,
                                     int           element,
                                     int           gaussPt );

    void extendAbaqusToVoigt( double*       stress6,
                              double*       stress,
                              double*       strain6,
                              const double* strain,
                              double*       dStrain6,
                              const double* dStrain,
                              int           nDirect,
                              int           nShear );

    void backToAbaqus( double* stress,
                       double* stress6,
                       double* dStressDDStrain,
                       double* dStressDDStrain66,
                       int     nDirect,
                       int     nShear );

    static const int sizeGeostaticDefinition = 6;
} // namespace userLibrary

namespace Abaqus {

    enum UelFlags1 {
        GeostaticStress = 61, // Geostatic stress field according to Abaqus Analysis User's Guide Tab. 5.1.2-1 Keys to
                              // procedure types.
    };
}
