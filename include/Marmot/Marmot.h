/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck,
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
 *
 * Matthias Neuner matthias.neuner@uibk.ac.at
 * Magdalena Schreter magdalena.schreter@uibk.ac.at
 *
 * This file is part of the MAteRialMOdellingToolbox (marmot).
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of marmot.
 * ---------------------------------------------------------------------
 */
#pragma once
#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotMaterial.h"
#include <functional>
#include <iostream>
#include <string>
#include <unordered_map>

namespace MarmotLibrary {

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
    OrthotropicLinearElastic           = 48,
    JointedHoekBrown                   = 49,
    ORDP                               = 50,
    B4                                 = 51,
    LinearElasticShrinkage             = 52,
    GCSCDP                             = 53,
    MenegottoPinto                     = 54,
    GCDP                               = 55,
    GosfordSandstone                   = 56,
    GMDamagedShearNeoHooke             = 57,
    SolidificationCDP                  = 58,
    SolidificationModLeon              = 59,
    ORDPNonLocal                       = 60,
    CosseratHoekBrown                  = 61,
    BulkMetallicGlass                  = 62,
    RedWildmoorSandstone               = 63,
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

    GCPS4  = 412,
    GCPS8  = 812,
    GCPS8R = 815,

    // Nonlocal mixed
    GCPS8RMixed  = 875,
    GCPE8RMixed  = 877,
    GC3D20RMixed = 2076,

    // Plane Stress - EAS
    GCPS4EAS2 = 41202,
    GCPS4EAS4 = 41204,
    GCPS4EAS5 = 41205,

    // Plane Strain
    GCPE4  = 417,
    GCPE4R = 418,
    GCPE8  = 817,
    GCPE8R = 818,

    // Plane Strain - EAS
    GCPE4EAS2 = 41702,
    GCPE4EAS4 = 41704,
    GCPE4EAS5 = 41705,

    GC3D8EAS3  = 81303,
    GC3D8EAS6b = 81361,
    GC3D8EAS9  = 81309,

    // Solid
    GC3D4   = 413,
    GC3D10  = 1013,
    GC3D8   = 813,
    GC3D8R  = 816,
    GC3D20  = 2013,
    GC3D20R = 2016,

    // Nonlocal, Updated Lagrange
    GC3D8UL  = 843,
    GCPE4UL  = 447,
    GCPE4RUL = 448,
    GCPE8RUL = 848,
    // FBar versions
    GC3D8ULFBar = 84301,
    GCPE4ULFBar = 44701,

    // Cosserat
    CCPE4   = 458,
    CCPE8R  = 858,
    CC3D20R = 2056,
    CC3D8   = 853,

    // Nonlocal Cosserat
    GCCPS4   = 465,
    GCCPE8R  = 868,
    GCCPS8R  = 865,
    GCC3D20R = 2066,
    GCC3D8   = 866,

    GMCPE8R  = 898,
    GMC3D20R = 2096,
    GMC3D8,
    GMCPE8RUL,
    GMC3D8UL,
    GMC3D20RUL,

    CHMCPE4UL,
    CHMCPE8RUL,
    CHMC3D20RUL,
  };

  // MaterialFactory
  //
  // - Allows materials to register themselve with their name and ID
  // - Allows the user to create instances of materials

  class MarmotMaterialFactory {
  public:
    using materialFactoryFunction = MarmotMaterial* (*)( const double* materialProperties,
                                                         int           nMaterialProperties,
                                                         int           materialNumber );
    MarmotMaterialFactory()       = delete;

    static MaterialCode getMaterialCodeFromName( const std::string& materialName );

    static MarmotMaterial* createMaterial( MaterialCode  material,
                                           const double* materialProperties,
                                           int           nMaterialProperties,
                                           int           materialNumber );

    static bool registerMaterial( MaterialCode            materialCode,
                                  const std::string&      materialName,
                                  materialFactoryFunction factoryFunction );

  private:
    static std::unordered_map< std::string, MaterialCode >             materialNameToCodeAssociation;
    static std::unordered_map< MaterialCode, materialFactoryFunction > materialFactoryFunctionByCode;
  };

  // ElementFactory
  //
  // - Allows elements to register themselve with their name and ID
  // - Allows the user to create instances of elements

  class MarmotElementFactory {
  public:
    using elementFactoryFunction = MarmotElement* (*)( int elementNumber );
    MarmotElementFactory()       = delete;

    static ElementCode getElementCodeFromName( const std::string& elementName );

    static MarmotElement* createElement( ElementCode elementCode, int elementNumber );

    static bool registerElement( const std::string&     elementName,
                                 ElementCode            elementCode,
                                 elementFactoryFunction factoryFunction );

  private:
    static std::unordered_map< std::string, ElementCode >            elementNameToCodeAssociation;
    static std::unordered_map< ElementCode, elementFactoryFunction > elementFactoryFunctionByCode;
  };

} // namespace MarmotLibrary
