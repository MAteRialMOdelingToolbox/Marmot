#pragma once 
#include <string>
#include "bftUel.h"
#include "bftMaterial.h"

namespace userLibrary{

    enum MaterialCode{
        LinearElastic=0,
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
        ShotLeonV2NonLocal=12,
        LinearElasticSolidificationCreep = 13,
        ModLeonAdaptive=14,
        ModLeonSemiExplicit=15,
        ModLeonSemiExplicitAdaptive=16,
    };

    MaterialCode getMaterialCodeFromName(const std::string& materialName);

    BftUel* UelFactory(
            int elementCode,
            const double* elementCoordinates, double* stateVars, int nStateVars, 
            const double* propertiesElement, int nPropertiesElement, int elementNumber, 
            MaterialCode material,
            int nStateVarsUmat, const double* propertiesUmat, int nPropertiesUmat);

    BftMaterial* bftMaterialFactory(
                                    MaterialCode material,
                                    double* stateVars,
                                    int nStateVars,
                                    const double* materialProperties, 
                                    int nMaterialProperties,
                                    int element, 
                                    int gaussPt
                                    );

    static const int sizeGeostaticDefinition = 6;
}

namespace Abaqus{
    
    enum UelFlags1{
        GeostaticStress=61,   // Geostatic stress field according to Abaqus Analysis User's Guide Tab. 5.1.2-1 Keys to procedure types.
    };
}


