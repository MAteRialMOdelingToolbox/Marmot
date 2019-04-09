#pragma once
#include "bftMaterial.h"

class BftMaterialHypoElasticNonLocal : public BftMaterial {

  public:
    using BftMaterial::BftMaterial;

    // Abstract methods
    virtual void computeStress( double*       stress,
                                double&       K_local,
                                double&       nonLocalRadius,
                                double*       dStressDDStrain,
                                double*       dK_localDDStrain,
                                double*       dStressDK,
                                const double* dStrain,
                                double        KOld,
                                double        dK,
                                const double* timeOld,
                                const double  dT,
                                double&       pNewDT ) = 0;

    virtual void computePlaneStress( double*       stress,
                                     double&       K_local,
                                     double&       nonLocalRadius,
                                     double*       dStressDDStrain,
                                     double*       dK_localDDStrain,
                                     double*       dStressDK,
                                     double*       dStrain,
                                     double        KOld,
                                     double        dK,
                                     const double* timeOld,
                                     const double  dT,
                                     double&       pNewDT );
};
