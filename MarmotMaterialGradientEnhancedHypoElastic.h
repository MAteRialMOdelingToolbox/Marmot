#pragma once
#include "Marmot/MarmotMaterialGradientEnhancedMechanical.h"

class MarmotMaterialGradientEnhancedHypoElastic : public MarmotMaterialGradientEnhancedMechanical {

  public:
    using MarmotMaterialGradientEnhancedMechanical::MarmotMaterialGradientEnhancedMechanical;
    
    // Abstract methods
    virtual void computeStress( double*       stress,
                                double&       K_local,
                                double&       nonLocalRadius,
                                double*       dStressDDDeformationGradient,
                                double*       dK_localDDeformationGradient,
                                double*       dStressDK,
                                const double* FOld,
                                const double* FNew,
                                const double  KOld,
                                const double  dK,
                                const double* timeOld,
                                const double  dT,
                                double&       pNewDT ) override;


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

    using MarmotMaterialGradientEnhancedMechanical::computePlaneStress;
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
