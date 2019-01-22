#pragma once
#include "bftMaterialMechanical.h"

class BftMaterialHypoElastic : public BftMaterialMechanical {

  public:
    using BftMaterialMechanical::BftMaterialMechanical;

    double characteristicElementLength;
    void   setCharacteristicElementLength( double length ) { characteristicElementLength = length; }

    // Abstract methods
    virtual void computeStress( double*       stress,
                                double*       dStressDDStrain,
                                const double* strainOld,
                                const double* dStrain,
                                const double* timeOld,
                                const double  dT,
                                double&       pNewDT ) = 0;

    virtual void computePlaneStress( double*       stress,
                                     double*       dStressDDStrain,
                                     const double* strainOld,
                                     double*       dStrain,
                                     const double* timeOld,
                                     const double  dT,
                                     double&       pNewDT );

    virtual void computeUniaxialStress( double*       stress,
                                        double*       dStressDDStrain,
                                        const double* strainOld,
                                        double*       dStrain,
                                        const double* timeOld,
                                        const double  dT,
                                        double&       pNewDT );
};
