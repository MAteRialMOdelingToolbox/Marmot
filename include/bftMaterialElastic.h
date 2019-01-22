#pragma once
#include "bftMaterial.h"

class BftMaterialMechanical : public BftMaterial{

  public:
    using BftMaterial::BftMaterial;

    // Abstract methods
    virtual void computeStress( double*       stress,
                                double*       dStressDDStrain, // d Chauchy d (Displacement Gradient) vs.  d Cauchy d Fnew ???
                                const double* FOld,
                                const double* FNew,
                                const double* timeOld,
                                const double  dT,
                                double&       pNewDT ) = 0;

    //virtual void computePlaneStress( double*       stress,
                                     //double*       dStressDDStrain,
                                     //const double* strainOld,
                                     //double*       dStrain,
                                     //const double* timeOld,
                                     //const double  dT,
                                     //double&       pNewDT );

    //virtual void computeUniaxialStress( double*       stress,
                                        //double*       dStressDDStrain,
                                        //const double* strainOld,
                                        //double*       dStrain,
                                        //const double* timeOld,
                                        //const double  dT,
                                        //double&       pNewDT );
};
