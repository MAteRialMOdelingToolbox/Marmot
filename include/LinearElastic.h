#pragma once
#include "MarmotMaterialHypoElastic.h"
#include <iostream>
#include <string>

class LinearElastic : public MarmotMaterialHypoElastic {
  protected:

    const double* materialProperties;

  public:
    using MarmotMaterialHypoElastic::MarmotMaterialHypoElastic;

    void computeStress( double*       stress,
                        double*       dStressDDStrain,
                        
                        const double* dStrain,
                        const double* timeOld,
                        const double  dT,
                        double&       pNewDT );

    void assignMaterialProperties ( const double* materialProperties, int nMaterialProperties );

    PermanentResultLocation getPermanentResultPointer( const std::string& result) 
    {
        return {nullptr, 0};
    };

    int getNumberOfRequiredStateVars() { return 0; }
};
