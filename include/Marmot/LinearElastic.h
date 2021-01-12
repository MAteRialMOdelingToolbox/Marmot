#pragma once
#include "Marmot/MarmotMaterialHypoElastic.h"
#include <iostream>
#include <string>

class LinearElastic : public MarmotMaterialHypoElastic {
  public:
    using MarmotMaterialHypoElastic::MarmotMaterialHypoElastic;

    LinearElastic( const double* materialProperties, int nMaterialProperties, int materialNumber );

  protected:
    const double& E;
    const double& nu;

    void computeStress( double* stress,
                        double* dStressDDStrain,

                        const double* dStrain,
                        const double* timeOld,
                        const double  dT,
                        double&       pNewDT );

    PermanentResultLocation getPermanentResultPointer( const std::string& result ) { return { nullptr, 0 }; };

    int getNumberOfRequiredStateVars() { return 0; }
};
