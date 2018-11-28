#pragma once
#include "bftMaterialHypoElastic.h"
#include <iostream>
#include <string>

class LinearElastic : public BftMaterialHypoElastic {
  public:
    using BftMaterialHypoElastic::BftMaterialHypoElastic;

    void computeStress( double*       stress,
                        double*       dStressDDStrain,
                        const double* strainOld,
                        const double* dStrain,
                        const double* timeOld,
                        const double  dT,
                        double&       pNewDT );

    double* getPermanentResultPointer( const std::string& result, int& resultLength ) { return nullptr; }

    int getNumberOfRequiredStateVars() { return 2; }

    void assignStateVars( double* stateVars, int nStateVars )
    {
        this->stateVars  = stateVars;
        this->nStateVars = nStateVars;
    };
};
