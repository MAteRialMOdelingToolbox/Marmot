#pragma once
#include "MarmotUtils.h"
#include <string>

class MarmotMaterial {

  protected:

    const double* materialProperties;
    const int nMaterialProperties;

    double* stateVars;
    int     nStateVars;

  public:
    const int materialNumber;

    MarmotMaterial( const double*materialProperties, int nMaterialProperties, int materialNumber );

    virtual ~MarmotMaterial();

    virtual int getNumberOfRequiredStateVars() = 0;

    virtual void assignStateVars( double* stateVars, int nStateVars );

    virtual PermanentResultLocation getPermanentResultPointer( const std::string& resultName ) = 0;

    const double* getAssignedStateVars();

    int getNumberOfAssignedStateVars();
};
