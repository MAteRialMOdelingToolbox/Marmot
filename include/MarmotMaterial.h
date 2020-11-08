#pragma once
#include "MarmotUtils.h"
#include <string>

class MarmotMaterial {

  protected:
    double* stateVars;
    int     nStateVars;

  public:
    const int materialNumber;

    MarmotMaterial( int materialNumber );

    virtual ~MarmotMaterial();

    virtual void assignMaterialProperties( const double* materialProperties, int nMaterialProperties ) = 0;

    virtual int getNumberOfRequiredStateVars() = 0;

    virtual void assignStateVars( double* stateVars, int nStateVars );

    virtual PermanentResultLocation getPermanentResultPointer( const std::string& resultName ) = 0;

    const double* getAssignedStateVars();

    int getNumberOfAssignedStateVars();
};
