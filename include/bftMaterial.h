#pragma once
#include <string>

class BftMaterial {

  public:
    double* stateVars;
    int     nStateVars;

    const double* materialProperties;
    const int     nMaterialProperties;

    const int elementID, gaussPt;

    BftMaterial( const double* materialProperties, int nMaterialProperties, int elementID, int gaussPt );

    virtual ~BftMaterial();

    virtual int getNumberOfRequiredStateVars() = 0;

    virtual void assignStateVars( double* stateVars, int nStateVars ) = 0;

    virtual double* getPermanentResultPointer( const std::string& resultName, int& resultLength ) = 0;

};
