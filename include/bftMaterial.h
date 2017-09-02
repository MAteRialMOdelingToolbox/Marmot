#pragma once
#include <string>

class BftMaterial{
    public:

        double* stateVars;
        const int nStateVars;

        const double* materialProperties;
        const int nMaterialProperties;

        const int elementID, gaussPt;

        BftMaterial(
                    double* stateVars, 
                    int nStateVars,
                    const double* materialProperties, 
                    int nMaterialProperties,
                    int elementID,
                    int gaussPt
                    ):
                    stateVars(stateVars),
                    nStateVars(nStateVars),
                    materialProperties(materialProperties),
                    nMaterialProperties(nMaterialProperties),
                    elementID(elementID),
                    gaussPt(gaussPt)
                    {};

        virtual ~BftMaterial(){};

        virtual double* getPermanentResultPointer(const std::string& resultName, int& resultLength) = 0;
};
