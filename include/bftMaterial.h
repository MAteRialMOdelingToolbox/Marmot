#pragma once
#include <string>

class BftMaterial{
    public:

        double* stateVars;
        int nStateVars;

        const double* materialProperties;
        const int nMaterialProperties;

        const int elementID, gaussPt;

        BftMaterial(
                    //double* stateVars, 
                    //int nStateVars,
                    const double* materialProperties, 
                    int nMaterialProperties,
                    int elementID,
                    int gaussPt
                    ):
                    //stateVars(stateVars),
                    //nStateVars(nStateVars),
                    stateVars(nullptr),
                    nStateVars(0),
                    materialProperties(materialProperties),
                    nMaterialProperties(nMaterialProperties),
                    elementID(elementID),
                    gaussPt(gaussPt)
                    {};

        virtual ~BftMaterial(){};

        virtual int getNumberOfRequiredStateVars() = 0;

        virtual void assignStateVars(double *stateVars, int nStateVars) { 
            this->stateVars = stateVars;
            this->nStateVars = nStateVars;};

        virtual double* getPermanentResultPointer(const std::string& resultName, int& resultLength) = 0;

};
