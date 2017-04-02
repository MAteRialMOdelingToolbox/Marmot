#pragma once
#include "bftTypedefs.h"

class BftUel{

    public:

    enum StateTypes{
        Sigma11,
        Sigma22,
        Sigma33,
        HydrostaticStress,
        GeostaticStress,
        UmatStateVars
    };


        //BftUel(                         const double* coordinates,
                                            //double* stateVarsTotal,
                                            //int nStateVarsTotal,
                                            //const double* propertiesElement,
                                            //int nPropertiesElement,
                                            //int noEl,
                                            //const bft::pUmatType umat,
                                            //int nStateVarsUmat, 
                                            //const double* propertiesUmat,
                                            //int nPropertiesUmat);
        virtual ~BftUel();
        //BftUel();
        virtual void computeYourself( const double* QTotal,
                                            const double* dQ,
                                            double* Pe,
                                            double* Ke,
                                            const double* time,
                                            double dT,
                                            double& pNewdT) = 0;

        virtual void setInitialConditions(StateTypes state, const double* values, int nValues)=0;
};

