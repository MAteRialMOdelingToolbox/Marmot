#pragma once

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

        enum DistributedLoadTypes{
            Pressure, 
        };

        BftUel(){};

        virtual ~BftUel(){};

        virtual void computeYourself( const double* QTotal,
                                            const double* dQ,
                                            double* Pe,
                                            double* Ke,
                                            const double* time,
                                            double dT,
                                            double& pNewdT) = 0;

        virtual void setInitialConditions(StateTypes state, const double* values){};

        virtual void computeDistributedLoad(
                                    DistributedLoadTypes loadType,
                                    double* P, 
                                    int elementFace, 
                                    const double* load,
                                    const double* time,
                                    double dT) {};
};

