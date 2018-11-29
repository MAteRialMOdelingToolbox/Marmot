#pragma once
#include "bftUelProperty.h"
#include <string>
#include <vector>

class BftUel {

  public:
    enum StateTypes { Sigma11, Sigma22, Sigma33, HydrostaticStress, GeostaticStress, UmatStateVars };

    enum DistributedLoadTypes {
        Pressure,
    };

    // todo: check if needed
    BftUel(){};

    virtual ~BftUel(){};

    virtual int getNumberOfRequiredStateVars() = 0;

    virtual std::vector<std::vector<std::string>> getNodeFields() = 0;

    virtual std::vector<int> getDofIndicesPermutationPattern() = 0;

    virtual int getNNodes() = 0;

    virtual int getNDofPerElement() = 0;

    virtual std::string getElementShape() = 0;

    virtual void assignStateVars( double* stateVars, int nStateVars ) = 0;

    virtual void assignProperty( const BftUelProperty& property ) = 0;

    virtual void initializeYourself( const double* coordinates ) = 0;

    virtual void computeYourself( const double* QTotal,
                                  const double* dQ,
                                  double*       Pint,
                                  double*       Ke,
                                  const double* time,
                                  double        dT,
                                  double&       pNewdT ) = 0;

    virtual void setInitialConditions( StateTypes state, const double* values ) = 0;

    virtual void computeDistributedLoad( DistributedLoadTypes loadType,
                                         double*              Pext,
                                         int                  elementFace,
                                         const double*        QTotal,
                                         const double*        load,
                                         const double*        time,
                                         double               dT ) = 0;

    virtual void computeBodyForce( double*       Pext,
                                   const double* load,
                                   const double* QTotal,
                                   const double* time,
                                   double        dT ) = 0;

    virtual double* getPermanentResultPointer( const std::string& resultName, int gaussPt, int& resultLength ) = 0;
};
