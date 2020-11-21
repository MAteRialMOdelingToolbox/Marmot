#pragma once
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotUtils.h"
#include <string>
#include <vector>

class MarmotElement {

  public:
    enum StateTypes { Sigma11, Sigma22, Sigma33, HydrostaticStress, GeostaticStress, MarmotMaterialStateVars, MarmotMaterialInitialization};

    enum DistributedLoadTypes {
        Pressure,
        SurfaceTorsion,
        SurfaceTraction,
    };

    virtual ~MarmotElement();

    virtual int getNumberOfRequiredStateVars() = 0;

    virtual std::vector<std::vector<std::string>> getNodeFields() = 0;

    virtual std::vector<int> getDofIndicesPermutationPattern() = 0;

    virtual int getNNodes() = 0;

    virtual int getNDofPerElement() = 0;

    virtual std::string getElementShape() = 0;

    virtual void assignStateVars( double* stateVars, int nStateVars ) = 0;

    virtual void assignProperty( const ElementProperties& property );

    virtual void assignProperty( const MarmotMaterialSection& property );

    virtual void initializeYourself( const double* coordinates ) = 0;

    virtual void setInitialConditions( StateTypes state, const double* values ) = 0;

    virtual void computeYourself( const double* QTotal,
                                  const double* dQ,
                                  double*       Pint,
                                  double*       K,
                                  const double* time,
                                  double        dT,
                                  double&       pNewdT ) = 0;

    virtual void computeDistributedLoad( DistributedLoadTypes loadType,
                                         double*              Pext,
                                         double*              K,
                                         int                  elementFace,
                                         const double*        load,
                                         const double*        QTotal,
                                         const double*        time,
                                         double               dT ) = 0;

    virtual void computeBodyForce( double*       Pext,
                                   double*       K,
                                   const double* load,
                                   const double* QTotal,
                                   const double* time,
                                   double        dT ) = 0;

    virtual PermanentResultLocation getPermanentResultPointer( const std::string& resultName, int gaussPt) = 0;
};
