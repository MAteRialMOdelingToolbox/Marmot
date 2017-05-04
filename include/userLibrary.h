#pragma once 
#include <string>

namespace bft{
    /*copy of bftTypedefs.h*/
        typedef void (*pUmatType)( double[],double[],double[],double&,double&,double&,double&,double[],double[],double&,const double[],const double[],const double[2],const double&,const double&,const double&,const double[],const double[],const char[80],const int&,const int&,const int&,const int&,const double[],const int&,const double[3],const double[9],double&,const double&,const double[9],const double[9],const int&,const int&,const int&,const int&,const int[4],const int&,const int);

        typedef void (*pSimpleUelWithUmatType)( double [], double [], double [], const int &, const double [], const int &, const double [], const double [], const double [], const double [2], const double &, const int& , double &, const int [], const int &, bft::pUmatType, int );

}


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
        virtual void computeYourself( const double* QTotal,
                                            const double* dQ,
                                            double* Pe,
                                            double* Ke,
                                            const double* time,
                                            double dT,
                                            double& pNewdT) = 0;

        virtual void setInitialConditions(StateTypes state, const double* values, int nValues)=0;
};

namespace userLibrary{

    bft::pUmatType getUmatById(int id);
    bft::pUmatType getUmatByName(const std::string& nameUpperCase);

    BftUel* UelFactory(int id, const double* elementCoordinates, double* stateVars, int nStateVars, 
            const double* propertiesElement, int nPropertiesElement, int elementNumber, 
            bft::pUmatType umat, int nStateVarsUmat, const double* propertiesUmat, int nPropertiesUmat);

    //bft::pSimpleUelWithUmatType getSimpleUelWithUmatById(int id);
    void umatPlaneStressWrapped(const bft::pUmatType, double[],double[],double[],double&,double&,double&,double&,double[],double[],double&,const double[],const double[],const double[2],const double&,const double&,const double&,const double[],const double[],const char[80],const int&, const int&,const int&,const int&,const double[],const int&,const double[3],const double[9],double&,const double&,const double[9],const double[9],const int&,const int&,const int&,const int&,const int[4],const int&,const int);


    static const int sizeGeostaticDefinition = 5;
}

namespace Abaqus{
    
    enum UelFlags1{
        Geostatic=61,
    };
}


