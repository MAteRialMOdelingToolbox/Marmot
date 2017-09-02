#pragma once 
#include <string>
#include "bftUel.h"
#include "bftMaterial.h"

//namespace bft{
    //[>copy of bftTypedefs.h<]
        //typedef void (*pUmatType)( double[],double[],double[],double&,double&,double&,double&,double[],double[],double&,const double[],const double[],const double[2],const double&,const double&,const double&,const double[],const double[],const char[80],const int&,const int&,const int&,const int&,const double[],const int&,const double[3],const double[9],double&,const double&,const double[9],const double[9],const int&,const int&,const int&,const int&,const int[4],const int&,const int);

        //typedef void (*pSimpleUelWithUmatType)( double [], double [], double [], const int &, const double [], const int &, const double [], const double [], const double [], const double [2], const double &, const int& , double &, const int [], const int &, bft::pUmatType, int );

//}

namespace userLibrary{

    //bft::pUmatType getUmatById(int id);
    //bft::pUmatType getUmatByName(const std::string& nameUpperCase);

    BftUel* UelFactory(int id, const double* elementCoordinates, double* stateVars, int nStateVars, 
            const double* propertiesElement, int nPropertiesElement, int elementNumber, 
            const std::string& bftMaterialName,
            int nStateVarsUmat, const double* propertiesUmat, int nPropertiesUmat);

    BftMaterial* bftMaterialFactory(const std::string& nameUpperCase ,
                                    double* stateVars,
                                    int nStateVars,
                                    const double* materialProperties, 
                                    int nMaterialProperties,
                                    int element, 
                                    int gaussPt
                                    );

    //void umatPlaneStressWrapped(const bft::pUmatType, double[],double[],double[],double&,double&,double&,double&,double[],double[],double&,const double[],const double[],const double[2],const double&,const double&,const double&,const double[],const double[],const char[80],const int&, const int&,const int&,const int&,const double[],const int&,const double[3],const double[9],double&,const double&,const double[9],const double[9],const int&,const int&,const int&,const int&,const int[4],const int&,const int);

    static const int sizeGeostaticDefinition = 6;
}

namespace Abaqus{
    
    enum UelFlags1{
        GeostaticStress=61,   // Geostatic stress field according to Abaqus Analysis User's Guide Tab. 5.1.2-1 Keys to procedure types.
    };
}


