#pragma once 
#include <string>

namespace bft{
    /*copy of bftTypedefs.h*/
        typedef void (*pUmatType)( double[],double[],double[],double&,double&,double&,double&,double[],double[],double&,const double[],const double[],const double[2],const double&,const double&,const double&,const double[],const double[],const char[80],const int&,const int&,const int&,const int&,const double[],const int&,const double[3],const double[9],double&,const double&,const double[9],const double[9],const int&,const int&,const int&,const int&,const int[4],const int&,const int);

        typedef void (*pSimpleUelWithUmatType)( double [], double [], double [], const int &, const double [], const int &, const double [], const double [], const double [], const double [2], const double &, const int& , double &, const int [], const int &, bft::pUmatType, int );

}

namespace userLibrary{

    bft::pUmatType getUmatById(int id);
    bft::pUmatType getUmatByName(const std::string& nameUpperCase);

    bft::pSimpleUelWithUmatType getSimpleUelWithUmatById(int id);
}

