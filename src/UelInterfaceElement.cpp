#include "UelInterfaceElement.h"

#include <iostream>
#include <string.h>

UelInterfaceElement::UelInterfaceElement(int elNumber, 
                                const double* coordinates, 
                                double* stateVars,
                                int nStateVars,
                                //double* stateVarsTemp,
                                const double* properties,
                                int nProperties,
                                const int* intProperties,
                                int nIntProperties,
                                pUmatType umat,
                                int nStateVarsUmat,
                                pSimpleUelWithUmatType uel):
    elNumber(elNumber),
    coordinates(coordinates),
    stateVars(stateVars),
    nStateVars(nStateVars),
    //stateVarsTemp(stateVarsTemp),
    properties(properties),
    nProperties(nProperties),
    intProperties(intProperties),
    nIntProperties(nIntProperties),
    umat(umat),
    nStateVarsUmat(nStateVarsUmat),
    uel(uel)
{
    stateVarsTemp = new double[nStateVars];
}

UelInterfaceElement::~UelInterfaceElement()
{
    delete stateVarsTemp;
}

void UelInterfaceElement::computeYourself(double* Pe, double* Ke, const double* UNew, const double* dU,  const double time[], double dTime, double &pNewDT )
{
    memcpy(stateVarsTemp, stateVars, nStateVars * sizeof(double));
    uel(Pe, Ke, stateVarsTemp, nStateVars, properties, nProperties, coordinates, UNew, dU, time, 
                     dTime, elNumber, pNewDT, intProperties, nIntProperties, umat, nStateVarsUmat);
    
}

void UelInterfaceElement::acceptLastState()
{
    memcpy(stateVars, stateVarsTemp, nStateVars* sizeof(double));
}

