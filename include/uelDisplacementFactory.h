#pragma once
#include "bftUel.h"
#include "bftTypedefs.h"

namespace UelDisplacementFactory{
    
    BftUel* generateUelCPS4(const double* coordinates, double* stateVars, int nStateVars, const double* propertiesElement, int nPropertiesElement,
                            int noEl, const bft::pUmatType umat, int nStateVarsUmat, const double* propertiesUmat, int nPropertiesUmat);

    BftUel* generateUelCPS8(const double* coordinates, double* stateVars, int nStateVars, const double* propertiesElement, int nPropertiesElement,
                            int noEl, const bft::pUmatType umat, int nStateVarsUmat, const double* propertiesUmat, int nPropertiesUmat);

    BftUel* generateUelCPE4(const double* coordinates, double* stateVars, int nStateVars, const double* propertiesElement, int nPropertiesElement,
                            int noEl, const bft::pUmatType umat, int nStateVarsUmat, const double* propertiesUmat, int nPropertiesUmat);

    BftUel* generateUelCPE8(const double* coordinates, double* stateVars, int nStateVars, const double* propertiesElement, int nPropertiesElement,
                            int noEl, const bft::pUmatType umat, int nStateVarsUmat, const double* propertiesUmat, int nPropertiesUmat);

    BftUel* generateUelCPS4R(const double* coordinates, double* stateVars, int nStateVars, const double* propertiesElement, int nPropertiesElement,
                            int noEl, const bft::pUmatType umat, int nStateVarsUmat, const double* propertiesUmat, int nPropertiesUmat);

    BftUel* generateUelCPS8R(const double* coordinates, double* stateVars, int nStateVars, const double* propertiesElement, int nPropertiesElement,
                            int noEl, const bft::pUmatType umat, int nStateVarsUmat, const double* propertiesUmat, int nPropertiesUmat);

    BftUel* generateUelCPE4R(const double* coordinates, double* stateVars, int nStateVars, const double* propertiesElement, int nPropertiesElement,
                            int noEl, const bft::pUmatType umat, int nStateVarsUmat, const double* propertiesUmat, int nPropertiesUmat);

    BftUel* generateUelCPE8R(const double* coordinates, double* stateVars, int nStateVars, const double* propertiesElement, int nPropertiesElement,
                            int noEl, const bft::pUmatType umat, int nStateVarsUmat, const double* propertiesUmat, int nPropertiesUmat);

}
