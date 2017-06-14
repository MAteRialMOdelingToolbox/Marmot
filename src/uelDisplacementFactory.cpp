#include "uelDisplacementFactory.h"
#include "uelDisplacementPlane.h"
#include "uelDisplacementSolid.h"

namespace UelDisplacementFactory{
    
    BftUel* generateUelCPS4(const double* coordinates, double* stateVars, int nStateVars, const double* propertiesElement, int nPropertiesElement,
                            int noEl, const bft::pUmatType umat, int nStateVarsUmat, const double* propertiesUmat, int nPropertiesUmat){
            return new UelDisplacementPlane< 4 >(coordinates, stateVars, nStateVars, propertiesElement, 
                                 nPropertiesElement, noEl, umat, nStateVarsUmat, propertiesUmat, 
                                 nPropertiesUmat,  bft::NumIntegration::IntegrationTypes::FullIntegration, 
                                 UelDisplacementPlane<4>::SectionType::PlaneStress); }

    BftUel* generateUelCPS8(const double* coordinates, double* stateVars, int nStateVars, const double* propertiesElement, int nPropertiesElement,
                            int noEl, const bft::pUmatType umat, int nStateVarsUmat, const double* propertiesUmat, int nPropertiesUmat){
            return new UelDisplacementPlane< 8 >(coordinates, stateVars, nStateVars, propertiesElement, 
                                 nPropertiesElement, noEl, umat, nStateVarsUmat, propertiesUmat, 
                                 nPropertiesUmat,  bft::NumIntegration::IntegrationTypes::FullIntegration, 
                                 UelDisplacementPlane<8>::SectionType::PlaneStress); }

    BftUel* generateUelCPE4(const double* coordinates, double* stateVars, int nStateVars, const double* propertiesElement, int nPropertiesElement,
                            int noEl, const bft::pUmatType umat, int nStateVarsUmat, const double* propertiesUmat, int nPropertiesUmat){
            return new UelDisplacementPlane< 4 >(coordinates, stateVars, nStateVars, propertiesElement, 
                                 nPropertiesElement, noEl, umat, nStateVarsUmat, propertiesUmat, 
                                 nPropertiesUmat,  bft::NumIntegration::IntegrationTypes::FullIntegration, 
                                 UelDisplacementPlane<4>::SectionType::PlaneStrain); }

    BftUel* generateUelCPE8(const double* coordinates, double* stateVars, int nStateVars, const double* propertiesElement, int nPropertiesElement,
                            int noEl, const bft::pUmatType umat, int nStateVarsUmat, const double* propertiesUmat, int nPropertiesUmat){
            return new UelDisplacementPlane< 8 >(coordinates, stateVars, nStateVars, propertiesElement, 
                                 nPropertiesElement, noEl, umat, nStateVarsUmat, propertiesUmat, 
                                 nPropertiesUmat,  bft::NumIntegration::IntegrationTypes::FullIntegration, 
                                 UelDisplacementPlane<8>::SectionType::PlaneStrain); }

    BftUel* generateUelCPS4R(const double* coordinates, double* stateVars, int nStateVars, const double* propertiesElement, int nPropertiesElement,
                            int noEl, const bft::pUmatType umat, int nStateVarsUmat, const double* propertiesUmat, int nPropertiesUmat){
            return new UelDisplacementPlane< 4 >(coordinates, stateVars, nStateVars, propertiesElement, 
                                 nPropertiesElement, noEl, umat, nStateVarsUmat, propertiesUmat, 
                                 nPropertiesUmat,  bft::NumIntegration::IntegrationTypes::ReducedIntegration, 
                                 UelDisplacementPlane<4>::SectionType::PlaneStress); }

    BftUel* generateUelCPS8R(const double* coordinates, double* stateVars, int nStateVars, const double* propertiesElement, int nPropertiesElement,
                            int noEl, const bft::pUmatType umat, int nStateVarsUmat, const double* propertiesUmat, int nPropertiesUmat){
            return new UelDisplacementPlane< 8 >(coordinates, stateVars, nStateVars, propertiesElement, 
                                 nPropertiesElement, noEl, umat, nStateVarsUmat, propertiesUmat, 
                                 nPropertiesUmat,  bft::NumIntegration::IntegrationTypes::ReducedIntegration, 
                                 UelDisplacementPlane<8>::SectionType::PlaneStress); }

    BftUel* generateUelCPE4R(const double* coordinates, double* stateVars, int nStateVars, const double* propertiesElement, int nPropertiesElement,
                            int noEl, const bft::pUmatType umat, int nStateVarsUmat, const double* propertiesUmat, int nPropertiesUmat){
            return new UelDisplacementPlane< 4 >(coordinates, stateVars, nStateVars, propertiesElement, 
                                 nPropertiesElement, noEl, umat, nStateVarsUmat, propertiesUmat, 
                                 nPropertiesUmat,  bft::NumIntegration::IntegrationTypes::ReducedIntegration, 
                                 UelDisplacementPlane<4>::SectionType::PlaneStrain); }

    BftUel* generateUelCPE8R(const double* coordinates, double* stateVars, int nStateVars, const double* propertiesElement, int nPropertiesElement,
                            int noEl, const bft::pUmatType umat, int nStateVarsUmat, const double* propertiesUmat, int nPropertiesUmat){
            return new UelDisplacementPlane< 8 >(coordinates, stateVars, nStateVars, propertiesElement, 
                                 nPropertiesElement, noEl, umat, nStateVarsUmat, propertiesUmat, 
                                 nPropertiesUmat,  bft::NumIntegration::IntegrationTypes::ReducedIntegration, 
                                 UelDisplacementPlane<8>::SectionType::PlaneStrain); }

    BftUel* generateUelC3D8(const double* coordinates, double* stateVars, int nStateVars, const double* propertiesElement, int nPropertiesElement,
                            int noEl, const bft::pUmatType umat, int nStateVarsUmat, const double* propertiesUmat, int nPropertiesUmat){
            return new UelDisplacementSolid< 8 >(coordinates, stateVars, nStateVars, propertiesElement, 
                                 nPropertiesElement, noEl, umat, nStateVarsUmat, propertiesUmat, 
                                 nPropertiesUmat,  bft::NumIntegration::IntegrationTypes::FullIntegration); }

    BftUel* generateUelC3D20(const double* coordinates, double* stateVars, int nStateVars, const double* propertiesElement, int nPropertiesElement,
                            int noEl, const bft::pUmatType umat, int nStateVarsUmat, const double* propertiesUmat, int nPropertiesUmat){
            return new UelDisplacementSolid< 20 >(coordinates, stateVars, nStateVars, propertiesElement, 
                                 nPropertiesElement, noEl, umat, nStateVarsUmat, propertiesUmat, 
                                 nPropertiesUmat,  bft::NumIntegration::IntegrationTypes::FullIntegration); }

    BftUel* generateUelC3D20R(const double* coordinates, double* stateVars, int nStateVars, const double* propertiesElement, int nPropertiesElement,
                            int noEl, const bft::pUmatType umat, int nStateVarsUmat, const double* propertiesUmat, int nPropertiesUmat){
            return new UelDisplacementSolid< 20 >(coordinates, stateVars, nStateVars, propertiesElement, 
                                 nPropertiesElement, noEl, umat, nStateVarsUmat, propertiesUmat, 
                                 nPropertiesUmat,  bft::NumIntegration::IntegrationTypes::ReducedIntegration); }



}
