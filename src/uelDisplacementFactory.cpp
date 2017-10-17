#include "uelDisplacementFactory.h"
#include "uelDisplacement.h"
//#include "uelDisplacement.h"

namespace UelDisplacementFactory{
    
    BftUel* generateUelCPS4(const double* coordinates, double* stateVars, int nStateVars, const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode materialCode, int nStateVarsUmat, const double* materialProperties, int nMaterialProperties){

            return new UelDisplacement<2, 4 >(coordinates, stateVars, nStateVars, elementProperties, 
                                 nElementProperties, noEl, materialCode , nStateVarsUmat, materialProperties, 
                                 nMaterialProperties,  bft::NumIntegration::IntegrationTypes::FullIntegration, 
                                 UelDisplacement<2, 4>::SectionType::PlaneStress); }

    BftUel* generateUelCPS8(const double* coordinates, double* stateVars, int nStateVars, const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode bftMaterialHypoElasticName, int nStateVarsUmat, const double* materialProperties, int nMaterialProperties){
            return new UelDisplacement<2, 8 >(coordinates, stateVars, nStateVars, elementProperties, 
                                 nElementProperties, noEl, bftMaterialHypoElasticName, nStateVarsUmat, materialProperties, 
                                 nMaterialProperties,  bft::NumIntegration::IntegrationTypes::FullIntegration, 
                                 UelDisplacement<2,8>::SectionType::PlaneStress); }

    BftUel* generateUelCPE4(const double* coordinates, double* stateVars, int nStateVars, const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode bftMaterialHypoElasticName, int nStateVarsUmat, const double* materialProperties, int nMaterialProperties){
            return new UelDisplacement<2, 4 >(coordinates, stateVars, nStateVars, elementProperties, 
                                 nElementProperties, noEl, bftMaterialHypoElasticName, nStateVarsUmat, materialProperties, 
                                 nMaterialProperties,  bft::NumIntegration::IntegrationTypes::FullIntegration, 
                                 UelDisplacement<2,4>::SectionType::PlaneStrain); }

    BftUel* generateUelCPE8(const double* coordinates, double* stateVars, int nStateVars, const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode bftMaterialHypoElasticName, int nStateVarsUmat, const double* materialProperties, int nMaterialProperties){
            return new UelDisplacement<2, 8 >(coordinates, stateVars, nStateVars, elementProperties, 
                                 nElementProperties, noEl, bftMaterialHypoElasticName, nStateVarsUmat, materialProperties, 
                                 nMaterialProperties,  bft::NumIntegration::IntegrationTypes::FullIntegration, 
                                 UelDisplacement<2,8>::SectionType::PlaneStrain); }

    BftUel* generateUelCPS4R(const double* coordinates, double* stateVars, int nStateVars, const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode bftMaterialHypoElasticName, int nStateVarsUmat, const double* materialProperties, int nMaterialProperties){
            return new UelDisplacement< 2,4 >(coordinates, stateVars, nStateVars, elementProperties, 
                                 nElementProperties, noEl, bftMaterialHypoElasticName, nStateVarsUmat, materialProperties, 
                                 nMaterialProperties,  bft::NumIntegration::IntegrationTypes::ReducedIntegration, 
                                 UelDisplacement<2,4>::SectionType::PlaneStress); }

    BftUel* generateUelCPS8R(const double* coordinates, double* stateVars, int nStateVars, const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode bftMaterialHypoElasticName, int nStateVarsUmat, const double* materialProperties, int nMaterialProperties){
            return new UelDisplacement< 2,8 >(coordinates, stateVars, nStateVars, elementProperties, 
                                 nElementProperties, noEl, bftMaterialHypoElasticName, nStateVarsUmat, materialProperties, 
                                 nMaterialProperties,  bft::NumIntegration::IntegrationTypes::ReducedIntegration, 
                                 UelDisplacement<2,8>::SectionType::PlaneStress); }

    BftUel* generateUelCPE4R(const double* coordinates, double* stateVars, int nStateVars, const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode bftMaterialHypoElasticName, int nStateVarsUmat, const double* materialProperties, int nMaterialProperties){
            return new UelDisplacement< 2,4 >(coordinates, stateVars, nStateVars, elementProperties, 
                                 nElementProperties, noEl, bftMaterialHypoElasticName, nStateVarsUmat, materialProperties, 
                                 nMaterialProperties,  bft::NumIntegration::IntegrationTypes::ReducedIntegration, 
                                 UelDisplacement<2,4>::SectionType::PlaneStrain); }

    BftUel* generateUelCPE8R(const double* coordinates, double* stateVars, int nStateVars, const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode bftMaterialHypoElasticName, int nStateVarsUmat, const double* materialProperties, int nMaterialProperties){
            return new UelDisplacement< 2,8 >(coordinates, stateVars, nStateVars, elementProperties, 
                                 nElementProperties, noEl, bftMaterialHypoElasticName, nStateVarsUmat, materialProperties, 
                                 nMaterialProperties,  bft::NumIntegration::IntegrationTypes::ReducedIntegration, 
                                 UelDisplacement<2,8>::SectionType::PlaneStrain); }

    BftUel* generateUelC3D8(const double* coordinates, double* stateVars, int nStateVars, const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode bftMaterialHypoElasticName, int nStateVarsUmat, const double* materialProperties, int nMaterialProperties){
            return new UelDisplacement< 3, 8 >(coordinates, stateVars, nStateVars, elementProperties, 
                                 nElementProperties, noEl, bftMaterialHypoElasticName, nStateVarsUmat, materialProperties, 
                                 nMaterialProperties,  bft::NumIntegration::IntegrationTypes::FullIntegration,
                                 UelDisplacement<3,8>::SectionType::Solid); }

    BftUel* generateUelC3D8R(const double* coordinates, double* stateVars, int nStateVars, const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode bftMaterialHypoElasticName, int nStateVarsUmat, const double* materialProperties, int nMaterialProperties){
            return new UelDisplacement< 3, 8 >(coordinates, stateVars, nStateVars, elementProperties, 
                                 nElementProperties, noEl, bftMaterialHypoElasticName, nStateVarsUmat, materialProperties, 
                                 nMaterialProperties,  bft::NumIntegration::IntegrationTypes::ReducedIntegration,
                                 UelDisplacement<3,8>::SectionType::Solid); }

    BftUel* generateUelC3D20(const double* coordinates, double* stateVars, int nStateVars, const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode bftMaterialHypoElasticName, int nStateVarsUmat, const double* materialProperties, int nMaterialProperties){
            return new UelDisplacement< 3, 20 >(coordinates, stateVars, nStateVars, elementProperties, 
                                 nElementProperties, noEl, bftMaterialHypoElasticName, nStateVarsUmat, materialProperties, 
                                 nMaterialProperties,  bft::NumIntegration::IntegrationTypes::FullIntegration,
                                 UelDisplacement<3,20>::SectionType::Solid); }

    BftUel* generateUelC3D20R(const double* coordinates, double* stateVars, int nStateVars, const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode bftMaterialHypoElasticName, int nStateVarsUmat, const double* materialProperties, int nMaterialProperties){
            return new UelDisplacement< 3, 20 >(coordinates, stateVars, nStateVars, elementProperties, 
                                 nElementProperties, noEl, bftMaterialHypoElasticName, nStateVarsUmat, materialProperties, 
                                 nMaterialProperties,  bft::NumIntegration::IntegrationTypes::ReducedIntegration,
                                 UelDisplacement<3,20>::SectionType::Solid); }



}
