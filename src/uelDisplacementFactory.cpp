#include "uelDisplacementFactory.h"
#include "uelDisplacementPlane.h"
#include "uelDisplacementSolid.h"

namespace UelDisplacementFactory{
    
    BftUel* generateUelCPS4(const double* coordinates, double* stateVars, int nStateVars, const double* elementProperties, int nElementProperties,
                            int noEl, const std::string& materialName, int nStateVarsUmat, const double* materialProperties, int nMaterialProperties){

            return new UelDisplacementPlane< 4 >(coordinates, stateVars, nStateVars, elementProperties, 
                                 nElementProperties, noEl, materialName , nStateVarsUmat, materialProperties, 
                                 nMaterialProperties,  bft::NumIntegration::IntegrationTypes::FullIntegration, 
                                 UelDisplacementPlane<4>::SectionType::PlaneStress); }

    //BftUel* generateUelCPS8(const double* coordinates, double* stateVars, int nStateVars, const double* elementProperties, int nElementProperties,
                            //int noEl, const std::string& bftMaterialHypoElasticName, int nStateVarsUmat, const double* materialProperties, int nMaterialProperties){
            //return new UelDisplacementPlane< 8 >(coordinates, stateVars, nStateVars, elementProperties, 
                                 //nElementProperties, noEl, bftMaterialHypoElasticName, nStateVarsUmat, materialProperties, 
                                 //nMaterialProperties,  bft::NumIntegration::IntegrationTypes::FullIntegration, 
                                 //UelDisplacementPlane<8>::SectionType::PlaneStress); }

    BftUel* generateUelCPE4(const double* coordinates, double* stateVars, int nStateVars, const double* elementProperties, int nElementProperties,
                            int noEl, const std::string& bftMaterialHypoElasticName, int nStateVarsUmat, const double* materialProperties, int nMaterialProperties){
            return new UelDisplacementPlane< 4 >(coordinates, stateVars, nStateVars, elementProperties, 
                                 nElementProperties, noEl, bftMaterialHypoElasticName, nStateVarsUmat, materialProperties, 
                                 nMaterialProperties,  bft::NumIntegration::IntegrationTypes::FullIntegration, 
                                 UelDisplacementPlane<4>::SectionType::PlaneStrain); }

    //BftUel* generateUelCPE8(const double* coordinates, double* stateVars, int nStateVars, const double* elementProperties, int nElementProperties,
                            //int noEl, const std::string& bftMaterialHypoElasticName, int nStateVarsUmat, const double* materialProperties, int nMaterialProperties){
            //return new UelDisplacementPlane< 8 >(coordinates, stateVars, nStateVars, elementProperties, 
                                 //nElementProperties, noEl, bftMaterialHypoElasticName, nStateVarsUmat, materialProperties, 
                                 //nMaterialProperties,  bft::NumIntegration::IntegrationTypes::FullIntegration, 
                                 //UelDisplacementPlane<8>::SectionType::PlaneStrain); }

    //BftUel* generateUelCPS4R(const double* coordinates, double* stateVars, int nStateVars, const double* elementProperties, int nElementProperties,
                            //int noEl, const std::string& bftMaterialHypoElasticName, int nStateVarsUmat, const double* materialProperties, int nMaterialProperties){
            //return new UelDisplacementPlane< 4 >(coordinates, stateVars, nStateVars, elementProperties, 
                                 //nElementProperties, noEl, bftMaterialHypoElasticName, nStateVarsUmat, materialProperties, 
                                 //nMaterialProperties,  bft::NumIntegration::IntegrationTypes::ReducedIntegration, 
                                 //UelDisplacementPlane<4>::SectionType::PlaneStress); }

    //BftUel* generateUelCPS8R(const double* coordinates, double* stateVars, int nStateVars, const double* elementProperties, int nElementProperties,
                            //int noEl, const std::string& bftMaterialHypoElasticName, int nStateVarsUmat, const double* materialProperties, int nMaterialProperties){
            //return new UelDisplacementPlane< 8 >(coordinates, stateVars, nStateVars, elementProperties, 
                                 //nElementProperties, noEl, bftMaterialHypoElasticName, nStateVarsUmat, materialProperties, 
                                 //nMaterialProperties,  bft::NumIntegration::IntegrationTypes::ReducedIntegration, 
                                 //UelDisplacementPlane<8>::SectionType::PlaneStress); }

    //BftUel* generateUelCPE4R(const double* coordinates, double* stateVars, int nStateVars, const double* elementProperties, int nElementProperties,
                            //int noEl, const std::string& bftMaterialHypoElasticName, int nStateVarsUmat, const double* materialProperties, int nMaterialProperties){
            //return new UelDisplacementPlane< 4 >(coordinates, stateVars, nStateVars, elementProperties, 
                                 //nElementProperties, noEl, bftMaterialHypoElasticName, nStateVarsUmat, materialProperties, 
                                 //nMaterialProperties,  bft::NumIntegration::IntegrationTypes::ReducedIntegration, 
                                 //UelDisplacementPlane<4>::SectionType::PlaneStrain); }

    //BftUel* generateUelCPE8R(const double* coordinates, double* stateVars, int nStateVars, const double* elementProperties, int nElementProperties,
                            //int noEl, const std::string& bftMaterialHypoElasticName, int nStateVarsUmat, const double* materialProperties, int nMaterialProperties){
            //return new UelDisplacementPlane< 8 >(coordinates, stateVars, nStateVars, elementProperties, 
                                 //nElementProperties, noEl, bftMaterialHypoElasticName, nStateVarsUmat, materialProperties, 
                                 //nMaterialProperties,  bft::NumIntegration::IntegrationTypes::ReducedIntegration, 
                                 //UelDisplacementPlane<8>::SectionType::PlaneStrain); }

    //BftUel* generateUelC3D8(const double* coordinates, double* stateVars, int nStateVars, const double* elementProperties, int nElementProperties,
                            //int noEl, const std::string& bftMaterialHypoElasticName, int nStateVarsUmat, const double* materialProperties, int nMaterialProperties){
            //return new UelDisplacementSolid< 8 >(coordinates, stateVars, nStateVars, elementProperties, 
                                 //nElementProperties, noEl, bftMaterialHypoElasticName, nStateVarsUmat, materialProperties, 
                                 //nMaterialProperties,  bft::NumIntegration::IntegrationTypes::FullIntegration); }

    //BftUel* generateUelC3D8R(const double* coordinates, double* stateVars, int nStateVars, const double* elementProperties, int nElementProperties,
                            //int noEl, const std::string& bftMaterialHypoElasticName, int nStateVarsUmat, const double* materialProperties, int nMaterialProperties){
            //return new UelDisplacementSolid< 8 >(coordinates, stateVars, nStateVars, elementProperties, 
                                 //nElementProperties, noEl, bftMaterialHypoElasticName, nStateVarsUmat, materialProperties, 
                                 //nMaterialProperties,  bft::NumIntegration::IntegrationTypes::ReducedIntegration); }

    //BftUel* generateUelC3D20(const double* coordinates, double* stateVars, int nStateVars, const double* elementProperties, int nElementProperties,
                            //int noEl, const std::string& bftMaterialHypoElasticName, int nStateVarsUmat, const double* materialProperties, int nMaterialProperties){
            //return new UelDisplacementSolid< 20 >(coordinates, stateVars, nStateVars, elementProperties, 
                                 //nElementProperties, noEl, bftMaterialHypoElasticName, nStateVarsUmat, materialProperties, 
                                 //nMaterialProperties,  bft::NumIntegration::IntegrationTypes::FullIntegration); }

    //BftUel* generateUelC3D20R(const double* coordinates, double* stateVars, int nStateVars, const double* elementProperties, int nElementProperties,
                            //int noEl, const std::string& bftMaterialHypoElasticName, int nStateVarsUmat, const double* materialProperties, int nMaterialProperties){
            //return new UelDisplacementSolid< 20 >(coordinates, stateVars, nStateVars, elementProperties, 
                                 //nElementProperties, noEl, bftMaterialHypoElasticName, nStateVarsUmat, materialProperties, 
                                 //nMaterialProperties,  bft::NumIntegration::IntegrationTypes::ReducedIntegration); }



}
