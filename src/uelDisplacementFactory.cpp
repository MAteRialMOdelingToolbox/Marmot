#include "uelDisplacementFactory.h"
#include "uelDisplacement.h"
#include "bftFiniteElementUelSpatialWrapper.h"

namespace UelDisplacementFactory{

    BftUel* generateUelT2D2( const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode materialCode, const double* materialProperties, int nMaterialProperties)
    {
            auto uelT2D2 =  std::unique_ptr<BftUel> ( new UelDisplacement<1, 2 >(elementProperties, 
                                 nElementProperties, noEl, materialCode ,  materialProperties, 
                                 nMaterialProperties,  bft::NumIntegration::IntegrationTypes::FullIntegration, 
                                 UelDisplacement<1, 2>::SectionType::UniaxialStress) );

            constexpr int indicesToBeWrapped[] = {0,1};
            constexpr int nIndicesToBeWrapped = 2;
            return new BftUelSpatialWrapper(2, 1, 2, 2, indicesToBeWrapped, nIndicesToBeWrapped, std::move( uelT2D2) );
    }
    
    BftUel* generateUelCPS4( const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode materialCode, const double* materialProperties, int nMaterialProperties){

            return new UelDisplacement<2, 4 >(elementProperties, 
                                 nElementProperties, noEl, materialCode, materialProperties, 
                                 nMaterialProperties,  bft::NumIntegration::IntegrationTypes::FullIntegration, 
                                 UelDisplacement<2, 4>::SectionType::PlaneStress); }

    BftUel* generateUelCPS8( const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode bftMaterialHypoElasticName, const double* materialProperties, int nMaterialProperties){
            return new UelDisplacement<2, 8 >(elementProperties, 
                                 nElementProperties, noEl, bftMaterialHypoElasticName, materialProperties, 
                                 nMaterialProperties,  bft::NumIntegration::IntegrationTypes::FullIntegration, 
                                 UelDisplacement<2,8>::SectionType::PlaneStress); }

    BftUel* generateUelCPE4( const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode bftMaterialHypoElasticName, const double* materialProperties, int nMaterialProperties){
            return new UelDisplacement<2, 4 >(elementProperties, 
                                 nElementProperties, noEl, bftMaterialHypoElasticName, materialProperties, 
                                 nMaterialProperties,  bft::NumIntegration::IntegrationTypes::FullIntegration, 
                                 UelDisplacement<2,4>::SectionType::PlaneStrain); }

    BftUel* generateUelCPE8( const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode bftMaterialHypoElasticName, const double* materialProperties, int nMaterialProperties){
            return new UelDisplacement<2, 8 >(elementProperties, 
                                 nElementProperties, noEl, bftMaterialHypoElasticName, materialProperties, 
                                 nMaterialProperties,  bft::NumIntegration::IntegrationTypes::FullIntegration, 
                                 UelDisplacement<2,8>::SectionType::PlaneStrain); }

    BftUel* generateUelCPS4R( const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode bftMaterialHypoElasticName, const double* materialProperties, int nMaterialProperties){
            return new UelDisplacement< 2,4 >(elementProperties, 
                                 nElementProperties, noEl, bftMaterialHypoElasticName, materialProperties, 
                                 nMaterialProperties,  bft::NumIntegration::IntegrationTypes::ReducedIntegration, 
                                 UelDisplacement<2,4>::SectionType::PlaneStress); }

    BftUel* generateUelCPS8R( const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode bftMaterialHypoElasticName, const double* materialProperties, int nMaterialProperties){
            return new UelDisplacement< 2,8 >(elementProperties, 
                                 nElementProperties, noEl, bftMaterialHypoElasticName, materialProperties, 
                                 nMaterialProperties,  bft::NumIntegration::IntegrationTypes::ReducedIntegration, 
                                 UelDisplacement<2,8>::SectionType::PlaneStress); }

    BftUel* generateUelCPE4R( const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode bftMaterialHypoElasticName, const double* materialProperties, int nMaterialProperties){
            return new UelDisplacement< 2,4 >(elementProperties, 
                                 nElementProperties, noEl, bftMaterialHypoElasticName, materialProperties, 
                                 nMaterialProperties,  bft::NumIntegration::IntegrationTypes::ReducedIntegration, 
                                 UelDisplacement<2,4>::SectionType::PlaneStrain); }

    BftUel* generateUelCPE8R( const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode bftMaterialHypoElasticName, const double* materialProperties, int nMaterialProperties){
            return new UelDisplacement< 2,8 >(elementProperties, 
                                 nElementProperties, noEl, bftMaterialHypoElasticName, materialProperties, 
                                 nMaterialProperties,  bft::NumIntegration::IntegrationTypes::ReducedIntegration, 
                                 UelDisplacement<2,8>::SectionType::PlaneStrain); }

    BftUel* generateUelC3D8( const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode bftMaterialHypoElasticName, const double* materialProperties, int nMaterialProperties){
            return new UelDisplacement< 3, 8 >(elementProperties, 
                                 nElementProperties, noEl, bftMaterialHypoElasticName, materialProperties, 
                                 nMaterialProperties,  bft::NumIntegration::IntegrationTypes::FullIntegration,
                                 UelDisplacement<3,8>::SectionType::Solid); }

    BftUel* generateUelC3D8R( const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode bftMaterialHypoElasticName, const double* materialProperties, int nMaterialProperties){
            return new UelDisplacement< 3, 8 >(elementProperties, 
                                 nElementProperties, noEl, bftMaterialHypoElasticName, materialProperties, 
                                 nMaterialProperties,  bft::NumIntegration::IntegrationTypes::ReducedIntegration,
                                 UelDisplacement<3,8>::SectionType::Solid); }

    BftUel* generateUelC3D20( const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode bftMaterialHypoElasticName, const double* materialProperties, int nMaterialProperties){
            return new UelDisplacement< 3, 20 >(elementProperties, 
                                 nElementProperties, noEl, bftMaterialHypoElasticName, materialProperties, 
                                 nMaterialProperties,  bft::NumIntegration::IntegrationTypes::FullIntegration,
                                 UelDisplacement<3,20>::SectionType::Solid); }

    BftUel* generateUelC3D20R( const double* elementProperties, int nElementProperties,
                            int noEl, userLibrary::MaterialCode bftMaterialHypoElasticName, const double* materialProperties, int nMaterialProperties){
            return new UelDisplacement< 3, 20 >(elementProperties, 
                                 nElementProperties, noEl, bftMaterialHypoElasticName, materialProperties, 
                                 nMaterialProperties,  bft::NumIntegration::IntegrationTypes::ReducedIntegration,
                                 UelDisplacement<3,20>::SectionType::Solid); }
}
