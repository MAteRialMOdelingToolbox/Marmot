#pragma once
#include "bftTypedefs.h"
#include "bftVoigt.h"
#include "bftConstants.h"
#include "bftUel.h"
#include "bftFiniteElement.h"
#include "bftGeometryElement.h"
#include "bftMaterialHypoElastic.h"
#include "userLibrary.h"
#include <iostream>
#include <vector>
#include <memory>

template <int nDim, int nNodes>
class UelDisplacement: public BftUel, public BftGeometryElement<nDim, nNodes>{
    public:

        typedef BftGeometryElement<nDim, nNodes> ParentGeometryElement;

        static constexpr int sizeLoadVector =   nNodes * nDim;
        static constexpr int nCoordinates =     nNodes * nDim;

        static constexpr int nStateVarsElement=     0;
        static constexpr int nStateVarsStressStrain = 12;

        typedef Matrix<double, sizeLoadVector, 1>                RhsSized;
        typedef Matrix<double, sizeLoadVector, sizeLoadVector>   KeSizedMatrix;

        Map<VectorXd>               stateVars;
        const int                   nStateVars;
        const Map<const VectorXd>   elementProperties;
        const Map<const VectorXd>   materialProperties;
        const int                   elLabel;
        //const bft::pUmatType        umat;
        const int                   nStateVarsMaterial;

        MatrixXd gaussPointList;
        VectorXd gaussWeights;

        std::vector< double >                                       detJAtGauss;        
        std::vector< typename ParentGeometryElement::JacobianSized> JAtGauss;
        std::vector< typename ParentGeometryElement::JacobianSized> JInvAtGauss;
        std::vector< typename ParentGeometryElement::dNdXiSized   > dNdXiAtGauss;
        std::vector< typename ParentGeometryElement::dNdXiSized   > dNdXAtGauss;
        std::vector< typename ParentGeometryElement::BSized       > BAtGauss;
        std::vector< std::unique_ptr< BftMaterialHypoElastic>     > materialAtGauss;

    public:

        UelDisplacement(const double* coordinates,
                double* stateVars,
                int nStateVars,
                const double* elementProperties,
                int nElementPropertiesElement,
                int noEl,
                const std::string& bftMaterialHypoElasticName,
                int nStateVarsMaterial,
                const double* materialProperties,
                int nMaterialProperties,
                bft::NumIntegration::IntegrationTypes integrationType
                );

        virtual void setInitialConditions(StateTypes state, const double* values);

        virtual void computeDistributedLoad( BftUel::DistributedLoadTypes loadType,
                double* P, 
                const int elementFace, 
                const double* load,
                const double* time,
                double dT)
        {}

        virtual void computeYourself( const double* QTotal,
                const double* dQ,
                double* Pe,
                double* Ke,
                const double* time,
                double dT,
                double& pNewdT) = 0;

        Ref<VectorXd> stateVarsMaterialAtGauss(int gaussPt)
        { 
            return stateVars.segment(gaussPt * (nStateVarsMaterial + nStateVarsStressStrain ), nStateVarsMaterial); 
        }

        Ref<bft::Vector6> stressAtGauss(int gaussPt)
        { 
            return stateVars.segment( nStateVarsMaterial + gaussPt * (nStateVarsMaterial + nStateVarsStressStrain ), 6); 
        }
        Ref<bft::Vector6>  strainAtGauss(int gaussPt)
        { 
            return stateVars.segment( nStateVarsMaterial + gaussPt * (nStateVarsMaterial + nStateVarsStressStrain ) + 6 ,6);
        }

        double* getPermanentResultPointer(const std::string& resultName, int gaussPt, int& resultLength)
        {
            if (resultName == "stress" ){
                resultLength = bft::Vgt::VoigtSize;
                return stressAtGauss(gaussPt).data();}
            else if (resultName == "strain" ){ 
                resultLength = bft::Vgt::VoigtSize;
                return strainAtGauss(gaussPt).data();}
            else if( resultName == "sdv"){
                resultLength = nStateVarsMaterial;
                return stateVarsMaterialAtGauss(gaussPt).data();}
            else
                return this->materialAtGauss[gaussPt]->getPermanentResultPointer(resultName, resultLength);
        }
};

template <int nDim, int nNodes>
UelDisplacement<nDim, nNodes>::UelDisplacement(const double* coords, 
        double* stateVars,
        int nStateVars,
        const double* properties,
        int nElementProperties,
        int noEl,
        //const bft::pUmatType umat,
        const std::string& materialName,
        int nStateVarsMaterial ,
        const double* materialProperties,
        int nMaterialProperties,
        bft::NumIntegration::IntegrationTypes integrationType
        ):
    ParentGeometryElement(coords),
    //stateVars(Map<VectorXd>(stateVars, nStateVars)),
    stateVars(stateVars, nStateVars),
    nStateVars(nStateVars),
    elementProperties(Map<const VectorXd>(properties, nElementProperties)),
    materialProperties(Map<const VectorXd>(materialProperties, nMaterialProperties)),
    elLabel(noEl),
    nStateVarsMaterial(nStateVarsMaterial)
{
    //std::cout << "Creation of Element" << std::endl;
    gaussPointList =    bft::NumIntegration::getGaussPointList(this->shape, integrationType);
    gaussWeights =      bft::NumIntegration::getGaussWeights( this->shape, integrationType);

    for(int i = 0; i < gaussPointList.rows(); i++){
        const typename ParentGeometryElement::XiSized& xiGauss = gaussPointList.row(i);

        dNdXiAtGauss    .push_back		(this->dNdXi(xiGauss));
        JAtGauss		.push_back		(this->Jacobian(dNdXiAtGauss [i] ));
        JInvAtGauss		.push_back		(JAtGauss [i] .inverse());
        dNdXAtGauss		.push_back		(this->dNdX(dNdXiAtGauss [i] , JInvAtGauss [i] ));
        BAtGauss		.push_back		(this->B(dNdXAtGauss [i] ));
        detJAtGauss		.push_back		(JAtGauss [i] .determinant());

        materialAtGauss .push_back      (std::unique_ptr<BftMaterialHypoElastic> (
                                            dynamic_cast<BftMaterialHypoElastic*>(
                                                userLibrary::bftMaterialFactory( 
                                                    materialName, 
                                                    stateVarsMaterialAtGauss(i).data(), 
                                                    nStateVarsMaterial, 
                                                    materialProperties, 
                                                    nMaterialProperties, 
                                                    elLabel, 
                                                    i))));
    }
    //std::cout << "End of Creation of Element" << std::endl;
}

template <int nDim, int nNodes>
void UelDisplacement<nDim, nNodes>::setInitialConditions(StateTypes state, const double* values)
{
    switch(state){

        case BftUel::GeostaticStress: 
            { 
                for(int i = 0; i < this->gaussPointList.rows(); i++) 
                {
                    typename ParentGeometryElement::XiSized coordAtGauss =  this->NB( this->N(this->gaussPointList.row(i))) * this->coordinates; 

                    Ref<bft::Vector6> stress( stressAtGauss(i) );

                    const double sigY1 = values[0];
                    const double sigY2 = values[2];
                    const double y1    = values[1];
                    const double y2    = values[3];

                    stress(1) = bft::Functions::linearInterpolation(coordAtGauss[1], y1, y2, sigY1, sigY2);  // sigma_y
                    stress(0) = values[4]*stress(1);  // sigma_x
                    stress(2) = values[5]*stress(1);}  // sigma_z

                    break; 
            }

        default: break;
    }
}

