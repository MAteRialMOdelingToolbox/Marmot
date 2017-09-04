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
        const int                   nStateVarsMaterial;

        struct GaussPt {
            Ref< bft::Vector6 > stress;
            Ref< bft::Vector6 > strain;
            Ref< VectorXd >     stateVarsMaterial;
            const double weight;
            typename ParentGeometryElement::XiSized       xi;
            typename ParentGeometryElement::JacobianSized J;
            typename ParentGeometryElement::JacobianSized JInv;
            double                                        detJ;
            typename ParentGeometryElement::dNdXiSized    dNdXi;
            typename ParentGeometryElement::dNdXiSized    dNdX;
            typename ParentGeometryElement::BSized        B;
            std::unique_ptr< BftMaterialHypoElastic>      material;

            GaussPt(Ref< bft::Vector6> stress,
                    Ref< bft::Vector6> strain,
                    Ref< VectorXd> stateVarsMaterial ,
                    double weight
                    ):
                stress(stress),
                strain(strain),
                stateVarsMaterial(stateVarsMaterial),
                weight(weight)
            {};
        };

        std::vector < GaussPt > gaussPts;

    public:

        UelDisplacement(const double* coordinates,
                double* stateVars,
                int nStateVars,
                const double* elementProperties,
                int nElementPropertiesElement,
                int noEl,
                userLibrary::MaterialCode material,
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

        Ref<VectorXd> stateVarsMaterialAtGauss(int gaussPt) { 
            return stateVars.segment(gaussPt * (nStateVarsMaterial + nStateVarsStressStrain ), nStateVarsMaterial); }

        Ref<bft::Vector6> stressAtGauss(int gaussPt) { 
            return stateVars.segment( nStateVarsMaterial + gaussPt * (nStateVarsMaterial + nStateVarsStressStrain ), 6); }

        Ref<bft::Vector6>  strainAtGauss(int gaussPt) { 
            return stateVars.segment( nStateVarsMaterial + gaussPt * (nStateVarsMaterial + nStateVarsStressStrain ) + 6 ,6); }

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
                return this->gaussPts[gaussPt].material->getPermanentResultPointer(resultName, resultLength);
        }
};

template <int nDim, int nNodes>
UelDisplacement<nDim, nNodes>::UelDisplacement(const double* coords, 
        double* stateVars,
        int nStateVars,
        const double* properties,
        int nElementProperties,
        int noEl,
        userLibrary::MaterialCode material,
        int nStateVarsMaterial ,
        const double* materialProperties,
        int nMaterialProperties,
        bft::NumIntegration::IntegrationTypes integrationType
        ):
    ParentGeometryElement(coords),
    stateVars(stateVars, nStateVars),
    nStateVars(nStateVars),
    elementProperties(Map<const VectorXd>(properties, nElementProperties)),
    materialProperties(Map<const VectorXd>(materialProperties, nMaterialProperties)),
    elLabel(noEl),
    nStateVarsMaterial(nStateVarsMaterial)
{
    MatrixXd gaussPointList =    bft::NumIntegration::getGaussPointList(this->shape, integrationType);
    VectorXd gaussWeights =      bft::NumIntegration::getGaussWeights( this->shape, integrationType);

    for(int i = 0; i < gaussPointList.rows(); i++){

        const typename ParentGeometryElement::XiSized& xi = gaussPointList.row(i);
        
        GaussPt gpt(stressAtGauss(i), 
                    strainAtGauss(i), 
                    stateVarsMaterialAtGauss(i), 
                    gaussWeights(i) );
        
        gpt.xi      =   xi;
        gpt.dNdXi   =   this->dNdXi(xi);
        gpt.J		=   this->Jacobian(gpt.dNdXi  );
        gpt.JInv    =   gpt.J.inverse();
        gpt.detJ	=   gpt.J.determinant();
        gpt.dNdX	=   this->dNdX(gpt.dNdXi, gpt.JInv );
        gpt.B		=   this->B(gpt.dNdX);
        gpt.material=   std::unique_ptr<BftMaterialHypoElastic>(
                                            dynamic_cast<BftMaterialHypoElastic*>(
                                                userLibrary::bftMaterialFactory( 
                                                    material,
                                                    stateVarsMaterialAtGauss(i).data(), 
                                                    nStateVarsMaterial, 
                                                    materialProperties, 
                                                    nMaterialProperties, 
                                                    elLabel, 
                                                    i)));
       gaussPts.push_back ( std::move (gpt) );
    }
}

template <int nDim, int nNodes>
void UelDisplacement<nDim, nNodes>::setInitialConditions(StateTypes state, const double* values)
{
    switch(state){

        case BftUel::GeostaticStress: 
            { 
                for(size_t i = 0; i < this->gaussPts.size(); i++) 
                {
                    GaussPt& gaussPt = this->gaussPts[i];
                    typename ParentGeometryElement::XiSized coordAtGauss =  this->NB( this->N(gaussPt.xi)) * this->coordinates; 

                    const double sigY1 = values[0];
                    const double sigY2 = values[2];
                    const double y1    = values[1];
                    const double y2    = values[3];

                    gaussPt.stress(1) = bft::Functions::linearInterpolation(coordAtGauss[1], y1, y2, sigY1, sigY2);  // sigma_y
                    gaussPt.stress(0) = values[4]*gaussPt.stress(1);  // sigma_x
                    gaussPt.stress(2) = values[5]*gaussPt.stress(1);}  // sigma_z

                    break; 
            }

        default: break;
    }
}

