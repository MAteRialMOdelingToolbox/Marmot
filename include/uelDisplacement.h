#pragma once
#include "bftTypedefs.h"
#include "bftVoigt.h"
#include "bftConstants.h"
#include "bftUel.h"
#include "bftFunctions.h"
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

    enum SectionType {
        PlaneStress,
        PlaneStrain,
        Solid,
    };

        typedef BftGeometryElement<nDim, nNodes> ParentGeometryElement;

        static constexpr int sizeLoadVector =   nNodes * nDim;
        static constexpr int nCoordinates =     nNodes * nDim;
        static constexpr int nStateVarsGPtAdditional = 12;

        typedef Matrix<double, sizeLoadVector, 1>                RhsSized;
        typedef Matrix<double, sizeLoadVector, sizeLoadVector>   KeSizedMatrix;

        typedef Matrix<double, ParentGeometryElement::VoigtSize, ParentGeometryElement::VoigtSize> CSized;
        typedef Matrix<double, ParentGeometryElement::VoigtSize, 1> Voigt;

        Map<VectorXd>               stateVars;
        const int                   nStateVars;
        const Map<const VectorXd>   elementProperties;
        const Map<const VectorXd>   materialProperties;
        const int                   elLabel;
        const int                   nStateVarsMaterial;
        SectionType                 sectionType;

        struct GaussPt {
            Ref< VectorXd >     stateVarsMaterial;
            Ref< bft::Vector6 > stress;
            Ref< bft::Vector6 > strain;
            const double weight;
            typename ParentGeometryElement::XiSized       xi;
            typename ParentGeometryElement::JacobianSized J;
            typename ParentGeometryElement::JacobianSized JInv;
            double                                        detJ;
            typename ParentGeometryElement::dNdXiSized    dNdXi;
            typename ParentGeometryElement::dNdXiSized    dNdX;
            typename ParentGeometryElement::BSized        B;
            std::unique_ptr< BftMaterialHypoElastic>      material;
            double intVol;

            GaussPt(Ref< VectorXd> stateVarsMaterial ,
                    Ref< bft::Vector6> stress,
                    Ref< bft::Vector6> strain,
                    double weight
                    ):
                stateVarsMaterial(stateVarsMaterial),
                stress(stress),
                strain(strain),
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
                bft::NumIntegration::IntegrationTypes integrationType,
                SectionType sectionType
                );

        virtual void setInitialConditions(StateTypes state, const double* values);

        virtual void computeDistributedLoad( BftUel::DistributedLoadTypes loadType,
                double* P, 
                const int elementFace, 
                const double* load,
                const double* time,
                double dT);

        virtual void computeYourself( const double* QTotal,
                const double* dQ,
                double* Pe,
                double* Ke,
                const double* time,
                double dT,
                double& pNewdT) ;

        double* getPermanentResultPointer(const std::string& resultName, int gaussPt, int& resultLength)
        {
            if (resultName == "stress" ){
                resultLength = bft::Vgt::VoigtSize;
                return gaussPts[gaussPt].stress.data(); }
            else if (resultName == "strain" ){ 
                resultLength = bft::Vgt::VoigtSize;
                return gaussPts[gaussPt].strain.data(); }
            else if( resultName == "sdv"){
                resultLength = nStateVarsMaterial;
                return gaussPts[gaussPt].stateVarsMaterial.data(); }
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
        bft::NumIntegration::IntegrationTypes integrationType,
        SectionType sectionType
        ):
    ParentGeometryElement(coords),
    stateVars(stateVars, nStateVars),
    nStateVars(nStateVars),
    elementProperties(Map<const VectorXd>(properties, nElementProperties)),
    materialProperties(Map<const VectorXd>(materialProperties, nMaterialProperties)),
    elLabel(noEl),
    nStateVarsMaterial(nStateVarsMaterial),
    sectionType(sectionType)
{

    MatrixXd gaussPointList =    bft::NumIntegration::getGaussPointList(this->shape, integrationType);
    VectorXd gaussWeights =      bft::NumIntegration::getGaussWeights( this->shape, integrationType);

    for(int i = 0; i < gaussPointList.rows(); i++){

        const typename ParentGeometryElement::XiSized& xi = gaussPointList.row(i);
        
        GaussPt gpt(this->stateVars.segment(                      i* (nStateVarsMaterial + nStateVarsGPtAdditional ), nStateVarsMaterial),
                    this->stateVars.segment( nStateVarsMaterial + i* (nStateVarsMaterial + nStateVarsGPtAdditional ), 6),
                    this->stateVars.segment( nStateVarsMaterial + i* (nStateVarsMaterial + nStateVarsGPtAdditional ) + 6 ,6),
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
                                                    gpt.stateVarsMaterial.data(), 
                                                    nStateVarsMaterial, 
                                                    materialProperties, 
                                                    nMaterialProperties, 
                                                    elLabel, 
                                                    i)));

        if( sectionType == SectionType::Solid){

            gpt.intVol = gpt.weight * gpt.detJ;
            gpt.material->setCharacteristicElementLength ( std::cbrt ( 8 * gpt.detJ )  ) ;}

        else if(sectionType == SectionType::PlaneStrain || 
                sectionType == SectionType::PlaneStress){

            const double& thickness = elementProperties[0];
            gpt.intVol = gpt.weight * gpt.detJ * thickness;
            gpt.material->setCharacteristicElementLength ( std::sqrt ( 4 * gpt.detJ )  ) ;}

       gaussPts.push_back ( std::move (gpt) );
    }

}
template <int nDim, int nNodes>
void UelDisplacement<nDim, nNodes>::computeYourself( const double* QTotal_,
        const double* dQ_,
        double* Pe_,
        double* Ke_,
        const double* time,
        double dT,
        double& pNewDT
        )
{

    using namespace bft;

    Map<const RhsSized>                  QTotal(QTotal_);				
    Map<const RhsSized>                  dQ(dQ_);			
    Map<KeSizedMatrix>                   Ke(Ke_);
    Map<RhsSized>                        Pe(Pe_);

    Ke.setZero();
    Pe.setZero();

    Voigt S, dE;
    CSized C;

    for(size_t i = 0; i < this->gaussPts.size(); i++){
        
        GaussPt& gaussPt = this->gaussPts[i];

        const typename ParentGeometryElement::BSized& B = gaussPt.B;
       
        dE = B * dQ;

        if constexpr (nDim == 2) {

            Vector6 dE6 = Vgt::planeVoigtToVoigt(  dE  ); 
            Matrix6 C66;

            if (sectionType == SectionType::PlaneStress) { 

                gaussPt.material->computePlaneStress(gaussPt.stress.data(), 
                                                     C66.data(),  
                                                     gaussPt.strain.data(),
                                                     dE6.data(), 
                                                     time, dT, pNewDT);
                
                C  = mechanics::getPlaneStressTangent(C66); }

            else if(sectionType == SectionType::PlaneStrain) {

                gaussPt.material->computeStress(gaussPt.stress.data(), 
                                                        C66.data(),  
                                                        gaussPt.strain.data(),
                                                        dE6.data(), 
                                                        time, dT, pNewDT);
                
                C =  mechanics::getPlaneStrainTangent(C66); 
            }

            S =  Vgt::voigtToPlaneVoigt(gaussPt.stress); 
            gaussPt.strain += dE6; 
        }

        else if constexpr (nDim==3) {

            if(sectionType == SectionType::Solid){

                gaussPt.material->computeStress(gaussPt.stress.data(), 
                                                C.data(),  
                                                gaussPt.strain.data(),
                                                dE.data(), 
                                                time, dT, pNewDT);
            }

            S = gaussPt.stress;
            gaussPt.strain += dE; 
        }

        if (pNewDT<1.0)
            return;

        Ke += B.transpose() * C * B * gaussPt.intVol;
        Pe -= B.transpose() * S * gaussPt.intVol;

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

template <int nDim, int nNodes>
void UelDisplacement<nDim, nNodes>::computeDistributedLoad( BftUel::DistributedLoadTypes loadType,
        double* P, 
        const int elementFace, 
        const double* load,
        const double* time,
        double dT){

    Map<RhsSized> fU(P);

    using namespace bft::FiniteElement::BoundaryElementFactory;
    VectorXd boundaryCoordIndices = getBoundaryNodeList(this->shape, elementFace);
    
    VectorXd boundaryCoordinates(boundaryCoordIndices.size());
    for(int i = 0; i < boundaryCoordIndices.size(); i++)
        boundaryCoordinates(i) = this->coordinates( boundaryCoordIndices(i) );

    switch(loadType){

        case BftUel::Pressure: { 
            const double p = load[0];
            
            if (std::abs(p)<bft::Constants::numZeroPos)
                return;

            VectorXd Pk = VectorXd::Zero(boundaryCoordIndices.size());
            MatrixXd gp =       getGaussPointList(this->shape); 
            VectorXd gpWeight = getGaussWeights(this->shape);  

            for(int i=0; i<gp.rows(); i++){
                MatrixXd xi = gp.row(i);        // necessary matrix mapping, as factory return type of gauss points is a matrix (with regard to future 3d elements) 
                VectorXd tractionVec = -p * getNormalVector(this->shape, boundaryCoordinates, xi);
                Pk += getIntVol(this->shape, boundaryCoordinates, xi)  * gpWeight.row(i) * tractionVec.transpose() * getNB(this->shape, xi);}
            
            if(nDim == 2)
                Pk *= elementProperties[0]; // thickness
            
            for(int i = 0; i < boundaryCoordIndices.size(); i++)
                fU( boundaryCoordIndices(i) ) +=  Pk(i);
            
            break;
        }
        default: {throw std::invalid_argument("Invalid Load Type specified");}
    }
}

