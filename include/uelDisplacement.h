#pragma once
#include "bftTypedefs.h"
#include "bftVoigt.h"
#include "bftConstants.h"
#include "bftUel.h"
#include "bftFunctions.h"
#include "bftFiniteElement.h"
#include "bftMath.h"
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
            UniaxialStress,
            PlaneStress,
            PlaneStrain,
            Solid,
        };

        typedef BftGeometryElement<nDim, nNodes> ParentGeometryElement;

        static constexpr int sizeLoadVector =   nNodes * nDim;
        static constexpr int nCoordinates =     nNodes * nDim;
        //static constexpr int nStateVarsGPtAdditional = 12;

        typedef Matrix<double, sizeLoadVector, 1>                RhsSized;
        typedef Matrix<double, sizeLoadVector, sizeLoadVector>   KeSizedMatrix;

        typedef Matrix<double, ParentGeometryElement::VoigtSize, ParentGeometryElement::VoigtSize> CSized;
        typedef Matrix<double, ParentGeometryElement::VoigtSize, 1> Voigt;

        //Map<VectorXd>               stateVars;
        //int                         nStateVars;

        const Map<const VectorXd>   elementProperties;
        const Map<const VectorXd>   materialProperties;
        const int                   elLabel;
        //int                         nStateVarsMaterial;
        const SectionType           sectionType;

        struct GaussPt {
            static constexpr int numberOfRequiredStateVars= 6 + 6;
            //Map< VectorXd >     stateVarsMaterial;
            Map< bft::Vector6 > stress;
            Map< bft::Vector6 > strain;
            const double weight;
            typename ParentGeometryElement::XiSized       xi;
            typename ParentGeometryElement::JacobianSized J;
            typename ParentGeometryElement::JacobianSized JInv;
            double                                        detJ;
            typename ParentGeometryElement::dNdXiSized    dNdXi;
            typename ParentGeometryElement::dNdXiSized    dNdX;
            typename ParentGeometryElement::BSized        B;
            std::unique_ptr< BftMaterialHypoElastic>      material;
            double                                        intVol;

            GaussPt(
                    //Ref< VectorXd> stateVarsMaterial ,
                    //Ref< bft::Vector6> stress,
                    //Ref< bft::Vector6> strain,
                    double weight
                   ):
                //stateVarsMaterial(stateVarsMaterial),
                //stress(stress),
                //strain(strain),
                //stateVarsMaterial(nullptr,0),
                stress(nullptr),
                strain(nullptr),
                weight(weight)
            {};
        };

        std::vector < GaussPt > gaussPts;

    public:

        UelDisplacement(const double* coordinates,
                //double* stateVars,
                //int nStateVars,
                const double* elementProperties,
                int nElementPropertiesElement,
                int noEl,
                userLibrary::MaterialCode material,
                //int nStateVarsMaterial,
                const double* materialProperties,
                int nMaterialProperties,
                bft::NumIntegration::IntegrationTypes integrationType,
                SectionType sectionType
                );

        virtual int getNumberOfRequiredStateVars();

        virtual void assignStateVars(double *stateVars, int nStateVars);

        virtual void initializeYourself();

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
                //resultLength = nStateVarsMaterial;
                //return gaussPts[gaussPt].stateVarsMaterial.data(); 
                resultLength = gaussPts[gaussPt].material->nStateVars;
                return gaussPts[gaussPt].material->stateVars; 
            }
            else
                return this->gaussPts[gaussPt].material->getPermanentResultPointer(resultName, resultLength);
        }
};

template <int nDim, int nNodes>
UelDisplacement<nDim, nNodes>::UelDisplacement(const double* coords, 
        //double* stateVars,
        //int nStateVars,
        const double* properties,
        int nElementProperties,
        int noEl,
        userLibrary::MaterialCode material,
        //int nStateVarsMaterial ,
        const double* materialProperties,
        int nMaterialProperties,
        bft::NumIntegration::IntegrationTypes integrationType,
        SectionType sectionType
        ):
    ParentGeometryElement(coords),
    //stateVars(nullptr, 0),
    //stateVars(stateVars, nStateVars),
    //nStateVars(nStateVars),
    elementProperties(Map<const VectorXd>(properties, nElementProperties)),
    materialProperties(Map<const VectorXd>(materialProperties, nMaterialProperties)),
    elLabel(noEl),
    //nStateVarsMaterial(nStateVarsMaterial),
    sectionType(sectionType)
{

    MatrixXd gaussPointList =    bft::NumIntegration::getGaussPointList(this->shape, integrationType);
    VectorXd gaussWeights =      bft::NumIntegration::getGaussWeights( this->shape, integrationType);

    for(int i = 0; i < gaussPointList.rows(); i++){

        const typename ParentGeometryElement::XiSized& xi = gaussPointList.row(i);

        GaussPt gpt( gaussWeights(i) );

        gpt.xi      =   xi;
        gpt.material =   std::unique_ptr<BftMaterialHypoElastic>(
                dynamic_cast<BftMaterialHypoElastic*>(
                    userLibrary::bftMaterialFactory( 
                        material,
                        //gpt.stateVarsMaterial.data(), 
                        //nStateVarsMaterial, 
                        materialProperties, 
                        nMaterialProperties, 
                        elLabel, 
                        i)));
        gaussPts.push_back ( std::move (gpt) );
    }
}

    template <int nDim, int nNodes>
int UelDisplacement<nDim, nNodes>::getNumberOfRequiredStateVars()
{
    return  ( gaussPts[0].material->getNumberOfRequiredStateVars() + GaussPt::numberOfRequiredStateVars )  * gaussPts.size();
}

    template <int nDim, int nNodes>
void UelDisplacement<nDim, nNodes>::assignStateVars(double *stateVars, int nStateVars)
{
    // assign the state vars using the 'placement new' operator
    //new (&this->stateVars) Map<VectorXd> (stateVars, nStateVars);
    //this->nStateVars = nStateVars;

    int nStateVarsMaterial =  nStateVars / gaussPts.size() - GaussPt::numberOfRequiredStateVars;

    for(size_t i = 0; i < gaussPts.size(); i++){

        GaussPt& gpt = gaussPts[i];
        //new (&gpt.stateVarsMaterial)  Map< VectorXd >      (this->stateVars.segment(                      i* (nStateVarsMaterial + nStateVarsGPtAdditional ), nStateVarsMaterial).data(), nStateVarsMaterial);
        //new (&gpt.stress)             Map< bft::Vector6 >  (this->stateVars.segment( nStateVarsMaterial + i* (nStateVarsMaterial + nStateVarsGPtAdditional ), 6).data());
        //new (&gpt.strain)             Map< bft::Vector6 >  (this->stateVars.segment( nStateVarsMaterial + i* (nStateVarsMaterial + nStateVarsGPtAdditional ) + 6 ,6).data());

        //new (&gpt.stateVarsMaterial)  Map< VectorXd >      (stateVars(                      i* (nStateVarsMaterial + nStateVarsGPtAdditional )), nStateVarsMaterial);
        gpt.material->assignStateVars( 
                stateVars +                      i * (nStateVarsMaterial + GaussPt::numberOfRequiredStateVars), nStateVarsMaterial);

        new (&gpt.stress) Map< bft::Vector6 > (
                stateVars + nStateVarsMaterial + i * (nStateVarsMaterial + GaussPt::numberOfRequiredStateVars), 6);

        new (&gpt.strain) Map< bft::Vector6 > (
                stateVars + nStateVarsMaterial + i * (nStateVarsMaterial + GaussPt::numberOfRequiredStateVars) + 6, 6);

        //gpt.material->assignStateVars( gpt.stateVarsMaterial.data(), nStateVarsMaterial );

    }

}

    template <int nDim, int nNodes>
void UelDisplacement<nDim, nNodes>::initializeYourself()
{
    //for(size_t i = 0; i < gaussPts.size(); i++){

    //GaussPt& gpt = gaussPts[i];

    for( GaussPt& gpt : gaussPts){

        //const typename ParentGeometryElement::XiSized& xi = gpt.xi;

        //gpt.xi      =   gpt.xi;
        gpt.dNdXi   =   this->dNdXi(gpt.xi);
        gpt.J		=   this->Jacobian(gpt.dNdXi  );
        gpt.JInv    =   gpt.J.inverse();
        gpt.detJ	=   gpt.J.determinant();
        gpt.dNdX	=   this->dNdX(gpt.dNdXi, gpt.JInv );
        gpt.B		=   this->B(gpt.dNdX);

        if( sectionType == SectionType::Solid){

            gpt.intVol = gpt.weight * gpt.detJ;
            gpt.material->setCharacteristicElementLength ( std::cbrt ( 8 * gpt.detJ )  ) ;}

        else if(sectionType == SectionType::PlaneStrain || 
                sectionType == SectionType::PlaneStress){

            const double& thickness = elementProperties[0];
            gpt.intVol = gpt.weight * gpt.detJ * thickness;
            gpt.material->setCharacteristicElementLength ( std::sqrt ( 4 * gpt.detJ )  ) ;}

        else if(sectionType == SectionType::UniaxialStress){

            const double& crossSection = elementProperties[0];
            gpt.intVol = gpt.weight * gpt.detJ * crossSection;
            gpt.material->setCharacteristicElementLength (  2 * gpt.detJ  ) ;}
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

    //for(size_t i = 0; i < this->gaussPts.size(); i++){

    //GaussPt& gaussPt = this->gaussPts[i];
    for(GaussPt& gaussPt : gaussPts){  

        const typename ParentGeometryElement::BSized& B = gaussPt.B;

        dE = B * dQ;

        if constexpr (nDim == 1)
        {
            Vector6 dE6; dE6 << dE, 0,0,0,0,0;
            Matrix6 C66;
            gaussPt.material->computeUniaxialStress(gaussPt.stress.data(), 
                    C66.data(),  
                    gaussPt.strain.data(),
                    dE6.data(), 
                    time, dT, pNewDT);

            C <<  mechanics::getUniaxialStressTangent(C66); 
        }

        else if constexpr (nDim == 2) {

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
        case BftUel::GeostaticStress: { 
                    //for(size_t i = 0; i < this->gaussPts.size(); i++) 
                    //{
                    //GaussPt& gaussPt = this->gaussPts[i];
                if constexpr (nDim>1) {
                    for(GaussPt& gaussPt : gaussPts){

                        typename ParentGeometryElement::XiSized coordAtGauss =  this->NB( this->N(gaussPt.xi)) * this->coordinates; 

                        const double sigY1 = values[0];
                        const double sigY2 = values[2];
                        const double y1    = values[1];
                        const double y2    = values[3];

                        gaussPt.stress(1) = bft::Math::linearInterpolation(coordAtGauss[1], y1, y2, sigY1, sigY2);  // sigma_y
                        gaussPt.stress(0) = values[4]*gaussPt.stress(1);  // sigma_x
                        gaussPt.stress(2) = values[5]*gaussPt.stress(1);
                    }  // sigma_z
                }
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

            switch(loadType){

                case BftUel::Pressure: { 
                                           const double p = load[0];

                                           if (std::abs(p)<bft::Constants::numZeroPos)
                                               return;

                                           bft::FiniteElement::BoundaryElement boundaryEl(this->shape, 
                                                   elementFace,
                                                   nDim, 
                                                   this->coordinates);

                                           VectorXd Pk =  - p * boundaryEl.expandBoundaryToParentVector( boundaryEl.computeNormalLoadVector() );

                                           if(nDim == 2)
                                               Pk *= elementProperties[0]; // thickness

                                           fU+= Pk;

                                           break;
                                       }
                default: {throw std::invalid_argument("Invalid Load Type specified");}
            }
        }

