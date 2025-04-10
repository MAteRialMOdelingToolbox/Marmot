/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck,
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
 *
 * Matthias Neuner matthias.neuner@uibk.ac.at
 * Alexander Dummer alexander.dummer@uibk.ac.at
 *
 * This file is part of the MAteRialMOdellingToolbox (marmot).
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of marmot.
 * ---------------------------------------------------------------------
 */
#pragma once

#include "Marmot/Marmot.h"
#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryElement.h"
#include "Marmot/MarmotGeostaticStress.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotTypedefs.h"
#include <iostream>
#include <memory>
#include <type_traits>
#include <vector>

namespace Marmot::Elements {

  template < int nDim, int nNodes >
  class DisplacementFiniteStrainULElement : public MarmotElement, public MarmotGeometryElement< nDim, nNodes > {

  public:
    enum SectionType {
      PlaneStress,
      PlaneStrain,
      Solid,
    };

    static constexpr int nDofPerNodeU = nDim; // Displacement   field U
    static constexpr int nCoordinates = nNodes * nDim;

    // block sizes

    static constexpr int bsU = nNodes * nDofPerNodeU;

    static constexpr int sizeLoadVector = bsU;

    static constexpr int idxU = 0;

    using ParentGeometryElement = MarmotGeometryElement< nDim, nNodes >;
    using Material              = MarmotMaterialFiniteStrain;

    using JacobianSized = typename ParentGeometryElement::JacobianSized;
    using NSized        = typename ParentGeometryElement::NSized;
    using dNdXiSized    = typename ParentGeometryElement::dNdXiSized;
    using BSized        = typename ParentGeometryElement::BSized;
    using XiSized       = typename ParentGeometryElement::XiSized;
    using RhsSized      = Eigen::Matrix< double, sizeLoadVector, 1 >;
    using KSizedMatrix  = Eigen::Matrix< double, sizeLoadVector, sizeLoadVector >;
    using USizedVector  = Eigen::Matrix< double, bsU, 1 >;

    using ForceSized = Eigen::Matrix< double, nDim, 1 >;

    Eigen::Map< const Eigen::VectorXd > elementProperties;
    const int                           elLabel;
    const SectionType                   sectionType;

    bool hasEigenDeformation;

    struct QuadraturePoint {

      const XiSized xi;
      const double  weight;

      dNdXiSized dNdX;
      double     J0xW;

      class QPStateVarManager : public MarmotStateVarVectorManager {

        inline const static auto layout = makeLayout( {
          { .name = "stress", .length = 9 },
          { .name = "F0 XX", .length = 1 },
          { .name = "F0 YY", .length = 1 },
          { .name = "F0 ZZ", .length = 1 },
          { .name = "begin of material state", .length = 0 },
        } );

      public:
        Eigen::Map< Marmot::Vector9d > stress;
        double&                        F0_XX;
        double&                        F0_YY;
        double&                        F0_ZZ;
        Eigen::Map< Eigen::VectorXd >  materialStateVars;

        static int getNumberOfRequiredStateVarsQuadraturePointOnly() { return layout.nRequiredStateVars; };

        QPStateVarManager( double* theStateVarVector, int nStateVars )
          : MarmotStateVarVectorManager( theStateVarVector, layout ),
            stress( &find( "stress" ) ),
            F0_XX( find( "F0 XX" ) ),
            F0_YY( find( "F0 YY" ) ),
            F0_ZZ( find( "F0 ZZ" ) ),
            materialStateVars( &find( "begin of material state" ),
                               nStateVars - getNumberOfRequiredStateVarsQuadraturePointOnly() ){};
      };

      std::unique_ptr< QPStateVarManager > managedStateVars;

      std::unique_ptr< Material > material;

      int getNumberOfRequiredStateVarsQuadraturePointOnly()
      {
        return QPStateVarManager::getNumberOfRequiredStateVarsQuadraturePointOnly();
      };

      int getNumberOfRequiredStateVars()
      {
        return getNumberOfRequiredStateVarsQuadraturePointOnly() + material->getNumberOfRequiredStateVars();
      };

      void assignStateVars( double* stateVars, int nStateVars )
      {
        managedStateVars = std::make_unique< QPStateVarManager >( stateVars, nStateVars );
        material->assignStateVars( managedStateVars->materialStateVars.data(),
                                   managedStateVars->materialStateVars.size() );
      }

      QuadraturePoint( XiSized xi, double weight )
        : xi( xi ), weight( weight ), dNdX( dNdXiSized::Zero() ), J0xW( 0.0 ){};
    };

    std::vector< QuadraturePoint > qps;

    DisplacementFiniteStrainULElement( int                                                 elementID,
                                       Marmot::FiniteElement::Quadrature::IntegrationTypes integrationType,
                                       SectionType                                         sectionType );

    int getNumberOfRequiredStateVars();

    std::vector< std::vector< std::string > > getNodeFields();

    std::vector< int > getDofIndicesPermutationPattern();

    int getNNodes() { return nNodes; }

    int getNSpatialDimensions() { return nDim; }

    int getNDofPerElement() { return sizeLoadVector; }

    std::string getElementShape() { return ParentGeometryElement::getElementShape(); }

    void assignStateVars( double* managedStateVars, int nStateVars );

    void assignProperty( const ElementProperties& MarmotElementProperty );

    void assignProperty( const MarmotMaterialSection& MarmotElementProperty );

    void assignNodeCoordinates( const double* coordinates );

    void initializeYourself();

    void setInitialConditions( StateTypes state, const double* values );

    void computeDistributedLoad( MarmotElement::DistributedLoadTypes loadType,
                                 double*                             P,
                                 double*                             K,
                                 const int                           elementFace,
                                 const double*                       load,
                                 const double*                       QTotal,
                                 const double*                       time,
                                 double                              dT );

    void computeBodyForce( double*       P,
                           double*       K,
                           const double* load,
                           const double* QTotal,
                           const double* time,
                           double        dT );

    void computeYourself( const double* QTotal,
                          const double* dQ,
                          double*       Pe,
                          double*       Ke,
                          const double* time,
                          double        dT,
                          double&       pNewdT );

    StateView getStateView( const std::string& stateName, int qpNumber );

    std::vector< double > getCoordinatesAtCenter();

    std::vector< std::vector< double > > getCoordinatesAtQuadraturePoints();

    int getNumberOfQuadraturePoints();
  };

  template < int nDim, int nNodes >
  StateView DisplacementFiniteStrainULElement< nDim, nNodes >::getStateView( const std::string& stateName,
                                                                             int                qpNumber )
  {
    const auto& qp = qps[qpNumber];
    if ( qp.managedStateVars->contains( stateName ) ) {
      return qp.managedStateVars->getStateView( stateName );
    }
    else {
      return qp.material->getStateView( stateName );
    }
  }

  template < int nDim, int nNodes >
  DisplacementFiniteStrainULElement< nDim, nNodes >::DisplacementFiniteStrainULElement(
    int                                                 elementID,
    Marmot::FiniteElement::Quadrature::IntegrationTypes integrationType,
    SectionType                                         sectionType )
    : ParentGeometryElement(),
      elementProperties( Eigen::Map< const Eigen::VectorXd >( nullptr, 0 ) ),
      elLabel( elementID ),
      sectionType( sectionType ),
      hasEigenDeformation( false )
  {
    for ( const auto& qpInfo : Marmot::FiniteElement::Quadrature::getGaussPointInfo( this->shape, integrationType ) ) {
      QuadraturePoint qp( qpInfo.xi, qpInfo.weight );
      qps.push_back( std::move( qp ) );
    }
  }

  template < int nDim, int nNodes >
  int DisplacementFiniteStrainULElement< nDim, nNodes >::getNumberOfRequiredStateVars()
  {
    return qps[0].getNumberOfRequiredStateVars() * qps.size();
  }

  template < int nDim, int nNodes >
  std::vector< std::vector< std::string > > DisplacementFiniteStrainULElement< nDim, nNodes >::getNodeFields()
  {
    using namespace std;

    static vector< vector< string > > nodeFields;
    if ( nodeFields.empty() )
      for ( int i = 0; i < nNodes; i++ ) {
        nodeFields.push_back( vector< string >() );
        nodeFields[i].push_back( "displacement" );
      }

    return nodeFields;
  }

  template < int nDim, int nNodes >
  std::vector< int > DisplacementFiniteStrainULElement< nDim, nNodes >::getDofIndicesPermutationPattern()
  {
    static std::vector< int > permutationPattern;
    if ( permutationPattern.empty() ) {
      for ( int i = 0; i < nNodes; i++ )
        for ( int j = 0; j < nDim; j++ )
          permutationPattern.push_back( i * nDim + j );
    }

    return permutationPattern;
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteStrainULElement< nDim, nNodes >::assignStateVars( double* managedStateVars, int nStateVars )
  {
    const int nQpStateVars = nStateVars / qps.size();

    for ( size_t i = 0; i < qps.size(); i++ ) {
      auto&   qp          = qps[i];
      double* qpStateVars = &managedStateVars[i * nQpStateVars];
      qp.assignStateVars( qpStateVars, nQpStateVars );
    }
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteStrainULElement< nDim, nNodes >::assignProperty(
    const ElementProperties& elementPropertiesInfo )
  {
    new ( &elementProperties ) Eigen::Map< const Eigen::VectorXd >( elementPropertiesInfo.elementProperties,
                                                                    elementPropertiesInfo.nElementProperties );
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteStrainULElement< nDim, nNodes >::assignProperty( const MarmotMaterialSection& section )
  {

    for ( auto& qp : qps ) {
      qp.material = std::unique_ptr< Material >(
        dynamic_cast< Material* >( MarmotLibrary::MarmotMaterialFactory::createMaterial( section.materialCode,
                                                                                         section.materialProperties,
                                                                                         section.nMaterialProperties,
                                                                                         elLabel ) ) );
    }
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteStrainULElement< nDim, nNodes >::assignNodeCoordinates( const double* coordinates )
  {
    ParentGeometryElement::assignNodeCoordinates( coordinates );
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteStrainULElement< nDim, nNodes >::initializeYourself()
  {
    for ( QuadraturePoint& qp : qps ) {

      const dNdXiSized dNdXi_ = this->dNdXi( qp.xi );

      const JacobianSized J    = this->Jacobian( dNdXi_ );
      const JacobianSized JInv = J.inverse();
      const double        detJ = J.determinant();

      qp.dNdX = this->dNdX( dNdXi_, JInv );

      if constexpr ( nDim == 3 ) {
        qp.J0xW = qp.weight * detJ;
      }
      if constexpr ( nDim == 2 ) {
        const double& thickness = elementProperties[0];
        qp.J0xW                 = qp.weight * detJ * thickness;
      }
      if constexpr ( nDim == 1 ) {
        const double& crossSection = elementProperties[0];
        qp.J0xW                    = qp.weight * detJ * crossSection;
      }
    }
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteStrainULElement< nDim, nNodes >::computeYourself( const double* qTotal,
                                                                           const double* dQ,
                                                                           double*       rightHandSide,
                                                                           double*       stiffnessMatrix,
                                                                           const double* time,
                                                                           double        dT,
                                                                           double&       pNewDT )
  {
    using namespace Fastor;

    const static Tensor< double, nDim, nDim > I(
      ( Eigen::Matrix< double, nDim, nDim >() << Eigen::Matrix< double, nDim, nDim >::Identity() ).finished().data() );

    // in  ...
    const auto qU_np = TensorMap< const double, nNodes, nDim >( qTotal );

    // ... and out: residuals and stiffness
    TensorMap< double, nNodes, nDim > r_U( rightHandSide );

    // temporary stiffness matrices, which are assembled into the large one at the end of the method
    Tensor< double, nDim, nNodes, nDim, nNodes > k_UU( 0.0 );

    Eigen::Map< Eigen::VectorXd > rhs( rightHandSide, sizeLoadVector );

    for ( auto& qp : qps ) {

      using namespace Marmot::FastorIndices;

      const auto& dNdX_ = qp.dNdX;

      const auto dNdX = Tensor< double, nDim, nNodes >( dNdX_.data(), ColumnMajor );

      const auto F_np = evaluate( einsum< Ai, jA >( qU_np, dNdX ) + I );

      const Material::Deformation< nDim > deformation = { F_np };

      const Material::TimeIncrement timeIncrement{ time[0], dT };

      Material::ConstitutiveResponse< nDim > response;
      Material::AlgorithmicModuli< nDim >    tangents;
      try {
        if constexpr ( nDim == 2 ) {

          if ( sectionType == SectionType::PlaneStrain ) {

            using namespace Marmot;

            Material::ConstitutiveResponse< 3 >
              response3D{ FastorStandardTensors::Tensor33d( qp.managedStateVars->stress.data(), Fastor::ColumnMajor ),
                          -1.0,
                          -1.0 };

            Material::AlgorithmicModuli< 3 > algorithmicModuli3D;

            Material::Deformation< 3 > deformation3D{
              expandTo3D( deformation.F ),
            };

            deformation3D.F( 2, 2 ) = 1.0;

            if ( hasEigenDeformation )
              qp.material->computePlaneStrain( response3D,
                                               algorithmicModuli3D,
                                               deformation3D,
                                               timeIncrement,
                                               { qp.managedStateVars->F0_XX,
                                                 qp.managedStateVars->F0_YY,
                                                 qp.managedStateVars->F0_ZZ } );
            else
              qp.material->computePlaneStrain( response3D, algorithmicModuli3D, deformation3D, timeIncrement );

            response = { reduceTo2D< U, U >( response3D.tau ), response3D.rho, response3D.elasticEnergyDensity };

            tangents = {
              reduceTo2D< U, U, U, U >( algorithmicModuli3D.dTau_dF ),
            };

            qp.managedStateVars->stress = Marmot::mapEigenToFastor( response3D.tau ).reshaped();
          }
        }
        else {
          response = { Marmot::FastorStandardTensors::Tensor33d( qp.managedStateVars->stress.data(), ColumnMajor ),
                       -1.0,
                       -1.0 };

          qp.material->computeStress( response, tangents, deformation, timeIncrement );

          // implicit conversion to col major
          qp.managedStateVars->stress = Marmot::mapEigenToFastor( response.tau ).reshaped();
        }
      }
      catch ( const std::runtime_error& ) {
        pNewDT = 0.25;
        return;
      }
      const auto dNdx = evaluate( einsum< ji, jA >( inv( F_np ), dNdX ) );

      const double& J0xW = qp.J0xW;

      const auto& tau = response.tau;

      const auto& t = tangents;

      // aux stiffness tensors
      const auto dTau_dqU = evaluate( +einsum< ijkl, lB >( t.dTau_dF, dNdX ) );

      // r[ node, dim ] (swap to abuse directly colmajor layout)
      // directly operate via TensorMap
      r_U -= ( +einsum< iA, ij >( dNdx, tau ) ) * J0xW;

      // K [dim, node, dim, node ]
      k_UU += ( +einsum< iA, ijkB, to_jAkB >( dNdx, dTau_dqU ) ) * J0xW;

      // geometric contribution
      k_UU += ( -einsum< kA, ij, iB, to_jAkB >( dNdx, tau, dNdx ) ) * J0xW;
    }
    // copy back to the subblocks using mighty Eigen block access,
    // note the layout swap rowmajor -> colmajor
    // Fastor offers similar functionality, but performancy is (slightly) inferior

    using namespace Eigen;

    Map< KSizedMatrix > K( stiffnessMatrix );

    K.template block< bsU, bsU >( idxU, idxU ) += Map< Matrix< double, bsU, bsU > >( torowmajor( k_UU ).data() );
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteStrainULElement< nDim, nNodes >::computeDistributedLoad(
    MarmotElement::DistributedLoadTypes loadType,
    double*                             rightHandSide,
    double*                             stiffnessMatrix,
    const int                           elementFace,
    const double*                       load,
    const double*                       QTotal_,
    const double*                       time,
    double                              dT )
  {

    Eigen::Map< USizedVector > r_U( rightHandSide );

    switch ( loadType ) {

    case MarmotElement::Pressure: {
      const double                                                          p = load[0];
      const Eigen::Map< const RhsSized >                                    QTotal( QTotal_ );
      const Eigen::Ref< const USizedVector >                                qU( QTotal.head( bsU ) );
      Eigen::Map< Eigen::Matrix< double, sizeLoadVector, sizeLoadVector > > K( stiffnessMatrix );
      Eigen::Ref< Eigen::Matrix< double, bsU, bsU > >                       kUU( K.topLeftCorner( bsU, bsU ) );

      const USizedVector             coordinates_np = this->coordinates + qU;
      FiniteElement::BoundaryElement boundaryEl( this->shape, elementFace, nDim, coordinates_np );

      Eigen::VectorXd Pb = -p * boundaryEl.computeSurfaceNormalVectorialLoadVector();
      Eigen::MatrixXd Kb = -p * boundaryEl.computeDSurfaceNormalVectorialLoadVector_dCoordinates();

      if ( nDim == 2 ) {
        Pb *= elementProperties[0]; // thickness
        Kb *= elementProperties[0];
      }

      boundaryEl.assembleIntoParentVectorial( Pb, r_U );
      boundaryEl.assembleIntoParentStiffnessVectorial( Kb, K );

      break;
    }
    case MarmotElement::SurfaceTraction: {

      FiniteElement::BoundaryElement boundaryEl( this->shape, elementFace, nDim, this->coordinates );

      const XiSized tractionVector( load );

      auto Pk = boundaryEl.computeVectorialLoadVector( tractionVector );
      if ( nDim == 2 )
        Pk *= elementProperties[0]; // thickness
      boundaryEl.assembleIntoParentVectorial( Pk, r_U );

      break;
    }
    default: {
      throw std::invalid_argument( "Invalid Load Type specified" );
    }
    }
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteStrainULElement< nDim, nNodes >::setInitialConditions(
    StateTypes    state,
    const double* initialConditionDefinition )
  {
    if constexpr ( nDim > 1 ) {
      switch ( state ) {

      case MarmotElement::MarmotMaterialInitialization: {
        for ( QuadraturePoint& qp : qps ) {

          qp.managedStateVars->F0_XX = 1.0;
          qp.managedStateVars->F0_YY = 1.0;
          qp.managedStateVars->F0_ZZ = 1.0;

          qp.material->initializeYourself();
        }
        break;
      }

      case MarmotElement::GeostaticStress: {

        for ( QuadraturePoint& qp : qps ) {

          XiSized coordAtGauss = this->NB( this->N( qp.xi ) ) * this->coordinates;

          const auto geostaticNormalStressComponents = Marmot::GeostaticStress::
            getGeostaticStressFromLinearDistribution( initialConditionDefinition, coordAtGauss[1] );

          const auto [F0_XX,
                      F0_YY,
                      F0_ZZ] = qp.material->findEigenDeformationForEigenStress( { qp.managedStateVars->F0_XX,
                                                                                  qp.managedStateVars->F0_YY,
                                                                                  qp.managedStateVars->F0_ZZ },
                                                                                geostaticNormalStressComponents );

          qp.managedStateVars->F0_XX = F0_XX;
          qp.managedStateVars->F0_YY = F0_YY;
          qp.managedStateVars->F0_ZZ = F0_ZZ;

          hasEigenDeformation = true;
        }
        break;
      }
      default: throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid initial condition" );
      }
    }
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteStrainULElement< nDim, nNodes >::computeBodyForce( double*       rightHandSide,
                                                                            double*       stiffnessMatrix,
                                                                            const double* load,

                                                                            const double* qTotal,
                                                                            const double* time,
                                                                            double        dT )
  {
    Eigen::Map< RhsSized >                                     r( rightHandSide );
    Eigen::Ref< USizedVector >                                 r_U( r.head( bsU ) );
    const Eigen::Map< const Eigen::Matrix< double, nDim, 1 > > f( load );

    for ( const auto& qp : qps )
      r_U += this->NB( this->N( qp.xi ) ).transpose() * f * qp.J0xW;
  }

  template < int nDim, int nNodes >
  std::vector< double > DisplacementFiniteStrainULElement< nDim, nNodes >::getCoordinatesAtCenter()
  {
    std::vector< double > coords( nDim );

    Eigen::Map< XiSized > coordsMap( &coords[0] );
    const auto            centerXi = XiSized::Zero();
    coordsMap                      = this->NB( this->N( centerXi ) ) * this->coordinates;
    return coords;
  }

  template < int nDim, int nNodes >
  std::vector< std::vector< double > > DisplacementFiniteStrainULElement< nDim,
                                                                          nNodes >::getCoordinatesAtQuadraturePoints()
  {
    std::vector< std::vector< double > > listedCoords;

    std::vector< double > coords( nDim );
    Eigen::Map< XiSized > coordsMap( &coords[0] );

    for ( const auto& qp : qps ) {
      coordsMap = this->NB( this->N( qp.xi ) ) * this->coordinates;
      listedCoords.push_back( coords );
    }

    return listedCoords;
  }

  template < int nDim, int nNodes >
  int DisplacementFiniteStrainULElement< nDim, nNodes >::getNumberOfQuadraturePoints()
  {
    return qps.size();
  }

  template < int nNodes >
  class AxiSymmetricDisplacementFiniteStrainULElement : public DisplacementFiniteStrainULElement< 2, nNodes > {

    using DisplacementFiniteStrainULElement< 2, nNodes >::DisplacementFiniteStrainULElement;

    void computeYourself( const double* QTotal,
                          const double* dQ,
                          double*       Pe,
                          double*       Ke,
                          const double* time,
                          double        dT,
                          double&       pNewdT );
  };

  template < int nNodes >
  void AxiSymmetricDisplacementFiniteStrainULElement< nNodes >::computeYourself( const double* qTotal,
                                                                                 const double* dQ,
                                                                                 double*       rightHandSide,
                                                                                 double*       stiffnessMatrix,
                                                                                 const double* time,
                                                                                 double        dT,
                                                                                 double&       pNewDT )
  {
    constexpr int nDim = 2;

    using Parent              = DisplacementFiniteStrainULElement< 2, nNodes >;
    const auto idxU           = Parent::idxU;
    const auto sizeLoadVector = DisplacementFiniteStrainULElement< 2, nNodes >::sizeLoadVector;
    using Material            = MarmotMaterialFiniteStrain;

    using namespace Fastor;

    const static Tensor< double, nDim, nDim > I(
      ( Eigen::Matrix< double, nDim, nDim >() << Eigen::Matrix< double, nDim, nDim >::Identity() ).finished().data() );

    // in  ...
    const auto qU_np = TensorMap< const double, nNodes, nDim >( qTotal );

    // ... and out: residuals and stiffness
    TensorMap< double, nNodes, nDim > r_U( rightHandSide );

    // temporary stiffness matrices, which are assembled into the large one at the end of the method
    Tensor< double, nDim, nNodes, nDim, nNodes > k_UU( 0.0 );

    Eigen::Map< Eigen::VectorXd > rhs( rightHandSide, sizeLoadVector );

    for ( auto& qp : Parent::qps ) {

      using namespace Marmot::FastorIndices;

      auto        N_    = this->N( qp.xi );
      const auto& dNdX_ = qp.dNdX;

      Eigen::Vector2d coords = this->NB( N_ ) * this->coordinates;
      const double    r      = coords[0];

      const auto N    = Tensor< double, nNodes >( N_.data() );
      const auto dNdX = Tensor< double, nDim, nNodes >( dNdX_.data(), ColumnMajor );

      const auto u_np = evaluate( einsum< A, Ai >( N, qU_np ) );

      const auto F_np = evaluate( einsum< Ai, jA >( qU_np, dNdX ) + I );

      const Material::Deformation< nDim > deformation = {
        F_np,
      };

      const Material ::TimeIncrement timeIncrement{ time[0], dT };

      Material::ConstitutiveResponse< nDim > response;
      Material::AlgorithmicModuli< nDim >    tangents;

      using namespace Marmot;
      Material::ConstitutiveResponse< 3 >
        response3D{ FastorStandardTensors::Tensor33d( qp.managedStateVars->stress.data(), Fastor::ColumnMajor ),
                    -1.0,
                    -1.0 };

      Material::AlgorithmicModuli< 3 > algorithmicModuli3D;

      Material::Deformation< 3 > deformation3D{ expandTo3D( deformation.F ) };

      deformation3D.F( 2, 2 ) = 1 + u_np[0] / r;

      try {
        qp.material->computePlaneStrain( response3D, algorithmicModuli3D, deformation3D, timeIncrement );
      }
      catch ( const std::runtime_error& ) {
        pNewDT = 0.25;
        return;
      }
      response = { reduceTo2D< U, U >( response3D.tau ), response3D.rho, response3D.elasticEnergyDensity };

      tangents = {
        reduceTo2D< U, U, U, U >( algorithmicModuli3D.dTau_dF ),
      };

      qp.managedStateVars->stress = Marmot::mapEigenToFastor( response3D.tau ).reshaped();

      const auto dNdx = evaluate( einsum< ji, jA >( inv( F_np ), dNdX ) );

      const double J0xWxRx2Pi = qp.J0xW * 2 * Constants::Pi * r;

      const auto& tau = response.tau;

      const auto& t = tangents;

      // aux stiffness tensors
      auto dTau_dqU = evaluate( +einsum< ijkl, lB >( t.dTau_dF, dNdX ) );

      Fastor::Tensor< double, 2, 2 > dTau33_dF_2D( 0.0 );
      for ( int i = 0; i < 2; i++ )
        for ( int j = 0; j < 2; j++ ) {
          dTau33_dF_2D( i, j ) = algorithmicModuli3D.dTau_dF( 2, 2, i, j );
        }
      auto dTau33_dqU = evaluate( +einsum< kl, lB >( dTau33_dF_2D, dNdX ) );

      // additional terms due to axisymmetry:

      for ( int B = 0; B < nNodes; B++ ) {
        for ( int i = 0; i < 2; i++ )
          for ( int j = 0; j < 2; j++ ) {
            dTau_dqU( i, j, 0, B ) += algorithmicModuli3D.dTau_dF( i, j, 2, 2 ) * N( B ) / r;
          }
        dTau33_dqU( 0, B ) += algorithmicModuli3D.dTau_dF( 2, 2, 2, 2 ) * N( B ) / r;
      }

      // r[ node, dim ] (swap to abuse directly colmajor layout)
      // directly operate via TensorMap
      r_U -= ( +einsum< iA, ij >( dNdx, tau ) ) * J0xWxRx2Pi;

      const double F33          = 1 + u_np[0] / r;
      const double invF33       = 1. / F33;
      const double dInvF33_dF33 = -1. / ( F33 * F33 );

      for ( int A = 0; A < nNodes; A++ ) {
        r_U( A, 0 ) -= ( +N( A ) * invF33 * response3D.tau( 2, 2 ) / r ) * J0xWxRx2Pi;

        for ( int B = 0; B < nNodes; B++ ) {
          for ( int j = 0; j < 2; j++ ) {
            k_UU( 0, A, j, B ) += ( +N( A ) * invF33 * dTau33_dqU( j, B ) / r ) * J0xWxRx2Pi;
          }
          const double dF33_dN_qU_0 = ( N( B ) / r );
          k_UU( 0, A, 0, B ) += ( +N( A ) * dInvF33_dF33 * dF33_dN_qU_0 * response3D.tau( 2, 2 ) / r ) * J0xWxRx2Pi;
        }
      }

      // K [dim, node, dim, node ]
      k_UU += ( +einsum< iA, ijkB, to_jAkB >( dNdx, dTau_dqU ) ) * J0xWxRx2Pi;
    }
    // copy back to the subblocks using mighty Eigen block access,
    // note the layout swap rowmajor -> colmajor
    // Fastor offers similar functionality, but performancy is (slightly) inferior

    using namespace Eigen;

    Map< typename Parent::KSizedMatrix > K( stiffnessMatrix );

    K.template block< Parent::bsU, Parent::bsU >( idxU, idxU ) += Map< Matrix< double, Parent::bsU, Parent::bsU > >(
      torowmajor( k_UU ).data() );
  }

} // namespace Marmot::Elements
