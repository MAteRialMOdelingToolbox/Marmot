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
 * Magdalena Schreter magdalena.schreter@uibk.ac.at
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
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryElement.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotLowerDimensionalStress.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include <iostream>
#include <memory>
#include <vector>

using namespace Marmot;
using namespace Eigen;

namespace Marmot::Elements {

  template < int nDim, int nNodes >
  class DisplacementFiniteElement : public MarmotElement, public MarmotGeometryElement< nDim, nNodes > {

  public:
    enum SectionType {
      UniaxialStress,
      PlaneStress,
      PlaneStrain,
      Solid,
    };

    static constexpr int sizeLoadVector = nNodes * nDim;
    static constexpr int nCoordinates   = nNodes * nDim;

    using ParentGeometryElement = MarmotGeometryElement< nDim, nNodes >;
    using JacobianSized         = typename ParentGeometryElement::JacobianSized;
    using dNdXiSized            = typename ParentGeometryElement::dNdXiSized;
    using BSized                = typename ParentGeometryElement::BSized;
    using XiSized               = typename ParentGeometryElement::XiSized;
    using RhsSized              = Matrix< double, sizeLoadVector, 1 >;
    using KeSizedMatrix         = Matrix< double, sizeLoadVector, sizeLoadVector >;
    using CSized                = Matrix< double, ParentGeometryElement::VoigtSize, ParentGeometryElement::VoigtSize >;
    using Voigt                 = Matrix< double, ParentGeometryElement::VoigtSize, 1 >;

    Map< const VectorXd > elementProperties;
    const int             elLabel;
    const SectionType     sectionType;

    struct QuadraturePoint {

      const XiSized xi;
      const double  weight;

      double detJ;
      double J0xW;
      BSized B;

      class QPStateVarManager : public MarmotStateVarVectorManager {

        inline const static auto layout = makeLayout( {
          { .name = "stress", .length = 6 },
          { .name = "strain", .length = 6 },
          { .name = "begin of material state", .length = 0 },
        } );

      public:
        mVector6d                     stress;
        mVector6d                     strain;
        Eigen::Map< Eigen::VectorXd > materialStateVars;

        static int getNumberOfRequiredStateVarsQuadraturePointOnly() { return layout.nRequiredStateVars; };

        QPStateVarManager( double* theStateVarVector, int nStateVars )
          : MarmotStateVarVectorManager( theStateVarVector, layout ),
            stress( &find( "stress" ) ),
            strain( &find( "strain" ) ),
            materialStateVars( &find( "begin of material state" ),
                               nStateVars - getNumberOfRequiredStateVarsQuadraturePointOnly() ){};
      };

      std::unique_ptr< QPStateVarManager > managedStateVars;

      std::unique_ptr< MarmotMaterialHypoElastic > material;

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
        : xi( xi ), weight( weight ), detJ( 0.0 ), J0xW( 0.0 ), B( BSized::Zero() ){};
    };

    std::vector< QuadraturePoint > qps;

    DisplacementFiniteElement( int                                         elementID,
                               FiniteElement::Quadrature::IntegrationTypes integrationType,
                               SectionType                                 sectionType );

    int getNumberOfRequiredStateVars();

    std::vector< std::vector< std::string > > getNodeFields();

    std::vector< int > getDofIndicesPermutationPattern();

    int getNNodes() { return nNodes; }

    int getNDofPerElement() { return sizeLoadVector; }

    std::string getElementShape() { return ParentGeometryElement::getElementShape(); }

    void assignStateVars( double* stateVars, int nStateVars );

    void assignProperty( const ElementProperties& marmotElementProperty );

    void assignProperty( const MarmotMaterialSection& marmotElementProperty );

    void initializeYourself( const double* coordinates );

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

    StateView getStateView( const std::string& stateName, int qpNumber )
    {
      const auto& qp = qps[qpNumber];

      if ( stateName == "sdv" ) {
        /* std::cout<< __PRETTY_FUNCTION__ << " on 'sdv' is discouraged and deprecated, please use precise state name"; */
        return { qp.managedStateVars->materialStateVars.data(), static_cast<int>(qp.managedStateVars->materialStateVars.size() ) };
      }

      if ( qp.managedStateVars->contains( stateName ) ) {
        return qp.managedStateVars->getStateView( stateName );
      }
      else {
        return qp.material->getStateView( stateName );
      }
    }
  };

  template < int nDim, int nNodes >
  DisplacementFiniteElement< nDim, nNodes >::DisplacementFiniteElement(
    int                                         elementID,
    FiniteElement::Quadrature::IntegrationTypes integrationType,
    SectionType                                 sectionType )
    : ParentGeometryElement(),
      elementProperties( Map< const VectorXd >( nullptr, 0 ) ),
      elLabel( elementID ),
      sectionType( sectionType )
  {
    for ( const auto& qpInfo : FiniteElement::Quadrature::getGaussPointInfo( this->shape, integrationType ) ) {
      QuadraturePoint qp( qpInfo.xi, qpInfo.weight );
      qps.push_back( std::move( qp ) );
    }
  }

  template < int nDim, int nNodes >
  int DisplacementFiniteElement< nDim, nNodes >::getNumberOfRequiredStateVars()
  {
    return qps[0].getNumberOfRequiredStateVars() * qps.size();
  }

  template < int nDim, int nNodes >
  std::vector< std::vector< std::string > > DisplacementFiniteElement< nDim, nNodes >::getNodeFields()
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
  std::vector< int > DisplacementFiniteElement< nDim, nNodes >::getDofIndicesPermutationPattern()
  {
    static std::vector< int > permutationPattern;
    if ( permutationPattern.empty() )
      for ( int i = 0; i < nNodes * nDim; i++ )
        permutationPattern.push_back( i );

    return permutationPattern;
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::assignStateVars( double* stateVars, int nStateVars )
  {
    const int nQpStateVars = nStateVars / qps.size();

    for ( size_t i = 0; i < qps.size(); i++ ) {
      auto&   qp          = qps[i];
      double* qpStateVars = stateVars + ( i * nQpStateVars );
      qp.assignStateVars( qpStateVars, nQpStateVars );
    }
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::assignProperty( const ElementProperties& elementPropertiesInfo )
  {
    new ( &elementProperties ) Eigen::Map< const Eigen::VectorXd >( elementPropertiesInfo.elementProperties,
                                                                    elementPropertiesInfo.nElementProperties );
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::assignProperty( const MarmotMaterialSection& section )
  {
    for ( auto& qp : qps )
      qp.material = std::unique_ptr< MarmotMaterialHypoElastic >( dynamic_cast< MarmotMaterialHypoElastic* >(
        MarmotLibrary::MarmotMaterialFactory::createMaterial( section.materialCode,
                                                              section.materialProperties,
                                                              section.nMaterialProperties,
                                                              elLabel ) ) );
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::initializeYourself( const double* coordinates )
  {
    ParentGeometryElement::initializeYourself( coordinates );

    for ( QuadraturePoint& qp : qps ) {

      const dNdXiSized    dNdXi = this->dNdXi( qp.xi );
      const JacobianSized J     = this->Jacobian( dNdXi );
      const JacobianSized JInv  = J.inverse();
      const dNdXiSized    dNdX  = this->dNdX( dNdXi, JInv );
      qp.detJ                   = J.determinant();
      qp.B                      = this->B( dNdX );

      if ( sectionType == SectionType::Solid ) {

        qp.J0xW = qp.weight * qp.detJ;
        qp.material->setCharacteristicElementLength( std::cbrt( 8 * qp.detJ ) );
      }

      else if ( sectionType == SectionType::PlaneStrain || sectionType == SectionType::PlaneStress ) {

        const double& thickness = elementProperties[0];
        qp.J0xW                 = qp.weight * qp.detJ * thickness;
        qp.material->setCharacteristicElementLength( std::sqrt( 4 * qp.detJ ) );
      }

      else if ( sectionType == SectionType::UniaxialStress ) {

        const double& crossSection = elementProperties[0];
        qp.J0xW                    = qp.weight * qp.detJ * crossSection;
        qp.material->setCharacteristicElementLength( 2 * qp.detJ );
      }
    }
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::computeYourself( const double* QTotal_,
                                                                   const double* dQ_,
                                                                   double*       Pe_,
                                                                   double*       Ke_,
                                                                   const double* time,
                                                                   double        dT,
                                                                   double&       pNewDT )
  {
    using namespace Marmot;

    Map< const RhsSized > QTotal( QTotal_ );
    Map< const RhsSized > dQ( dQ_ );
    Map< KeSizedMatrix >  Ke( Ke_ );
    Map< RhsSized >       Pe( Pe_ );

    Voigt  S, dE;
    CSized C;

    for ( QuadraturePoint& qp : qps ) {

      const BSized& B = qp.B;

      dE = B * dQ;

      if constexpr ( nDim == 1 ) {
        Vector6d dE6 = ( Vector6d() << dE, 0, 0, 0, 0, 0 ).finished();
        Matrix6d C66;
        qp.material
          ->computeUniaxialStress( qp.managedStateVars->stress.data(), C66.data(), dE6.data(), time, dT, pNewDT );

        C << ContinuumMechanics::UniaxialStress::getUniaxialStressTangent( C66 );
        S( 0 ) = qp.managedStateVars->stress( 0 );
      }

      else if constexpr ( nDim == 2 ) {

        Vector6d dE6 = ContinuumMechanics::VoigtNotation::planeVoigtToVoigt( dE );
        Matrix6d C66;

        if ( sectionType == SectionType::PlaneStress ) {

          qp.material
            ->computePlaneStress( qp.managedStateVars->stress.data(), C66.data(), dE6.data(), time, dT, pNewDT );

          C = ContinuumMechanics::PlaneStress::getPlaneStressTangent( C66 );
        }

        else if ( sectionType == SectionType::PlaneStrain ) {

          qp.material->computeStress( qp.managedStateVars->stress.data(), C66.data(), dE6.data(), time, dT, pNewDT );

          C = ContinuumMechanics::PlaneStrain::getPlaneStrainTangent( C66 );
        }

        S = ContinuumMechanics::VoigtNotation::voigtToPlaneVoigt( qp.managedStateVars->stress );
        qp.managedStateVars->strain += dE6;
      }

      else if constexpr ( nDim == 3 ) {

        if ( sectionType == SectionType::Solid ) {
          qp.material->computeStress( qp.managedStateVars->stress.data(), C.data(), dE.data(), time, dT, pNewDT );
        }

        S = qp.managedStateVars->stress;
        qp.managedStateVars->strain += dE;
      }

      if ( pNewDT < 1.0 )
        return;

      Ke += B.transpose() * C * B * qp.J0xW;
      Pe -= B.transpose() * S * qp.J0xW;
    }
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::setInitialConditions( StateTypes state, const double* values )
  {
    switch ( state ) {
    case MarmotElement::GeostaticStress: {
      if constexpr ( nDim == 3 ) {
        for ( QuadraturePoint& qp : qps ) {

          XiSized coordAtGauss = this->NB( this->N( qp.xi ) ) * this->coordinates;

          const double sigY1 = values[0];
          const double sigY2 = values[2];
          const double y1    = values[1];
          const double y2    = values[3];

          using namespace Math;
          qp.managedStateVars->stress( 1 ) = linearInterpolation( coordAtGauss[1], y1, y2, sigY1, sigY2 ); // sigma_y
          qp.managedStateVars->stress( 0 ) = values[4] * qp.managedStateVars->stress( 1 );                 // sigma_x
          qp.managedStateVars->stress( 2 ) = values[5] * qp.managedStateVars->stress( 1 );
        }
      }
      break;
    }
    case MarmotElement::MarmotMaterialStateVars: {
      throw std::invalid_argument( "Please use initializeStateVars directly on material" );
    }
    default: throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid initial condition" );
    }
    /* } */
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::computeDistributedLoad( MarmotElement::DistributedLoadTypes loadType,
                                                                          double*                             P,
                                                                          double*                             K,
                                                                          const int     elementFace,
                                                                          const double* load,
                                                                          const double* QTotal,
                                                                          const double* time,
                                                                          double        dT )
  {
    Map< RhsSized > fU( P );

    switch ( loadType ) {

    case MarmotElement::Pressure: {
      const double p = load[0];

      FiniteElement::BoundaryElement boundaryEl( this->shape, elementFace, nDim, this->coordinates );

      VectorXd Pk = -p * boundaryEl.computeNormalLoadVector();

      if ( nDim == 2 )
        Pk *= elementProperties[0]; // thickness

      boundaryEl.assembleIntoParentVector( Pk, fU );

      break;
    }
    default: {
      throw std::invalid_argument( "Invalid Load Type specified" );
    }
    }
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::computeBodyForce( double*       P_,
                                                                    double*       K,
                                                                    const double* load,
                                                                    const double* QTotal,
                                                                    const double* time,
                                                                    double        dT )
  {
    Map< RhsSized >                              Pe( P_ );
    const Map< const Matrix< double, nDim, 1 > > f( load );

    for ( const auto& qp : qps )
      Pe += this->NB( this->N( qp.xi ) ).transpose() * f * qp.J0xW;
  }
} // namespace Marmot::Elements
