/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck.
 *
 * Thomas Mader thomas.mader@boku.ac.at
 *
 * This file is part of the MAteRialMOdellingToolbox (marmot).
 * LGPL v2.1+, see LICENSE.md at the top level directory of marmot.
 * ---------------------------------------------------------------------
 *
 * Gradient-enhanced (implicit-gradient / nonlocal) finite-strain
 * displacement element.  It is the non-micropolar sibling of
 * GradientEnhancedMicropolarULFiniteElement: the micro-rotation field W is
 * removed, leaving the displacement field U (nDim DOFs/node) coupled to a
 * single scalar nonlocal-damage field N (1 DOF/node) that regularises
 * softening.  Updated-Lagrange kinematics as in DisplacementFiniteStrainULElement.
 * Consumes a MarmotMaterialGradientEnhancedFiniteStrain.
 */
#pragma once

#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryElement.h"
#include "Marmot/MarmotGeostaticStress.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrain.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrainFactory.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotTypedefs.h"
#include <memory>
#include <vector>

namespace Marmot::Elements {

  template < int nDim, int nNodes >
  class GradientEnhancedFiniteStrainDisplacementElement : public MarmotElement,
                                                          public MarmotGeometryElement< nDim, nNodes > {

  public:
    enum SectionType {
      PlaneStress,
      PlaneStrain,
      Solid,
    };

    static constexpr int nDofPerNodeU = nDim; // Displacement    field U
    static constexpr int nDofPerNodeN = 1;    // Nonlocal damage  field N

    static constexpr int nCoordinates = nNodes * nDim;

    static constexpr int bsU = nNodes * nDofPerNodeU;
    static constexpr int bsN = nNodes * nDofPerNodeN;

    static constexpr int sizeLoadVector = bsU + bsN;

    static constexpr int idxU = 0;
    static constexpr int idxN = idxU + bsU;

    using ParentGeometryElement = MarmotGeometryElement< nDim, nNodes >;
    using Material              = MarmotMaterialGradientEnhancedFiniteStrain;

    using JacobianSized = typename ParentGeometryElement::JacobianSized;
    using NSized        = typename ParentGeometryElement::NSized;
    using dNdXiSized    = typename ParentGeometryElement::dNdXiSized;
    using XiSized       = typename ParentGeometryElement::XiSized;
    using RhsSized      = Eigen::Matrix< double, sizeLoadVector, 1 >;
    using KSizedMatrix  = Eigen::Matrix< double, sizeLoadVector, sizeLoadVector >;
    using USizedVector  = Eigen::Matrix< double, bsU, 1 >;

    Eigen::Map< const Eigen::VectorXd > elementProperties;
    const int                           elLabel;
    const SectionType                   sectionType;
    bool                                hasEigenDeformation;

    struct QuadraturePoint {

      const XiSized xi;
      const double  weight;

      dNdXiSized dNdX;
      double     detJ;
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
      std::unique_ptr< Material >          material;

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
      }

      QuadraturePoint( XiSized xi, double weight )
        : xi( xi ), weight( weight ), dNdX( dNdXiSized::Zero() ), detJ( 0.0 ), J0xW( 0.0 ){};
    };

    std::vector< QuadraturePoint > qps;

    GradientEnhancedFiniteStrainDisplacementElement( int                                                 elementID,
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
                                 double                              time,
                                 double                              dT );

    void computeBodyForce( double* P, double* K, const double* load, const double* QTotal, double time, double dT );

    void computeKernels( const double* QTotal, const double* dQ, double* Pe, double* Ke, double time, double dT );

    void computeKernelsExplicit( const double* QTotal, const double* dQ, double* Pe, double time, double dT );

    StateView getStateView( const std::string& stateName, int qpNumber );

    std::vector< double > getCoordinatesAtCenter();

    std::vector< std::vector< double > > getCoordinatesAtQuadraturePoints();

    int getNumberOfQuadraturePoints();
  };

  template < int nDim, int nNodes >
  StateView GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::getStateView( const std::string& stateName,
                                                                                           int qpNumber )
  {
    const auto& qp = qps[qpNumber];
    if ( qp.managedStateVars->contains( stateName ) ) {
      return qp.managedStateVars->getStateView( stateName );
    }
    else {
      return qp.material->getStateView( stateName, qp.managedStateVars->materialStateVars.data() );
    }
  }

  template < int nDim, int nNodes >
  GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::GradientEnhancedFiniteStrainDisplacementElement(
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
  int GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::getNumberOfRequiredStateVars()
  {
    return qps[0].getNumberOfRequiredStateVars() * qps.size();
  }

  template < int nDim, int nNodes >
  std::vector< std::vector< std::string > > GradientEnhancedFiniteStrainDisplacementElement< nDim,
                                                                                             nNodes >::getNodeFields()
  {
    using namespace std;

    static vector< vector< string > > nodeFields;
    if ( nodeFields.empty() )
      for ( int i = 0; i < nNodes; i++ ) {
        nodeFields.push_back( vector< string >() );
        nodeFields[i].push_back( "displacement" );
        nodeFields[i].push_back( "nonlocal damage" );
      }

    return nodeFields;
  }

  template < int nDim, int nNodes >
  std::vector< int > GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::getDofIndicesPermutationPattern()
  {
    static std::vector< int > permutationPattern;
    if ( permutationPattern.empty() ) {
      for ( int i = 0; i < nNodes; i++ )
        for ( int j = 0; j < nDim; j++ )
          permutationPattern.push_back( i * ( nDim + 1 ) + j );
      for ( int i = 0; i < nNodes; i++ )
        permutationPattern.push_back( i * ( nDim + 1 ) + nDim );
    }

    return permutationPattern;
  }

  template < int nDim, int nNodes >
  void GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::assignStateVars( double* managedStateVars,
                                                                                         int     nStateVars )
  {
    const int nQpStateVars = nStateVars / qps.size();

    for ( size_t i = 0; i < qps.size(); i++ ) {
      auto&   qp          = qps[i];
      double* qpStateVars = &managedStateVars[i * nQpStateVars];
      qp.assignStateVars( qpStateVars, nQpStateVars );
    }
  }

  template < int nDim, int nNodes >
  void GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::assignProperty(
    const ElementProperties& elementPropertiesInfo )
  {
    new ( &elementProperties ) Eigen::Map< const Eigen::VectorXd >( elementPropertiesInfo.elementProperties,
                                                                    elementPropertiesInfo.nElementProperties );
  }

  template < int nDim, int nNodes >
  void GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::assignProperty(
    const MarmotMaterialSection& section )
  {
    for ( auto& qp : qps ) {
      qp.material = std::unique_ptr< Material >(
        MarmotLibrary::MarmotMaterialGradientEnhancedFiniteStrainFactory::createMaterial( section.materialName,
                                                                                          section.materialProperties,
                                                                                          section.nMaterialProperties,
                                                                                          elLabel ) );
    }
  }

  template < int nDim, int nNodes >
  void GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::assignNodeCoordinates(
    const double* coordinates )
  {
    ParentGeometryElement::assignNodeCoordinates( coordinates );
  }

  template < int nDim, int nNodes >
  void GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::initializeYourself()
  {
    for ( QuadraturePoint& qp : qps ) {

      const dNdXiSized dNdXi_ = this->dNdXi( qp.xi );

      const JacobianSized J    = this->Jacobian( dNdXi_ );
      const JacobianSized JInv = J.inverse();
      const double        detJ = J.determinant();
      qp.detJ                  = detJ;

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
  void GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::computeKernels( const double* qTotal,
                                                                                        const double* dQ,
                                                                                        double*       rightHandSide,
                                                                                        double*       stiffnessMatrix,
                                                                                        double        time,
                                                                                        double        dT )
  {
    using namespace Fastor;

    static const Tensor< double, nDim, nDim > I(
      ( Eigen::Matrix< double, nDim, nDim >() << Eigen::Matrix< double, nDim, nDim >::Identity() ).finished().data() );

    // in  ...
    const auto qU_np = TensorMap< const double, nNodes, nDim >( qTotal );
    const auto qN_np = TensorMap< const double, nNodes >( qTotal + idxN );

    // ... and out: residuals
    TensorMap< double, nNodes, nDim > r_U( rightHandSide );
    TensorMap< double, nNodes >       r_N( rightHandSide + idxN );

    // temporary stiffness sub-blocks
    Tensor< double, nDim, nNodes, nDim, nNodes > k_UU( 0.0 );
    Tensor< double, nDim, nNodes, nNodes >       k_UN( 0.0 );
    Tensor< double, nNodes, nDim, nNodes >       k_NU( 0.0 );
    Tensor< double, nNodes, nNodes >             k_NN( 0.0 );

    Eigen::Map< Eigen::VectorXd > rhs( rightHandSide, sizeLoadVector );

    for ( auto& qp : qps ) {

      using namespace Marmot::FastorIndices;

      auto        N_    = this->N( qp.xi );
      const auto& dNdX_ = qp.dNdX;

      const auto N    = Tensor< double, nNodes >( N_.data() );
      const auto dNdX = Tensor< double, nDim, nNodes >( dNdX_.data(), ColumnMajor );

      const auto F_np = evaluate( einsum< Ai, jA >( qU_np, dNdX ) + I );

      const double nonlocalField = inner( N, qN_np );

      const typename Material::Deformation< nDim > deformation = { F_np, nonlocalField };

      const typename Material::TimeIncrement timeIncrement{ time, dT };

      typename Material::ConstitutiveResponse< nDim >
        response( Tensor< double, nDim, nDim >( qp.managedStateVars->stress.data(), ColumnMajor ),
                  0.0,
                  0.0,
                  0.0,
                  0.0,
                  qp.managedStateVars->materialStateVars.data() );
      typename Material::AlgorithmicModuli< nDim > tangents;

      if constexpr ( nDim == 2 ) {

        if ( sectionType == SectionType::PlaneStrain ) {

          using namespace Marmot;

          typename Material::ConstitutiveResponse< 3 >
            response3D( FastorStandardTensors::Tensor33d( qp.managedStateVars->stress.data(), Fastor::ColumnMajor ),
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        qp.managedStateVars->materialStateVars.data() );

          typename Material::AlgorithmicModuli< 3 > algorithmicModuli3D;

          typename Material::Deformation< 3 > deformation3D{ expandTo3D( deformation.F ), deformation.N };
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

          response.tau            = reduceTo2D< U, U >( response3D.tau );
          response.L              = response3D.L;
          response.nonLocalRadius = response3D.nonLocalRadius;

          tangents.dTau_dF = reduceTo2D< U, U, U, U >( algorithmicModuli3D.dTau_dF );
          tangents.dTau_dN = reduceTo2D< U, U >( algorithmicModuli3D.dTau_dN );
          tangents.dL_dF   = reduceTo2D< U, U >( algorithmicModuli3D.dL_dF );
          tangents.dL_dN   = algorithmicModuli3D.dL_dN;

          qp.managedStateVars->stress = Marmot::mapEigenToFastor( response3D.tau ).reshaped();
        }
        else {
          throw std::runtime_error( "Plane stress is not implemented for gradient-enhanced finite strain materials." );
        }
      }
      else {
        qp.material->computeStress( response, tangents, deformation, timeIncrement );
        qp.managedStateVars->stress = Marmot::mapEigenToFastor( response.tau ).reshaped();
      }

      const auto dNdx = evaluate( einsum< ji, jA >( inv( F_np ), dNdX ) );

      const double& J0xW = qp.J0xW;

      const auto&  S          = response.tau;
      const double localField = response.L;
      const double c          = response.nonLocalRadius * response.nonLocalRadius;

      const auto& t = tangents;

      // aux stiffness tensors
      const auto dS_dqU = evaluate( +einsum< ijkl, lB >( t.dTau_dF, dNdX ) );
      const auto dS_dqN = evaluate( +einsum< ij, B >( t.dTau_dN, N ) );
      const auto dL_dqU = evaluate( +einsum< kl, lB >( t.dL_dF, dNdX ) );

      // residuals (Pe = +internal force, matching DisplacementFiniteStrainULElement)
      r_U += ( +einsum< iA, ij >( dNdx, S ) ) * J0xW;
      r_N += ( N * nonlocalField + c * einsum< iA, iB, B >( dNdX, dNdX, qN_np ) - N * localField ) * J0xW;

      // stiffness
      k_UU += ( +einsum< iA, ijkB, to_jAkB >( dNdx, dS_dqU ) ) * J0xW;
      k_UU += ( -einsum< kA, ij, iB, to_jAkB >( dNdx, S, dNdx ) ) * J0xW; // geometric
      k_UN += ( +einsum< iA, ijB, to_jAB >( dNdx, dS_dqN ) ) * J0xW;
      k_NU += ( -einsum< A, kB >( N, dL_dqU ) ) * J0xW;
      k_NN += ( +einsum< A, B >( N, N ) + einsum< iA, iB >( dNdX, dNdX ) * c ) * J0xW;
      k_NN += ( -einsum< A, B >( N, N ) * t.dL_dN ) * J0xW;
    }

    using namespace Eigen;

    Map< KSizedMatrix > K( stiffnessMatrix );

    K.template block< bsU, bsU >( idxU, idxU ) += Map< Matrix< double, bsU, bsU > >( torowmajor( k_UU ).data() );
    K.template block< bsU, bsN >( idxU, idxN ) += Map< Matrix< double, bsU, bsN > >( torowmajor( k_UN ).data() );
    K.template block< bsN, bsU >( idxN, idxU ) += Map< Matrix< double, bsN, bsU > >( torowmajor( k_NU ).data() );
    K.template block< bsN, bsN >( idxN, idxN ) += Map< Matrix< double, bsN, bsN > >( torowmajor( k_NN ).data() );
  }

  template < int nDim, int nNodes >
  void GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::computeKernelsExplicit( const double*,
                                                                                                const double*,
                                                                                                double*,
                                                                                                double,
                                                                                                double )
  {
    throw std::runtime_error(
      "Explicit kernels are not implemented for the gradient-enhanced finite strain displacement element." );
  }

  template < int nDim, int nNodes >
  void GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::computeDistributedLoad(
    MarmotElement::DistributedLoadTypes loadType,
    double*                             rightHandSide,
    double*                             stiffnessMatrix,
    const int                           elementFace,
    const double*                       load,
    const double*                       QTotal_,
    double                              time,
    double                              dT )
  {
    Eigen::Map< USizedVector > r_U( rightHandSide );

    switch ( loadType ) {

    case MarmotElement::Pressure: {
      const double                                                          p = load[0];
      const Eigen::Map< const RhsSized >                                    QTotal( QTotal_ );
      const Eigen::Ref< const USizedVector >                                qU( QTotal.head( bsU ) );
      Eigen::Map< Eigen::Matrix< double, sizeLoadVector, sizeLoadVector > > K( stiffnessMatrix );

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
  void GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::setInitialConditions(
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

          qp.material->initializeYourself( qp.managedStateVars->materialStateVars.data(),
                                           qp.managedStateVars->materialStateVars.size() );
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
                      F0_ZZ] = qp.material
                                 ->findEigenDeformationForEigenStress( { qp.managedStateVars->F0_XX,
                                                                         qp.managedStateVars->F0_YY,
                                                                         qp.managedStateVars->F0_ZZ },
                                                                       geostaticNormalStressComponents,
                                                                       qp.managedStateVars->materialStateVars.data() );

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
  void GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::computeBodyForce( double*       rightHandSide,
                                                                                          double*       stiffnessMatrix,
                                                                                          const double* load,
                                                                                          const double* qTotal,
                                                                                          double        time,
                                                                                          double        dT )
  {
    Eigen::Map< RhsSized >                                     r( rightHandSide );
    Eigen::Ref< USizedVector >                                 r_U( r.head( bsU ) );
    const Eigen::Map< const Eigen::Matrix< double, nDim, 1 > > f( load );

    for ( const auto& qp : qps )
      r_U += this->NB( this->N( qp.xi ) ).transpose() * f * qp.J0xW;
  }

  template < int nDim, int nNodes >
  std::vector< double > GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::getCoordinatesAtCenter()
  {
    std::vector< double > coords( nDim );

    Eigen::Map< XiSized > coordsMap( &coords[0] );
    const auto            centerXi = XiSized::Zero();
    coordsMap                      = this->NB( this->N( centerXi ) ) * this->coordinates;
    return coords;
  }

  template < int nDim, int nNodes >
  std::vector< std::vector< double > > GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::
    getCoordinatesAtQuadraturePoints()
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
  int GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::getNumberOfQuadraturePoints()
  {
    return qps.size();
  }

} // namespace Marmot::Elements
