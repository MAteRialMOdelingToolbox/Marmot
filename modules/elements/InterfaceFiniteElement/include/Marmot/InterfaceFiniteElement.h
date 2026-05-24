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
 * Alexandros Stathas alexandros.stathas@boku.ac.at
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

/**
 * @file InterfaceFiniteElement.h
 * @brief Interface finite element formulation for displacement-jump based interface mechanics.
 *
 * This file defines the templated class `Marmot::Elements::InterfaceFiniteElement`
 * and its full in-header implementation. The element evaluates interface traction-like
 * quantities from displacement jumps and average surface kinematics and assembles
 * residual/tangent contributions at quadrature points.
 */
#pragma once

#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryInterfaceElement.h"
#include "Marmot/MarmotInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotTypedefs.h"

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

namespace Marmot::Elements {

  /**
   * @class InterfaceFiniteElement
   * @tparam nDim Spatial embedding dimension.
   * @tparam nNodes Number of element nodes.
   * @brief Interface finite element with displacement-jump kinematics.
   *
   * The element combines geometric interface operators from
   * `MarmotGeometryInterfaceElement<nDim, nNodes>` with an interface-material
   * update (`MarmotInterfaceMaterialHypoElastic`) at each quadrature point.
   * It stores quadrature-point state variables and assembles:
   * - element residual vector,
   * - algorithmic tangent matrix,
   * - zero inertia terms (current formulation).
   */
  template < int nDim, int nNodes >
  class InterfaceFiniteElement : public MarmotElement, public MarmotGeometryInterfaceElement< nDim, nNodes > {

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /**
     * @brief Section type selector for the interface element.
     */
    enum SectionType {
      Interface,
    };

    static constexpr int nDofPerNodeU    = nDim;
    static constexpr int nInterfaceNodes = nNodes / 2;

    static constexpr int sizeLoadVector = nNodes * nDim;
    static constexpr int nCoordinates   = nNodes * nDim;

    static constexpr int nSideDofs = nInterfaceNodes * nDim;
    static constexpr int nTensor   = nDim * nDim;

    using ParentGeometryElement = MarmotGeometryInterfaceElement< nDim, nNodes >;

    using XiSized              = typename ParentGeometryElement::XiSized;
    using NSized               = typename ParentGeometryElement::NSized;
    using dNdXiSized           = typename ParentGeometryElement::dNdXiSized;
    using SurfaceJacobianSized = typename ParentGeometryElement::SurfaceJacobianSized;
    using MetricSized          = typename ParentGeometryElement::MetricSized;
    using GradSized            = typename ParentGeometryElement::GradSized;

    using VectorDim = typename ParentGeometryElement::VectorDim;
    using TensorDim = typename ParentGeometryElement::TensorDim;

    using NMatrixSized     = typename ParentGeometryElement::NMatrixSized;
    using NJumpMatrixSized = typename ParentGeometryElement::NJumpMatrixSized;

    using BSurfaceSized    = typename ParentGeometryElement::BSurfaceSized;
    using BAvgSurfaceSized = typename ParentGeometryElement::BAvgSurfaceSized;

    using RhsSized      = Eigen::Matrix< double, sizeLoadVector, 1 >;
    using KeSizedMatrix = Eigen::Matrix< double, sizeLoadVector, sizeLoadVector >;

    using ForceSized                = Eigen::Matrix< double, nDim, 1 >;
    using SurfaceStressSized        = Eigen::Matrix< double, nTensor, 1 >;
    using InterfaceDisplSized       = Eigen::Matrix< double, 2 * nDim, 1 >;
    using InterfaceSurfaceGradSized = Eigen::Matrix< double, 2 * nTensor, 1 >;

    using QMatrixSized = Eigen::Matrix< double, nDim, nDim, Eigen::RowMajor >;
    using ZMatrixSized = Eigen::Matrix< double, nTensor, nTensor, Eigen::RowMajor >;
    using HMatrixSized = Eigen::Matrix< double, nDim, nTensor, Eigen::RowMajor >;
    using YMatrixSized = Eigen::Matrix< double, nTensor, nTensor, Eigen::RowMajor >;

    using Material = MarmotInterfaceMaterialHypoElastic;

    Eigen::Map< const Eigen::VectorXd > elementProperties;

    const int         elLabel;
    const SectionType sectionType;

    /**
     * @brief Container for all per-quadrature-point data.
     */
    struct QuadraturePoint {
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW

      const XiSized xi;
      const double  weight;

      double detJ;
      double sqrtDetG;
      double J0xW;

      NSized               N;
      dNdXiSized           dNdXi;
      SurfaceJacobianSized J;
      MetricSized          G;
      GradSized            gradN;

      VectorDim normal;
      TensorDim normalProjection;
      TensorDim tangentProjection;

      /*
       * One-side operators:
       *
       * NmatSide:
       *   maps one side's nodal displacement vector to u at the interface qp.
       *
       * BmatSide:
       *   maps one side's nodal displacement vector to the surface-gradient
       *   quantity passed to the material.
       */
      NMatrixSized  NmatSide;
      BSurfaceSized BmatSide;

      /*
       * Whole-element operators:
       *
       * NmatJump:
       *   jump u = u_top - u_bottom
       *   NmatJump = [ -Nside , +Nside ]
       *
       * BmatAverage:
       *   grad_s u_avg = 0.5 * (grad_s u_bottom + grad_s u_top)
       *   BmatAverage = 0.5 * [ Bside , Bside ]
       */
      NJumpMatrixSized NmatJump;
      BAvgSurfaceSized BmatAverage;

      /**
       * @brief Named state-variable manager for interface quadrature points.
       */
      class QPStateVarManager : public MarmotStateVarVectorManager {

        /*
         * Persistent state layout for accumulated force, surface stress,
         * displacement, surface strain, and material state variables.
         *
         * The displacement and surface strain entries store the accumulated
         * top/bottom quantities:
         *
         *   displacement   = [u_top, u_bottom]
         *   surface strain = [grad_s u_top, grad_s u_bottom]
         */
        inline const static auto layout = makeLayout( {
          { .name = "force", .length = nDim },
          { .name = "alignment padding", .length = nDim % 2 },
          { .name = "surface stress", .length = nDim * nDim },
          { .name = "displacement", .length = 2 * nDim },
          { .name = "surface strain", .length = 2 * nDim * nDim },
          { .name   = "state block alignment padding",
            .length = ( 4 - ( ( nDim + ( nDim % 2 ) + nDim * nDim + 2 * nDim + 2 * nDim * nDim ) % 4 ) ) % 4 },
          { .name = "begin of material state", .length = 0 },
        } );

      public:
        Eigen::Map< ForceSized >                force;
        Eigen::Map< SurfaceStressSized >        surfaceStress;
        Eigen::Map< InterfaceDisplSized >       displacement;
        Eigen::Map< InterfaceSurfaceGradSized > surfaceStrain;

        Eigen::Map< Eigen::VectorXd > materialStateVars;

        static int getNumberOfRequiredStateVarsQuadraturePointOnly() { return layout.nRequiredStateVars; }

        QPStateVarManager( double* theStateVarVector, int nStateVars )
          : MarmotStateVarVectorManager( theStateVarVector, layout ),
            force( &find( "force" ) ),
            surfaceStress( &find( "surface stress" ) ),
            displacement( &find( "displacement" ) ),
            surfaceStrain( &find( "surface strain" ) ),
            materialStateVars( &find( "begin of material state" ),
                               nStateVars - getNumberOfRequiredStateVarsQuadraturePointOnly() )
        {
        }
      };

      std::unique_ptr< QPStateVarManager > managedStateVars;
      std::unique_ptr< Material >          material;

      /**
       * @brief Number of non-material state variables stored at this quadrature point.
       */
      int getNumberOfRequiredStateVarsQuadraturePointOnly()
      {
        return QPStateVarManager::getNumberOfRequiredStateVarsQuadraturePointOnly();
      }

      /**
       * @brief Total number of state variables including material state.
       */
      int getNumberOfRequiredStateVars()
      {
        return getNumberOfRequiredStateVarsQuadraturePointOnly() + material->getNumberOfRequiredStateVars();
      }

      /**
       * @brief Assign state-variable memory to this quadrature point.
       */
      void assignStateVars( double* stateVars, int nStateVars )
      {
        managedStateVars = std::make_unique< QPStateVarManager >( stateVars, nStateVars );

        material->initializeYourself( managedStateVars->materialStateVars.data(),
                                      managedStateVars->materialStateVars.size() );
      }

      /**
       * @brief Construct a quadrature-point data container.
       */
      QuadraturePoint( XiSized xi, double weight )
        : xi( xi ),
          weight( weight ),
          detJ( 0.0 ),
          sqrtDetG( 0.0 ),
          J0xW( 0.0 ),
          N( NSized::Zero() ),
          dNdXi( dNdXiSized::Zero() ),
          J( SurfaceJacobianSized::Zero() ),
          G( MetricSized::Zero() ),
          gradN( GradSized::Zero() ),
          normal( VectorDim::Zero() ),
          normalProjection( TensorDim::Zero() ),
          tangentProjection( TensorDim::Zero() ),
          NmatSide( NMatrixSized::Zero() ),
          BmatSide( BSurfaceSized::Zero() ),
          NmatJump( NJumpMatrixSized::Zero() ),
          BmatAverage( BAvgSurfaceSized::Zero() )
      {
      }
    };

    std::vector< QuadraturePoint, Eigen::aligned_allocator< QuadraturePoint > > qps;

    /**
     * @brief Construct an interface finite element.
     * @param[in] elementID Element label.
     * @param[in] integrationType Quadrature rule type.
     * @param[in] sectionType Interface section type.
     */
    InterfaceFiniteElement( int                                         elementID,
                            FiniteElement::Quadrature::IntegrationTypes integrationType,
                            SectionType                                 sectionType = SectionType::Interface );

    /**
     * @brief Return number of required state variables per element.
     */
    int getNumberOfRequiredStateVars();

    /**
     * @brief Return nodal primary fields associated with this element.
     */
    std::vector< std::vector< std::string > > getNodeFields();

    /**
     * @brief Return DOF permutation pattern.
     */
    std::vector< int > getDofIndicesPermutationPattern();

    int getNNodes() { return nNodes; }

    int getNSpatialDimensions() { return nDim; }

    int getNDofPerElement() { return sizeLoadVector; }

    /*
     * Important:
     * ParentGeometryElement::getElementShape() returns the computational
     * interface shape, e.g. "iquad4".
     *
     * EnSight/ParaView do not understand "iquad4" as a geometry keyword.
     * Therefore, for output/visualization we map interface elements to valid
     * EnSight element names.
     *
     * The computational element remains an interface element. This only affects
     * the geometry keyword written to result files.
     */
    std::string getElementShape()
    {
      if constexpr ( nDim == 3 && nNodes == 8 ) {
        return "hexa8";
      }
      else if constexpr ( nDim == 2 && nNodes == 4 ) {
        return "line2";
      }
      else {
        return ParentGeometryElement::getElementShape();
      }
    }

    /**
     * @brief Assign element state-variable memory to quadrature points.
     */
    void assignStateVars( double* stateVars, int nStateVars );

    /**
     * @brief Assign element-level properties (e.g., interface thickness).
     */
    void assignProperty( const ElementProperties& marmotElementProperty );

    /**
     * @brief Assign a material section to all quadrature points.
     */
    void assignProperty( const MarmotMaterialSection& marmotElementProperty );

    /**
     * @brief Assign a material by name and property array to all quadrature points.
     */
    void assignMaterial( const std::string& materialName, const double* materialProperties, int nMaterialProperties );

    /**
     * @brief Assign nodal coordinates.
     */
    void assignNodeCoordinates( const double* coordinates );

    /**
     * @brief Initialize element geometric operators and quadrature-point geometry data.
     */
    void initializeYourself();

    /**
     * @brief Set initial conditions for supported state categories.
     */
    void setInitialConditions( StateTypes state, const double* values );

    /**
     * @brief Distributed-load routine (currently not implemented).
     */
    void computeDistributedLoad( MarmotElement::DistributedLoadTypes loadType,
                                 double*                             P,
                                 double*                             K,
                                 const int                           elementFace,
                                 const double*                       load,
                                 const double*                       QTotal,
                                 const double*                       time,
                                 double                              dT );

    /**
     * @brief Body-force routine (currently not implemented).
     */
    void computeBodyForce( double*       P,
                           double*       K,
                           const double* load,
                           const double* QTotal,
                           const double* time,
                           double        dT );

    /**
     * @brief Assemble residual and tangent for one increment.
     */
    void computeYourself( const double* QTotal,
                          const double* dQ,
                          double*       Pe,
                          double*       Ke,
                          const double* time,
                          double        dT,
                          double&       pNewdT );

    /**
     * @brief Compute consistent mass contribution (currently zero matrix).
     */
    void computeConsistentInertia( double* M );

    /**
     * @brief Compute lumped mass contribution (currently zero vector).
     */
    void computeLumpedInertia( double* M );

    /**
     * @brief Access a named state view at a quadrature point.
     */
    StateView getStateView( const std::string& stateName, int qpNumber )
    {
      const auto& qp = qps[qpNumber];

      if ( qp.managedStateVars->contains( stateName ) ) {
        return qp.managedStateVars->getStateView( stateName );
      }

      if ( stateName == "sdv" ) {
        std::cout << __PRETTY_FUNCTION__ << " on 'sdv' is discouraged and deprecated, please use precise state name";
        return { qp.managedStateVars->materialStateVars.data(),
                 static_cast< int >( qp.managedStateVars->materialStateVars.size() ) };
      }

      return qp.material->getStateView( stateName, qp.managedStateVars->materialStateVars.data() );
    }

    /**
     * @brief Coordinates of the element center (on reference interface side).
     */
    std::vector< double > getCoordinatesAtCenter();

    /**
     * @brief Coordinates of all quadrature points (on reference interface side).
     */
    std::vector< std::vector< double > > getCoordinatesAtQuadraturePoints();

    /**
     * @brief Number of quadrature points.
     */
    int getNumberOfQuadraturePoints();
  };

  /**
   * @name Template implementation (header-only)
   * @brief In-header method definitions for `InterfaceFiniteElement`.
   */
  ///@{

  template < int nDim, int nNodes >
  InterfaceFiniteElement< nDim, nNodes >::InterfaceFiniteElement(
    int                                         elementID,
    FiniteElement::Quadrature::IntegrationTypes integrationType,
    SectionType                                 sectionType )
    : ParentGeometryElement(),
      elementProperties( Eigen::Map< const Eigen::VectorXd >( nullptr, 0 ) ),
      elLabel( elementID ),
      sectionType( sectionType )
  {
    const auto qpInfos = FiniteElement::Quadrature::getGaussPointInfo( this->shape, integrationType );

    for ( const auto& qpInfo : qpInfos ) {
      QuadraturePoint qp( qpInfo.xi, qpInfo.weight );
      qps.push_back( std::move( qp ) );
    }
  }

  template < int nDim, int nNodes >
  int InterfaceFiniteElement< nDim, nNodes >::getNumberOfRequiredStateVars()
  {
    return qps[0].getNumberOfRequiredStateVars() * qps.size();
  }

  template < int nDim, int nNodes >
  std::vector< std::vector< std::string > > InterfaceFiniteElement< nDim, nNodes >::getNodeFields()
  {
    using namespace std;

    static vector< vector< string > > nodeFields;

    if ( nodeFields.empty() ) {
      for ( int i = 0; i < nNodes; i++ ) {
        nodeFields.push_back( vector< string >() );
        nodeFields[i].push_back( "displacement" );
      }
    }

    return nodeFields;
  }

  template < int nDim, int nNodes >
  std::vector< int > InterfaceFiniteElement< nDim, nNodes >::getDofIndicesPermutationPattern()
  {
    static std::vector< int > permutationPattern;

    if ( permutationPattern.empty() ) {
      for ( int i = 0; i < nNodes * nDim; i++ )
        permutationPattern.push_back( i );
    }

    return permutationPattern;
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::assignStateVars( double* stateVars, int nStateVars )
  {
    const int nQpStateVars = nStateVars / qps.size();

    for ( size_t i = 0; i < qps.size(); i++ ) {
      auto&   qp          = qps[i];
      double* qpStateVars = stateVars + ( i * nQpStateVars );
      qp.assignStateVars( qpStateVars, nQpStateVars );
    }
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::assignProperty( const ElementProperties& elementPropertiesInfo )
  {
    new ( &elementProperties ) Eigen::Map< const Eigen::VectorXd >( elementPropertiesInfo.elementProperties,
                                                                    elementPropertiesInfo.nElementProperties );
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::assignProperty( const MarmotMaterialSection& section )
  {
    for ( auto& qp : qps ) {
      qp.material = std::unique_ptr< Material >( dynamic_cast< Material* >(
        MarmotLibrary::MarmotInterfaceMaterialHypoElasticFactory::createMaterial( section.materialName,
                                                                                  section.materialProperties,
                                                                                  section.nMaterialProperties,
                                                                                  elLabel ) ) );

      if ( !qp.material ) {
        throw std::invalid_argument(
          MakeString() << __PRETTY_FUNCTION__
                       << ": invalid material assigned; cannot cast to MarmotInterfaceMaterialHypoElastic!" );
      }
    }
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::assignMaterial( const std::string& materialName,
                                                               const double*      materialProperties,
                                                               int                nMaterialProperties )
  {
    for ( auto& qp : qps ) {
      qp.material = std::unique_ptr< Material >( dynamic_cast< Material* >(
        MarmotLibrary::MarmotInterfaceMaterialHypoElasticFactory::createMaterial( materialName,
                                                                                  materialProperties,
                                                                                  nMaterialProperties,
                                                                                  elLabel ) ) );

      if ( !qp.material ) {
        throw std::invalid_argument(
          MakeString() << __PRETTY_FUNCTION__
                       << ": invalid material assigned; cannot cast to MarmotInterfaceMaterialHypoElastic!" );
      }
    }
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::assignNodeCoordinates( const double* coordinates )
  {
    ParentGeometryElement::assignNodeCoordinates( coordinates );
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::initializeYourself()
  {
    const double thickness = elementProperties.size() > 0 ? elementProperties[0] : 1.0;

    for ( QuadraturePoint& qp : qps ) {
      const bool fullyProjectedB = ( nDim == 3 );
      const auto geom            = this->evaluateAt( qp.xi, 0, fullyProjectedB );

      qp.N                 = geom.N;
      qp.dNdXi             = geom.dNdXi;
      qp.J                 = geom.J;
      qp.G                 = geom.G;
      qp.sqrtDetG          = geom.sqrtDetG;
      qp.detJ              = geom.sqrtDetG;
      qp.gradN             = geom.gradN;
      qp.normal            = geom.n;
      qp.normalProjection  = geom.normalProjection;
      qp.tangentProjection = geom.tangentProjection;

      qp.NmatSide    = geom.NmatSide;
      qp.BmatSide    = geom.BmatSide;
      qp.NmatJump    = geom.NmatJump;
      qp.BmatAverage = geom.BmatAverage;

      qp.J0xW = qp.weight * qp.sqrtDetG * thickness;

      if ( qp.material ) {
        if constexpr ( nDim == 3 ) {
          qp.material->setCharacteristicElementLength( std::sqrt( qp.sqrtDetG ) );
        }
        else if constexpr ( nDim == 2 ) {
          qp.material->setCharacteristicElementLength( qp.sqrtDetG );
        }
      }
    }
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::computeYourself( const double* QTotal_,
                                                                const double* dQ_,
                                                                double*       Pe_,
                                                                double*       Ke_,
                                                                const double* time,
                                                                double        dT,
                                                                double&       pNewDT )
  {
    (void)QTotal_;

    Eigen::Map< const RhsSized > dQ( dQ_ );
    Eigen::Map< KeSizedMatrix >  Ke( Ke_ );
    Eigen::Map< RhsSized >       Pe( Pe_ );

    constexpr int halfSize = nNodes * nDim / 2;

    for ( QuadraturePoint& qp : qps ) {
      const auto& Nside = qp.NmatSide;
      const auto& Bside = qp.BmatSide;
      const auto& Njump = qp.NmatJump;
      const auto& Bavg  = qp.BmatAverage;

      const auto dQBottom = dQ.template segment< halfSize >( 0 );
      const auto dQTop    = dQ.template segment< halfSize >( halfSize );

      InterfaceDisplSized dU_GPs;
      dU_GPs.template segment< nDim >( 0 )    = Nside * dQTop;
      dU_GPs.template segment< nDim >( nDim ) = Nside * dQBottom;

      InterfaceSurfaceGradSized dSurface_strain_GPs;
      dSurface_strain_GPs.template segment< nTensor >( 0 )       = Bside * dQTop;
      dSurface_strain_GPs.template segment< nTensor >( nTensor ) = Bside * dQBottom;

      ForceSized         force          = qp.managedStateVars->force;
      SurfaceStressSized surface_stress = qp.managedStateVars->surfaceStress;

      QMatrixSized Q_ij;
      ZMatrixSized Z_ijkl;
      HMatrixSized H_ijk;
      YMatrixSized Y_ijkl;

      Q_ij.setZero();
      Z_ijkl.setZero();
      H_ijk.setZero();
      Y_ijkl.setZero();

      Material::State         materialState{ force.data(),
                                     surface_stress.data(),
                                     qp.managedStateVars->materialStateVars.data() };
      Material::Tangents      materialTangents{ Q_ij.data(), Z_ijkl.data(), H_ijk.data(), Y_ijkl.data() };
      Material::Deformation   materialDeformation{ dU_GPs.data(), dSurface_strain_GPs.data(), qp.normal.data() };
      Material::TimeIncrement materialTimeIncrement{ time, dT, pNewDT };

      qp.material->computeStress( materialState, materialTangents, materialDeformation, materialTimeIncrement );

      if ( pNewDT < 1.0 )
        return;

      qp.managedStateVars->force         = force;
      qp.managedStateVars->surfaceStress = surface_stress;
      qp.managedStateVars->displacement += dU_GPs;
      qp.managedStateVars->surfaceStrain += dSurface_strain_GPs;

      Pe -= Njump.transpose() * force * qp.J0xW;
      Pe -= Bavg.transpose() * surface_stress * qp.J0xW;

      Ke += ( Njump.transpose() * Q_ij * Njump + Bavg.transpose() * Z_ijkl * Bavg + Bavg.transpose() * Y_ijkl * Bavg +
              Njump.transpose() * H_ijk * Bavg + Bavg.transpose() * H_ijk.transpose() * Njump ) *
            qp.J0xW;
    }
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::setInitialConditions( StateTypes state, const double* values )
  {
    switch ( state ) {
    case MarmotElement::MarmotMaterialInitialization: {
      /* Interface material state is managed through qp.managedStateVars in v26.05 style. */
      break;
    }

    case MarmotElement::MarmotMaterialStateVars: {
      throw std::invalid_argument( "Please use initializeStateVars directly on material" );
    }

    default:
      throw std::invalid_argument( MakeString()
                                   << __PRETTY_FUNCTION__ << ": invalid initial condition for InterfaceFiniteElement" );
    }
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::computeDistributedLoad( MarmotElement::DistributedLoadTypes loadType,
                                                                       double*                             P,
                                                                       double*                             K,
                                                                       const int                           elementFace,
                                                                       const double*                       load,
                                                                       const double*                       QTotal,
                                                                       const double*                       time,
                                                                       double                              dT )
  {
    throw std::invalid_argument(
      MakeString() << __PRETTY_FUNCTION__ << ": distributed loads are not implemented for InterfaceFiniteElement." );
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::computeBodyForce( double*       P,
                                                                 double*       K,
                                                                 const double* load,
                                                                 const double* QTotal,
                                                                 const double* time,
                                                                 double        dT )
  {
    throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__
                                              << ": body forces are not implemented for InterfaceFiniteElement." );
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::computeConsistentInertia( double* M )
  {
    Eigen::Map< KeSizedMatrix > Me( M );
    Me.setZero();
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::computeLumpedInertia( double* M )
  {
    Eigen::Map< RhsSized > Me( M );
    Me.setZero();
  }

  template < int nDim, int nNodes >
  std::vector< double > InterfaceFiniteElement< nDim, nNodes >::getCoordinatesAtCenter()
  {
    std::vector< double > coords( nDim );

    Eigen::Map< VectorDim > coordsMap( coords.data() );

    const auto centerXi = XiSized::Zero();
    const auto Ncenter  = this->N( centerXi );
    const auto Nmat     = this->NMatrix( Ncenter );

    const auto xSide = this->getSideCoordinates( 0 );

    coordsMap = Nmat * xSide;

    return coords;
  }

  template < int nDim, int nNodes >
  std::vector< std::vector< double > > InterfaceFiniteElement< nDim, nNodes >::getCoordinatesAtQuadraturePoints()
  {
    std::vector< std::vector< double > > listedCoords;

    for ( const auto& qp : qps ) {
      std::vector< double > coords( nDim );

      Eigen::Map< VectorDim > coordsMap( coords.data() );

      const auto xSide = this->getSideCoordinates( 0 );
      coordsMap        = qp.NmatSide * xSide;

      listedCoords.push_back( coords );
    }

    return listedCoords;
  }

  template < int nDim, int nNodes >
  int InterfaceFiniteElement< nDim, nNodes >::getNumberOfQuadraturePoints()
  {
    return static_cast< int >( qps.size() );
  }

  ///@}

} // namespace Marmot::Elements
