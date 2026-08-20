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

#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryElement.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotTypedefs.h"
#include <memory>
#include <vector>

namespace Marmot::Elements {

  /**
   * @class Marmot::Elements::DisplacementPressureFiniteStrainElement
   * @tparam nNodes         Number of nodes carrying the displacement field (a fully-integrated
   *                        quadratic serendipity shape: 10 for Tetra10, 20 for Hexa20).
   * @tparam nPressureNodes Number of nodes carrying the (scalar) pressure field. Must be the
   *                        corner-node subset of the displacement nodes (4 for Tetra10, 8 for
   *                        Hexa20), listed first in the element's node ordering.
   *
   * @brief Mixed displacement-pressure finite-strain solid element for (nearly) incompressible
   * hyperelasticity, following the classical Simo-Taylor-Pister two-field formulation.
   *
   * @details The functional
   * \f[
   *   \Pi(\boldsymbol u, p) = \int_{\Omega_0} \Psi_\mathrm{iso}(\boldsymbol F) \, \mathrm dV_0
   *   + \int_{\Omega_0} p\,(J-1) \, \mathrm dV_0
   * \f]
   * is stationary w.r.t. both the displacement field \f$\boldsymbol u\f$ and an independently
   * interpolated pressure field \f$p\f$ (the Lagrange multiplier enforcing the incompressibility
   * constraint \f$J=\det\boldsymbol F=1\f$ weakly). The total Kirchhoff stress is
   * \f$\boldsymbol\tau = \boldsymbol\tau_\mathrm{iso}(\boldsymbol F) + p\,J\,\boldsymbol I\f$, where
   * \f$\boldsymbol\tau_\mathrm{iso}\f$ is supplied by the attached (purely isochoric)
   * MarmotMaterialFiniteStrain instance (e.g. Marmot::Materials::Ogden). Any material implementing
   * that interface can be used, as long as it carries no volumetric energy term of its own -
   * otherwise the volumetric response would be counted twice.
   *
   * The displacement field is interpolated with the full quadratic serendipity shape
   * (\p nNodes), the pressure field with the continuous linear shape on the corner-node subset
   * (\p nPressureNodes). This P2/P1-type pairing is inf-sup (LBB) stable, so no pressure
   * stabilization or selective/reduced integration is required - both fields are integrated with
   * the same full Gauss rule.
   *
   * @note Restricted to 3D solid elements (Tetra10, Hexa20). A 2D plane-strain variant would
   * additionally require the finite-strain plane-strain 3D-embedding machinery used by
   * DisplacementFiniteStrainULElement and is not implemented here.
   */
  template < int nNodes, int nPressureNodes >
  class DisplacementPressureFiniteStrainElement : public MarmotElement {

  public:
    /// @brief Number of spatial dimensions (fixed to 3; see class-level note).
    static constexpr int nDim = 3;

    /// @brief Displacement degrees of freedom per node.
    static constexpr int nDofPerNodeU = nDim;
    /// @brief Pressure degrees of freedom per node (scalar).
    static constexpr int nDofPerNodeP = 1;

    /// @brief Block size for the displacement field.
    static constexpr int sizeDoFU = nNodes * nDofPerNodeU;
    /// @brief Block size for the pressure field.
    static constexpr int sizeDoFP = nPressureNodes * nDofPerNodeP;
    /// @brief Total number of degrees of freedom of the element.
    static constexpr int sizeLoadVector = sizeDoFU + sizeDoFP;

    /// @brief Geometry helper for the displacement field.
    using DisplacementGeometryElement = MarmotGeometryElement< nDim, nNodes >;
    /// @brief Geometry helper for the pressure field (corner-node subset).
    using PressureGeometryElement = MarmotGeometryElement< nDim, nPressureNodes >;
    /// @brief Material interface expected by this element.
    using Material = MarmotMaterialFiniteStrain;

    using JacobianSized = typename DisplacementGeometryElement::JacobianSized;
    using dNdXiSizedU   = typename DisplacementGeometryElement::dNdXiSized;
    using NSizedP       = typename PressureGeometryElement::NSized;
    using XiSized       = typename DisplacementGeometryElement::XiSized;

    using RhsSized      = Eigen::Matrix< double, sizeLoadVector, 1 >;
    using KeSizedMatrix = Eigen::Matrix< double, sizeLoadVector, sizeLoadVector >;
    using USizedVector  = Eigen::Matrix< double, sizeDoFU, 1 >;
    using PSizedVector  = Eigen::Matrix< double, sizeDoFP, 1 >;

    /// @brief Geometry element for the displacement field.
    DisplacementGeometryElement displacementGeometryElement;
    /// @brief Geometry element for the pressure field.
    PressureGeometryElement pressureGeometryElement;

    /// @brief Element label (ID).
    const int elLabel;

    /// @struct QuadraturePoint
    /// @brief Data and state associated with a quadrature point.
    struct QuadraturePoint {

      const XiSized xi;
      const double  weight;

      double      detJ;
      double      J0xW;
      dNdXiSizedU dNdX; ///< Physical (reference) shape function derivatives of the displacement field.
      NSizedP     Np;   ///< Pressure shape function values.

      /// @class QPStateVarManager
      /// @brief Named state variable layout at the quadrature point.
      class QPStateVarManager : public MarmotStateVarVectorManager {

        /// \hideinitializer
        inline const static auto layout = makeLayout( {
          { .name = "stress", .length = 9 },
          { .name = "total strain energy", .length = 1 },
          { .name = "elastic energy", .length = 1 },
          { .name = "dissipation", .length = 1 },
          { .name = "begin of material state", .length = 0 },
        } );

      public:
        Eigen::Map< Marmot::Vector9d > stress;
        double&                        totalStrainEnergy;
        double&                        elasticEnergy;
        double&                        dissipation;
        Eigen::Map< Eigen::VectorXd >  materialStateVars;

        static int getNumberOfRequiredStateVarsQuadraturePointOnly() { return layout.nRequiredStateVars; };

        QPStateVarManager( double* theStateVarVector, int nStateVars )
          : MarmotStateVarVectorManager( theStateVarVector, layout ),
            stress( &find( "stress" ) ),
            totalStrainEnergy( find( "total strain energy" ) ),
            elasticEnergy( find( "elastic energy" ) ),
            dissipation( find( "dissipation" ) ),
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
        : xi( xi ), weight( weight ), detJ( 0.0 ), J0xW( 0.0 ), dNdX( dNdXiSizedU::Zero() ), Np( NSizedP::Zero() ){};
    };

    /// @brief Quadrature points owned by the element.
    std::vector< QuadraturePoint > qps;

    /**
     * @brief Construct the element.
     * @param elementID       Element label.
     * @param integrationType Full or reduced integration; the same rule is shared by both fields
     * (see class-level note on inf-sup stability).
     */
    DisplacementPressureFiniteStrainElement( int                                                 elementID,
                                             Marmot::FiniteElement::Quadrature::IntegrationTypes integrationType );

    int getNumberOfRequiredStateVars();

    /** @brief Node fields: "displacement" for every node, additionally "pressure" for the first
     * \p nPressureNodes nodes. */
    std::vector< std::vector< std::string > > getNodeFields();

    std::vector< int > getDofIndicesPermutationPattern();

    int getNNodes() { return nNodes; }

    int getNSpatialDimensions() { return nDim; }

    int getNDofPerElement() { return sizeLoadVector; }

    std::string getElementShape() { return displacementGeometryElement.getElementShape(); }

    void assignStateVars( double* stateVars, int nStateVars );

    /// @brief Creates the per-quadrature-point isochoric finite-strain material instance.
    void assignProperty( const MarmotMaterialSection& section );

    /** @brief Assigns nodal coordinates to both the displacement and pressure geometry elements.
     * @note Assumes that the corner (pressure) nodes are listed before the higher-order
     * (displacement-only) nodes. */
    void assignNodeCoordinates( const double* coordinates );

    void initializeYourself();

    /** @brief Initializes material state.
     * @details Only Marmot::MarmotMaterialInitialization is supported. */
    void setInitialConditions( StateTypes state, const double* values );

    void computeDistributedLoad( MarmotElement::DistributedLoadTypes loadType,
                                 double*                             Pext,
                                 double*                             K,
                                 int                                 elementFace,
                                 const double*                       load,
                                 const double*                       QTotal,
                                 double                              time,
                                 double                              dT );

    void computeBodyForce( double* Pext, double* K, const double* load, const double* QTotal, double time, double dT );

    /**
     * @brief Assembles the element residual and tangent stiffness.
     * @details Per quadrature point: computes \f$\boldsymbol F\f$ from the total displacement
     * dofs, evaluates the attached material for \f$\boldsymbol\tau_\mathrm{iso}\f$ and
     * \f$\partial\boldsymbol\tau_\mathrm{iso}/\partial\boldsymbol F\f$, forms the total stress
     * \f$\boldsymbol\tau=\boldsymbol\tau_\mathrm{iso}+p J\boldsymbol I\f$ and its derivative, then
     * assembles
     * \f[
     *   \boldsymbol R_u = \sum_{qp} (\nabla_x \boldsymbol N_u)^{\!\top} \boldsymbol\tau \, J_0 w,
     *   \qquad
     *   R_p = \sum_{qp} \boldsymbol N_p^{\!\top} (J-1) \, J_0 w,
     * \f]
     * with tangent blocks \f$K_{uu}\f$ (material + geometric stiffness, exactly as for a
     * single-field UL finite-strain element, but using the total stress/tangent),
     * \f$K_{up}=K_{pu}^{\!\top} = \sum_{qp} J\,(\nabla_x \boldsymbol N_u) \otimes \boldsymbol N_p \,
     * J_0 w\f$, and \f$K_{pp}=0\f$ (the functional is linear in \f$p\f$, which is exactly why the
     * inf-sup stable P2/P1 pairing is required - no pressure stabilization compensates for a
     * vanishing \f$K_{pp}\f$ block).
     */
    void computeKernels( const double* QTotal, const double* dQ, double* Pint, double* K, double time, double dT );

    StateView getStateView( const std::string& stateName, int qpNumber );

    void computeInternalEnergy( double& internalEnergy );

    std::vector< double > getCoordinatesAtCenter();

    std::vector< std::vector< double > > getCoordinatesAtQuadraturePoints();

    int getNumberOfQuadraturePoints();
  };

  template < int nNodes, int nPressureNodes >
  DisplacementPressureFiniteStrainElement< nNodes, nPressureNodes >::DisplacementPressureFiniteStrainElement(
    int                                                 elementID,
    Marmot::FiniteElement::Quadrature::IntegrationTypes integrationType )
    : elLabel( elementID )
  {
    for ( const auto& qpInfo :
          Marmot::FiniteElement::Quadrature::getGaussPointInfo( displacementGeometryElement.shape, integrationType ) )
      qps.emplace_back( qpInfo.xi, qpInfo.weight );
  }

  template < int nNodes, int nPressureNodes >
  int DisplacementPressureFiniteStrainElement< nNodes, nPressureNodes >::getNumberOfRequiredStateVars()
  {
    return qps[0].getNumberOfRequiredStateVars() * qps.size();
  }

  template < int nNodes, int nPressureNodes >
  std::vector< std::vector< std::string > > DisplacementPressureFiniteStrainElement< nNodes,
                                                                                     nPressureNodes >::getNodeFields()
  {
    static const std::vector< std::vector< std::string > > nodeFields = [] {
      std::vector< std::vector< std::string > > nodeFields;
      for ( int i = 0; i < nNodes; i++ ) {
        nodeFields.push_back( { "displacement" } );
        if ( i < nPressureNodes )
          nodeFields[i].push_back( "pressure" );
      }
      return nodeFields;
    }();
    return nodeFields;
  }

  template < int nNodes, int nPressureNodes >
  std::vector< int > DisplacementPressureFiniteStrainElement< nNodes,
                                                              nPressureNodes >::getDofIndicesPermutationPattern()
  {
    static const std::vector< int > permutationPattern = [] {
      std::vector< int > permutationPattern;
      for ( int i = 0; i < nNodes; i++ )
        for ( int j = 0; j < nDim; j++ )
          permutationPattern.push_back( i * nDim + ( i < nPressureNodes ? i : nPressureNodes ) + j );
      for ( int i = 0; i < nPressureNodes; i++ )
        permutationPattern.push_back( i * ( nDim + nDofPerNodeP ) + nDim );
      return permutationPattern;
    }();
    return permutationPattern;
  }

  template < int nNodes, int nPressureNodes >
  void DisplacementPressureFiniteStrainElement< nNodes, nPressureNodes >::assignStateVars( double* stateVars,
                                                                                           int     nStateVars )
  {
    const int nQpStateVars = nStateVars / qps.size();
    for ( size_t i = 0; i < qps.size(); i++ ) {
      auto&   qp          = qps[i];
      double* qpStateVars = &stateVars[i * nQpStateVars];
      qp.assignStateVars( qpStateVars, nQpStateVars );
    }
  }

  template < int nNodes, int nPressureNodes >
  void DisplacementPressureFiniteStrainElement< nNodes, nPressureNodes >::assignProperty(
    const MarmotMaterialSection& section )
  {
    for ( auto& qp : qps ) {
      qp.material = std::unique_ptr< Material >(
        MarmotLibrary::MarmotMaterialFiniteStrainFactory::createMaterial( section.materialName,
                                                                          section.materialProperties,
                                                                          section.nMaterialProperties,
                                                                          elLabel ) );
    }
  }

  template < int nNodes, int nPressureNodes >
  void DisplacementPressureFiniteStrainElement< nNodes, nPressureNodes >::assignNodeCoordinates(
    const double* coordinates )
  {
    // Assumes that the corner nodes (pressure nodes) are listed before the higher-order nodes.
    displacementGeometryElement.assignNodeCoordinates( coordinates );
    pressureGeometryElement.assignNodeCoordinates( coordinates );
  }

  template < int nNodes, int nPressureNodes >
  void DisplacementPressureFiniteStrainElement< nNodes, nPressureNodes >::initializeYourself()
  {
    for ( QuadraturePoint& qp : qps ) {
      const dNdXiSizedU   dNdXi_ = displacementGeometryElement.dNdXi( qp.xi );
      const JacobianSized J      = displacementGeometryElement.Jacobian( dNdXi_ );
      const JacobianSized JInv   = J.inverse();

      qp.detJ = J.determinant();
      qp.dNdX = displacementGeometryElement.dNdX( dNdXi_, JInv );
      qp.Np   = pressureGeometryElement.N( qp.xi );
      qp.J0xW = qp.weight * qp.detJ;
    }
  }

  template < int nNodes, int nPressureNodes >
  void DisplacementPressureFiniteStrainElement< nNodes, nPressureNodes >::setInitialConditions( StateTypes    state,
                                                                                                const double* values )
  {
    switch ( state ) {
    case MarmotElement::MarmotMaterialInitialization: {
      for ( QuadraturePoint& qp : qps )
        qp.material->initializeYourself( qp.managedStateVars->materialStateVars.data(),
                                         qp.managedStateVars->materialStateVars.size() );
      break;
    }
    default:
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__
                                                << ": invalid or unsupported initial condition for this element" );
    }
  }

  template < int nNodes, int nPressureNodes >
  void DisplacementPressureFiniteStrainElement< nNodes, nPressureNodes >::computeBodyForce( double*       Pext_,
                                                                                            double*       K,
                                                                                            const double* load,
                                                                                            const double* QTotal,
                                                                                            double        time,
                                                                                            double        dT )
  {
    Eigen::Map< RhsSized >                                     Pext( Pext_ );
    Eigen::Ref< USizedVector >                                 fU( Pext.head( sizeDoFU ) );
    const Eigen::Map< const Eigen::Matrix< double, nDim, 1 > > f( load );

    for ( const auto& qp : qps )
      fU += displacementGeometryElement.NB( displacementGeometryElement.N( qp.xi ) ).transpose() * f * qp.J0xW;
  }

  template < int nNodes, int nPressureNodes >
  void DisplacementPressureFiniteStrainElement< nNodes, nPressureNodes >::computeDistributedLoad(
    MarmotElement::DistributedLoadTypes loadType,
    double*                             Pext_,
    double*                             K,
    int                                 elementFace,
    const double*                       load,
    const double*                       QTotal_,
    double                              time,
    double                              dT )
  {
    Eigen::Map< RhsSized >     Pext( Pext_ );
    Eigen::Ref< USizedVector > fU( Pext.head( sizeDoFU ) );

    switch ( loadType ) {
    case MarmotElement::Pressure: {
      const double p = load[0];

      Marmot::FiniteElement::BoundaryElement boundaryEl( displacementGeometryElement.shape,
                                                         elementFace,
                                                         nDim,
                                                         displacementGeometryElement.coordinates );

      const Eigen::VectorXd Pk = -p * boundaryEl.computeSurfaceNormalVectorialLoadVector();
      boundaryEl.assembleIntoParentVectorial( Pk, fU );
      break;
    }
    case MarmotElement::SurfaceTraction: {
      Marmot::FiniteElement::BoundaryElement boundaryEl( displacementGeometryElement.shape,
                                                         elementFace,
                                                         nDim,
                                                         displacementGeometryElement.coordinates );

      const XiSized         tractionVector( load );
      const Eigen::VectorXd Pk = boundaryEl.computeVectorialLoadVector( tractionVector );
      boundaryEl.assembleIntoParentVectorial( Pk, fU );
      break;
    }
    default: throw std::invalid_argument( "Invalid Load Type specified" );
    }
  }

  template < int nNodes, int nPressureNodes >
  void DisplacementPressureFiniteStrainElement< nNodes, nPressureNodes >::computeKernels( const double* QTotal_,
                                                                                          const double* dQ_,
                                                                                          double*       Pint_,
                                                                                          double*       K_,
                                                                                          double        time,
                                                                                          double        dT )
  {
    using namespace Fastor;
    using namespace Marmot;
    using namespace Marmot::FastorIndices;
    using namespace Marmot::FastorStandardTensors;

    const auto                             qU_np = TensorMap< const double, nNodes, nDim >( QTotal_ );
    const Eigen::Map< const PSizedVector > qP( QTotal_ + sizeDoFU );

    Eigen::Map< KeSizedMatrix >       Ke( K_ );
    TensorMap< double, nNodes, nDim > r_U( Pint_ );
    Eigen::Map< PSizedVector >        r_P( Pint_ + sizeDoFU );

    Eigen::Matrix< double, sizeDoFU, sizeDoFP > kUP_local = Eigen::Matrix< double, sizeDoFU, sizeDoFP >::Zero();

    Tensor< double, nDim, nNodes, nDim, nNodes > k_UU( 0.0 );

    for ( auto& qp : qps ) {

      const auto& dNdX_ = qp.dNdX;
      const auto  dNdX  = Tensor< double, nDim, nNodes >( dNdX_.data(), ColumnMajor );

      const auto F_np = evaluate( einsum< Ai, jA >( qU_np, dNdX ) + Spatial3D::I );

      const double p = ( qp.Np * qP )( 0, 0 );

      const auto [J, dJ_dF] = ContinuumMechanics::DeformationMeasures::FirstOrderDerived::volumeRatio( F_np );

      const Material::Deformation< 3 > deformation{ F_np };
      const Material::TimeIncrement    timeIncrement{ time, dT };

      Material::ConstitutiveResponse< 3 > response( Tensor33d( qp.managedStateVars->stress.data(), ColumnMajor ),
                                                    qp.managedStateVars->elasticEnergy / qp.J0xW,
                                                    qp.managedStateVars->dissipation / qp.J0xW,
                                                    qp.managedStateVars->materialStateVars.data() );
      Material::AlgorithmicModuli< 3 >    tangents;

      qp.material->computeStress( response, tangents, deformation, timeIncrement );

      qp.managedStateVars->stress            = Marmot::mapEigenToFastor( response.tau ).reshaped();
      qp.managedStateVars->elasticEnergy     = response.elasticEnergyDensity * qp.J0xW;
      qp.managedStateVars->dissipation       = response.dissipation * qp.J0xW;
      qp.managedStateVars->totalStrainEnergy = ( response.elasticEnergyDensity + response.dissipation ) * qp.J0xW;

      // total Kirchhoff stress: isochoric material response + volumetric pressure term p*J*I
      const Tensor33d tauTotal = response.tau + ( p * J ) * Spatial3D::I;

      Tensor3333d dTauTotal_dF = tangents.dTau_dF;
      for ( int i = 0; i < 3; ++i )
        for ( int j = 0; j < 3; ++j )
          for ( int k = 0; k < 3; ++k )
            for ( int K = 0; K < 3; ++K )
              dTauTotal_dF( i, j, k, K ) += p * Spatial3D::I( i, j ) * dJ_dF( k, K );

      const auto dNdx = evaluate( einsum< ji, jA >( inv( F_np ), dNdX ) );

      const double& J0xW = qp.J0xW;

      r_U += ( +einsum< iA, ij >( dNdx, tauTotal ) ) * J0xW;

      const auto dTau_dqU = evaluate( +einsum< ijkl, lB >( dTauTotal_dF, dNdX ) );
      k_UU += ( +einsum< iA, ijkB, to_jAkB >( dNdx, dTau_dqU ) ) * J0xW;
      k_UU += ( -einsum< kA, ij, iB, to_jAkB >( dNdx, tauTotal, dNdx ) ) * J0xW;

      r_P += qp.Np.transpose() * ( J - 1.0 ) * J0xW;

      for ( int A = 0; A < nNodes; ++A )
        for ( int j = 0; j < nDim; ++j )
          for ( int C = 0; C < nPressureNodes; ++C )
            kUP_local( A * nDim + j, C ) += dNdx( j, A ) * J * qp.Np( C ) * J0xW;
    }

    using namespace Eigen;

    Ke.template block< sizeDoFU, sizeDoFU >( 0, 0 ) += Map< Matrix< double, sizeDoFU, sizeDoFU > >(
      torowmajor( k_UU ).data() );
    Ke.template block< sizeDoFU, sizeDoFP >( 0, sizeDoFU ) += kUP_local;
    Ke.template block< sizeDoFP, sizeDoFU >( sizeDoFU, 0 ) += kUP_local.transpose();
  }

  template < int nNodes, int nPressureNodes >
  StateView DisplacementPressureFiniteStrainElement< nNodes, nPressureNodes >::getStateView(
    const std::string& stateName,
    int                qpNumber )
  {
    const auto& qp = qps[qpNumber];
    if ( qp.managedStateVars->contains( stateName ) )
      return qp.managedStateVars->getStateView( stateName );
    else
      return qp.material->getStateView( stateName, qp.managedStateVars->materialStateVars.data() );
  }

  template < int nNodes, int nPressureNodes >
  void DisplacementPressureFiniteStrainElement< nNodes, nPressureNodes >::computeInternalEnergy(
    double& internalEnergy )
  {
    internalEnergy = 0.0;
    for ( const auto& qp : qps )
      internalEnergy += qp.managedStateVars->totalStrainEnergy;
  }

  template < int nNodes, int nPressureNodes >
  std::vector< double > DisplacementPressureFiniteStrainElement< nNodes, nPressureNodes >::getCoordinatesAtCenter()
  {
    std::vector< double > coords( nDim );
    Eigen::Map< XiSized > coordsMap( &coords[0] );
    const auto            centerXi = XiSized::Zero();
    coordsMap                      = displacementGeometryElement.NB( displacementGeometryElement.N( centerXi ) ) *
                displacementGeometryElement.coordinates;
    return coords;
  }

  template < int nNodes, int nPressureNodes >
  std::vector< std::vector< double > > DisplacementPressureFiniteStrainElement< nNodes, nPressureNodes >::
    getCoordinatesAtQuadraturePoints()
  {
    std::vector< std::vector< double > > listedCoords;
    std::vector< double >                coords( nDim );
    Eigen::Map< XiSized >                coordsMap( &coords[0] );

    for ( const auto& qp : qps ) {
      coordsMap = displacementGeometryElement.NB( displacementGeometryElement.N( qp.xi ) ) *
                  displacementGeometryElement.coordinates;
      listedCoords.push_back( coords );
    }
    return listedCoords;
  }

  template < int nNodes, int nPressureNodes >
  int DisplacementPressureFiniteStrainElement< nNodes, nPressureNodes >::getNumberOfQuadraturePoints()
  {
    return qps.size();
  }

} // namespace Marmot::Elements
