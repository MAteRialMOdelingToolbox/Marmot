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
#include "Marmot/EmbeddedBondKinematics.h"
#include "Marmot/MarmotBondSlipLaw.h"
#include "Marmot/MarmotBondSlipLawFactory.h"
#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryElement.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace Marmot::Discretizations::Elements {

  /**
   * @class Marmot::Discretizations::Elements::EmbeddedBondElement
   * @brief Bond-slip coupling of an embedded reinforcement bar (truss) to the host continuum element it lies in.
   *
   * @tparam nDim        number of spatial dimensions, 2 or 3
   * @tparam nBarNodes number of nodes of the bar element, 2 or 3 (end, end, mid)
   * @tparam nHostNodes  number of nodes of the host element (Tri3, Quad4, Quad8, Tetra4, Hexa8, Hexa20; the hosts with
   * a continuum element in Marmot)
   *
   * The nodes of the element are the nodes of the reinforcement (truss) element, followed by the nodes of the host
   * element (bar first: the element is drawn as its bar; the large-slip element puts the host first, to be drawn as the
   * host on which its channel lies). The reinforcement keeps its own, independent displacement field; the element
   * couples it to the host by the bond stress, which a pluggable bond-slip law (Marmot::Materials::BondSlipLaw)
   * computes from the slip, i.e., the displacement of the bar relative to the host material point it lies at: \f[
   *   \boldsymbol{\Delta}(\eta) = \sum_i M_i(\eta)\, \boldsymbol{u}_i - \sum_a N_a(\boldsymbol{\xi}(\eta))\,
   *   \boldsymbol{u}_a , \qquad \boldsymbol{s} = \boldsymbol{R}^T \boldsymbol{\Delta},
   * \f]
   * with the shape functions \f$M_i\f$ of the bar, \f$N_a\f$ of the host, the host parent coordinates
   * \f$\boldsymbol{\xi}(\eta)\f$ of the bar point \f$\eta\f$ (found by inverse mapping at initialization), and the
   * local frame \f$\boldsymbol{R} = [\boldsymbol{t}, \boldsymbol{n}_1, \boldsymbol{n}_2]\f$ of the bar axis. The
   * internal forces and the tangent are
   * \f[
   *   \boldsymbol{f} = \int \boldsymbol{L}^T \boldsymbol{R}\, \boldsymbol{\tau}\; p\, \mathrm{d}S, \qquad
   *   \boldsymbol{K} = \int \boldsymbol{L}^T \boldsymbol{R}\, \frac{\partial \boldsymbol{\tau}}{\partial
   *   \boldsymbol{s}}\, \boldsymbol{R}^T \boldsymbol{L}\; p\, \mathrm{d}S
   * \f]
   * with the interpolation operator \f$\boldsymbol{L}\f$ of \f$\boldsymbol{\Delta}\f$ and the bar perimeter \f$p\f$.
   *
   * One element covers the part \f$[\eta_s, \eta_e]\f$ of the bar element that lies inside the host element; a bar
   * element crossing several host elements is coupled by several bond elements. A preprocessor (e.g., the embedded
   * reinforcement generator of EdelweissFE) computes these parts.
   *
   * The formulation is geometrically linear: the frame is the one of the reference configuration, which is
   * appropriate for small rotations of the bar.
   *
   * Element properties: [perimeter \f$p\f$, \f$\eta_s\f$, \f$\eta_e\f$, (number of integration points)]. The bond
   * stress is integrated with a Gauss-Lobatto rule (2 to 5 points; default: 2 for linear and 3 for quadratic bars),
   * which places integration points at the ends of the bar part and avoids the traction oscillations of Gauss rules
   * for stiff bond.
   *
   * Material: a bond-slip law registered in the Marmot::Registry::BondSlipLawFactory.
   */
  template < int nDim, int nBarNodes, int nHostNodes >
  class EmbeddedBondElement : public Marmot::Discretizations::Elements::Element {

    static_assert( nDim == 2 || nDim == 3, "EmbeddedBondElement: nDim must be 2 or 3" );
    static_assert( nBarNodes == 2 || nBarNodes == 3, "EmbeddedBondElement: nBarNodes must be 2 or 3" );

  public:
    /// number of nodes: the bar nodes, followed by the host nodes
    static constexpr int nNodes = nBarNodes + nHostNodes;
    /// number of element dofs (displacements only)
    static constexpr int sizeLoadVector = nNodes * nDim;

    using HostGeometry = Marmot::Discretizations::Elements::GeometryElement< nDim, nHostNodes >; ///< host shape
    using VectorDim    = Eigen::Matrix< double, nDim, 1 >;                                       ///< spatial vector
    using MatrixDim    = Eigen::Matrix< double, nDim, nDim >;                                    ///< spatial tensor
    using RhsSized     = Eigen::Matrix< double, sizeLoadVector, 1 >;              ///< element force vector
    using KSized       = Eigen::Matrix< double, sizeLoadVector, sizeLoadVector >; ///< element stiffness matrix
    using LSized       = Eigen::Matrix< double, nDim, sizeLoadVector >;           ///< relative displacement operator
    using BarN         = Eigen::Matrix< double, 1, nBarNodes >;                   ///< bar shape functions

    /// Element label
    const int elLabel;

    /// Element properties: [perimeter, eta start, eta end, (number of integration points)]
    Eigen::Map< const Eigen::VectorXd > elementProperties;

    /// Reference coordinates of the reinforcement nodes, [dim, node]
    Eigen::Map< const Eigen::Matrix< double, nDim, nBarNodes > > barCoordinates;

    /// The geometry of the host element
    HostGeometry hostGeometry;

    /// The state of an integration point
    class QPStateVarManager : public Marmot::Utilities::StateVarVectorManager {

      /// \hideinitializer
      inline const static auto layout = makeLayout( {
        { .name = "slip", .length = 3 },
        { .name = "bond stress", .length = 3 },
        { .name = "elastic energy", .length = 1 },
        { .name = "dissipation", .length = 1 },
        { .name = "begin of material state", .length = 0 },
      } );

    public:
      Eigen::Map< Eigen::Vector3d > slip;
      Eigen::Map< Eigen::Vector3d > bondStress;
      double&                       elasticEnergy;
      double&                       dissipation;
      Eigen::Map< Eigen::VectorXd > materialStateVars;

      static int getNumberOfRequiredStateVarsQuadraturePointOnly() { return layout.nRequiredStateVars; };

      QPStateVarManager( double* theStateVarVector, int nStateVars )
        : Marmot::Utilities::StateVarVectorManager( theStateVarVector, layout ),
          slip( &find( "slip" ) ),
          bondStress( &find( "bond stress" ) ),
          elasticEnergy( find( "elastic energy" ) ),
          dissipation( find( "dissipation" ) ),
          materialStateVars( &find( "begin of material state" ),
                             nStateVars - getNumberOfRequiredStateVarsQuadraturePointOnly() ){};
    };

    /// An integration point along the bar
    struct QuadraturePoint {
      double                           eta    = 0; ///< parent coordinate of the bar element
      double                           weight = 0; ///< weight of the Gauss-Lobatto rule on [-1, 1]
      Eigen::Matrix< double, nDim, 1 > hostXi;     ///< parent coordinates in the host element
      LSized                           L;          ///< interpolation operator of the relative displacement
      MatrixDim                        R;          ///< local frame [t, n1, (n2)] as columns
      double perimeterTimesDS = 0;                 ///< perimeter x reference length of the point (integration weight)
      std::unique_ptr< QPStateVarManager >              managedStateVars;
      std::unique_ptr< Marmot::Materials::BondSlipLaw > bondSlipLaw;
    };

    /// the integration points
    std::vector< QuadraturePoint > qps;

    /**
     * @brief Construct the element with the default number of integration points.
     * @param elementID the element label
     */
    explicit EmbeddedBondElement( int elementID )
      : elLabel( elementID ), elementProperties( nullptr, 0 ), barCoordinates( nullptr )
    {
      qps.resize( defaultNumberOfIntegrationPoints() );
    }

    /// 2 Gauss-Lobatto points for linear bars, 3 for quadratic bars (exact for the stiffness of a linear law)
    static constexpr int defaultNumberOfIntegrationPoints() { return nBarNodes == 2 ? 2 : 3; }

    /// shape functions of the bar element at the bar parent coordinate @p eta
    static BarN barN( double eta )
    {
      if constexpr ( nBarNodes == 2 )
        return Marmot::Discretizations::Elements::Spatial1D::Bar2::N( eta );
      else
        return Marmot::Discretizations::Elements::Spatial1D::Bar3::N( eta );
    }

    /// derivatives of the bar shape functions with respect to @p eta
    static BarN barDNdEta( double eta )
    {
      if constexpr ( nBarNodes == 2 )
        return Marmot::Discretizations::Elements::Spatial1D::Bar2::dNdXi( eta );
      else
        return Marmot::Discretizations::Elements::Spatial1D::Bar3::dNdXi( eta );
    }

    /// The tolerance of "inside the host element" in parent coordinates: looser than the one EdelweissFE uses to
    /// place the ends of the bond parts on host boundaries, so that such points are always accepted
    static constexpr double insideTolerance = 1e-6;

    /**
     * @brief Parent coordinates of a point in the host element: the inverse isoparametric mapping of the host
     * geometry (GeometryElement::parentCoordinates).
     * @param[in] X The point.
     * @return The parent coordinates.
     * @throws std::invalid_argument if the point lies outside of the host element
     */
    VectorDim findHostParentCoordinates( const VectorDim& X ) const
    {
      VectorDim xi;
      try {
        xi = hostGeometry.parentCoordinates( X );
      }
      catch ( const std::runtime_error& ) {
        // no convergence: only possible far outside of the element
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << "EmbeddedBondElement " << elLabel
                                     << ": a bond integration point lies outside of the host element, "
                                     << X.transpose() );
      }
      if ( !HostGeometry::isInsideParentDomain( xi, insideTolerance ) )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << "EmbeddedBondElement " << elLabel
                                     << ": a bond integration point lies outside of the host element, "
                                        "parent coordinates "
                                     << xi.transpose() );
      return xi;
    }

    int getNumberOfRequiredStateVars() override
    {
      if ( !qps[0].bondSlipLaw )
        throw std::logic_error( Marmot::Utilities::MakeString()
                                << "EmbeddedBondElement " << elLabel
                                << ": the bond-slip law must be assigned before the state variables" );
      return static_cast< int >( qps.size() ) * ( QPStateVarManager::getNumberOfRequiredStateVarsQuadraturePointOnly() +
                                                  qps[0].bondSlipLaw->getNumberOfRequiredStateVars() );
    }

    std::vector< std::vector< std::string > > getNodeFields() override
    {
      return std::vector< std::vector< std::string > >( nNodes, { "displacement" } );
    }

    std::vector< int > getDofIndicesPermutationPattern() override
    {
      std::vector< int > pattern( sizeLoadVector );
      for ( int i = 0; i < sizeLoadVector; i++ )
        pattern[i] = i;
      return pattern;
    }

    int getNNodes() override { return nNodes; }

    int getNSpatialDimensions() override { return nDim; }

    int getNDofPerElement() override { return sizeLoadVector; }

    std::string getElementShape() override { return nBarNodes == 2 ? "bar2" : "bar3"; }

    void assignStateVars( double* stateVars, int nStateVars ) override
    {
      const int nQpStateVars = nStateVars / static_cast< int >( qps.size() );
      for ( size_t i = 0; i < qps.size(); i++ )
        qps[i].managedStateVars = std::make_unique< QPStateVarManager >( stateVars + i * nQpStateVars, nQpStateVars );
    }

    void assignProperty( const ElementProperties& property ) override
    {
      if ( property.nElementProperties < 3 )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__
                                     << ": requires the properties [perimeter, eta start, eta end, "
                                        "(number of integration points)]" );
      new ( &elementProperties )
        Eigen::Map< const Eigen::VectorXd >( property.elementProperties, property.nElementProperties );

      const double perimeter = elementProperties[0], etaStart = elementProperties[1], etaEnd = elementProperties[2];
      if ( !( perimeter > 0 ) )
        throw std::invalid_argument( Marmot::Utilities::MakeString() << "EmbeddedBondElement " << elLabel
                                                                     << ": requires perimeter > 0, got " << perimeter );
      if ( !( -1 <= etaStart && etaStart < etaEnd && etaEnd <= 1 ) )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << "EmbeddedBondElement " << elLabel
                                     << ": requires -1 <= eta start < eta end <= 1, got [" << etaStart << ", " << etaEnd
                                     << "]" );

      const int nQps = property.nElementProperties > 3 ? static_cast< int >( std::lround( elementProperties[3] ) )
                                                       : defaultNumberOfIntegrationPoints();
      EmbeddedBondKinematics::gaussLobatto( nQps ); // validates
      // resizing the integration points would destroy the bond-slip laws assigned with them
      if ( qps[0].bondSlipLaw && nQps != static_cast< int >( qps.size() ) )
        throw std::logic_error( Marmot::Utilities::MakeString()
                                << "EmbeddedBondElement " << elLabel
                                << ": the element properties (number of integration points) must be assigned before "
                                   "the bond-slip law" );
      qps.resize( nQps );
    }

    void assignProperty( const Marmot::Discretizations::MaterialSection& section ) override
    {
      for ( auto& qp : qps )
        qp.bondSlipLaw = std::unique_ptr< Marmot::Materials::BondSlipLaw >(
          Marmot::Registry::BondSlipLawFactory::createBondSlipLaw( section.materialName,
                                                                   section.materialProperties,
                                                                   section.nMaterialProperties,
                                                                   elLabel ) );
    }

    void assignNodeCoordinates( const double* coordinates ) override
    {
      new ( &barCoordinates ) Eigen::Map< const Eigen::Matrix< double, nDim, nBarNodes > >( coordinates );
      hostGeometry.assignNodeCoordinates( coordinates + nBarNodes * nDim );
    }

    void initializeYourself() override
    {
      const double perimeter = elementProperties[0];
      const double etaStart  = elementProperties[1];
      const double etaEnd    = elementProperties[2];
      const auto&  rule      = EmbeddedBondKinematics::gaussLobatto( static_cast< int >( qps.size() ) );

      for ( size_t q = 0; q < qps.size(); q++ ) {
        auto& qp  = qps[q];
        qp.eta    = etaStart + 0.5 * ( etaEnd - etaStart ) * ( rule[q].xi( 0 ) + 1 );
        qp.weight = rule[q].weight;

        const BarN      M       = barN( qp.eta );
        const VectorDim X       = barCoordinates * M.transpose();
        const VectorDim dX_dEta = barCoordinates * barDNdEta( qp.eta ).transpose();
        const double    dS      = dX_dEta.norm() * 0.5 * ( etaEnd - etaStart ) * qp.weight;

        qp.perimeterTimesDS = perimeter * dS;
        qp.R                = EmbeddedBondKinematics::localFrame< nDim >( dX_dEta / dX_dEta.norm() );
        qp.hostXi           = findHostParentCoordinates( X );

        const auto N = hostGeometry.N( qp.hostXi );
        qp.L.setZero();
        for ( int i = 0; i < nBarNodes; i++ )
          qp.L.template block< nDim, nDim >( 0, i * nDim ) = M( i ) * MatrixDim::Identity();
        for ( int a = 0; a < nHostNodes; a++ )
          qp.L.template block< nDim, nDim >( 0, ( nBarNodes + a ) * nDim ) = -N( a ) * MatrixDim::Identity();
      }
    }

    void setInitialConditions( StateTypes state, const double* ) override
    {
      switch ( state ) {
      case Marmot::Discretizations::Elements::Element::MarmotMaterialInitialization: {
        for ( auto& qp : qps )
          qp.bondSlipLaw->initializeYourself( qp.managedStateVars->materialStateVars.data(),
                                              static_cast< int >( qp.managedStateVars->materialStateVars.size() ) );
        break;
      }
      default:
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__ << ": invalid initial condition" );
      }
    }

    void computeKernels( const double* QTotal_, const double* dQ_, double* Pe_, double* Ke_, double time, double dT )
      override
    {
      Eigen::Map< const RhsSized > QTotal( QTotal_ );
      Eigen::Map< const RhsSized > dQ( dQ_ );
      Eigen::Map< RhsSized >       Pe( Pe_ );
      Eigen::Map< KSized >         Ke( Ke_ );

      for ( auto& qp : qps ) {
        auto& sv = *qp.managedStateVars;

        Eigen::Vector3d slip = Eigen::Vector3d::Zero(), dSlip = Eigen::Vector3d::Zero();
        slip.template head< nDim >()  = qp.R.transpose() * ( qp.L * QTotal );
        dSlip.template head< nDim >() = qp.R.transpose() * ( qp.L * dQ );

        Marmot::Materials::BondSlipLaw::State state{ sv.bondStress,
                                                     sv.elasticEnergy / qp.perimeterTimesDS,
                                                     sv.dissipation / qp.perimeterTimesDS,
                                                     sv.materialStateVars.data() };
        Eigen::Matrix3d                       dBondStress_dSlip;
        qp.bondSlipLaw->computeBondStress( state, dBondStress_dSlip, slip, dSlip, { time, dT } );

        sv.slip          = slip;
        sv.bondStress    = state.bondStress;
        sv.elasticEnergy = state.elasticEnergyDensity * qp.perimeterTimesDS;
        sv.dissipation   = state.dissipation * qp.perimeterTimesDS;

        const LSized dSlip_dU = qp.R.transpose() * qp.L;
        Pe += dSlip_dU.transpose() * state.bondStress.template head< nDim >() * qp.perimeterTimesDS;
        Ke += dSlip_dU.transpose() * dBondStress_dSlip.template topLeftCorner< nDim, nDim >() * dSlip_dU *
              qp.perimeterTimesDS;
      }
    }

    void computeKernelsExplicit( const double* QTotal, const double* dQ, double* Pe, double time, double dT ) override
    {
      KSized K = KSized::Zero();
      computeKernels( QTotal, dQ, Pe, K.data(), time, dT );
    }

    void computeDistributedLoad( DistributedLoadTypes,
                                 double*,
                                 double*,
                                 int,
                                 const double*,
                                 const double*,
                                 double,
                                 double ) override
    {
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__ << ": distributed loads are not supported by bond elements" );
    }

    /// the bond has no volume: no body forces
    void computeBodyForce( double*, double*, const double*, const double*, double, double ) override {}

    /// the bond has no mass
    void computeLumpedInertia( double* ) override {}

    /// the bond has no mass
    void computeConsistentInertia( double* ) override {}

    void computeInternalEnergy( double& internalEnergy ) override
    {
      for ( const auto& qp : qps )
        internalEnergy += qp.managedStateVars->elasticEnergy;
    }

    Marmot::Utilities::StateView getStateView( const std::string& stateName, int qpNumber ) override
    {
      const auto& qp = qps[qpNumber];
      if ( qp.managedStateVars->contains( stateName ) )
        return qp.managedStateVars->getStateView( stateName );
      return qp.bondSlipLaw->getStateView( stateName, qp.managedStateVars->materialStateVars.data() );
    }

    std::vector< double > getCoordinatesAtCenter() override
    {
      const VectorDim X = barCoordinates * barN( 0.5 * ( elementProperties[1] + elementProperties[2] ) ).transpose();
      return std::vector< double >( X.data(), X.data() + nDim );
    }

    std::vector< std::vector< double > > getCoordinatesAtQuadraturePoints() override
    {
      std::vector< std::vector< double > > result;
      for ( const auto& qp : qps ) {
        const VectorDim X = barCoordinates * barN( qp.eta ).transpose();
        result.emplace_back( X.data(), X.data() + nDim );
      }
      return result;
    }

    int getNumberOfQuadraturePoints() override { return static_cast< int >( qps.size() ); }
  };

} // namespace Marmot::Discretizations::Elements
