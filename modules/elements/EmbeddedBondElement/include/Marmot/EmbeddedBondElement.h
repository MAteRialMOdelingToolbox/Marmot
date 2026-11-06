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
#include <string>
#include <vector>

namespace Marmot::Discretizations::Elements {

  /**
   * @class Marmot::Discretizations::Elements::EmbeddedBondElement
   * @brief Bond-slip coupling of an embedded reinforcement bar (truss) to the host continuum element it lies in.
   *
   * @tparam nDim        number of spatial dimensions, 2 or 3
   * @tparam nRebarNodes number of nodes of the reinforcement element, 2 or 3 (end, end, mid)
   * @tparam nHostNodes  number of nodes of the host element (Quad4, Quad8, Tetra4, Tetra10, Hexa8, Hexa20)
   *
   * The nodes of the element are the nodes of the reinforcement (truss) element, followed by the nodes of the host
   * element. The reinforcement keeps its own, independent displacement field; the element couples it to the host by
   * the bond stress, which a pluggable bond-slip law (Marmot::Materials::BondSlipLaw) computes from the slip, i.e., the
   * displacement of the bar relative to the host material point it lies at:
   * \f[
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
  template < int nDim, int nRebarNodes, int nHostNodes >
  class EmbeddedBondElement : public Marmot::Discretizations::Elements::Element {

    static_assert( nDim == 2 || nDim == 3, "EmbeddedBondElement: nDim must be 2 or 3" );
    static_assert( nRebarNodes == 2 || nRebarNodes == 3, "EmbeddedBondElement: nRebarNodes must be 2 or 3" );

  public:
    static constexpr int nNodes         = nRebarNodes + nHostNodes;
    static constexpr int sizeLoadVector = nNodes * nDim;

    using HostGeometry = Marmot::Discretizations::Elements::GeometryElement< nDim, nHostNodes >;
    using VectorDim    = Eigen::Matrix< double, nDim, 1 >;
    using MatrixDim    = Eigen::Matrix< double, nDim, nDim >;
    using RhsSized     = Eigen::Matrix< double, sizeLoadVector, 1 >;
    using KSized       = Eigen::Matrix< double, sizeLoadVector, sizeLoadVector >;
    using LSized       = Eigen::Matrix< double, nDim, sizeLoadVector >;
    using RebarN       = Eigen::Matrix< double, 1, nRebarNodes >;

    /// Element label
    const int elLabel;

    /// Element properties: [perimeter, eta start, eta end, (number of integration points)]
    Eigen::Map< const Eigen::VectorXd > elementProperties;

    /// Reference coordinates of the reinforcement nodes, [dim, node]
    Eigen::Map< const Eigen::Matrix< double, nDim, nRebarNodes > > rebarCoordinates;

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
      double                                            eta    = 0; ///< parent coordinate of the bar element
      double                                            weight = 0; ///< weight of the Gauss-Lobatto rule on [-1, 1]
      Eigen::Matrix< double, nDim, 1 >                  hostXi;     ///< parent coordinates in the host element
      LSized                                            L; ///< interpolation operator of the relative displacement
      MatrixDim                                         R; ///< local frame [t, n1, (n2)] as columns
      double                                            pdS = 0; ///< perimeter x reference length of the point
      std::unique_ptr< QPStateVarManager >              managedStateVars;
      std::unique_ptr< Marmot::Materials::BondSlipLaw > bondSlipLaw;
    };

    /// the integration points
    std::vector< QuadraturePoint > qps;

    explicit EmbeddedBondElement( int elementID )
      : elLabel( elementID ), elementProperties( nullptr, 0 ), rebarCoordinates( nullptr )
    {
      qps.resize( defaultNumberOfIntegrationPoints() );
    }

    static constexpr int defaultNumberOfIntegrationPoints() { return nRebarNodes == 2 ? 2 : 3; }

    /// Gauss-Lobatto points and weights on [-1, 1]
    static std::vector< std::pair< double, double > > gaussLobatto( int n )
    {
      switch ( n ) {
      case 2: return { { -1., 1. }, { 1., 1. } };
      case 3: return { { -1., 1. / 3 }, { 0., 4. / 3 }, { 1., 1. / 3 } };
      case 4: {
        const double a = std::sqrt( 1. / 5 );
        return { { -1., 1. / 6 }, { -a, 5. / 6 }, { a, 5. / 6 }, { 1., 1. / 6 } };
      }
      case 5: {
        const double a = std::sqrt( 3. / 7 );
        return { { -1., 0.1 }, { -a, 49. / 90 }, { 0., 32. / 45 }, { a, 49. / 90 }, { 1., 0.1 } };
      }
      default:
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << "EmbeddedBondElement: 2 to 5 integration points are supported, " << n
                                     << " were requested" );
      }
    }

    static RebarN rebarN( double eta )
    {
      if constexpr ( nRebarNodes == 2 )
        return Marmot::Discretizations::Elements::Spatial1D::Bar2::N( eta );
      else
        return Marmot::Discretizations::Elements::Spatial1D::Bar3::N( eta );
    }

    static RebarN rebarDNdEta( double eta )
    {
      if constexpr ( nRebarNodes == 2 )
        return Marmot::Discretizations::Elements::Spatial1D::Bar2::dNdXi( eta );
      else
        return Marmot::Discretizations::Elements::Spatial1D::Bar3::dNdXi( eta );
    }

    /// local frame of the unit tangent t, columns [t, n1, (n2)]
    static MatrixDim localFrame( const VectorDim& t )
    {
      MatrixDim R;
      if constexpr ( nDim == 2 ) {
        R.col( 0 ) = t;
        R.col( 1 ) = VectorDim( -t( 1 ), t( 0 ) );
      }
      else {
        // the global axis least aligned with t gives a well conditioned first normal
        int axis;
        t.cwiseAbs().minCoeff( &axis );
        const Eigen::Vector3d e  = Eigen::Vector3d::Unit( axis );
        const Eigen::Vector3d n1 = t.cross( e ).normalized();
        R.col( 0 )               = t;
        R.col( 1 )               = n1;
        R.col( 2 )               = t.cross( n1 );
      }
      return R;
    }

    /**
     * @brief Parent coordinates of a point in the host element (inverse isoparametric mapping).
     * @throws std::invalid_argument if the point lies outside of the host element
     */
    VectorDim findHostParentCoordinates( const VectorDim& X ) const
    {
      VectorDim xi = VectorDim::Zero();
      if ( hostGeometry.shape == Marmot::Discretizations::Elements::ElementShapes::Tetra4 ||
           hostGeometry.shape == Marmot::Discretizations::Elements::ElementShapes::Tetra10 )
        xi.setConstant( 0.25 );

      for ( int iteration = 0;; iteration++ ) {
        const VectorDim residual = hostGeometry.NB( hostGeometry.N( xi ) ) * hostGeometry.coordinates - X;
        const MatrixDim J        = hostGeometry.Jacobian( hostGeometry.dNdXi( xi ) );
        const VectorDim dXi      = -J.inverse() * residual;
        xi += dXi;
        if ( dXi.norm() < 1e-13 )
          break;
        if ( iteration > 50 )
          throw std::invalid_argument( Marmot::Utilities::MakeString()
                                       << "EmbeddedBondElement " << elLabel
                                       << ": inverse mapping into the host element did not converge" );
      }

      constexpr double tol = 1e-6;
      bool             inside;
      if ( hostGeometry.shape == Marmot::Discretizations::Elements::ElementShapes::Tetra4 ||
           hostGeometry.shape == Marmot::Discretizations::Elements::ElementShapes::Tetra10 )
        inside = xi.minCoeff() >= -tol && xi.sum() <= 1 + tol;
      else
        inside = xi.cwiseAbs().maxCoeff() <= 1 + tol;

      if ( !inside )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << "EmbeddedBondElement " << elLabel
                                     << ": a bond integration point lies outside of the host element, "
                                        "parent coordinates "
                                     << xi.transpose() );
      return xi;
    }

    int getNumberOfRequiredStateVars() override
    {
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

    std::string getElementShape() override { return nRebarNodes == 2 ? "bar2" : "bar3"; }

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

      const int nQps = property.nElementProperties > 3 ? static_cast< int >( std::lround( elementProperties[3] ) )
                                                       : defaultNumberOfIntegrationPoints();
      gaussLobatto( nQps ); // validates
      qps.clear();
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
      new ( &rebarCoordinates ) Eigen::Map< const Eigen::Matrix< double, nDim, nRebarNodes > >( coordinates );
      hostGeometry.assignNodeCoordinates( coordinates + nRebarNodes * nDim );
    }

    void initializeYourself() override
    {
      const double perimeter = elementProperties[0];
      const double etaStart  = elementProperties[1];
      const double etaEnd    = elementProperties[2];
      const auto   rule      = gaussLobatto( static_cast< int >( qps.size() ) );

      for ( size_t q = 0; q < qps.size(); q++ ) {
        auto& qp  = qps[q];
        qp.eta    = etaStart + 0.5 * ( etaEnd - etaStart ) * ( rule[q].first + 1 );
        qp.weight = rule[q].second;

        const RebarN    M  = rebarN( qp.eta );
        const VectorDim X  = rebarCoordinates * M.transpose();
        const VectorDim G  = rebarCoordinates * rebarDNdEta( qp.eta ).transpose();
        const double    dS = G.norm() * 0.5 * std::abs( etaEnd - etaStart ) * qp.weight;

        qp.pdS    = perimeter * dS;
        qp.R      = localFrame( G / G.norm() );
        qp.hostXi = findHostParentCoordinates( X );

        const auto N = hostGeometry.N( qp.hostXi );
        qp.L.setZero();
        for ( int i = 0; i < nRebarNodes; i++ )
          qp.L.template block< nDim, nDim >( 0, i * nDim ) = M( i ) * MatrixDim::Identity();
        for ( int a = 0; a < nHostNodes; a++ )
          qp.L.template block< nDim, nDim >( 0, ( nRebarNodes + a ) * nDim ) = -N( a ) * MatrixDim::Identity();
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
                                                     sv.elasticEnergy / qp.pdS,
                                                     sv.dissipation / qp.pdS,
                                                     sv.materialStateVars.data() };
        Eigen::Matrix3d                       C;
        qp.bondSlipLaw->computeBondStress( state, C, slip, dSlip, { time, dT } );

        sv.slip          = slip;
        sv.bondStress    = state.bondStress;
        sv.elasticEnergy = state.elasticEnergyDensity * qp.pdS;
        sv.dissipation   = state.dissipation * qp.pdS;

        const LSized RtL = qp.R.transpose() * qp.L;
        Pe += RtL.transpose() * state.bondStress.template head< nDim >() * qp.pdS;
        Ke += RtL.transpose() * C.template topLeftCorner< nDim, nDim >() * RtL * qp.pdS;
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
      const VectorDim X = rebarCoordinates *
                          rebarN( 0.5 * ( elementProperties[1] + elementProperties[2] ) ).transpose();
      return std::vector< double >( X.data(), X.data() + nDim );
    }

    std::vector< std::vector< double > > getCoordinatesAtQuadraturePoints() override
    {
      std::vector< std::vector< double > > result;
      for ( const auto& qp : qps ) {
        const VectorDim X = rebarCoordinates * rebarN( qp.eta ).transpose();
        result.emplace_back( X.data(), X.data() + nDim );
      }
      return result;
    }

    int getNumberOfQuadraturePoints() override { return static_cast< int >( qps.size() ); }
  };

} // namespace Marmot::Discretizations::Elements
