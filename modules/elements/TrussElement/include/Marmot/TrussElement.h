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
#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include <Eigen/Dense>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

namespace Marmot::Discretizations::Elements {

  /**
   * @class Marmot::Discretizations::Elements::TrussElement
   * @brief Truss (bar) element in two or three spatial dimensions, for small and finite strains.
   *
   * @tparam nDim   number of spatial dimensions, 2 or 3
   * @tparam nNodes number of nodes, 2 (linear) or 3 (quadratic, node order: end, end, mid)
   *
   * The truss carries only an axial force. The material sees a uniaxial stress state along the truss axis, which is
   * found by the uniaxial stress reductions of the material base classes, so any 3D material of Marmot can be used:
   *
   * - @b SmallStrain: hypoelastic materials (Marmot::Materials::MaterialHypoElastic). The axial strain is
   *   \f$\varepsilon = \boldsymbol{T} \cdot \frac{\partial \boldsymbol{u}}{\partial S}\f$ with the reference axis
   *   \f$\boldsymbol{T}\f$ and the reference arc length \f$S\f$ (geometrically linear).
   * - @b FiniteStrain: finite strain materials (Marmot::Materials::MaterialFiniteStrain). The material is evaluated in
   * the co-rotated frame of the truss, \f$\boldsymbol{F} = \mathrm{diag}( \lambda, \lambda_2, \lambda_3 )\f$ with the
   *   axial stretch \f$\lambda = |\boldsymbol{g}| / |\boldsymbol{G}|\f$ of the current and the reference tangent
   *   vectors, and the lateral stretches \f$\lambda_2, \lambda_3\f$ found such that the lateral stresses vanish. Rigid
   *   body rotations are therefore exactly stress free. The internal forces follow from the nominal axial stress
   *   \f$P = \tau_{11} / \lambda\f$ and the current axis \f$\boldsymbol{n}\f$,
   *   \f$\boldsymbol{f}_a = \int P \, \boldsymbol{n} \frac{\partial N_a}{\partial S} A_0 \, \mathrm{d}S\f$, with the
   *   consistent tangent
   *   \f$\boldsymbol{K}_{ab} = \int \left[ \frac{\mathrm{d} P}{\mathrm{d}\lambda} \boldsymbol{n} \otimes
   *   \boldsymbol{n} + \frac{P}{\lambda} ( \boldsymbol{I} - \boldsymbol{n} \otimes \boldsymbol{n} ) \right]
   *   \frac{\partial N_a}{\partial S} \frac{\partial N_b}{\partial S} A_0 \, \mathrm{d}S\f$.
   *
   * Element property: the (reference) cross section area \f$A_0\f$.
   *
   * Quadrature point states: @b stress (axial Cauchy stress), @b strain (small strain: axial strain; finite strain:
   * logarithmic axial strain \f$\ln\lambda\f$), @b normal force, @b elastic energy, @b dissipation (both integrated
   * over the quadrature point's volume), @b lateral stretches and @b kirchhoff stress (finite strain only), followed
   * by the material state.
   */
  template < int nDim, int nNodes >
  class TrussElement : public Marmot::Discretizations::Elements::Element {

    static_assert( nDim == 2 || nDim == 3, "TrussElement: nDim must be 2 or 3" );
    static_assert( nNodes == 2 || nNodes == 3, "TrussElement: nNodes must be 2 or 3" );

  public:
    /// Kinematic formulation
    enum Kinematics {
      SmallStrain,
      FiniteStrain,
    };

    static constexpr int sizeLoadVector = nNodes * nDim;

    using NSized          = Eigen::Matrix< double, 1, nNodes >;
    using CoordinateArray = Eigen::Matrix< double, nDim, nNodes >;
    using VectorDim       = Eigen::Matrix< double, nDim, 1 >;
    using MatrixDim       = Eigen::Matrix< double, nDim, nDim >;
    using RhsSized        = Eigen::Matrix< double, sizeLoadVector, 1 >;
    using KSized          = Eigen::Matrix< double, sizeLoadVector, sizeLoadVector >;
    using BSized          = Eigen::Matrix< double, 1, sizeLoadVector >;

    /// Element label
    const int elLabel;
    /// Kinematic formulation of this instance
    const Kinematics kinematics;

    /// Element properties: [cross section area]
    Eigen::Map< const Eigen::VectorXd > elementProperties;
    /// Reference nodal coordinates, [dim, node]
    Eigen::Map< const CoordinateArray > coordinates;

    /// The state of a quadrature point, see the class documentation for the layout
    class QPStateVarManager : public Marmot::Utilities::StateVarVectorManager {

      /// \hideinitializer
      inline const static auto layout = makeLayout( {
        { .name = "stress", .length = 1 },
        { .name = "strain", .length = 1 },
        { .name = "normal force", .length = 1 },
        { .name = "elastic energy", .length = 1 },
        { .name = "dissipation", .length = 1 },
        { .name = "lateral stretches", .length = 2 },
        { .name = "kirchhoff stress", .length = 9 },
        { .name = "begin of material state", .length = 0 },
      } );

    public:
      double&                       stress;
      double&                       strain;
      double&                       normalForce;
      double&                       elasticEnergy;
      double&                       dissipation;
      Eigen::Map< Eigen::Vector2d > lateralStretches;
      Eigen::Map< Eigen::Matrix3d > kirchhoffStress;
      Eigen::Map< Eigen::VectorXd > materialStateVars;

      static int getNumberOfRequiredStateVarsQuadraturePointOnly() { return layout.nRequiredStateVars; };

      QPStateVarManager( double* theStateVarVector, int nStateVars )
        : Marmot::Utilities::StateVarVectorManager( theStateVarVector, layout ),
          stress( find( "stress" ) ),
          strain( find( "strain" ) ),
          normalForce( find( "normal force" ) ),
          elasticEnergy( find( "elastic energy" ) ),
          dissipation( find( "dissipation" ) ),
          lateralStretches( &find( "lateral stretches" ) ),
          kirchhoffStress( &find( "kirchhoff stress" ) ),
          materialStateVars( &find( "begin of material state" ),
                             nStateVars - getNumberOfRequiredStateVarsQuadraturePointOnly() ){};
    };

    /// A quadrature point along the truss axis
    struct QuadraturePoint {
      const double xi;       ///< parent coordinate
      const double weight;   ///< quadrature weight

      NSized    N;           ///< shape functions
      NSized    dNdS;        ///< derivatives of the shape functions w.r.t. the reference arc length
      VectorDim T;           ///< reference unit tangent
      double    dSdXi = 0.0; ///< reference length per parent length
      double    A0xW  = 0.0; ///< cross section x dS/dxi x weight, i.e., the reference volume of the point
      BSized    B;           ///< small strain axial strain operator

      std::unique_ptr< QPStateVarManager >                       managedStateVars;
      std::unique_ptr< Marmot::Materials::MaterialHypoElastic >  materialHypoElastic;
      std::unique_ptr< Marmot::Materials::MaterialFiniteStrain > materialFiniteStrain;

      QuadraturePoint( double xi, double weight ) : xi( xi ), weight( weight ) {}

      int getNumberOfRequiredMaterialStateVars() const
      {
        return materialHypoElastic ? materialHypoElastic->getNumberOfRequiredStateVars()
                                   : materialFiniteStrain->getNumberOfRequiredStateVars();
      }

      double getDensity() const
      {
        const double* sv = managedStateVars->materialStateVars.data();
        return materialHypoElastic ? materialHypoElastic->getDensity( sv ) : materialFiniteStrain->getDensity( sv );
      }
    };

    /// the quadrature points
    std::vector< QuadraturePoint > qps;

    /**
     * @brief Construct a truss element.
     * @param elementID  element label
     * @param kinematics small or finite strain formulation
     */
    TrussElement( int elementID, Kinematics kinematics )
      : elLabel( elementID ), kinematics( kinematics ), elementProperties( nullptr, 0 ), coordinates( nullptr )
    {
      // exact for straight trusses of constant cross section and linear elastic material
      for ( const auto& qpInfo : Marmot::Discretizations::Elements::Quadrature::
              getGaussPointInfo( shape(),
                                 Marmot::Discretizations::Elements::Quadrature::IntegrationTypes::ReducedIntegration ) )
        qps.emplace_back( qpInfo.xi( 0 ), qpInfo.weight );
    }

    static constexpr Marmot::Discretizations::Elements::ElementShapes shape()
    {
      return nNodes == 2 ? Marmot::Discretizations::Elements::ElementShapes::Bar2
                         : Marmot::Discretizations::Elements::ElementShapes::Bar3;
    }

    static NSized N( double xi )
    {
      if constexpr ( nNodes == 2 )
        return Marmot::Discretizations::Elements::Spatial1D::Bar2::N( xi );
      else
        return Marmot::Discretizations::Elements::Spatial1D::Bar3::N( xi );
    }

    static NSized dNdXi( double xi )
    {
      if constexpr ( nNodes == 2 )
        return Marmot::Discretizations::Elements::Spatial1D::Bar2::dNdXi( xi );
      else
        return Marmot::Discretizations::Elements::Spatial1D::Bar3::dNdXi( xi );
    }

    int getNumberOfRequiredStateVars() override
    {
      return static_cast< int >( qps.size() ) * ( QPStateVarManager::getNumberOfRequiredStateVarsQuadraturePointOnly() +
                                                  qps[0].getNumberOfRequiredMaterialStateVars() );
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

    std::string getElementShape() override { return nNodes == 2 ? "bar2" : "bar3"; }

    void assignStateVars( double* stateVars, int nStateVars ) override
    {
      const int nQpStateVars = nStateVars / static_cast< int >( qps.size() );
      for ( size_t i = 0; i < qps.size(); i++ )
        qps[i].managedStateVars = std::make_unique< QPStateVarManager >( stateVars + i * nQpStateVars, nQpStateVars );
    }

    void assignProperty( const ElementProperties& property ) override
    {
      if ( property.nElementProperties < 1 )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__
                                     << ": a truss element requires its cross section area as property" );
      new ( &elementProperties )
        Eigen::Map< const Eigen::VectorXd >( property.elementProperties, property.nElementProperties );
    }

    void assignProperty( const Marmot::Discretizations::MaterialSection& section ) override
    {
      for ( auto& qp : qps ) {
        if ( kinematics == SmallStrain ) {
          qp.materialHypoElastic = std::unique_ptr< Marmot::Materials::MaterialHypoElastic >(
            Marmot::Registry::MaterialHypoElasticFactory::createMaterial( section.materialName,
                                                                          section.materialProperties,
                                                                          section.nMaterialProperties,
                                                                          elLabel ) );
          qp.materialHypoElastic->setCharacteristicElementLength( 2. * qp.dSdXi );
        }
        else
          qp.materialFiniteStrain = std::unique_ptr< Marmot::Materials::MaterialFiniteStrain >(
            Marmot::Registry::MaterialFiniteStrainFactory::createMaterial( section.materialName,
                                                                           section.materialProperties,
                                                                           section.nMaterialProperties,
                                                                           elLabel ) );
      }
    }

    void assignNodeCoordinates( const double* coords ) override
    {
      new ( &coordinates ) Eigen::Map< const CoordinateArray >( coords );
    }

    void initializeYourself() override
    {
      for ( auto& qp : qps ) {
        qp.N                   = N( qp.xi );
        const NSized    dNdXi_ = dNdXi( qp.xi );
        const VectorDim G      = coordinates * dNdXi_.transpose();
        qp.dSdXi               = G.norm();
        qp.T                   = G / qp.dSdXi;
        qp.dNdS                = dNdXi_ / qp.dSdXi;
        qp.A0xW                = elementProperties[0] * qp.dSdXi * qp.weight;
        for ( int a = 0; a < nNodes; a++ )
          qp.B.template segment< nDim >( a * nDim ) = qp.dNdS( a ) * qp.T.transpose();
      }
    }

    void setInitialConditions( StateTypes state, const double* values ) override
    {
      switch ( state ) {
      case Marmot::Discretizations::Elements::Element::MarmotMaterialInitialization: {
        for ( auto& qp : qps ) {
          auto& sv = *qp.managedStateVars;
          sv.lateralStretches.setOnes();
          if ( qp.materialHypoElastic )
            qp.materialHypoElastic->initializeYourself( sv.materialStateVars.data(), sv.materialStateVars.size() );
          else
            qp.materialFiniteStrain->initializeYourself( sv.materialStateVars.data(), sv.materialStateVars.size() );
        }
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

        if ( kinematics == SmallStrain ) {
          Marmot::Materials::MaterialHypoElastic::state1D state;
          state.stress               = sv.stress;
          state.elasticEnergyDensity = sv.elasticEnergy / qp.A0xW;
          state.dissipation          = sv.dissipation / qp.A0xW;
          state.stateVars            = sv.materialStateVars.data();

          double C = 0;
          qp.materialHypoElastic->computeUniaxialStress( state, C, qp.B * dQ, { time, dT } );

          sv.stress        = state.stress;
          sv.strain        = qp.B * QTotal;
          sv.normalForce   = state.stress * elementProperties[0];
          sv.elasticEnergy = state.elasticEnergyDensity * qp.A0xW;
          sv.dissipation   = state.dissipation * qp.A0xW;

          Pe += qp.B.transpose() * state.stress * qp.A0xW;
          Ke += qp.B.transpose() * C * qp.B * qp.A0xW;
        }

        else {
          using Material  = Marmot::Materials::MaterialFiniteStrain;
          using Tensor33d = Fastor::Tensor< double, 3, 3 >;

          const Eigen::Map< const CoordinateArray > U( QTotal_ );
          // the stretch from the reference tangent G and its change w = dU/dS, with lambda - 1 free of cancellation:
          // exactly zero without displacement, for any orientation of the truss
          const VectorDim G              = qp.T;
          const VectorDim w              = U * qp.dNdS.transpose();
          const VectorDim g              = G + w;
          const double    gNorm          = g.norm();
          const double    lambdaMinusOne = ( 2.0 * G.dot( w ) + w.dot( w ) ) / ( 1.0 + gNorm );
          const double    lambda         = 1.0 + lambdaMinusOne;
          const VectorDim n              = g / gNorm;

          Material::Deformation< 3 > deformation{ Tensor33d( 0.0 ) };
          deformation.F( 0, 0 ) = lambda;
          // the converged lateral stretches of the last update are the initial guess
          deformation.F( 1, 1 ) = sv.lateralStretches( 0 ) > 0 ? sv.lateralStretches( 0 ) : 1.0;
          deformation.F( 2, 2 ) = sv.lateralStretches( 1 ) > 0 ? sv.lateralStretches( 1 ) : 1.0;

          Tensor33d tauOld;
          Eigen::Map< Eigen::Matrix3d >( tauOld.data() ) = sv.kirchhoffStress; // symmetric: storage order irrelevant

          Material::ConstitutiveResponse< 3 > response( tauOld,
                                                        sv.elasticEnergy / qp.A0xW,
                                                        sv.dissipation / qp.A0xW,
                                                        sv.materialStateVars.data() );
          double                              dTau_dLambda = 0;
          qp.materialFiniteStrain->computeUniaxialStress( response, dTau_dLambda, deformation, { time, dT } );

          const double tau = response.tau( 0, 0 );
          const double J   = lambda * deformation.F( 1, 1 ) * deformation.F( 2, 2 );
          const double P   = tau / lambda;
          const double dP  = dTau_dLambda / lambda - tau / ( lambda * lambda );

          sv.lateralStretches << deformation.F( 1, 1 ), deformation.F( 2, 2 );
          sv.kirchhoffStress = Eigen::Map< const Eigen::Matrix3d >( response.tau.data() );
          sv.stress          = tau / J;
          sv.strain          = std::log1p( lambdaMinusOne );
          sv.normalForce     = P * elementProperties[0];
          sv.elasticEnergy   = response.elasticEnergyDensity * qp.A0xW;
          sv.dissipation     = response.dissipation * qp.A0xW;

          const MatrixDim nn       = n * n.transpose();
          const MatrixDim kTangent = dP * nn + P / lambda * ( MatrixDim::Identity() - nn );

          for ( int a = 0; a < nNodes; a++ ) {
            Pe.template segment< nDim >( a * nDim ) += P * n * qp.dNdS( a ) * qp.A0xW;
            for ( int b = 0; b < nNodes; b++ )
              Ke.template block< nDim, nDim >( a * nDim, b * nDim ) += kTangent * qp.dNdS( a ) * qp.dNdS( b ) * qp.A0xW;
          }
        }
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
                                   << __PRETTY_FUNCTION__
                                   << ": distributed loads are not supported by truss elements" );
    }

    void computeBodyForce( double* P_, double*, const double* load, const double*, double, double ) override
    {
      Eigen::Map< RhsSized >        P( P_ );
      Eigen::Map< const VectorDim > b( load );
      for ( const auto& qp : qps )
        for ( int a = 0; a < nNodes; a++ )
          P.template segment< nDim >( a * nDim ) += qp.N( a ) * b * qp.A0xW;
    }

    void computeConsistentInertia( double* M_ ) override
    {
      Eigen::Map< KSized > M( M_ );
      // the mass matrix needs a rule one order higher than the stiffness
      const double rho = qps[0].getDensity();
      for ( const auto& qpInfo : Marmot::Discretizations::Elements::Quadrature::
              getGaussPointInfo( shape(),
                                 Marmot::Discretizations::Elements::Quadrature::IntegrationTypes::FullIntegration ) ) {
        const double    xi   = qpInfo.xi( 0 );
        const NSized    N_   = N( xi );
        const VectorDim G    = coordinates * dNdXi( xi ).transpose();
        const double    mass = rho * elementProperties[0] * G.norm() * qpInfo.weight;
        for ( int a = 0; a < nNodes; a++ )
          for ( int b = 0; b < nNodes; b++ )
            M.template block< nDim, nDim >( a * nDim, b * nDim ) += N_( a ) * N_( b ) * mass * MatrixDim::Identity();
      }
    }

    void computeLumpedInertia( double* M_ ) override
    {
      // row sum lumping, which is positive for both the linear and the quadratic bar (1/6, 1/6, 2/3)
      KSized M = KSized::Zero();
      computeConsistentInertia( M.data() );
      Eigen::Map< RhsSized > m( M_ );
      m += M.rowwise().sum();
    }

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

      double* materialStateVars = qp.managedStateVars->materialStateVars.data();
      return qp.materialHypoElastic ? qp.materialHypoElastic->getStateView( stateName, materialStateVars )
                                    : qp.materialFiniteStrain->getStateView( stateName, materialStateVars );
    }

    std::vector< double > getCoordinatesAtCenter() override
    {
      const VectorDim x = coordinates * N( 0.0 ).transpose();
      return std::vector< double >( x.data(), x.data() + nDim );
    }

    std::vector< std::vector< double > > getCoordinatesAtQuadraturePoints() override
    {
      std::vector< std::vector< double > > result;
      for ( const auto& qp : qps ) {
        const VectorDim x = coordinates * N( qp.xi ).transpose();
        result.emplace_back( x.data(), x.data() + nDim );
      }
      return result;
    }

    int getNumberOfQuadraturePoints() override { return static_cast< int >( qps.size() ); }
  };

} // namespace Marmot::Discretizations::Elements
