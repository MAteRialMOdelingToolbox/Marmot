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
#include "Marmot/MarmotConsistentMass.h"
#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryElement.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
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
   * The truss carries only an axial force. The material sees a uniaxial stress state along the truss axis, found by
   * the uniaxial stress reductions of the material base classes. Element and material go in pairs:
   *
   * - @b T<d>D<n> (SmallStrain, geometrically linear): hypoelastic materials
   *   (Marmot::Materials::MaterialHypoElastic), e.g. LINEARELASTIC, VONMISES;
   * - @b T<d>D<n>CR (FiniteStrain): finite strain materials (Marmot::Materials::MaterialFiniteStrain) whose response
   *   to a diagonal deformation gradient is diagonal, e.g. the isotropic COMPRESSIBLENEOHOOKE,
   *   FINITESTRAINJ2PLASTICITY, VONMISES/HUGHES-WINGET.
   *
   * - @b SmallStrain: the axial strain is \f$\varepsilon = \boldsymbol{T} \cdot \frac{\partial \boldsymbol{u}}
   *   {\partial S}\f$ with the reference unit tangent \f$\boldsymbol{T}\f$ and the reference arc length \f$S\f$. As
   *   any geometrically linear element (unlike T2D2 in Abaqus with NLGEOM), a finite rotation strains the truss: use
   *   the CR variant whenever rotations are not small.
   * - @b FiniteStrain (total Lagrangian): the material is evaluated in the frame co-rotated with the bar's current
   *   tangent at the quadrature point, \f$\boldsymbol{F} = \mathrm{diag}( \lambda, \lambda_2, \lambda_3 )\f$ with the
   *   axial stretch \f$\lambda = |\boldsymbol{g}|\f$ of the current tangent \f$\boldsymbol{g} = \partial \boldsymbol{x}
   *   / \partial S\f$, and the lateral stretches \f$\lambda_2, \lambda_3\f$ found such that the lateral stresses
   *   vanish. Rigid body rotations are therefore exactly stress free. For 2 nodes the frame is the chord and the
   *   element is the co-rotational truss (Crisfield, Non-linear Finite Element Analysis of Solids and Structures,
   *   Vol. 1, ch. 3), hence "CR"; for 3 nodes the frame follows the curved bar point by point. The internal forces
   *   follow from the nominal axial stress \f$P = \tau_{11} / \lambda\f$ (the first Piola-Kirchhoff stress) and the
   *   current axis \f$\boldsymbol{n}\f$, \f$\boldsymbol{f}_a = \int P \, \boldsymbol{n} \frac{\partial N_a}
   *   {\partial S} A_0 \, \mathrm{d}S\f$, with the consistent tangent
   *   \f$\boldsymbol{K}_{ab} = \int \left[ \frac{\mathrm{d} P}{\mathrm{d}\lambda} \boldsymbol{n} \otimes
   *   \boldsymbol{n} + \frac{P}{\lambda} ( \boldsymbol{I} - \boldsymbol{n} \otimes \boldsymbol{n} ) \right]
   *   \frac{\partial N_a}{\partial S} \frac{\partial N_b}{\partial S} A_0 \, \mathrm{d}S\f$: the first term is the
   *   material stiffness along the bar, the second the geometric (initial stress) stiffness transverse to it, with
   *   \f$P / \lambda = S_{11}\f$ the second Piola-Kirchhoff stress -- the term behind the snap-through of a shallow
   *   truss.
   *
   * Element property: the (reference) cross section area \f$A_0\f$.
   *
   * Quadrature: Gauss rule exact to degree \f$2p - 1\f$ (1 point for 2 nodes, 2 points for 3 nodes), exact for the
   * stiffness of straight trusses with a linear elastic material.
   *
   * Limitations: a truss has no lateral stiffness. A node connected only to collinear small-strain trusses (or to CR
   * trusses without tension) is free laterally and the stiffness matrix is singular; a quadratic truss additionally
   * has internal zero-energy modes of its mid node transverse to the axis (one in 2D, two in 3D). Restrain them by
   * embedding, constraints, boundary conditions or pre-tension.
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
      SmallStrain,  ///< geometrically linear, hypoelastic materials
      FiniteStrain, ///< total Lagrangian in the co-rotated frame, finite strain materials
    };

    /// The number of degrees of freedom: nDim displacements per node
    static constexpr int sizeLoadVector = nNodes * nDim;

    using NSized          = Eigen::Matrix< double, 1, nNodes >;                              ///< shape functions (row)
    using CoordinateArray = Eigen::Matrix< double, nDim, nNodes >;                           ///< nodal coordinates
    using VectorDim       = Eigen::Matrix< double, nDim, 1 >;                                ///< a spatial vector
    using MatrixDim       = Eigen::Matrix< double, nDim, nDim >;                             ///< a spatial tensor
    using RhsSized        = Eigen::Matrix< double, sizeLoadVector, 1 >;                      ///< element force vector
    using KSized          = Eigen::Matrix< double, sizeLoadVector, sizeLoadVector >;         ///< element matrix
    using BSized          = Eigen::Matrix< double, 1, sizeLoadVector >;                      ///< axial strain operator
    using ParentGeometry  = Marmot::Discretizations::Elements::GeometryElement< 1, nNodes >; ///< the parent line

    /// Element label
    const int elLabel;
    /// Kinematic formulation of this instance
    const Kinematics kinematics;

    /// Element properties: [cross section area]
    Eigen::Map< const Eigen::VectorXd > elementProperties;
    /// Reference nodal coordinates, [dim, node]
    Eigen::Map< const CoordinateArray > coordinates;

    /// The state of a quadrature point, see the class documentation for the layout; the lateral stretches and the
    /// Kirchhoff stress exist only for the finite strain formulation
    class QPStateVarManager : public Marmot::Utilities::StateVarVectorManager {

      /// \hideinitializer
      inline const static auto smallStrainLayout = makeLayout( {
        { .name = "stress", .length = 1 },
        { .name = "strain", .length = 1 },
        { .name = "normal force", .length = 1 },
        { .name = "elastic energy", .length = 1 },
        { .name = "dissipation", .length = 1 },
        { .name = "begin of material state", .length = 0 },
      } );

      /// \hideinitializer
      inline const static auto finiteStrainLayout = makeLayout( {
        { .name = "stress", .length = 1 },
        { .name = "strain", .length = 1 },
        { .name = "normal force", .length = 1 },
        { .name = "elastic energy", .length = 1 },
        { .name = "dissipation", .length = 1 },
        { .name = "lateral stretches", .length = 2 },
        { .name = "kirchhoff stress", .length = 9 },
        { .name = "begin of material state", .length = 0 },
      } );

      static const StateVarVectorLayout& layout( Kinematics kinematics )
      {
        return kinematics == SmallStrain ? smallStrainLayout : finiteStrainLayout;
      }

    public:
      double&                       stress;            ///< axial Cauchy stress
      double&                       strain;            ///< axial strain (finite strain: logarithmic)
      double&                       normalForce;       ///< normal force
      double&                       elasticEnergy;     ///< elastic energy of the point's volume
      double&                       dissipation;       ///< dissipation of the point's volume
      Eigen::Map< Eigen::Vector2d > lateralStretches;  ///< finite strain: the converged lateral stretches
      Eigen::Map< Eigen::Matrix3d > kirchhoffStress;   ///< finite strain: the Kirchhoff stress in the co-rotated frame
      Eigen::Map< Eigen::VectorXd > materialStateVars; ///< the state of the material

      /// The number of state variables of a point without the material's
      /// @param[in] kinematics The kinematic formulation.
      /// @return The number of state variables.
      static int getNumberOfRequiredStateVarsQuadraturePointOnly( Kinematics kinematics )
      {
        return layout( kinematics ).nRequiredStateVars;
      }

      /// @param[in] theStateVarVector The state of the point.
      /// @param[in] nStateVars        Its length.
      /// @param[in] kinematics        The kinematic formulation.
      QPStateVarManager( double* theStateVarVector, int nStateVars, Kinematics kinematics )
        : Marmot::Utilities::StateVarVectorManager( theStateVarVector, layout( kinematics ) ),
          stress( find( "stress" ) ),
          strain( find( "strain" ) ),
          normalForce( find( "normal force" ) ),
          elasticEnergy( find( "elastic energy" ) ),
          dissipation( find( "dissipation" ) ),
          lateralStretches( kinematics == FiniteStrain ? &find( "lateral stretches" ) : nullptr ),
          kirchhoffStress( kinematics == FiniteStrain ? &find( "kirchhoff stress" ) : nullptr ),
          materialStateVars( &find( "begin of material state" ),
                             nStateVars - getNumberOfRequiredStateVarsQuadraturePointOnly( kinematics ) ){};
    };

    /// A quadrature point along the truss axis
    struct QuadraturePoint {
      const double xi;        ///< parent coordinate
      const double weight;    ///< quadrature weight

      NSized    N;            ///< shape functions
      NSized    dN_dS;        ///< derivatives of the shape functions w.r.t. the reference arc length
      VectorDim T;            ///< reference unit tangent dX/dS
      double    dS_dXi = 0.0; ///< reference length per parent length, |dX/dxi|
      double    dV0    = 0.0; ///< reference volume of the point: A0 dS/dxi weight
      BSized    dStrain_dU;   ///< small strain: axial strain operator, eps = dStrain_dU . U

      std::unique_ptr< QPStateVarManager >                       managedStateVars;
      std::unique_ptr< Marmot::Materials::MaterialHypoElastic >  materialHypoElastic;
      std::unique_ptr< Marmot::Materials::MaterialFiniteStrain > materialFiniteStrain;

      /// @param[in] xi     Parent coordinate.
      /// @param[in] weight Quadrature weight.
      QuadraturePoint( double xi, double weight ) : xi( xi ), weight( weight ) {}

      /// @return Whether the material has been assigned.
      bool hasMaterial() const { return materialHypoElastic || materialFiniteStrain; }

      /// @return The number of state variables of the material.
      int getNumberOfRequiredMaterialStateVars() const
      {
        return materialHypoElastic ? materialHypoElastic->getNumberOfRequiredStateVars()
                                   : materialFiniteStrain->getNumberOfRequiredStateVars();
      }

      /// @return The mass density of the material.
      double getDensity() const
      {
        const double* sv = managedStateVars->materialStateVars.data();
        return materialHypoElastic ? materialHypoElastic->getDensity( sv ) : materialFiniteStrain->getDensity( sv );
      }

      /// @param[in] length The characteristic element length (regularisation of softening materials).
      void setCharacteristicElementLength( double length )
      {
        if ( materialHypoElastic )
          materialHypoElastic->setCharacteristicElementLength( length );
        else
          materialFiniteStrain->setCharacteristicElementLength( length );
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
      // degree 2p - 1: exact for the stiffness of straight trusses of constant cross section and linear elastic
      // material (1 point for 2 nodes, 2 points for 3 nodes)
      using namespace Marmot::Discretizations::Elements::Quadrature;
      for ( const auto& point : getQuadraturePoints( shape(), gaussExactToDegree( 2 * ( nNodes - 1 ) - 1 ) ) )
        qps.emplace_back( point.xi( 0 ), point.weight );
    }

    /// @return The shape of the element: Bar2 or Bar3.
    static constexpr Marmot::Discretizations::Elements::ElementShapes shape()
    {
      return nNodes == 2 ? Marmot::Discretizations::Elements::ElementShapes::Bar2
                         : Marmot::Discretizations::Elements::ElementShapes::Bar3;
    }

    /// @param[in] xi Parent coordinate.
    /// @return The shape functions (node order: end, end, mid).
    static NSized N( double xi )
    {
      if constexpr ( nNodes == 2 )
        return Marmot::Discretizations::Elements::Spatial1D::Bar2::N( xi );
      else
        return Marmot::Discretizations::Elements::Spatial1D::Bar3::N( xi );
    }

    /// @param[in] xi Parent coordinate.
    /// @return The derivatives of the shape functions w.r.t. the parent coordinate.
    static NSized dNdXi( double xi )
    {
      if constexpr ( nNodes == 2 )
        return Marmot::Discretizations::Elements::Spatial1D::Bar2::dNdXi( xi );
      else
        return Marmot::Discretizations::Elements::Spatial1D::Bar3::dNdXi( xi );
    }

    /// The number of state variables; known once the material section is assigned (before or after
    /// initializeYourself)
    int getNumberOfRequiredStateVars() override
    {
      if ( !qps[0].hasMaterial() )
        throw std::logic_error( Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << ": element " << elLabel
                                                                << ": assign the material section first" );
      return static_cast< int >( qps.size() ) *
             ( QPStateVarManager::getNumberOfRequiredStateVarsQuadraturePointOnly( kinematics ) +
               qps[0].getNumberOfRequiredMaterialStateVars() );
    }

    /// @return One node field per node: displacement.
    std::vector< std::vector< std::string > > getNodeFields() override
    {
      return std::vector< std::vector< std::string > >( nNodes, { "displacement" } );
    }

    /// @return The identity: the dofs are ordered node by node.
    std::vector< int > getDofIndicesPermutationPattern() override
    {
      std::vector< int > pattern( sizeLoadVector );
      for ( int i = 0; i < sizeLoadVector; i++ )
        pattern[i] = i;
      return pattern;
    }

    /// @return The number of nodes.
    int getNNodes() override { return nNodes; }

    /// @return The number of spatial dimensions.
    int getNSpatialDimensions() override { return nDim; }

    /// @return The number of degrees of freedom.
    int getNDofPerElement() override { return sizeLoadVector; }

    /// @return The EnSight shape: bar2 or bar3.
    std::string getElementShape() override { return nNodes == 2 ? "bar2" : "bar3"; }

    /// Assign the state variables, nQp equal blocks.
    /// @param[in] stateVars  The state variables of the element.
    /// @param[in] nStateVars Their number (at least getNumberOfRequiredStateVars()).
    void assignStateVars( double* stateVars, int nStateVars ) override
    {
      if ( nStateVars < getNumberOfRequiredStateVars() )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__ << ": element " << elLabel << ": too few state variables" );
      const int nQpStateVars = nStateVars / static_cast< int >( qps.size() );
      for ( size_t i = 0; i < qps.size(); i++ )
        qps[i].managedStateVars = std::make_unique< QPStateVarManager >( stateVars + i * nQpStateVars,
                                                                         nQpStateVars,
                                                                         kinematics );
    }

    /// Assign the element properties: [cross section area].
    /// @param[in] property The element properties.
    void assignProperty( const ElementProperties& property ) override
    {
      if ( property.nElementProperties < 1 )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__
                                     << ": a truss element requires its cross section area as property" );
      new ( &elementProperties )
        Eigen::Map< const Eigen::VectorXd >( property.elementProperties, property.nElementProperties );
    }

    /// Create one material per quadrature point: a hypoelastic one for T<d>D<n>, a finite strain one for
    /// T<d>D<n>CR.
    /// @param[in] section The material section.
    void assignProperty( const Marmot::Discretizations::MaterialSection& section ) override
    {
      for ( auto& qp : qps ) {
        try {
          if ( kinematics == SmallStrain )
            qp.materialHypoElastic = std::unique_ptr< Marmot::Materials::MaterialHypoElastic >(
              Marmot::Registry::MaterialHypoElasticFactory::createMaterial( section.materialName,
                                                                            section.materialProperties,
                                                                            section.nMaterialProperties,
                                                                            elLabel ) );
          else
            qp.materialFiniteStrain = std::unique_ptr< Marmot::Materials::MaterialFiniteStrain >(
              Marmot::Registry::MaterialFiniteStrainFactory::createMaterial( section.materialName,
                                                                             section.materialProperties,
                                                                             section.nMaterialProperties,
                                                                             elLabel ) );
        }
        catch ( const std::invalid_argument& e ) {
          throw std::invalid_argument(
            Marmot::Utilities::MakeString()
            << e.what() << "\nelement " << elLabel << ": "
            << ( kinematics == SmallStrain
                   ? "the geometrically linear truss T<d>D<n> takes hypoelastic materials; for finite strain "
                     "materials use T<d>D<n>CR"
                   : "the finite strain truss T<d>D<n>CR takes finite strain materials (e.g., COMPRESSIBLENEOHOOKE, "
                     "VONMISES/HUGHES-WINGET); for hypoelastic materials use T<d>D<n>" ) );
        }
      }
      setCharacteristicElementLength();
    }

    /// The characteristic element length of the materials, the element length (volume-equivalent length of the
    /// parent line), as soon as both the materials and the geometry (initializeYourself) are known
    void setCharacteristicElementLength()
    {
      if ( !qps[0].hasMaterial() || !( qps[0].dS_dXi > 0 ) )
        return;
      for ( auto& qp : qps )
        qp.setCharacteristicElementLength( ParentGeometry::volumeEquivalentLength( qp.dS_dXi ) );
    }

    /// @param[in] coords The reference nodal coordinates, node by node.
    void assignNodeCoordinates( const double* coords ) override
    {
      new ( &coordinates ) Eigen::Map< const CoordinateArray >( coords );
    }

    /// The geometry at the quadrature points: tangent, shape function derivatives w.r.t. the arc length, volume and
    /// the small strain operator.
    void initializeYourself() override
    {
      if ( elementProperties.size() == 0 )
        throw std::logic_error( Marmot::Utilities::MakeString()
                                << __PRETTY_FUNCTION__ << ": element " << elLabel
                                << ": assign the element properties before initializing" );
      for ( auto& qp : qps ) {
        qp.N                   = N( qp.xi );
        const NSized    dN_dXi = dNdXi( qp.xi );
        const VectorDim dX_dXi = coordinates * dN_dXi.transpose();
        qp.dS_dXi              = dX_dXi.norm();
        qp.T                   = dX_dXi / qp.dS_dXi;
        qp.dN_dS               = dN_dXi / qp.dS_dXi;
        qp.dV0                 = elementProperties[0] * qp.dS_dXi * qp.weight;
        for ( int a = 0; a < nNodes; a++ )
          qp.dStrain_dU.template segment< nDim >( a * nDim ) = qp.dN_dS( a ) * qp.T.transpose();
      }
      setCharacteristicElementLength();
    }

    /// Initialize the material states.
    /// @param[in] state Only MarmotMaterialInitialization.
    void setInitialConditions( StateTypes state, const double* ) override
    {
      switch ( state ) {
      case Marmot::Discretizations::Elements::Element::MarmotMaterialInitialization: {
        for ( auto& qp : qps ) {
          auto& sv = *qp.managedStateVars;
          if ( qp.materialHypoElastic )
            qp.materialHypoElastic->initializeYourself( sv.materialStateVars.data(),
                                                        static_cast< int >( sv.materialStateVars.size() ) );
          else {
            sv.lateralStretches.setOnes();
            qp.materialFiniteStrain->initializeYourself( sv.materialStateVars.data(),
                                                         static_cast< int >( sv.materialStateVars.size() ) );
          }
        }
        break;
      }
      default:
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__ << ": invalid initial condition" );
      }
    }

    /// Internal forces and tangent, see the class documentation.
    /// @param[in]  QTotal_ The total nodal displacements.
    /// @param[in]  dQ_     Their increment.
    /// @param[out] Pe_     The internal forces (added).
    /// @param[out] Ke_     The tangent (added).
    /// @param[in]  time    (Pseudo-)time.
    /// @param[in]  dT      (Pseudo-)time increment.
    void computeKernels( const double* QTotal_, const double* dQ_, double* Pe_, double* Ke_, double time, double dT )
      override
    {
      if ( kinematics == SmallStrain )
        computeKernelsSmallStrain( QTotal_, dQ_, Pe_, Ke_, time, dT );
      else
        computeKernelsFiniteStrain( QTotal_, Pe_, Ke_, time, dT );
    }

    /// Geometrically linear: eps = T . du/dS, the material in uniaxial stress.
    /// @param[in]  QTotal_ The total nodal displacements.
    /// @param[in]  dQ_     Their increment.
    /// @param[out] Pe_     The internal forces (added).
    /// @param[out] Ke_     The tangent (added).
    /// @param[in]  time    (Pseudo-)time.
    /// @param[in]  dT      (Pseudo-)time increment.
    void computeKernelsSmallStrain( const double* QTotal_,
                                    const double* dQ_,
                                    double*       Pe_,
                                    double*       Ke_,
                                    double        time,
                                    double        dT )
    {
      Eigen::Map< const RhsSized > QTotal( QTotal_ );
      Eigen::Map< const RhsSized > dQ( dQ_ );
      Eigen::Map< RhsSized >       Pe( Pe_ );
      Eigen::Map< KSized >         Ke( Ke_ );

      for ( auto& qp : qps ) {
        auto&                                           sv = *qp.managedStateVars;
        Marmot::Materials::MaterialHypoElastic::state1D state;
        state.stress               = sv.stress;
        state.elasticEnergyDensity = sv.elasticEnergy / qp.dV0;
        state.dissipation          = sv.dissipation / qp.dV0;
        state.stateVars            = sv.materialStateVars.data();

        double dStress_dStrain = 0;
        qp.materialHypoElastic->computeUniaxialStress( state, dStress_dStrain, qp.dStrain_dU * dQ, { time, dT } );

        sv.stress        = state.stress;
        sv.strain        = qp.dStrain_dU * QTotal;
        sv.normalForce   = state.stress * elementProperties[0];
        sv.elasticEnergy = state.elasticEnergyDensity * qp.dV0;
        sv.dissipation   = state.dissipation * qp.dV0;

        Pe += qp.dStrain_dU.transpose() * state.stress * qp.dV0;
        Ke += qp.dStrain_dU.transpose() * dStress_dStrain * qp.dStrain_dU * qp.dV0;
      }
    }

    /// Finite strain: the material in the frame co-rotated with the current tangent, see the class documentation.
    /// @param[in]  QTotal_ The total nodal displacements.
    /// @param[out] Pe_     The internal forces (added).
    /// @param[out] Ke_     The tangent (added).
    /// @param[in]  time    (Pseudo-)time.
    /// @param[in]  dT      (Pseudo-)time increment.
    void computeKernelsFiniteStrain( const double* QTotal_, double* Pe_, double* Ke_, double time, double dT )
    {
      using Material  = Marmot::Materials::MaterialFiniteStrain;
      using Tensor33d = Fastor::Tensor< double, 3, 3 >;

      Eigen::Map< RhsSized >                    Pe( Pe_ );
      Eigen::Map< KSized >                      Ke( Ke_ );
      const Eigen::Map< const CoordinateArray > U( QTotal_ );

      for ( auto& qp : qps ) {
        auto& sv = *qp.managedStateVars;

        // the stretch from the reference unit tangent T and the displacement gradient du_dS, with lambda - 1 free of
        // cancellation: exactly zero without displacement, for any orientation of the truss
        const VectorDim du_dS          = U * qp.dN_dS.transpose();
        const VectorDim dx_dS          = qp.T + du_dS;
        const double    dx_dSNorm      = dx_dS.norm();
        const double    lambdaMinusOne = ( 2.0 * qp.T.dot( du_dS ) + du_dS.dot( du_dS ) ) / ( 1.0 + dx_dSNorm );
        const double    lambda         = 1.0 + lambdaMinusOne;
        const VectorDim n              = dx_dS / dx_dSNorm;

        Material::Deformation< 3 > deformation{ Tensor33d( 0.0 ) };
        deformation.F( 0, 0 ) = lambda;
        // the converged lateral stretches of the last update are the initial guess
        deformation.F( 1, 1 ) = sv.lateralStretches( 0 ) > 0 ? sv.lateralStretches( 0 ) : 1.0;
        deformation.F( 2, 2 ) = sv.lateralStretches( 1 ) > 0 ? sv.lateralStretches( 1 ) : 1.0;

        Tensor33d tauOld;
        Eigen::Map< Eigen::Matrix3d >( tauOld.data() ) = sv.kirchhoffStress; // symmetric: storage order irrelevant

        Material::ConstitutiveResponse< 3 > response( tauOld,
                                                      sv.elasticEnergy / qp.dV0,
                                                      sv.dissipation / qp.dV0,
                                                      sv.materialStateVars.data() );
        double                              dTau_dLambda = 0;
        qp.materialFiniteStrain->computeUniaxialStress( response, dTau_dLambda, deformation, { time, dT } );

        const double tau        = response.tau( 0, 0 );
        const double J          = lambda * deformation.F( 1, 1 ) * deformation.F( 2, 2 );
        const double P          = tau / lambda; // first Piola-Kirchhoff stress
        const double dP_dLambda = dTau_dLambda / lambda - tau / ( lambda * lambda );

        sv.lateralStretches << deformation.F( 1, 1 ), deformation.F( 2, 2 );
        sv.kirchhoffStress = Eigen::Map< const Eigen::Matrix3d >( response.tau.data() );
        sv.stress          = tau / J;
        sv.strain          = std::log1p( lambdaMinusOne );
        sv.normalForce     = P * elementProperties[0];
        sv.elasticEnergy   = response.elasticEnergyDensity * qp.dV0;
        sv.dissipation     = response.dissipation * qp.dV0;

        // material stiffness along the bar, geometric (initial stress) stiffness transverse to it, P / lambda = S11
        const MatrixDim nn         = n * n.transpose();
        const MatrixDim kMaterial  = dP_dLambda * nn;
        const MatrixDim kGeometric = P / lambda * ( MatrixDim::Identity() - nn );
        const MatrixDim kTangent   = kMaterial + kGeometric;

        for ( int a = 0; a < nNodes; a++ ) {
          Pe.template segment< nDim >( a * nDim ) += P * n * qp.dN_dS( a ) * qp.dV0;
          for ( int b = 0; b < nNodes; b++ )
            Ke.template block< nDim, nDim >( a * nDim, b * nDim ) += kTangent * qp.dN_dS( a ) * qp.dN_dS( b ) * qp.dV0;
        }
      }
    }

    /// The internal forces only (explicit dynamics).
    /// @param[in]  QTotal The total nodal displacements.
    /// @param[in]  dQ     Their increment.
    /// @param[out] Pe     The internal forces (added).
    /// @param[in]  time   (Pseudo-)time.
    /// @param[in]  dT     (Pseudo-)time increment.
    void computeKernelsExplicit( const double* QTotal, const double* dQ, double* Pe, double time, double dT ) override
    {
      KSized K = KSized::Zero();
      computeKernels( QTotal, dQ, Pe, K.data(), time, dT );
    }

    /// Not supported: a truss has no surfaces.
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

    /// The consistent nodal forces of a body force (per volume), in the reference configuration.
    /// @param[out] P_   The nodal forces (added).
    /// @param[in]  load The body force vector.
    void computeBodyForce( double* P_, double*, const double* load, const double*, double, double ) override
    {
      Eigen::Map< RhsSized >        P( P_ );
      Eigen::Map< const VectorDim > b( load );
      for ( const auto& qp : qps )
        for ( int a = 0; a < nNodes; a++ )
          P.template segment< nDim >( a * nDim ) += qp.N( a ) * b * qp.dV0;
    }

    /// The consistent mass matrix, integrated with ConsistentMass::integrationRule(); the density of the nearest
    /// quadrature point of the element.
    /// @param[out] M_ The mass matrix (added).
    void computeConsistentInertia( double* M_ ) override
    {
      Eigen::Map< KSized > M( M_ );
      struct PointXi {
        Eigen::VectorXd xi;
      };
      std::vector< PointXi > elementPoints;
      for ( const auto& qp : qps )
        elementPoints.push_back( { Eigen::VectorXd::Constant( 1, qp.xi ) } );

      for ( const auto& point : ConsistentMass::integrationRule( shape() ) ) {
        const double    xi     = point.xi( 0 );
        const double    rho    = qps[ConsistentMass::nearestQuadraturePoint( point.xi, elementPoints )].getDensity();
        const NSized    N_     = N( xi );
        const VectorDim dX_dXi = coordinates * dNdXi( xi ).transpose();
        const double    mass   = rho * elementProperties[0] * dX_dXi.norm() * point.weight;
        for ( int a = 0; a < nNodes; a++ )
          for ( int b = 0; b < nNodes; b++ )
            M.template block< nDim, nDim >( a * nDim, b * nDim ) += N_( a ) * N_( b ) * mass * MatrixDim::Identity();
      }
    }

    /**
     * @brief The critical time step of explicit dynamics: the stiffest axial mode of the bar with the row-sum lumped
     * mass (computeLumpedInertia).
     *
     * @details With the reference length \f$L\f$ and the wave speed \f$c\f$ of the material: 2 nodes, masses
     * \f$\rho A L (1/2, 1/2)\f$, \f$\omega_{\max}^2 = 4 c^2 / L^2\f$, \f$\Delta t = 2 / \omega_{\max} = L / c\f$;
     * 3 nodes, masses \f$\rho A L (1/6, 1/6, 2/3)\f$, \f$\omega_{\max}^2 = 24 c^2 / L^2\f$, \f$\Delta t = L / (\sqrt{6}
     * c)\f$. The wave speed is the materials' maximum (P-wave) one, larger than the bar wave speed
     * \f$\sqrt{E/\rho}\f$: the estimate is on the safe side.
     * @param[out] criticalTimeStep The critical time step.
     * @param[in]  QTotal           The total nodal displacements (the finite strain materials are evaluated at the
     * current stretch).
     */
    void computeCriticalTimeStepForExplicitDynamics( double& criticalTimeStep, const double* QTotal ) override
    {
      const double timeStepFactor = nNodes == 2 ? 1.0 : 1.0 / std::sqrt( 6.0 );
      double       length         = 0.0;
      for ( const auto& qp : qps )
        length += qp.dS_dXi * qp.weight;

      const Eigen::Map< const CoordinateArray > U( QTotal );
      criticalTimeStep = std::numeric_limits< double >::max();
      for ( const auto& qp : qps ) {
        auto&  sv = *qp.managedStateVars; // the material needs mutable state pointers
        double c  = 0.0;
        if ( qp.materialHypoElastic ) {
          Marmot::Math::Vector6d stress = Marmot::Math::Vector6d::Zero();
          stress( 0 )                   = sv.stress;
          c                             = qp.materialHypoElastic->getMaximumWaveSpeed(
            Marmot::Materials::MaterialHypoElastic::state3D( stress,
                                                             sv.elasticEnergy / qp.dV0,
                                                             sv.dissipation / qp.dV0,
                                                             sv.materialStateVars.data() ) );
        }
        else {
          const VectorDim                dx_dS = qp.T + U * qp.dN_dS.transpose();
          Fastor::Tensor< double, 3, 3 > F( 0.0 );
          F( 0, 0 ) = dx_dS.norm();
          F( 1, 1 ) = sv.lateralStretches( 0 ) > 0 ? sv.lateralStretches( 0 ) : 1.0;
          F( 2, 2 ) = sv.lateralStretches( 1 ) > 0 ? sv.lateralStretches( 1 ) : 1.0;
          c         = qp.materialFiniteStrain->getMaximumWaveSpeed( sv.materialStateVars.data(), F );
        }
        if ( !( c > 0.0 ) )
          throw std::runtime_error( Marmot::Utilities::MakeString()
                                    << __PRETTY_FUNCTION__ << ": element " << elLabel << ": non-positive wave speed" );
        criticalTimeStep = std::min( criticalTimeStep, timeStepFactor * length / c );
      }
    }

    /// Row sum lumping, which is positive for both the linear and the quadratic bar (1/6, 1/6, 2/3).
    /// @param[out] M_ The lumped masses per dof (added).
    void computeLumpedInertia( double* M_ ) override
    {
      KSized M = KSized::Zero();
      computeConsistentInertia( M.data() );
      Eigen::Map< RhsSized > m( M_ );
      m += M.rowwise().sum();
    }

    /// @param[out] internalEnergy The elastic energy of the element.
    void computeInternalEnergy( double& internalEnergy ) override
    {
      internalEnergy = 0.0;
      for ( const auto& qp : qps )
        internalEnergy += qp.managedStateVars->elasticEnergy;
    }

  protected:
    /// @param[in] stateName The name of a state of the point or of its material.
    /// @param[in] qpNumber  The quadrature point.
    /// @return A view of the state.
    Marmot::Utilities::StateView getStateViewAtQuadraturePoint( const std::string& stateName, int qpNumber ) override
    {
      const auto& qp = qps[qpNumber];
      if ( qp.managedStateVars->contains( stateName ) )
        return qp.managedStateVars->getStateView( stateName );

      double* materialStateVars = qp.managedStateVars->materialStateVars.data();
      return qp.materialHypoElastic ? qp.materialHypoElastic->getStateView( stateName, materialStateVars )
                                    : qp.materialFiniteStrain->getStateView( stateName, materialStateVars );
    }

  public:
    /// @return The coordinates of the center (parent coordinate 0).
    std::vector< double > getCoordinatesAtCenter() override
    {
      const VectorDim x = coordinates * N( 0.0 ).transpose();
      return std::vector< double >( x.data(), x.data() + nDim );
    }

    /// @return The reference coordinates of the quadrature points.
    std::vector< std::vector< double > > getCoordinatesAtQuadraturePoints() override
    {
      std::vector< std::vector< double > > result;
      for ( const auto& qp : qps ) {
        const VectorDim x = coordinates * N( qp.xi ).transpose();
        result.emplace_back( x.data(), x.data() + nDim );
      }
      return result;
    }

    /// @return The number of quadrature points.
    int getNumberOfQuadraturePoints() override { return static_cast< int >( qps.size() ); }
  };

} // namespace Marmot::Discretizations::Elements
