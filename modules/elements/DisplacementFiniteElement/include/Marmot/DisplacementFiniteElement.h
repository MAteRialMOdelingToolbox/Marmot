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
#include "Marmot/MarmotGeometryElement.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotLowerDimensionalStress.h"
#include "Marmot/MarmotMassLumping.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

using namespace Marmot;
using namespace Eigen;

namespace Marmot::Elements {

  /**
   * @class Marmot::Elements::DisplacementFiniteElement
   * @tparam nDim Number of spatial dimensions (1, 2, or 3).
   * @tparam nNodes Number of element nodes.
   * @brief Displacement-based finite element template.
   * @details Uses linearized kinematics (small strains) and supports 1D, 2D, and 3D
   * formulations via a section assumption. Holds quadrature-point state and
   * delegates constitutive updates to Marmot materials while assembling element
   * residuals, tangents, and mass matrices.
   */
  template < int nDim, int nNodes >
  class DisplacementFiniteElement : public MarmotElement, public MarmotGeometryElement< nDim, nNodes > {

  public:
    /**
     * @brief Kinematic section assumption used by the element.
     * @details Controls which constitutive call is performed at quadrature points.
     * - UniaxialStress
     * - PlaneStress
     * - PlaneStrain
     * - Solid
     */
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
    using CSized                = Matrix< double, ParentGeometryElement::voigtSize, ParentGeometryElement::voigtSize >;
    using Voigt                 = Matrix< double, ParentGeometryElement::voigtSize, 1 >;

    /** Element-level properties (e.g., thickness for 2D, area for 1D). */
    Map< const VectorXd > elementProperties;
    /** Element label (ID) used for logging and material creation. */
    const int elLabel;
    /** Section assumption applied by this element instance. */
    const SectionType sectionType;

    /**
     * @brief Data and state associated with a quadrature point.
     * @details Holds parent coordinates, integration weight, Jacobian determinant,
     *          kinematic (strain-displacement) B-matrix, and a material instance with
     *          managed state variables.
     */
    struct QuadraturePoint {

      const XiSized xi;
      const double  weight;

      double detJ;
      double J0xW;
      BSized B;

      /**
       * @brief Manager for per-quadrature-point state variables.
       * @details Provides named accessors to stress \f$\sig\f$, strain \f$\eps\f$
       * and the material state vector. The layout is [stress(6), strain(6), begin of material state(...)]
       * in 3D Voigt notation.
       */
      class QPStateVarManager : public MarmotStateVarVectorManager {

        /// \hideinitializer
        inline const static auto layout = makeLayout( {
          { .name = "stress", .length = 6 },
          { .name = "strain", .length = 6 },
          { .name = "total strain energy", .length = 1 },
          { .name = "elastic strain energy", .length = 1 },
          { .name = "dissipation", .length = 1 },
          { .name = "begin of material state", .length = 0 },
        } );

      public:
        mVector6d                     stress;
        mVector6d                     strain;
        double&                       totalStrainEnergy;
        double&                       elasticStrainEnergy;
        double&                       dissipation;
        Eigen::Map< Eigen::VectorXd > materialStateVars;

        static int getNumberOfRequiredStateVarsQuadraturePointOnly() { return layout.nRequiredStateVars; };

        QPStateVarManager( double* theStateVarVector, int nStateVars )
          : MarmotStateVarVectorManager( theStateVarVector, layout ),
            stress( &find( "stress" ) ),
            strain( &find( "strain" ) ),
            totalStrainEnergy( find( "total strain energy" ) ),
            elasticStrainEnergy( find( "elastic strain energy" ) ),
            dissipation( find( "dissipation" ) ),
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
      }

      QuadraturePoint( XiSized xi, double weight )
        : xi( xi ), weight( weight ), detJ( 0.0 ), J0xW( 0.0 ), B( BSized::Zero() ){};
    };

    /// Quadrature points owned by the element (one per integration point).
    std::vector< QuadraturePoint > qps;

    /**
     * @brief Construct element with ID, quadrature rule and section assumption.
     * @param elementID Unique element label.
     * @param integrationType Integration (quadrature) rule.
     * @param sectionType Section assumption (1D/2D/3D).
     */
    DisplacementFiniteElement( int                                         elementID,
                               FiniteElement::Quadrature::IntegrationTypes integrationType,
                               SectionType                                 sectionType );

    /** @brief Total number of required state variables for this element (sum over all quadrature points). */
    int getNumberOfRequiredStateVars();

    /** @brief Node-level fields exposed by the element. Returns ["displacement"] for each node. */
    std::vector< std::vector< std::string > > getNodeFields();

    /** @brief Permutation pattern from local DOF ordering to solver ordering (identity by default). */
    std::vector< int > getDofIndicesPermutationPattern();

    /** @brief Number of nodes of this element type. */
    int getNNodes() { return nNodes; }

    /** @brief Number of spatial dimensions. */
    int getNSpatialDimensions() { return nDim; }

    /** @brief Number of degrees of freedom per element (nNodes * nDim). */
    int getNDofPerElement() { return sizeLoadVector; }

    /** @brief Geometric shape of the element (as reported by the parent geometry element). */
    std::string getElementShape() { return ParentGeometryElement::getElementShape(); }

    /** @brief Map the provided element state vector to all quadrature points. */
    void assignStateVars( double* stateVars, int nStateVars );

    /** @brief Assign element properties (e.g., thickness in 2D, area in 1D). */
    void assignProperty( const ElementProperties& marmotElementProperty );

    /** @brief Assign material section and instantiate per-quadrature-point materials. */
    void assignProperty( const MarmotMaterialSection& marmotElementProperty );

    /** @brief Provide nodal coordinates to the parent geometry element. */
    void assignNodeCoordinates( const double* coordinates );

    /** @brief Precompute geometry-related quantities at quadrature points (B, detJ, J0xW). */
    void initializeYourself();

    /**
     * @brief Initialize state or materials.
     * @param state MarmotMaterialInitialization, GeostaticStress or MarmotMaterialStateVars.
     * @param values For GeostaticStress: [sigmaY(z1), y1, sigmaY(z2), y2, kx, kz].
     */
    void setInitialConditions( StateTypes state, const double* values );

    /**
     * @brief Assemble distributed surface loads on a boundary face.
     * @details Pressure and traction contributions are integrated on the boundary \f$\Gamma_e\f$:
     * \f[
     * \mathbf{P}_e^{(p)} = - \int_{\Gamma_e} p\, \mathbf{N}^\mathsf{T} \mathbf{n}\, \mathrm{d}\Gamma,\qquad
     * \mathbf{P}_e^{(t)} = \int_{\Gamma_e} \mathbf{N}^\mathsf{T} \mathbf{t}\, \mathrm{d}\Gamma.
     * \f]
     * @param loadType Pressure or SurfaceTraction.
     * @param P Element RHS contribution (accumulated).
     * @param K Optional stiffness contribution (unused).
     * @param elementFace Boundary face index.
     * @param load Pressure magnitude or traction vector (size nDim).
     * @param QTotal Total DOF vector (unused).
     * @param time Current time data forwarded to materials.
     * @param dT Time increment.
     */
    void computeDistributedLoad( MarmotElement::DistributedLoadTypes loadType,
                                 double*                             P,
                                 double*                             K,
                                 const int                           elementFace,
                                 const double*                       load,
                                 const double*                       QTotal,
                                 double                              time,
                                 double                              dT );

    /**
     * @brief Assemble body force contribution.
     * @details Integrates \f$\mathbf{P}_e^{(b)} = \int_{\Omega_e} \mathbf{N}^\mathsf{T} \mathbf{f}\,
     * \mathrm{d}\Omega\f$.
     */
    void computeBodyForce( double* P, double* K, const double* load, const double* QTotal, double time, double dT );

    /**
     * @brief Compute internal force and consistent tangent stiffness.
     * @details Uses the small-strain relation \f$\Delta\boldsymbol{\varepsilon}=\mathbf{B}\,\Delta\mathbf{u}\f$ and
     * integrates
     * \f[
     * \mathbf{K}_e = \sum_{qp} \mathbf{B}^\mathsf{T} \mathbf{C} \mathbf{B}\, J_0 w,\qquad
     * \mathbf{P}_e = \sum_{qp} \mathbf{B}^\mathsf{T} \boldsymbol{\sigma}\, J_0 w.
     * \f]
     * If pNewdT<1, the routine returns early to signal time step reduction.
     * @param QTotal Total displacement vector.
     * @param dQ Incremental displacement.
     * @param Pe Internal force vector (accumulated).
     * @param Ke Tangent stiffness matrix (accumulated).
     * @param time Time data forwarded to materials.
     * @param dT Time increment.
     */
    void computeKernels( const double* QTotal, const double* dQ, double* Pe, double* Ke, double time, double dT );

    /**
     * @brief Compute internal force only (no tangent stiffness).
     * @details Uses the small-strain relation \f$\Delta\boldsymbol{\varepsilon}=\mathbf{B}\,\Delta\mathbf{u}\f$ and
     * integrates
     * \f[
     * \mathbf{P}_e = \sum_{qp} \mathbf{B}^\mathsf{T} \boldsymbol{\sigma}\, J_0 w.
     * \f]
     * @param QTotal Total displacement vector.
     * @param dQ Incremental displacement.
     * @param Pe Internal force vector (accumulated).
     * @param time Time data forwarded to materials.
     * @param dT Time increment.
     */
    void computeKernelsExplicit( const double* QTotal, const double* dQ, double* Pe, double time, double dT );
    /**
     * @brief Compute consistent mass matrix using material density.
     * @details \f$\mathbf{M}_e = \sum_{qp} \rho\, \mathbf{N}^\mathsf{T}\mathbf{N}\, J_0 w\f$.
     */
    void computeConsistentInertia( double* M );

    /**
     * @brief Compute the lumped (diagonal) mass matrix.
     * @details Uses the manifold-based lumping scheme according to
     * Yang, Zheng & Sivaselvan (2017) "A rigorous and unified mass lumping scheme for higher-order
     * elements", CMAME 319, 491-514. The hexa20 weight the derivation below yields is the split
     * assessed by Duczek & Gravenkamp (2019) "Critical assessment of different mass lumping schemes
     * for higher order serendipity finite elements", CMAME 350, 836-897.
     * The lumped mass entries are computed using a weighted shape function
     * \f$\hat{N} = w\,N + (1-w)\,N_\mathrm{lin}\f$,
     * where \f$N\f$ is the high-order shape function and \f$N_\mathrm{lin}\f$ is the corresponding
     * linear (corner-node) shape function on the same element.
     *
     * The blend weight \f$w\f$ cannot be a constant, which is the subtlety here. A corner node's
     * lumped mass is \f$w S^{N}_i + (1-w) S^{\mathrm{lin}}_i\f$, and for a serendipity element
     * \f$S^{N}_i < 0 < S^{\mathrm{lin}}_i\f$, so positivity requires
     * \f[ w < w_\mathrm{max} = \min_i \frac{S^{\mathrm{lin}}_i}{S^{\mathrm{lin}}_i - S^{N}_i}. \f]
     * That limit is element-dependent: \f$0.75\f$ for a quad8, but exactly \f$0.50\f$ for a
     * hexa20, where a hard-coded \f$\tfrac{1}{2}\f$ therefore sits precisely on the boundary and
     * yields an exactly zero corner mass for any regular (affinely-mapped) element.
     *
     * The weight is therefore derived per element by
     * Marmot::FiniteElement::MassLumping::manifoldBlendWeight(), which returns
     * \f$w = \min(\tfrac{1}{2}, \tfrac{2}{3}\,w_\mathrm{max})\f$: exactly \f$\tfrac{1}{2}\f$
     * wherever that is safe -- every 2D serendipity element, and every linear element, where the
     * result does not depend on \f$w\f$ at all -- and \f$\tfrac{1}{3}\f$ for a hexa20. It
     * reproduces -- analytically, and to rounding in floating point -- the values a
     * per-element-type special case would give, without needing one, and
     * the critical time step reads its mass distribution from the same helper so the two cannot
     * disagree.
     *
     * @note The element total is independent of \f$w\f$: the blend only moves mass between the
     * corner and the remaining nodes. An incorrect weight therefore leaves the element mass, and
     * hence the model mass, perfectly correct -- and is invisible to any check on totals.
     */
    void computeLumpedInertia( double* M );

    /**
     * @brief Compute the critical time step for explicit dynamics.
     * @param criticalTimeStep Output parameter for the computed critical time step.
     * @details The estimate is \f$l / c\f$, scaled by the factor
     * Marmot::FiniteElement::MassLumping::timeStepFactorFromMassDistribution() derives from the
     * same lumped mass fractions computeLumpedInertia() assembles: \f$l/c\f$ is the stable
     * increment for an element whose mass is spread uniformly over its nodes, which lumping does
     * not do, and the lightest node sets the highest frequency. \f$l\f$ is twice the smallest
     * singular value of the Jacobian, i.e. the element's smallest physical extent, so a sliver is
     * not mistaken for its volume-equivalent cube. The minimum over all quadrature points is
     * returned.
     *
     * @warning This corrects the mass-distribution and element-distortion parts of the estimate
     * only. It does NOT correct for polynomial order, and that residue is large: \f$l/c\f$ is a
     * linear-element formula, while a quadratic element's highest free eigenfrequency lies well
     * above what it predicts. A convergence study on a 20-node bar settles only around a courant
     * number of 0.1-0.2, i.e. roughly a further factor of five is unaccounted for. Closing that
     * properly wants an eigenvalue-based estimate rather than another factor.
     */
    void computeCriticalTimeStepForExplicitDynamics( double& criticalTimeStep, const double* QTotal );

    /**
     * @brief Compute the internal energy of the element.
     * @param internalEnergy Output parameter for the computed internal energy.
     */
    void computeInternalEnergy( double& internalEnergy );

    /**
     * @brief Access a named state view at a quadrature point.
     * @note Using "sdv" returns the raw material state vector and is deprecated.
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

      else {
        return qp.material->getStateView( stateName, qp.managedStateVars->materialStateVars.data() );
      }
    }

    /** @brief Get physical coordinates at the element center. */
    std::vector< double > getCoordinatesAtCenter();

    /** @brief Get physical coordinates at each quadrature point. */
    std::vector< std::vector< double > > getCoordinatesAtQuadraturePoints();

    /** @brief Number of quadrature points of this element. */
    int getNumberOfQuadraturePoints();
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

    static const vector< vector< string > > nodeFields = [] {
      vector< vector< string > > nodeFields;
      for ( int i = 0; i < nNodes; i++ ) {
        nodeFields.push_back( vector< string >() );
        nodeFields[i].push_back( "displacement" );
      }
      return nodeFields;
    }();

    return nodeFields;
  }

  template < int nDim, int nNodes >
  std::vector< int > DisplacementFiniteElement< nDim, nNodes >::getDofIndicesPermutationPattern()
  {
    static const std::vector< int > permutationPattern = [] {
      std::vector< int > permutationPattern;
      for ( int i = 0; i < nNodes * nDim; i++ )
        permutationPattern.push_back( i );
      return permutationPattern;
    }();

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
    for ( auto& qp : qps ) {
      qp.material = std::unique_ptr< MarmotMaterialHypoElastic >(
        MarmotLibrary::MarmotMaterialHypoElasticFactory::createMaterial( section.materialName,
                                                                         section.materialProperties,
                                                                         section.nMaterialProperties,
                                                                         elLabel ) );

      if ( !qp.material )
        throw std::invalid_argument( MakeString()
                                     << __PRETTY_FUNCTION__
                                     << ": invalid material assigned; cannot cast to MarmotMaterialHypoElastic!" );

      /* Deliberately VOLUME based, and deliberately not the length
       * computeCriticalTimeStepForExplicitDynamics uses. The material's characteristic length is a
       * regularisation length -- it sets the width over which a softening law dissipates its
       * fracture energy -- and the volume-equivalent length is the right measure for that. The
       * stability estimate needs the opposite: the element's SMALLEST extent, because that is what
       * bounds the highest frequency. The two definitions are not interchangeable; do not unify
       * them.
       */
      if constexpr ( nDim == 3 )
        qp.material->setCharacteristicElementLength( std::cbrt( 8 * qp.detJ ) );
      if constexpr ( nDim == 2 )
        qp.material->setCharacteristicElementLength( std::sqrt( 4 * qp.detJ ) );
      if constexpr ( nDim == 1 )
        qp.material->setCharacteristicElementLength( 2 * qp.detJ );
    }
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::assignNodeCoordinates( const double* coordinates )
  {
    ParentGeometryElement::assignNodeCoordinates( coordinates );
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::initializeYourself()
  {
    for ( QuadraturePoint& qp : qps ) {
      const dNdXiSized    dNdXi = this->dNdXi( qp.xi );
      const JacobianSized J     = this->Jacobian( dNdXi );
      const JacobianSized JInv  = J.inverse();
      const dNdXiSized    dNdX  = this->dNdX( dNdXi, JInv );
      qp.detJ                   = J.determinant();
      qp.B                      = this->B( dNdX );

      if constexpr ( nDim == 3 ) {
        qp.J0xW = qp.weight * qp.detJ;
      }
      if constexpr ( nDim == 2 ) {
        const double& thickness = elementProperties[0];
        qp.J0xW                 = qp.weight * qp.detJ * thickness;
      }
      if constexpr ( nDim == 1 ) {
        const double& crossSection = elementProperties[0];
        qp.J0xW                    = qp.weight * qp.detJ * crossSection;
      }
    }
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::computeKernels( const double* QTotal_,
                                                                  const double* dQ_,
                                                                  double*       Pe_,
                                                                  double*       Ke_,
                                                                  double        time,
                                                                  double        dT )
  {
    using namespace Marmot;
    using namespace ContinuumMechanics::VoigtNotation;

    Map< const RhsSized > QTotal( QTotal_ );
    Map< const RhsSized > dQ( dQ_ );
    Map< KeSizedMatrix >  Ke( Ke_ );
    Map< RhsSized >       Pe( Pe_ );

    Voigt  S  = Voigt::Zero();
    Voigt  dE = Voigt::Zero();
    CSized C  = CSized::Zero();

    for ( QuadraturePoint& qp : qps ) {

      const BSized& B             = qp.B;
      dE                          = B * dQ;
      double elasticEnergyDensity = qp.J0xW > 0.0 ? qp.managedStateVars->elasticStrainEnergy / qp.J0xW : 0.0;
      double dissipation          = qp.J0xW > 0.0 ? qp.managedStateVars->dissipation / qp.J0xW : 0.0;

      if constexpr ( nDim == 1 ) {

        MarmotMaterialHypoElastic::state1D  state;
        MarmotMaterialHypoElastic::timeInfo timeInfo;

        // set state info
        state.stress = reduce3DVoigt< ParentGeometryElement::voigtSize >( qp.managedStateVars->stress )( 0 );
        state.elasticEnergyDensity = elasticEnergyDensity;
        state.dissipation          = dissipation;
        state.stateVars            = qp.managedStateVars->materialStateVars.data();

        // set time info
        timeInfo.time = time;
        timeInfo.dT   = dT;
        qp.material->computeUniaxialStress( state, C[0], dE[0], timeInfo );
        Eigen::VectorXd stress1D( 1 );
        stress1D( 0 )               = state.stress;
        qp.managedStateVars->stress = make3DVoigt< ParentGeometryElement::voigtSize >( stress1D );
        elasticEnergyDensity        = state.elasticEnergyDensity;
        dissipation                 = state.dissipation;
      }

      else if constexpr ( nDim == 2 ) {

        if ( sectionType == SectionType::PlaneStress ) {

          MarmotMaterialHypoElastic::state2D  state;
          MarmotMaterialHypoElastic::timeInfo timeInfo;

          // set state info
          state.stress               = reduce3DVoigt< ParentGeometryElement::voigtSize >( qp.managedStateVars->stress );
          state.elasticEnergyDensity = elasticEnergyDensity;
          state.dissipation          = dissipation;
          state.stateVars            = qp.managedStateVars->materialStateVars.data();

          // set time info
          timeInfo.time = time;
          timeInfo.dT   = dT;
          qp.material->computePlaneStress( state, C, dE, timeInfo );
          qp.managedStateVars->stress = make3DVoigt< ParentGeometryElement::voigtSize >( state.stress );
          S                           = state.stress;
          elasticEnergyDensity        = state.elasticEnergyDensity;
          dissipation                 = state.dissipation;
        }

        else if ( sectionType == SectionType::PlaneStrain ) {

          Vector6d dE6 = planeVoigtToVoigt( dE );
          Matrix6d C66;

          // Vector6d S6 = qp.managedStateVars->stress;
          MarmotMaterialHypoElastic::state3D  state;
          MarmotMaterialHypoElastic::timeInfo timeInfo;

          // set state info
          state.stress               = qp.managedStateVars->stress;
          state.elasticEnergyDensity = elasticEnergyDensity;
          state.dissipation          = dissipation;
          state.stateVars            = qp.managedStateVars->materialStateVars.data();

          // set time info
          timeInfo.time = time;
          timeInfo.dT   = dT;
          qp.material->computeStress( state, C66, dE6, timeInfo );
          qp.managedStateVars->stress = state.stress;

          S                    = reduce3DVoigt< ParentGeometryElement::voigtSize >( state.stress );
          C                    = ContinuumMechanics::PlaneStrain::getPlaneStrainTangent( C66 );
          elasticEnergyDensity = state.elasticEnergyDensity;
          dissipation          = state.dissipation;
        }
      }

      else if constexpr ( nDim == 3 ) {
        if ( sectionType == SectionType::Solid ) {

          MarmotMaterialHypoElastic::state3D  state;
          MarmotMaterialHypoElastic::timeInfo timeInfo;

          // set state info
          state.stress               = qp.managedStateVars->stress;
          state.elasticEnergyDensity = elasticEnergyDensity;
          state.dissipation          = dissipation;
          state.stateVars            = qp.managedStateVars->materialStateVars.data();

          // set time info
          timeInfo.time = time;
          timeInfo.dT   = dT;
          qp.material->computeStress( state, C, dE, timeInfo );
          qp.managedStateVars->stress = state.stress;
          S                           = state.stress;
          elasticEnergyDensity        = state.elasticEnergyDensity;
          dissipation                 = state.dissipation;
        }
      }

      qp.managedStateVars->elasticStrainEnergy = elasticEnergyDensity * qp.J0xW;
      qp.managedStateVars->dissipation         = dissipation * qp.J0xW;
      qp.managedStateVars->totalStrainEnergy   = ( elasticEnergyDensity + dissipation ) * qp.J0xW;
      qp.managedStateVars->strain += make3DVoigt< ParentGeometryElement::voigtSize >( dE );

      Ke += B.transpose() * C * B * qp.J0xW;
      Pe += B.transpose() * S * qp.J0xW;
    }
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::computeKernelsExplicit( const double* QTotal_,
                                                                          const double* dQ_,
                                                                          double*       Pe_,
                                                                          double        time,
                                                                          double        dT )
  {
    using namespace Marmot;
    using namespace ContinuumMechanics::VoigtNotation;

    Map< const RhsSized > QTotal( QTotal_ );
    Map< const RhsSized > dQ( dQ_ );
    Map< RhsSized >       Pe( Pe_ );

    Voigt S  = Voigt::Zero();
    Voigt dE = Voigt::Zero();

    for ( QuadraturePoint& qp : qps ) {

      const BSized& B             = qp.B;
      dE                          = B * dQ;
      double elasticEnergyDensity = qp.J0xW > 0.0 ? qp.managedStateVars->elasticStrainEnergy / qp.J0xW : 0.0;
      double dissipation          = qp.J0xW > 0.0 ? qp.managedStateVars->dissipation / qp.J0xW : 0.0;

      if constexpr ( nDim == 1 ) {

        throw std::runtime_error( "Explicit uniaxial stress not implemented yet" );
      }

      else if constexpr ( nDim == 2 ) {

        if ( sectionType == SectionType::PlaneStress ) {

          // TODO. implement explicit plane stress update
          MarmotMaterialHypoElastic::state2D  state;
          MarmotMaterialHypoElastic::timeInfo timeInfo;

          // set state info
          state.stress               = reduce3DVoigt< ParentGeometryElement::voigtSize >( qp.managedStateVars->stress );
          state.elasticEnergyDensity = elasticEnergyDensity;
          state.dissipation          = dissipation;
          state.stateVars            = qp.managedStateVars->materialStateVars.data();

          // set time info
          timeInfo.time = time;
          timeInfo.dT   = dT;
          Matrix3d C    = Matrix3d::Zero();
          qp.material->computePlaneStress( state, C, dE, timeInfo );
          qp.managedStateVars->stress = make3DVoigt< ParentGeometryElement::voigtSize >( state.stress );
          S                           = state.stress;
          elasticEnergyDensity        = state.elasticEnergyDensity;
          dissipation                 = state.dissipation;
        }

        else if ( sectionType == SectionType::PlaneStrain ) {

          Vector6d dE6 = planeVoigtToVoigt( dE );

          // Vector6d S6 = qp.managedStateVars->stress;
          MarmotMaterialHypoElastic::state3D  state;
          MarmotMaterialHypoElastic::timeInfo timeInfo;

          // set state info
          state.stress               = qp.managedStateVars->stress;
          state.elasticEnergyDensity = elasticEnergyDensity;
          state.dissipation          = dissipation;
          state.stateVars            = qp.managedStateVars->materialStateVars.data();

          // set time info
          timeInfo.time = time;
          timeInfo.dT   = dT;
          qp.material->computeStressExplicit( state, dE6, timeInfo );
          qp.managedStateVars->stress = state.stress;

          S                    = reduce3DVoigt< ParentGeometryElement::voigtSize >( state.stress );
          elasticEnergyDensity = state.elasticEnergyDensity;
          dissipation          = state.dissipation;
        }
      }

      else if constexpr ( nDim == 3 ) {
        if ( sectionType == SectionType::Solid ) {

          MarmotMaterialHypoElastic::state3D  state;
          MarmotMaterialHypoElastic::timeInfo timeInfo;

          // set state info
          state.stress               = qp.managedStateVars->stress;
          state.elasticEnergyDensity = elasticEnergyDensity;
          state.dissipation          = dissipation;
          state.stateVars            = qp.managedStateVars->materialStateVars.data();

          // set time info
          timeInfo.time = time;
          timeInfo.dT   = dT;
          qp.material->computeStressExplicit( state, dE, timeInfo );
          qp.managedStateVars->stress = state.stress;
          S                           = state.stress;
          elasticEnergyDensity        = state.elasticEnergyDensity;
          dissipation                 = state.dissipation;
        }
      }

      qp.managedStateVars->elasticStrainEnergy = elasticEnergyDensity * qp.J0xW;
      qp.managedStateVars->dissipation         = dissipation * qp.J0xW;
      qp.managedStateVars->totalStrainEnergy   = ( elasticEnergyDensity + dissipation ) * qp.J0xW;
      qp.managedStateVars->strain += make3DVoigt< ParentGeometryElement::voigtSize >( dE );

      Pe += B.transpose() * S * qp.J0xW;
    }
  }
  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::setInitialConditions( StateTypes state, const double* values )
  {
    switch ( state ) {
    case MarmotElement::MarmotMaterialInitialization: {
      for ( QuadraturePoint& qp : qps ) {
        qp.material->initializeYourself( qp.managedStateVars->materialStateVars.data(),
                                         qp.managedStateVars->materialStateVars.size() );
      }
      break;
    }
    case MarmotElement::GeostaticStress: {
      if ( nDim >= 2 )
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
      break;
    }
    case MarmotElement::MarmotMaterialStateVars: {
      throw std::invalid_argument( "Please use initializeStateVars directly on material" );
    }
    default: throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid initial condition" );
    }
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::computeDistributedLoad( MarmotElement::DistributedLoadTypes loadType,
                                                                          double*                             P,
                                                                          double*                             K,
                                                                          const int     elementFace,
                                                                          const double* load,
                                                                          const double* QTotal,
                                                                          double        time,
                                                                          double        dT )
  {
    Map< RhsSized > fU( P );

    switch ( loadType ) {

    case MarmotElement::Pressure: {
      const double p = load[0];

      FiniteElement::BoundaryElement boundaryEl( this->shape, elementFace, nDim, this->coordinates );

      VectorXd Pk = -p * boundaryEl.computeSurfaceNormalVectorialLoadVector();

      if ( nDim == 2 )
        Pk *= elementProperties[0]; // thickness

      boundaryEl.assembleIntoParentVectorial( Pk, fU );

      break;
    }
    case MarmotElement::SurfaceTraction: {

      FiniteElement::BoundaryElement boundaryEl( this->shape, elementFace, nDim, this->coordinates );

      const XiSized tractionVector( load );

      auto Pk = boundaryEl.computeVectorialLoadVector( tractionVector );
      if ( nDim == 2 )
        Pk *= elementProperties[0]; // thickness
      boundaryEl.assembleIntoParentVectorial( Pk, fU );

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
                                                                    double        time,
                                                                    double        dT )
  {
    Map< RhsSized >                              Pe( P_ );
    const Map< const Matrix< double, nDim, 1 > > f( load );

    for ( const auto& qp : qps )
      Pe += this->NB( this->N( qp.xi ) ).transpose() * f * qp.J0xW;
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::computeConsistentInertia( double* M )
  {
    Map< KeSizedMatrix > Me( M );
    Me.setZero();

    for ( const auto& qp : qps ) {
      const auto   N_  = this->NB( this->N( qp.xi ) );
      const double rho = qp.material->getDensity( qp.managedStateVars->materialStateVars.data() );
      Me += N_.transpose() * N_ * qp.J0xW * rho;
    }
  }
  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::computeLumpedInertia( double* M )
  {
    Map< RhsSized > LMM( M );
    LMM.setZero();

    /* Row sums of the consistent mass matrix, from this element's own shape functions and from the
     * linear (corner-node) shape functions of the same element. Density-free: the blend weight and
     * the mass distribution are properties of the element geometry alone, and deriving them in one
     * place is what keeps this function and computeCriticalTimeStepForExplicitDynamics consistent
     * with each other -- a time step derived from a different mass distribution than the one
     * actually assembled is exactly the kind of inconsistency that surfaces as an unexplained
     * instability rather than as a clean failure.
     */
    constexpr int   nNodesLinear     = ( 1 << nDim );
    auto            linGeometryEl    = MarmotGeometryElement< nDim, nNodesLinear >();
    Eigen::VectorXd rowSumsHighOrder = Eigen::VectorXd::Zero( nNodes );
    Eigen::VectorXd rowSumsLinear    = Eigen::VectorXd::Zero( nNodesLinear );
    for ( const auto& qp : qps ) {
      rowSumsHighOrder += Eigen::VectorXd( this->N( qp.xi ) ) * qp.J0xW;
      rowSumsLinear += Eigen::VectorXd( linGeometryEl.N( qp.xi ) ) * qp.J0xW;
    }
    const double weight = FiniteElement::MassLumping::manifoldBlendWeight( rowSumsHighOrder, rowSumsLinear );

    for ( const auto& qp : qps ) {
      const auto N_    = this->N( qp.xi );
      const auto N_lin = linGeometryEl.N( qp.xi );

      VectorXd N_weighted = weight * ( N_ );
      N_weighted.head( nNodesLinear ) += ( 1.0 - weight ) * N_lin;

      const double rho = qp.material->getDensity( qp.managedStateVars->materialStateVars.data() );
      VectorXd     m_  = N_weighted * qp.J0xW * rho;
      for ( int i = 0; i < nNodes; i++ ) {
        for ( int d = 0; d < nDim; d++ )
          LMM( i * nDim + d ) += m_( i );
      }
    }
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::computeCriticalTimeStepForExplicitDynamics( double& criticalTimeStep,
                                                                                              const double* QTotal )
  {

    /* The l / c estimate below assumes the element's mass is spread UNIFORMLY over its nodes, each
     * carrying 1 / nNodes of it. The lumping scheme computeLumpedInertia applies does not do that:
     * under the manifold-based blend a hexa20 corner node carries less than the uniform share, and
     * the lightest node sets the highest frequency (omega = sqrt( k / m )), so the stable increment
     * scales with sqrt( m_min / m_uniform ).
     *
     * Ignoring it puts the default courant number of 0.8 ABOVE the true limit for every 20-node
     * element. That does not present as a marginally noisy run but as violent exponential
     * divergence -- and only where the material provides no damping, so a 20-node run on a
     * viscously regularised material can sit just above the limit and look perfectly healthy. That
     * is worse than a clean failure, because it makes the fault look model-specific.
     *
     * The fractions come from the same helper computeLumpedInertia uses, so the two cannot drift
     * apart. Exactly 1 for a linear element on a regular mesh, so nothing changes there, and
     * slightly below 1 for a distorted element -- which is correct, a distorted element does have a
     * tighter limit than its volume alone suggests.
     */
    constexpr int   nNodesLinear     = ( 1 << nDim );
    auto            linGeometryEl    = MarmotGeometryElement< nDim, nNodesLinear >();
    Eigen::VectorXd rowSumsHighOrder = Eigen::VectorXd::Zero( nNodes );
    Eigen::VectorXd rowSumsLinear    = Eigen::VectorXd::Zero( nNodesLinear );
    for ( const auto& qp : qps ) {
      rowSumsHighOrder += Eigen::VectorXd( this->N( qp.xi ) ) * qp.J0xW;
      rowSumsLinear += Eigen::VectorXd( linGeometryEl.N( qp.xi ) ) * qp.J0xW;
    }
    const double weight = FiniteElement::MassLumping::manifoldBlendWeight( rowSumsHighOrder, rowSumsLinear );
    const double lumpedMassTimeStepFactor = FiniteElement::MassLumping::timeStepFactorFromMassDistribution(
      FiniteElement::MassLumping::manifoldMassFractions( rowSumsHighOrder, rowSumsLinear, weight ) );

    criticalTimeStep = std::numeric_limits< double >::max();
    for ( const auto& qp : qps ) {
      /* The characteristic length has to be the element's SMALLEST physical extent, not a
       * volume-averaged one. cbrt( 8 * detJ ) and its lower-dimensional analogues are volume
       * based: for a sliver -- thin in one direction but not the others -- the volume stays
       * moderate while the thin dimension collapses, so they OVERESTIMATE the length and hence
       * the stable time step. Refining a distorted parent element is precisely how slivers are
       * produced, so an h-adaptive explicit run integrates its most distorted elements above
       * their stability limit and diverges thousands of increments later.
       *
       * The Jacobian maps the natural cube [-1,1]^nDim onto the element, so twice its smallest
       * singular value IS that smallest physical extent. For a well-shaped element this
       * reproduces the previous expressions exactly -- a cube of side h gives h either way --
       * so the estimate is tightened only where it was previously wrong.
       */
      const JacobianSized J_                          = this->Jacobian( this->dNdXi( qp.xi ) );
      const double        characteristicElementLength = 2.0 *
                                                 Eigen::JacobiSVD< JacobianSized >( J_ ).singularValues().minCoeff();

      MarmotMaterialHypoElastic::state3D state( qp.managedStateVars->stress,
                                                qp.managedStateVars->elasticStrainEnergy / qp.J0xW,
                                                qp.managedStateVars->dissipation / qp.J0xW,
                                                qp.managedStateVars->materialStateVars.data() );

      const double c = qp.material->getMaximumWaveSpeed( state );
      if ( c <= 0.0 )
        throw std::runtime_error( "Non-positive wave speed encountered in computeCriticalTimeStepForExplicitDynamics" );
      const double dt = lumpedMassTimeStepFactor * characteristicElementLength / c;
      if ( dt < criticalTimeStep )
        criticalTimeStep = dt;
    }
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::computeInternalEnergy( double& internalEnergy )
  {
    internalEnergy = 0.0;
    for ( const auto& qp : qps ) {
      internalEnergy += qp.managedStateVars->totalStrainEnergy;
    }
  }

  template < int nDim, int nNodes >
  std::vector< double > DisplacementFiniteElement< nDim, nNodes >::getCoordinatesAtCenter()
  {
    std::vector< double > coords( nDim );

    Eigen::Map< XiSized > coordsMap( &coords[0] );
    const auto            centerXi = XiSized::Zero();
    coordsMap                      = this->NB( this->N( centerXi ) ) * this->coordinates;
    return coords;
  }

  template < int nDim, int nNodes >
  std::vector< std::vector< double > > DisplacementFiniteElement< nDim, nNodes >::getCoordinatesAtQuadraturePoints()
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
  int DisplacementFiniteElement< nDim, nNodes >::getNumberOfQuadraturePoints()
  {
    return qps.size();
  }
} // namespace Marmot::Elements
