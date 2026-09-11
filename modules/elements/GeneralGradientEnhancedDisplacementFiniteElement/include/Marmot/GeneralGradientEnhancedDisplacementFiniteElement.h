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
#include "Marmot/MarmotMaterialGeneralGradientEnhancedHypoElastic.h"
#include "Marmot/MarmotMaterialGeneralGradientEnhancedHypoElasticFactory.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include <memory>
#include <vector>

using namespace Marmot;
using namespace Eigen;

namespace Marmot::Elements {

  /**
   * @class Marmot::Elements::GeneralGradientEnhancedDisplacementFiniteElement
   * @tparam nDim Number of spatial dimensions (1, 2, or 3).
   * @tparam nNodes Number of element nodes for the displacement field.
   * @tparam nNonlocalVariables Number of nonlocal variables.
   * @tparam nNonLocalNodes Number of nonlocal element nodes. Usually the order of the nonlocal interpolation is similar
   * or one order lower than the local interpolation, so nNonLocalNodes <= nNodes.
   * @brief General Gradient Enhanced Displacement-based finite element template.
   * @details Uses linearized kinematics (small strains) and supports 1D, 2D, and 3D
   * formulations via a section assumption. Holds quadrature-point state and
   * delegates constitutive updates to Marmot materials while assembling element
   * residuals and tangents, taking non-local effects into account.
   */
  template < int nDim, int nNodes, int nNonlocalVariables = 1, int nNonLocalNodes = nNodes >
  class GeneralGradientEnhancedDisplacementFiniteElement : public MarmotElement {

  public:
    /**
     * @brief Kinematic section assumption used by the element.
     * @details Controls which constitutive call is performed at quadrature points.
     * - PlaneStress
     * - PlaneStrain
     * - Solid
     */
    enum SectionType {
      PlaneStress,
      PlaneStrain,
      Solid,
    };

    static constexpr int nDofPerNodeU = nDim;
    static constexpr int nDofPerNodeK = nNonlocalVariables;

    static constexpr int sizeLoadVector = nNodes * nDofPerNodeU + nNonLocalNodes * nDofPerNodeK;
    static constexpr int nCoordinates   = nNodes * nDim;

    static constexpr int sizeDoFU = nNodes * nDofPerNodeU;
    static constexpr int sizeDoFK = nNonLocalNodes * nDofPerNodeK;

    MarmotGeometryElement< nDim, nNodes >         localGeometryElement;
    MarmotGeometryElement< nDim, nNonLocalNodes > nonLocalGeometryElement;

    using ParentGeometryElement = MarmotGeometryElement< nDim, nNodes >;
    using JacobianSized         = typename ParentGeometryElement::JacobianSized;
    using NSized                = typename ParentGeometryElement::NSized;
    using dNdXiSized            = typename ParentGeometryElement::dNdXiSized;
    using BSized                = typename ParentGeometryElement::BSized;
    using XiSized               = typename ParentGeometryElement::XiSized;
    using RhsSized              = Matrix< double, sizeLoadVector, 1 >;
    using KeSizedMatrix         = Matrix< double, sizeLoadVector, sizeLoadVector >;
    using USizedVector          = Matrix< double, sizeDoFU, 1 >;
    using KSizedVector          = Matrix< double, sizeDoFK, 1 >;
    using CSized                = Matrix< double, ParentGeometryElement::voigtSize, ParentGeometryElement::voigtSize >;
    using Voigt                 = Matrix< double, ParentGeometryElement::voigtSize, 1 >;

    using ParentGeometryElementK = MarmotGeometryElement< nDim, nNonLocalNodes >;
    using JacobianSizedK         = typename ParentGeometryElementK::JacobianSized;
    using NSizedK                = typename ParentGeometryElementK::NSized;
    using dNdXiSizedK            = typename ParentGeometryElementK::dNdXiSized;
    using BSizedK                = typename ParentGeometryElementK::BSized;

    /** Element-level properties (e.g., thickness for 2D, area for 1D). */
    Map< const VectorXd > elementProperties;
    /** Element label (ID) used for logging and material creation. */
    const int elLabel;

    /** Section assumption applied by this element instance. */
    const SectionType sectionType;

    /**
     * @brief Data and state associated with a quadrature point.
     * @details Holds parent coordinates, integration weight, Jacobian determinant, shape functions N and their
     * gradients dNdX, kinematic (strain-displacement) B-matrices for both local and non-local evaluations, and a
     * material instance with managed state variables.
     */
    struct QuadraturePoint {

      const XiSized xi;
      const double  weight;

      double      detJ;
      double      J0xW;
      NSized      N;
      dNdXiSized  dNdX;
      BSized      B;
      NSizedK     N_K;
      dNdXiSizedK dNdX_K;
      BSizedK     B_K;

      /**
       * @brief Manager for per-quadrature-point state variables.
       * @details Provides named accessors to stress \f$\sigma\f$, strain \f$\varepsilon\f$
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

      std::unique_ptr< MarmotMaterialGeneralGradientEnhancedHypoElastic< nNonlocalVariables > > material;

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
        /* material->assignStateVars( managedStateVars->materialStateVars.data(), */
        /*                            managedStateVars->materialStateVars.size() ); */
      }

      QuadraturePoint( XiSized xi, double weight )
        : xi( xi ),
          weight( weight ),
          detJ( 0.0 ),
          J0xW( 0.0 ),
          N( NSized::Zero() ),
          dNdX( dNdXiSized::Zero() ),
          B( BSized::Zero() ),
          N_K( NSizedK::Zero() ),
          dNdX_K( dNdXiSizedK::Zero() ),
          B_K( BSizedK::Zero() ){};
    };

    /// Quadrature points owned by the element (one per integration point).
    std::vector< QuadraturePoint > qps;

    /**
     * @brief Construct element with ID, quadrature rule and section assumption.
     * @param elementID Unique element label.
     * @param integrationType Integration (quadrature) rule.
     * @param sectionType Section assumption (1D/2D/3D).
     */
    GeneralGradientEnhancedDisplacementFiniteElement( int                                         elementID,
                                                      FiniteElement::Quadrature::IntegrationTypes integrationType,
                                                      SectionType                                 sectionType );

    /** @brief Total number of required state variables for this element (sum over all quadrature points). */
    int getNumberOfRequiredStateVars();

    /** @brief Node-level fields exposed by the element. Returns ["displacement", "nonlocal field", "strain symmetric"]
     * for the respective nodes. */
    std::vector< std::vector< std::string > > getNodeFields();

    /** @brief The permutation pattern for the residual vector and the stiffness matrix to aggregate all entries in
     * order to resemble the defined fields nodewise.  */
    std::vector< int > getDofIndicesPermutationPattern();

    /** @brief Number of nodes of this element type. */
    int getNNodes() { return nNodes; }

    /** @brief Number of spatial dimensions. */
    int getNSpatialDimensions() { return nDim; }

    /** @brief Number of degrees of freedom per element. */
    int getNDofPerElement() { return sizeLoadVector; }

    /** @brief Geometric shape of the element (as reported by the local geometry element). */
    std::string getElementShape() { return localGeometryElement.getElementShape(); }

    /** @brief Map the provided element state vector to all quadrature points. */
    void assignStateVars( double* stateVars, int nStateVars );

    /** @brief Assign element properties (e.g., thickness in 2D, area in 1D). */
    void assignProperty( const ElementProperties& marmotElementProperty );

    /** @brief Assign material section and instantiate per-quadrature-point materials. */
    void assignProperty( const MarmotMaterialSection& marmotElementProperty );

    /** @brief Provide nodal coordinates to local and non-local geometry elements. */
    void assignNodeCoordinates( const double* coordinates );

    /** @brief Precompute geometry-related quantities at quadrature points (B, detJ, J0xW, N_K, dNdX_K). */
    void initializeYourself();

    /**
     * @brief Initialize state or materials.
     * @param state MarmotMaterialInitialization, GeostaticStress or MarmotMaterialStateVars.
     * @param values For GeostaticStress: [sigmaY(z1), y1, sigmaY(z2), y2, kx, kz].
     */
    void setInitialConditions( StateTypes state, const double* values );

    /**
     * @brief Assemble distributed surface loads on a boundary face.
     * @details Pressure contributions are integrated on the boundary \f$\Gamma_e\f$:
     * \f[
     * \mathbf{\fextp} = - \int_{\Gamma_e} p\, \mathbf{N}^\mathsf{T} \mathbf{n}\, \mathrm{d}\Gamma\, .
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
     * @details Integrates \f$\mathbf{\fextb} = \int_{\Omega_e} \mathbf{N}^\mathsf{T} \mathbf{f}\,
     * \mathrm{d}\Omega\f$.
     */
    void computeBodyForce( double* P,
                           double* K,

                           const double* load,
                           const double* QTotal,
                           double        time,
                           double        dT );

    /**
     * @brief Compute internal force and consistent tangent stiffness.
     * @details Uses the small-strain relation \f$\Delta \eps = \mathbf{B}\, \Delta \mathbf{\qu}\f$,
     * interpolates the non-local variables as \f$ \knl = \Nk\, \mathbf{\qk}\f$. Internal forces and consistent tangents
     * are evaluated by Gauss quadrature:
     * \f[
     * \mathbf{K}_e = \begin{bmatrix} \mathbf{K}_{uu} & \mathbf{K}_{uk} \\ \mathbf{K}_{ku} & \mathbf{K}_{kk}
     * \end{bmatrix},\qquad
     * \mathbf{\fuint} = \sum_{qp} \mathbf{B}^\mathsf{T}\, \sig\, J_0\, w_{qp}\, .
     * \f]
     * \f[
     * \mathbf{\fk} = \sum_{qp} \left (\mathbf{\Nk}^\mathsf{T}\, \knl\, + c\, \partial_\mathbf{x}
     * \mathbf{\Nk}^\mathsf{T}\,\partial_\mathbf{x} \mathbf{\Nk}\, \mathbf{\qk} - \mathbf{\Nk}^\mathsf{T}\, \kl \right )
     * J_0\, w_{qp} \, .
     * \f]
     * The stiffness submatrices are evaluated using the following expressions:
     * \f[
     * \mathbf{K}_e = \begin{bmatrix} \mathbf{K}_{uu} & \mathbf{K}_{uk} \\ \mathbf{K}_{ku} & \mathbf{K}_{kk}
     * \end{bmatrix},\qquad
     * \mathbf{\fuint} = \sum_{qp} \mathbf{B}^\mathsf{T}\, \sig\, J_0\, w_{qp}\, .
     * \f]
     * \f[
     * \mathbf{\fk} = \sum_{qp} \left (\mathbf{\Nk}^\mathsf{T}\, \knl\, + c\, \partial_\mathbf{x}
     * \mathbf{\Nk}^\mathsf{T}\,\partial_\mathbf{x} \mathbf{\Nk}\, \mathbf{\qk} - \mathbf{\Nk}^\mathsf{T}\, \kl \right )
     * J_0\, w_{qp} \, .
     * \f]
     * If pNewdT<1, the routine returns early to signal time step reduction.
     * @param QTotal Total displacement vector in field-wise format: \f$\mathbf{q} = [\mathbf{\qu},
     * \mathbf{\qk}]^\mathsf{T}\f$.
     * @param dQ Incremental displacement.
     * @param Pe Internal force vector (accumulated).
     * @param Ke Tangent stiffness matrix (accumulated).
     * @param time Time data forwarded to materials.
     * @param dT Time increment.
     */
    void computeKernels( const double* QTotal, const double* dQ, double* Pe, double* Ke, double time, double dT );

    /**
     * @brief Compute internal force without tangents.
     * @details Uses the small-strain relation \f$\Delta \eps = \mathbf{B}\, \Delta \mathbf{\qu}\f$,
     * interpolates the non-local variables as \f$ \knl = \Nk\, \mathbf{\qk}\f$. Internal forces
     * are evaluated by Gauss quadrature:
     * \mathbf{\fk} = \sum_{qp} \left (\mathbf{\Nk}^\mathsf{T}\, \knl\, + c\, \partial_\mathbf{x}
     * \mathbf{\Nk}^\mathsf{T}\,\partial_\mathbf{x} \mathbf{\Nk}\, \mathbf{\qk} - \mathbf{\Nk}^\mathsf{T}\, \kl \right )
     * J_0\, w_{qp} \, .
     * \f]
     * If pNewdT<1, the routine returns early to signal time step reduction.
     * @param QTotal Total displacement vector in field-wise format: \f$\mathbf{q} = [\mathbf{\qu},
     * \mathbf{\qk}]^\mathsf{T}\f$.
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
     * @brief Compute the critical time step for explicit dynamics based on
     * the dilatational wave speed and the element size.
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
     * @brief Compute the internal energy of the element by summing the strain energy contributions from all quadrature
     * points.
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

      if ( stateName == "sdv" ) {
        MarmotJournal::warningToMSG( MakeString()
                                     << __PRETTY_FUNCTION__
                                     << " on 'sdv' is discouraged and deprecated, please use precise state name" );
        return { qp.managedStateVars->materialStateVars.data(),
                 static_cast< int >( qp.managedStateVars->materialStateVars.size() ) };
      }

      if ( qp.managedStateVars->contains( stateName ) ) {
        return qp.managedStateVars->getStateView( stateName );
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

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVariables, nNonLocalNodes >::
    GeneralGradientEnhancedDisplacementFiniteElement( int                                         elementID,
                                                      FiniteElement::Quadrature::IntegrationTypes integrationType,
                                                      SectionType                                 sectionType )
    : elementProperties( nullptr, 0 ), elLabel( elementID ), sectionType( sectionType )
  {
    for ( const auto& qpInfo :
          FiniteElement::Quadrature::getGaussPointInfo( localGeometryElement.shape, integrationType ) ) {
      qps.emplace_back( qpInfo.xi, qpInfo.weight );
    }
  }

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  int GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVariables, nNonLocalNodes >::
    getNumberOfRequiredStateVars()
  {
    return qps[0].getNumberOfRequiredStateVars() * qps.size();
  }

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  void GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVariables, nNonLocalNodes >::
    assignStateVars( double* stateVars, int nStateVars )
  {
    const int nQpStateVars = nStateVars / qps.size();

    for ( size_t i = 0; i < qps.size(); i++ ) {
      auto&   qp          = qps[i];
      double* qpStateVars = stateVars + ( i * nQpStateVars );
      qp.assignStateVars( qpStateVars, nQpStateVars );
    }
  }

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  void GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVariables, nNonLocalNodes >::
    assignProperty( const ElementProperties& elementPropertiesInfo )
  {
    new ( &elementProperties )
      Map< const VectorXd >( elementPropertiesInfo.elementProperties, elementPropertiesInfo.nElementProperties );
  }

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  void GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVariables, nNonLocalNodes >::
    assignProperty( const MarmotMaterialSection& section )
  {
    for ( auto& qp : qps ) {
      qp.material = std::unique_ptr< MarmotMaterialGeneralGradientEnhancedHypoElastic< nNonlocalVariables > >(
        MarmotLibrary::MarmotMaterialGeneralGradientEnhancedHypoElasticFactory< nNonlocalVariables >::
          createMaterial( section.materialName, section.materialProperties, section.nMaterialProperties, elLabel ) );
    }
  }

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  std::vector< std::vector< std::string > > GeneralGradientEnhancedDisplacementFiniteElement<
    nDim,
    nNodes,
    nNonlocalVariables,
    nNonLocalNodes >::getNodeFields()
  {
    using namespace std;

    static const vector< vector< string > > nodeFields = [] {
      vector< vector< string > > nodeFields;
      for ( int i = 0; i < nNodes; i++ ) {
        nodeFields.push_back( vector< string >() );
        nodeFields[i].push_back( "displacement" );
        if ( i < nNonLocalNodes ) {
          if constexpr ( nNonlocalVariables == 6 )
            nodeFields[i].push_back( "strain symmetric" );
          else {
            nodeFields[i].push_back( "nonlocal damage" );
            for ( int j = 1; j < nNonlocalVariables; j++ )
              nodeFields[i].push_back( "nonlocal damage " + to_string( j + 1 ) );
          }
        }
      }
      return nodeFields;
    }();
    return nodeFields;
  }

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  std::vector< int > GeneralGradientEnhancedDisplacementFiniteElement<
    nDim,
    nNodes,
    nNonlocalVariables,
    nNonLocalNodes >::getDofIndicesPermutationPattern()
  {
    static const std::vector< int > permutationPattern = [] {
      std::vector< int > permutationPattern;
      for ( int i = 0; i < nNodes; i++ )
        for ( int j = 0; j < nDim; j++ )
          permutationPattern.push_back( i * nDim + nNonlocalVariables * ( i < nNonLocalNodes ? i : nNonLocalNodes ) +
                                        j );
      for ( int j = 0; j < nNonlocalVariables; j++ )
        for ( int i = 0; i < nNonLocalNodes; i++ )
          permutationPattern.push_back( i * ( nDim + nNonlocalVariables ) + nDim + j );
      return permutationPattern;
    }();
    return permutationPattern;
  }

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  void GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVariables, nNonLocalNodes >::
    assignNodeCoordinates( const double* coordinates )
  {
    localGeometryElement.assignNodeCoordinates( coordinates );
    nonLocalGeometryElement.assignNodeCoordinates( coordinates );
    // This assumes that the corner nodes (vertices) are listed before the mid-edge nodes!
  }

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  void GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVariables, nNonLocalNodes >::
    initializeYourself()
  {
    for ( QuadraturePoint& qp : qps ) {
      const auto          dNdXi = localGeometryElement.dNdXi( qp.xi );
      const JacobianSized J     = localGeometryElement.Jacobian( dNdXi );
      const JacobianSized JInv  = J.inverse();
      qp.detJ                   = J.determinant();
      qp.N                      = localGeometryElement.N( qp.xi );
      qp.dNdX                   = localGeometryElement.dNdX( dNdXi, JInv );
      qp.B                      = localGeometryElement.B( qp.dNdX );

      const auto           dNdXi_K = nonLocalGeometryElement.dNdXi( qp.xi );
      const JacobianSizedK J_K     = nonLocalGeometryElement.Jacobian( dNdXi_K );
      const JacobianSizedK JInv_K  = J_K.inverse();
      qp.N_K                       = nonLocalGeometryElement.N( qp.xi );
      qp.dNdX_K                    = nonLocalGeometryElement.dNdX( dNdXi_K, JInv_K );
      qp.B_K                       = nonLocalGeometryElement.B( qp.dNdX_K );

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

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  void GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVariables, nNonLocalNodes >::
    computeKernels( const double* QTotal_, const double* dQ_, double* Pe_, double* Ke_, double time, double dT )
  {

    Map< const RhsSized > QTotal( QTotal_ );
    Map< const RhsSized > dQ( dQ_ );
    Map< KeSizedMatrix >  Ke( Ke_ );
    Map< RhsSized >       Pe( Pe_ );

    // incremental nodal displacements Q and nodal internal parameters qK
    const Ref< const USizedVector > dQU( dQ.head( sizeDoFU ) );
    const Ref< const KSizedVector > dQK( dQ.tail( sizeDoFK ) );
    const Ref< const USizedVector > qU( QTotal.head( sizeDoFU ) );
    const Ref< const KSizedVector > qK( QTotal.tail( sizeDoFK ) );

    // substiffness matrices
    // [k_qq    k_qQK
    //  k_QKQ   k_QKQK ]
    Ref< Matrix< double, sizeDoFU, sizeDoFU > > kUU( Ke.topLeftCorner( sizeDoFU, sizeDoFU ) );
    Ref< Matrix< double, sizeDoFU, sizeDoFK > > kUK( Ke.topRightCorner( sizeDoFU, sizeDoFK ) );
    Ref< Matrix< double, sizeDoFK, sizeDoFU > > kKU( Ke.bottomLeftCorner( sizeDoFK, sizeDoFU ) );
    Ref< Matrix< double, sizeDoFK, sizeDoFK > > kKK( Ke.bottomRightCorner( sizeDoFK, sizeDoFK ) );

    // Righthandside
    Ref< USizedVector > fU( Pe.head( sizeDoFU ) );
    Ref< KSizedVector > fK( Pe.tail( sizeDoFK ) );

    using namespace Marmot;
    using namespace ContinuumMechanics::VoigtNotation;

    using response  = typename MarmotMaterialGeneralGradientEnhancedHypoElastic< nNonlocalVariables >::response;
    using tangents  = typename MarmotMaterialGeneralGradientEnhancedHypoElastic< nNonlocalVariables >::tangents;
    using increment = typename MarmotMaterialGeneralGradientEnhancedHypoElastic< nNonlocalVariables >::increment;

    for ( size_t i = 0; i < this->qps.size(); i++ ) {

      QuadraturePoint& qp = this->qps[i];

      const BSized&      B      = qp.B;
      const NSizedK&     N_K    = qp.N_K;
      const dNdXiSizedK& dNdX_K = qp.dNdX_K;

      Voigt dE = B * dQU;

      // delta of _K at Gausspoint
      Vector< double, nNonlocalVariables > K;
      Vector< double, nNonlocalVariables > dK;

      for ( size_t n = 0; n < nNonlocalVariables; n++ ) {
        Eigen::Index idx = static_cast< Eigen::Index >( n ) * static_cast< Eigen::Index >( nNonLocalNodes );
        K( n )           = N_K * qK.segment( idx, nNonLocalNodes );
        dK( n )          = N_K * dQK.segment( idx, nNonLocalNodes );
      }

      response  res;
      tangents  tan;
      increment inc;
      if constexpr ( nDim == 2 ) {
        Vector6d dE6             = ContinuumMechanics::VoigtNotation::planeVoigtToVoigt( dE );
        res.stress               = qp.managedStateVars->stress;
        res.elasticEnergyDensity = qp.managedStateVars->elasticStrainEnergy / qp.J0xW;
        res.dissipation          = qp.managedStateVars->dissipation / qp.J0xW;
        res.stateVars            = qp.managedStateVars->materialStateVars.data();
        inc                      = { dE6, K, dK, time, dT };
        CSized C                 = CSized::Zero();
        Voigt  S                 = Voigt::Zero();

        if ( sectionType == SectionType::PlaneStress ) {
          qp.material->computePlaneStress( res, tan, inc );
          S = ContinuumMechanics::VoigtNotation::voigtToPlaneVoigt( res.stress );
          C = ContinuumMechanics::PlaneStress::getPlaneStressTangent( tan.dStressddStrain );
        }
        else if ( sectionType == SectionType::PlaneStrain ) {
          qp.material->computeStress( res, tan, inc );
          S = ContinuumMechanics::VoigtNotation::voigtToPlaneVoigt( res.stress );
          C = ContinuumMechanics::PlaneStrain::getPlaneStrainTangent( tan.dStressddStrain );
        }
        else {
          throw std::invalid_argument( "Invalid section type for 2D element, expected PlaneStress or PlaneStrain" );
        }

        fU += B.transpose() * S * qp.J0xW;
        kUU += B.transpose() * C * B * qp.J0xW;

        for ( int n = 0; n < nNonlocalVariables; n++ ) {
          Eigen::Index idx = n * nNonLocalNodes;
          fK.segment( idx, nNonLocalNodes ) += ( N_K.transpose() * K( n ) +
                                                 res.c( n ) * dNdX_K.transpose() * dNdX_K *
                                                   qK.segment( idx, nNonLocalNodes ) -
                                                 N_K.transpose() * res.KLocal( n ) ) *
                                               qp.J0xW;
          const auto dSdK         = ContinuumMechanics::VoigtNotation::voigtToPlaneVoigt( tan.dStressddK.col( n ) );
          const auto dK_Local_dDE = ContinuumMechanics::VoigtNotation::voigtToPlaneVoigt(
            tan.dKLocalddStrain.row( n ).transpose() );

          kUK.block( 0, idx, sizeDoFU, nNonLocalNodes ) += B.transpose() * dSdK * N_K * qp.J0xW;
          kKU.block( idx, 0, nNonLocalNodes, sizeDoFU ) += N_K.transpose() * -dK_Local_dDE.transpose() * B * qp.J0xW;
          kKK.block( idx, idx, nNonLocalNodes, nNonLocalNodes ) += ( N_K.transpose() * N_K +
                                                                     res.c( n ) * dNdX_K.transpose() * dNdX_K +
                                                                     tan.dcddK( n ) * dNdX_K.transpose() * dNdX_K *
                                                                       qK.segment( idx, nNonLocalNodes ) * N_K -
                                                                     N_K.transpose() * tan.dKLocalddK( n, n ) * N_K ) *
                                                                   qp.J0xW;
        }
      }

      else if ( nDim == 3 ) {
        if ( sectionType == Solid ) {

          res.stress               = qp.managedStateVars->stress;
          res.elasticEnergyDensity = qp.managedStateVars->elasticStrainEnergy / qp.J0xW;
          res.dissipation          = qp.managedStateVars->dissipation / qp.J0xW;
          res.stateVars            = qp.managedStateVars->materialStateVars.data();
          inc                      = { dE, K, dK, time, dT };
          qp.material->computeStress( res, tan, inc );

          fU += B.transpose() * res.stress * qp.J0xW;
          kUU += B.transpose() * tan.dStressddStrain * B * qp.J0xW;

          for ( int n = 0; n < nNonlocalVariables; n++ ) {
            Eigen::Index idx = n * nNonLocalNodes;
            fK.segment( idx, nNonLocalNodes ) += ( N_K.transpose() * K( n ) +
                                                   res.c( n ) * dNdX_K.transpose() * dNdX_K *
                                                     qK.segment( idx, nNonLocalNodes ) -
                                                   N_K.transpose() * res.KLocal( n ) ) *
                                                 qp.J0xW;

            kUK.block( 0, idx, sizeDoFU, nNonLocalNodes ) += B.transpose() * tan.dStressddK.col( n ) * N_K * qp.J0xW;
            kKU.block( idx, 0, nNonLocalNodes, sizeDoFU ) += N_K.transpose() * -tan.dKLocalddStrain.row( n ) * B *
                                                             qp.J0xW;

            kKK.block( idx, idx, nNonLocalNodes, nNonLocalNodes ) += ( N_K.transpose() * N_K +
                                                                       res.c( n ) * dNdX_K.transpose() * dNdX_K +
                                                                       tan.dcddK( n ) * dNdX_K.transpose() * dNdX_K *
                                                                         qK.segment( idx, nNonLocalNodes ) * N_K -
                                                                       N_K.transpose() * tan.dKLocalddK( n, n ) *
                                                                         N_K ) *
                                                                     qp.J0xW;
          }
        }
        else
          throw std::invalid_argument( "Invalid section type for 3D element! Must be Solid" );
      }
      qp.managedStateVars->stress              = res.stress;
      qp.managedStateVars->elasticStrainEnergy = res.elasticEnergyDensity * qp.J0xW;
      qp.managedStateVars->dissipation         = res.dissipation * qp.J0xW;
      qp.managedStateVars->totalStrainEnergy   = ( res.elasticEnergyDensity + res.dissipation ) * qp.J0xW;
      qp.managedStateVars->strain += make3DVoigt< ParentGeometryElement::voigtSize >( dE );
    }
  }

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  void GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVariables, nNonLocalNodes >::
    computeKernelsExplicit( const double* QTotal_, const double* dQ_, double* Pe_, double time, double dT )
  {

    Map< const RhsSized > QTotal( QTotal_ );
    Map< const RhsSized > dQ( dQ_ );
    Map< RhsSized >       Pe( Pe_ );

    // incremental nodal displacements Q and nodal internal parameters qK
    const Ref< const USizedVector > dQU( dQ.head( sizeDoFU ) );
    const Ref< const KSizedVector > dQK( dQ.tail( sizeDoFK ) );
    const Ref< const USizedVector > qU( QTotal.head( sizeDoFU ) );
    const Ref< const KSizedVector > qK( QTotal.tail( sizeDoFK ) );

    // Righthandside
    Ref< USizedVector > fU( Pe.head( sizeDoFU ) );
    Ref< KSizedVector > fK( Pe.tail( sizeDoFK ) );

    using namespace Marmot;
    using namespace ContinuumMechanics::VoigtNotation;

    using response  = typename MarmotMaterialGeneralGradientEnhancedHypoElastic< nNonlocalVariables >::response;
    using increment = typename MarmotMaterialGeneralGradientEnhancedHypoElastic< nNonlocalVariables >::increment;

    for ( size_t i = 0; i < this->qps.size(); i++ ) {

      QuadraturePoint& qp = this->qps[i];

      const BSized&      B      = qp.B;
      const NSizedK&     N_K    = qp.N_K;
      const dNdXiSizedK& dNdX_K = qp.dNdX_K;

      Voigt dE = B * dQU;

      // delta of _K at Gausspoint
      Vector< double, nNonlocalVariables > K;
      Vector< double, nNonlocalVariables > dK;

      for ( size_t n = 0; n < nNonlocalVariables; n++ ) {
        Eigen::Index idx = static_cast< Eigen::Index >( n ) * static_cast< Eigen::Index >( nNonLocalNodes );
        K( n )           = N_K * qK.segment( idx, nNonLocalNodes );
        dK( n )          = N_K * dQK.segment( idx, nNonLocalNodes );
      }

      response  res;
      increment inc;
      if constexpr ( nDim == 2 ) {
        Vector6d dE6             = ContinuumMechanics::VoigtNotation::planeVoigtToVoigt( dE );
        res.stress               = qp.managedStateVars->stress;
        res.elasticEnergyDensity = qp.managedStateVars->elasticStrainEnergy / qp.J0xW;
        res.dissipation          = qp.managedStateVars->dissipation / qp.J0xW;
        res.stateVars            = qp.managedStateVars->materialStateVars.data();
        inc                      = { dE6, K, dK, time, dT };
        Voigt S                  = Voigt::Zero();

        if ( sectionType == SectionType::PlaneStress ) {
          qp.material->computePlaneStressExplicit( res, inc );
          S = ContinuumMechanics::VoigtNotation::voigtToPlaneVoigt( res.stress );
        }
        else if ( sectionType == SectionType::PlaneStrain ) {
          qp.material->computeStressExplicit( res, inc );
          S = ContinuumMechanics::VoigtNotation::voigtToPlaneVoigt( res.stress );
        }
        else {
          throw std::invalid_argument( "Invalid section type for 2D element, expected PlaneStress or PlaneStrain" );
        }

        fU += B.transpose() * S * qp.J0xW;

        for ( int n = 0; n < nNonlocalVariables; n++ ) {
          Eigen::Index idx = n * nNonLocalNodes;
          fK.segment( idx, nNonLocalNodes ) += ( N_K.transpose() * K( n ) +
                                                 res.c( n ) * dNdX_K.transpose() * dNdX_K *
                                                   qK.segment( idx, nNonLocalNodes ) -
                                                 N_K.transpose() * res.KLocal( n ) ) *
                                               qp.J0xW;
        }
      }

      else if ( nDim == 3 ) {
        if ( sectionType == Solid ) {

          res.stress               = qp.managedStateVars->stress;
          res.elasticEnergyDensity = qp.managedStateVars->elasticStrainEnergy / qp.J0xW;
          res.dissipation          = qp.managedStateVars->dissipation / qp.J0xW;
          res.stateVars            = qp.managedStateVars->materialStateVars.data();
          inc                      = { dE, K, dK, time, dT };
          qp.material->computeStressExplicit( res, inc );

          fU += B.transpose() * res.stress * qp.J0xW;

          for ( int n = 0; n < nNonlocalVariables; n++ ) {
            Eigen::Index idx = n * nNonLocalNodes;
            fK.segment( idx, nNonLocalNodes ) += ( N_K.transpose() * K( n ) +
                                                   res.c( n ) * dNdX_K.transpose() * dNdX_K *
                                                     qK.segment( idx, nNonLocalNodes ) -
                                                   N_K.transpose() * res.KLocal( n ) ) *
                                                 qp.J0xW;
          }
        }
        else
          throw std::invalid_argument( "Invalid section type for 3D element! Must be Solid" );
      }
      qp.managedStateVars->stress = res.stress;
      qp.managedStateVars->strain += make3DVoigt< ParentGeometryElement::voigtSize >( dE );
      qp.managedStateVars->elasticStrainEnergy = res.elasticEnergyDensity * qp.J0xW;
      qp.managedStateVars->dissipation         = res.dissipation * qp.J0xW;
      qp.managedStateVars->totalStrainEnergy   = ( res.elasticEnergyDensity + res.dissipation ) * qp.J0xW;
    }
  }

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  void GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVariables, nNonLocalNodes >::
    computeConsistentInertia( double* M )
  {
    Map< KeSizedMatrix > Me( M );
    Me.setZero();

    for ( const auto& qp : qps ) {
      const auto                  N_  = localGeometryElement.NB( localGeometryElement.N( qp.xi ) );
      const NSizedK&              N_K = qp.N_K;
      const double                rho = qp.material->getDensity( qp.managedStateVars->materialStateVars.data() );
      const std::vector< double > eta = qp.material->getNonlocalViscosity(
        qp.managedStateVars->materialStateVars.data() );
      Me.topLeftCorner( sizeDoFU, sizeDoFU ) += N_.transpose() * N_ * qp.J0xW * rho;
      for ( int n = 0; n < nNonlocalVariables; n++ ) {
        Eigen::Index idx = n * nNonLocalNodes;
        Me.bottomRightCorner( sizeDoFK, sizeDoFK )
          .block( idx, idx, nNonLocalNodes, nNonLocalNodes ) += N_K.transpose() * N_K * qp.J0xW * eta[n];
      }
    }
  }

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  void GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVariables, nNonLocalNodes >::
    computeLumpedInertia( double* M )
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
    constexpr int nNodesLinear = ( 1 << nDim );

    // computeLumpedInertia() only knows how to weight the non-local block against either the
    // displacement interpolation itself (equal-order, nNonLocalNodes == nNodes) or the linear
    // corner-node interpolation (reduced-order, nNonLocalNodes == nNodesLinear); assigning N_lin
    // to a differently-sized non-local weight vector below would silently misbehave for any other
    // nNonLocalNodes.
    static_assert( nNodes == nNonLocalNodes || nNonLocalNodes == nNodesLinear,
                   "GeneralGradientEnhancedDisplacementFiniteElement::computeLumpedInertia requires the "
                   "non-local field to use either the displacement interpolation order (nNonLocalNodes == "
                   "nNodes) or linear corner-node interpolation (nNonLocalNodes == 2^nDim)." );

    auto            linGeometryEl    = MarmotGeometryElement< nDim, nNodesLinear >();
    Eigen::VectorXd rowSumsHighOrder = Eigen::VectorXd::Zero( nNodes );
    Eigen::VectorXd rowSumsLinear    = Eigen::VectorXd::Zero( nNodesLinear );
    for ( const auto& qp : qps ) {
      rowSumsHighOrder += Eigen::VectorXd( localGeometryElement.N( qp.xi ) ) * qp.J0xW;
      rowSumsLinear += Eigen::VectorXd( linGeometryEl.N( qp.xi ) ) * qp.J0xW;
    }
    const double weight = FiniteElement::MassLumping::manifoldBlendWeight( rowSumsHighOrder, rowSumsLinear );

    for ( const auto& qp : qps ) {
      const auto N_    = localGeometryElement.N( qp.xi );
      const auto N_lin = linGeometryEl.N( qp.xi );

      VectorXd N_weighted = weight * ( N_ );
      N_weighted.head( nNodesLinear ) += ( 1.0 - weight ) * N_lin;

      // when nNodes == nNonlocalNodes
      VectorXd N_weighted_nonlocal;
      if ( nNodes != nNonLocalNodes ) {
        N_weighted_nonlocal = VectorXd::Zero( nNonLocalNodes );
        N_weighted_nonlocal = N_lin;
      }
      else {
        N_weighted_nonlocal = N_weighted;
      }

      const double                rho = qp.material->getDensity( qp.managedStateVars->materialStateVars.data() );
      const std::vector< double > eta = qp.material->getNonlocalViscosity(
        qp.managedStateVars->materialStateVars.data() );
      VectorXd m_ = N_weighted * qp.J0xW * rho;
      for ( int i = 0; i < nNodes; i++ ) {
        for ( int d = 0; d < nDim; d++ )
          LMM( i * nDim + d ) += m_( i );
      }
      for ( int n = 0; n < nNonlocalVariables; n++ ) {
        Eigen::Index idx = n * nNonLocalNodes;
        VectorXd     mK  = N_weighted_nonlocal * qp.J0xW * eta[n];
        for ( int i = 0; i < nNonLocalNodes; i++ )
          LMM( sizeDoFU + idx + i ) += mK( i );
      }
    }
  }

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  void GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVariables, nNonLocalNodes >::
    computeCriticalTimeStepForExplicitDynamics( double& criticalTimeStep, const double* QTotal )
  {
    using response = typename MarmotMaterialGeneralGradientEnhancedHypoElastic< nNonlocalVariables >::response;

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
     *
     * Note this covers the DISPLACEMENT field only. The nonlocal field's own forward-Euler limit,
     * set by its viscosity against the Helmholtz operator, is still not checked anywhere -- see the
     * TODO below.
     */
    constexpr int   nNodesLinear     = ( 1 << nDim );
    auto            linGeometryEl    = MarmotGeometryElement< nDim, nNodesLinear >();
    Eigen::VectorXd rowSumsHighOrder = Eigen::VectorXd::Zero( nNodes );
    Eigen::VectorXd rowSumsLinear    = Eigen::VectorXd::Zero( nNodesLinear );
    for ( const auto& qp : qps ) {
      rowSumsHighOrder += Eigen::VectorXd( localGeometryElement.N( qp.xi ) ) * qp.J0xW;
      rowSumsLinear += Eigen::VectorXd( linGeometryEl.N( qp.xi ) ) * qp.J0xW;
    }
    const double weight = FiniteElement::MassLumping::manifoldBlendWeight( rowSumsHighOrder, rowSumsLinear );
    const double lumpedMassTimeStepFactor = FiniteElement::MassLumping::timeStepFactorFromMassDistribution(
      FiniteElement::MassLumping::manifoldMassFractions( rowSumsHighOrder, rowSumsLinear, weight ) );

    // TODO: current implementation ignores nonlocal variables
    criticalTimeStep = std::numeric_limits< double >::max();
    for ( const auto& qp : qps ) {
      /* The characteristic length has to be the element's SMALLEST physical extent, not a
       * volume-averaged one. cbrt( 8 * detJ ) and its lower-dimensional analogues are volume
       * based: for a sliver -- thin in one direction but not the others -- the volume stays
       * moderate while the thin dimension collapses, so they OVERESTIMATE the length and hence
       * the stable time step. Refining a distorted parent element is precisely how slivers are
       * produced, so an h-adaptive explicit run integrates its most distorted elements above
       * their stability limit and diverges thousands of increments later, with nothing in the
       * log connecting the divergence to the refinement that caused it.
       *
       * The Jacobian maps the natural cube [-1,1]^nDim onto the element, so twice its smallest
       * singular value IS that smallest physical extent. For a well-shaped element this
       * reproduces the previous expressions exactly -- a cube of side h gives h either way --
       * so the estimate is tightened only where it was previously wrong.
       */
      const JacobianSized J_ = localGeometryElement.Jacobian( localGeometryElement.dNdXi( qp.xi ) );
      const double        characteristicElementLength = 2.0 *
                                                 Eigen::JacobiSVD< JacobianSized >( J_ ).singularValues().minCoeff();
      response waveSpeedResponse;
      waveSpeedResponse.stress = qp.managedStateVars->stress;
      waveSpeedResponse.KLocal.setZero();
      waveSpeedResponse.c.setZero();
      waveSpeedResponse.stateVars            = qp.managedStateVars->materialStateVars.data();
      waveSpeedResponse.elasticEnergyDensity = qp.managedStateVars->elasticStrainEnergy / qp.J0xW;
      waveSpeedResponse.dissipation          = qp.managedStateVars->dissipation / qp.J0xW;

      const double c = qp.material->getMaximumWaveSpeed( waveSpeedResponse );
      if ( c <= 0.0 )
        throw std::runtime_error( "Material returned non-positive wave speed, cannot compute critical time step" );
      const double& l  = characteristicElementLength;
      double        dt = lumpedMassTimeStepFactor * l / c;
      if ( dt < criticalTimeStep )
        criticalTimeStep = dt;
    }
  }

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  void GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVariables, nNonLocalNodes >::
    computeInternalEnergy( double& internalEnergy )
  {
    internalEnergy = 0.0;
    for ( const auto& qp : qps ) {
      internalEnergy += qp.managedStateVars->totalStrainEnergy;
    }
  }

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  void GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVariables, nNonLocalNodes >::
    computeDistributedLoad( MarmotElement::DistributedLoadTypes loadType,
                            double*                             P,
                            double*                             K,
                            const int                           elementFace,
                            const double*                       load,
                            const double*                       QTotal,
                            double                              time,
                            double                              dT )
  {

    Map< RhsSized >     P_( P );
    Ref< USizedVector > fU( P_.head( sizeDoFU ) );

    switch ( loadType ) {

    case MarmotElement::Pressure: {
      const double p = load[0];

      FiniteElement::BoundaryElement boundaryEl( localGeometryElement.shape,
                                                 elementFace,
                                                 nDim,
                                                 localGeometryElement.coordinates );

      VectorXd Pk = -p * boundaryEl.computeSurfaceNormalVectorialLoadVector();

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

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  void GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVariables, nNonLocalNodes >::
    setInitialConditions( StateTypes state, const double* values )
  {
    if constexpr ( nDim > 1 ) {
      switch ( state ) {
      case MarmotElement::GeostaticStress: {
        for ( QuadraturePoint& qp : qps ) {

          XiSized coordAtGauss = localGeometryElement.NB( localGeometryElement.N( qp.xi ) ) *
                                 localGeometryElement.coordinates;

          const double sigY1 = values[0];
          const double sigY2 = values[2];
          const double y1    = values[1];
          const double y2    = values[3];

          using namespace Math;
          qp.managedStateVars->stress( 1 ) = linearInterpolation( coordAtGauss[1], y1, y2, sigY1, sigY2 ); // sigma_y
          qp.managedStateVars->stress( 0 ) = values[4] * qp.managedStateVars->stress( 1 );                 // sigma_x
          qp.managedStateVars->stress( 2 ) = values[5] * qp.managedStateVars->stress( 1 );                 // sigma_z
        }
        break;
      }
      case MarmotElement::MarmotMaterialStateVars: {
        throw std::invalid_argument( "Please use initializeStateVars directly on material" );
      }
      case MarmotElement::MarmotMaterialInitialization: {
        for ( QuadraturePoint& qp : qps ) {
          qp.material->initializeYourself( qp.managedStateVars->materialStateVars.data(),
                                           qp.managedStateVars->materialStateVars.size() );
        }
        break;
      }
      default: throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid initial condition" );
      }
    }
  }

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  void GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVariables, nNonLocalNodes >::
    computeBodyForce( double*       P_,
                      double*       K,
                      const double* load,

                      const double* QTotal,
                      double        time,
                      double        dT )
  {
    Map< RhsSized >                              P( P_ );
    Ref< USizedVector >                          Pe( P.head( sizeDoFU ) );
    const Map< const Matrix< double, nDim, 1 > > f( load );

    for ( const auto& qp : qps )
      Pe += localGeometryElement.NB( localGeometryElement.N( qp.xi ) ).transpose() * f * qp.J0xW;
  }

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  std::vector< double > GeneralGradientEnhancedDisplacementFiniteElement< nDim,
                                                                          nNodes,
                                                                          nNonlocalVariables,
                                                                          nNonLocalNodes >::getCoordinatesAtCenter()
  {
    std::vector< double > coords( nDim );

    Eigen::Map< XiSized > coordsMap( &coords[0] );
    const auto            centerXi = XiSized::Zero();
    coordsMap = localGeometryElement.NB( localGeometryElement.N( centerXi ) ) * localGeometryElement.coordinates;
    return coords;
  }

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  std::vector< std::vector< double > > GeneralGradientEnhancedDisplacementFiniteElement<
    nDim,
    nNodes,
    nNonlocalVariables,
    nNonLocalNodes >::getCoordinatesAtQuadraturePoints()
  {
    std::vector< std::vector< double > > listedCoords;

    std::vector< double > coords( nDim );
    Eigen::Map< XiSized > coordsMap( &coords[0] );

    for ( const auto& qp : qps ) {
      coordsMap = localGeometryElement.NB( localGeometryElement.N( qp.xi ) ) * localGeometryElement.coordinates;
      listedCoords.push_back( coords );
    }

    return listedCoords;
  }

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  int GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVariables, nNonLocalNodes >::
    getNumberOfQuadraturePoints()
  {
    return qps.size();
  }
} // namespace Marmot::Elements
