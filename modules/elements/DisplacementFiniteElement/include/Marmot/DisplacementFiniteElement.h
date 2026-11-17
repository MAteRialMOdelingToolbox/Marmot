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
#include "Marmot/MarmotBulkViscosity.h"
#include "Marmot/MarmotConsistentMass.h"
#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotEnergyGradientProvider.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryElement.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotLowerDimensionalStress.h"
#include "Marmot/MarmotMassLumping.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotMaterialHypoElasticEnergyTangents.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <vector>

namespace Marmot::Discretizations::Elements {

  /**
   * @class Marmot::Discretizations::Elements::DisplacementFiniteElement
   * @tparam nDim Number of spatial dimensions (1, 2, or 3).
   * @tparam nNodes Number of element nodes.
   * @brief Displacement-based finite element template.
   * @details Uses linearized kinematics (small strains) and supports 1D, 2D (plane stress, plane strain and
   * axisymmetric), and 3D formulations via a section assumption. Holds quadrature-point state and
   * delegates constitutive updates to Marmot materials while assembling element
   * residuals, tangents, and mass matrices.
   */
  template < int nDim, int nNodes >
  class DisplacementFiniteElement : public Marmot::Discretizations::Elements::Element,
                                    public Marmot::Discretizations::Elements::GeometryElement< nDim, nNodes >,
                                    public Marmot::Discretizations::Elements::EnergyGradientProvider {

  public:
    /**
     * @brief Kinematic section assumption used by the element.
     * @details Controls the kinematics and which constitutive call is performed at quadrature points.
     * - UniaxialStress
     * - PlaneStress
     * - PlaneStrain
     * - Solid
     * - Axisymmetric: 2D only, \f$r = x_1 \geq 0\f$ and \f$z = x_2\f$, with the hoop strain
     *   \f$\varepsilon_{\theta\theta} = u_r / r\f$; integrated over the full ring, \f$\mathrm{d}V = 2\pi r\,
     *   \mathrm{d}A\f$.
     */
    enum SectionType {
      UniaxialStress,
      PlaneStress,
      PlaneStrain,
      Solid,
      Axisymmetric,
    };

    static constexpr int sizeLoadVector = nNodes * nDim;
    static constexpr int nCoordinates   = nNodes * nDim;

    using ParentGeometryElement = Marmot::Discretizations::Elements::GeometryElement< nDim, nNodes >;
    using JacobianSized         = typename ParentGeometryElement::JacobianSized;
    using dNdXiSized            = typename ParentGeometryElement::dNdXiSized;
    using BSized                = typename ParentGeometryElement::BSized;
    using BSizedAxisymmetric    = typename ParentGeometryElement::BSizedAxisymmetric;
    using XiSized               = typename ParentGeometryElement::XiSized;
    using NSized                = typename ParentGeometryElement::NSized;
    using RhsSized              = Eigen::Matrix< double, sizeLoadVector, 1 >;
    using KeSizedMatrix         = Eigen::Matrix< double, sizeLoadVector, sizeLoadVector >;
    using CSized = Eigen::Matrix< double, ParentGeometryElement::voigtSize, ParentGeometryElement::voigtSize >;
    using Voigt  = Eigen::Matrix< double, ParentGeometryElement::voigtSize, 1 >;

    /** Element-level properties (e.g., thickness for plane stress and plane strain, area for 1D). */
    Eigen::Map< const Eigen::VectorXd > elementProperties;
    /** Element label (ID) used for logging and material creation. */
    const int elLabel;
    /** Section assumption applied by this element instance. */
    const SectionType sectionType;
    /**
     * @brief Coefficients of the artificial bulk viscosity, assigned via the named property
     * "bulk viscosity". They default to zero, which is inactive: unless that property is
     * assigned, no viscous stress is formed and the element integrates exactly what it
     * integrated before the device existed.
     */
    Discretizations::Elements::BulkViscosity::Coefficients bulkViscosityCoefficients;

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
       * @brief The characteristic length for the critical time step at this quadrature point,
       * GeometryElement::characteristicLengthForCriticalTimeStepAt: the smallest physical extent of a quadrilateral or
       * hexahedron, 1 / |grad N|_F of a linear simplex.
       * @details Cached because the artificial bulk viscosity needs it on every explicit increment,
       * whereas the stable time increment that shares the definition is asked for rarely. The
       * element's nodal coordinates do not change, so the cached value cannot go stale.
       */
      double characteristicLengthForCriticalTimeStep = 0.0;

      /**
       * @brief Wave speed the artificial bulk viscosity is scaled with, cached on first use.
       * @details A current wave speed costs a full constitutive evaluation, affordable for the
       * stable increment (asked for rarely) and not per quadrature point per explicit increment.
       * What is cached is the UNDAMAGED speed: the term damps the highest frequency the MESH can
       * carry, which the undamaged material sets, and holding it fixed as the material softens
       * damps slightly harder -- the safe direction here. Empty means "not yet computed"; a
       * reference of zero is cached like any other value. It is
       * also the reference the optional degradation is measured against; see
       * Marmot::Discretizations::Elements::BulkViscosity::degradationFactor.
       *
       * @warning Captured on the first explicit increment, from the state the element has THEN. A
       * hypoelastic material carries its damage in its state variables, so a run beginning from an
       * already-damaged state -- a restart, or an explicit step after an implicit one -- takes that
       * degraded speed as its reference and measures no degradation afterwards. The viscosity is
       * then simply not degraded, which removes more energy rather than less, but it is not what
       * was asked for. The gradient-enhanced element has no such ambiguity: its reference is the
       * speed at a ZERO non-local field.
       */
      std::optional< double > referenceWaveSpeed;

      /**
       * @brief Manager for per-quadrature-point state variables.
       * @details Provides named accessors to stress \f$\sig\f$, strain \f$\eps\f$
       * and the material state vector. The layout is [stress(6), strain(6), begin of material state(...)]
       * in 3D Voigt notation.
       */
      class QPStateVarManager : public Marmot::Utilities::StateVarVectorManager {

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
        Marmot::Math::mVector6d       stress;
        Marmot::Math::mVector6d       strain;
        double&                       totalStrainEnergy;
        double&                       elasticStrainEnergy;
        double&                       dissipation;
        Eigen::Map< Eigen::VectorXd > materialStateVars;

        static int getNumberOfRequiredStateVarsQuadraturePointOnly() { return layout.nRequiredStateVars; };

        QPStateVarManager( double* theStateVarVector, int nStateVars )
          : Marmot::Utilities::StateVarVectorManager( theStateVarVector, layout ),
            stress( &find( "stress" ) ),
            strain( &find( "strain" ) ),
            totalStrainEnergy( find( "total strain energy" ) ),
            elasticStrainEnergy( find( "elastic strain energy" ) ),
            dissipation( find( "dissipation" ) ),
            materialStateVars( &find( "begin of material state" ),
                               nStateVars - getNumberOfRequiredStateVarsQuadraturePointOnly() ){};
      };

      std::unique_ptr< QPStateVarManager > managedStateVars;

      std::unique_ptr< Marmot::Materials::MaterialHypoElastic > material;

      /// The material's energy tangent capability, nullptr if it has none
      const Marmot::Materials::HypoElasticEnergyTangentProvider* energyTangentProvider = nullptr;

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
     * @param quadratureRule Quadrature rule: family and degree of exactness.
     * @param sectionType Section assumption (1D/2D/3D).
     * @throws std::invalid_argument for an axisymmetric section with nDim other than 2.
     */
    DisplacementFiniteElement( int                                                   elementID,
                               Discretizations::Elements::Quadrature::QuadratureRule quadratureRule,
                               SectionType                                           sectionType );

    /** @brief Total number of required state variables for this element (sum over all quadrature points). */
    int getNumberOfRequiredStateVars() override;

    /** @brief Node-level fields exposed by the element. Returns ["displacement"] for each node. */
    std::vector< std::vector< std::string > > getNodeFields() override;

    /** @brief Permutation pattern from local DOF ordering to solver ordering (identity by default). */
    std::vector< int > getDofIndicesPermutationPattern() override;

    /** @brief Number of nodes of this element type. */
    int getNNodes() override { return nNodes; }

    /** @brief Number of spatial dimensions. */
    int getNSpatialDimensions() override { return nDim; }

    /** @brief Number of degrees of freedom per element (nNodes * nDim). */
    int getNDofPerElement() override { return sizeLoadVector; }

    /** @brief Geometric shape of the element (as reported by the parent geometry element). */
    std::string getElementShape() override { return ParentGeometryElement::getElementShape(); }

    /** @brief Map the provided element state vector to all quadrature points. */
    void assignStateVars( double* stateVars, int nStateVars ) override;

    /** @brief Assign element properties (e.g., thickness in 2D, area in 1D). */
    void assignProperty( const Marmot::Discretizations::ElementProperties& marmotElementProperty ) override;

    /** @brief Assign material section and instantiate per-quadrature-point materials. */
    void assignProperty( const Marmot::Discretizations::MaterialSection& marmotElementProperty ) override;

    /**
     * @brief Assign a named element property.
     * @param propertyName One of "bulk viscosity" or "bulk viscosity damage degradation".
     * @param properties For "bulk viscosity": the two dimensionless coefficients \f$b_1\f$
     *        (linear) and \f$b_2\f$ (quadratic), in that order. For "bulk viscosity damage
     *        degradation": the single exponent \f$n\f$ of the optional degradation with the
     *        material's loss of stiffness.
     * @param nProperties Number of values behind that pointer: 2 for "bulk viscosity", 1 for
     *        "bulk viscosity damage degradation".
     * @throws std::invalid_argument if the name is not understood, if the count does not match
     *         what the property expects, or if a coefficient is negative.
     */
    void assignProperty( const std::string& propertyName, const double* properties, int nProperties ) override;

    /** @brief The named properties this element understands. */
    std::vector< std::string > getPropertyNames() const override;

    /** @brief Provide nodal coordinates to the parent geometry element. */
    void assignNodeCoordinates( const double* coordinates ) override;

    /**
     * @brief Precompute geometry-related quantities at quadrature points (B, detJ, J0xW).
     * @details \f$J_0 w = w \det\mathbf{J}\, s(\boldsymbol{\xi})\f$ with the section factor of sectionFactorAt().
     * @throws std::invalid_argument for an axisymmetric section with a node at \f$r < -10^{-12} h\f$, \f$h\f$ the
     * diagonal of the element's bounding box, or a quadrature point at \f$r \leq 0\f$.
     */
    void initializeYourself() override;

    /**
     * @brief Factor between the volume and the reference area or length of the element at @p xi.
     * @param xi Parametric coordinates.
     * @return The cross section (1D), the thickness (plane stress and plane strain), \f$2\pi r(\boldsymbol{\xi})\f$
     * (axisymmetric), or 1 (solid).
     */
    double sectionFactorAt( const XiSized& xi ) const;

    /**
     * @brief Initialize state or materials.
     * @param state MarmotMaterialInitialization, GeostaticStress or MarmotMaterialStateVars.
     * @param values For GeostaticStress: [sigmaY(z1), y1, sigmaY(z2), y2, kx, kz].
     */
    void setInitialConditions( StateTypes state, const double* values ) override;

    /**
     * @brief Assemble distributed surface loads on a boundary face.
     * @details Pressure and traction contributions are integrated on the boundary \f$\Gamma_e\f$:
     * \f[
     * \mathbf{P}_e^{(p)} = - \int_{\Gamma_e} p\, \mathbf{N}^\mathsf{T} \mathbf{n}\, \mathrm{d}\Gamma,\qquad
     * \mathbf{P}_e^{(t)} = \int_{\Gamma_e} \mathbf{N}^\mathsf{T} \mathbf{t}\, \mathrm{d}\Gamma.
     * \f]
     * Axisymmetric elements integrate over the ring, \f$\mathrm{d}\Gamma = 2\pi r\, \mathrm{d}s\f$.
     * @param loadType Pressure or SurfaceTraction.
     * @param P Element RHS contribution (accumulated).
     * @param K Optional stiffness contribution (unused).
     * @param elementFace Boundary face index.
     * @param load Pressure magnitude or traction vector (size nDim).
     * @param QTotal Total DOF vector (unused).
     * @param time Current time data forwarded to materials.
     * @param dT Time increment.
     */
    void computeDistributedLoad( Marmot::Discretizations::Elements::Element::DistributedLoadTypes loadType,
                                 double*                                                          P,
                                 double*                                                          K,
                                 int                                                              elementFace,
                                 const double*                                                    load,
                                 const double*                                                    QTotal,
                                 double                                                           time,
                                 double                                                           dT ) override;

    /**
     * @brief Assemble body force contribution.
     * @details Integrates \f$\mathbf{P}_e^{(b)} = \int_{\Omega_e} \mathbf{N}^\mathsf{T} \mathbf{f}\,
     * \mathrm{d}\Omega\f$.
     */
    void computeBodyForce( double* P, double* K, const double* load, const double* QTotal, double time, double dT )
      override;

    /**
     * @brief Compute internal force and consistent tangent stiffness.
     * @details Uses the small-strain relation \f$\Delta\boldsymbol{\varepsilon}=\mathbf{B}\,\Delta\mathbf{u}\f$ and
     * integrates
     * \f[
     * \mathbf{K}_e = \sum_{qp} \mathbf{B}^\mathsf{T} \mathbf{C} \mathbf{B}\, J_0 w,\qquad
     * \mathbf{P}_e = \sum_{qp} \mathbf{B}^\mathsf{T} \boldsymbol{\sigma}\, J_0 w.
     * \f]
     * Axisymmetric elements use the 4-row operator (rr, zz, \f$\theta\theta\f$, rz) with the hoop row
     * \f$N_i / r\f$.
     * @param QTotal Total displacement vector.
     * @param dQ Incremental displacement.
     * @param Pe Internal force vector (accumulated).
     * @param Ke Tangent stiffness matrix (accumulated).
     * @param time Time data forwarded to materials.
     * @param dT Time increment.
     */
    void computeKernels( const double* QTotal, const double* dQ, double* Pe, double* Ke, double time, double dT )
      override;

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
    void computeKernelsExplicit( const double* QTotal, const double* dQ, double* Pe, double time, double dT ) override;
    /**
     * @brief Compute consistent mass matrix using material density.
     * @details \f$\mathbf{M}_e = \sum_{p} \rho\, \mathbf{N}^\mathsf{T}\mathbf{N}\, J_0 w\f$ over the
     * points \f$p\f$ of the full Gauss rule of the element's shape, also for a reduced-integration
     * element, whose own rule would leave the mass rank-deficient; the section factor (sectionFactorAt()) is
     * evaluated at \f$p\f$, the density is taken from the nearest quadrature point of the element. See
     * Marmot::Discretizations::Elements::ConsistentMass.
     */
    void computeConsistentInertia( double* M ) override;

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
     * Marmot::Discretizations::Elements::MassLumping::manifoldBlendWeight(), which returns
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
    void computeLumpedInertia( double* M ) override;

    /**
     * @brief Compute the critical time step for explicit dynamics.
     * @param[out] criticalTimeStep Computed critical time step.
     * @param[in] QTotal Total (current) degrees of freedom of the element.
     * @details The estimate is \f$l / c\f$, scaled by the factor
     * Marmot::Discretizations::Elements::MassLumping::timeStepFactorFromMassDistribution() derives from the
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
    void computeCriticalTimeStepForExplicitDynamics( double& criticalTimeStep, const double* QTotal ) override;

    /**
     * @brief Compute the internal energy of the element.
     * @param internalEnergy Output parameter for the computed internal energy.
     */
    void computeInternalEnergy( double& internalEnergy ) override;

    /**
     * @brief Assign the arrays computeKernels() writes the energy gradients to, see EnergyGradientProvider.
     * @param dElasticEnergy_dQ Gradient of the elastic energy; nullptr, together with @p dDissipation_dQ, stops it.
     * @param dDissipation_dQ Gradient of the dissipation.
     * @details Supported for the solid, plane strain and axisymmetric sections with a material that is a
     * Marmot::Materials::HypoElasticEnergyTangentProvider.
     */
    void assignEnergyGradientOutput( double* dElasticEnergy_dQ, double* dDissipation_dQ ) override;

    /**
     * @brief Access a named state view at a quadrature point.
     * @note Using "sdv" returns the raw material state vector and is deprecated.
     */
    Marmot::Utilities::StateView getStateView( const std::string& stateName, int qpNumber ) override
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

    /** @brief Get physical coordinates at the element center (the centroid of the reference domain). */
    std::vector< double > getCoordinatesAtCenter() override;

    /** @brief Get physical coordinates at each quadrature point. */
    std::vector< std::vector< double > > getCoordinatesAtQuadraturePoints() override;

    /** @brief Number of quadrature points of this element. */
    int getNumberOfQuadraturePoints() override;

  private:
    /// Outputs of the energy gradients, see assignEnergyGradientOutput()
    double* dElasticEnergy_dQOutput = nullptr;
    double* dDissipation_dQOutput   = nullptr;

    /**
     * @brief Throws std::invalid_argument unless the section and the materials of all quadrature points provide the
     * energy gradients: solid, plane strain and axisymmetric sections, materials that are a
     * Marmot::Materials::HypoElasticEnergyTangentProvider.
     */
    void checkEnergyGradientSupport() const;

    /**
     * @brief Internal force and, unless @p Ke is nullptr (explicit), the tangent; shared by computeKernels() and
     * computeKernelsExplicit().
     * @param dQ Incremental displacement.
     * @param Pe Internal force vector (accumulated).
     * @param Ke Tangent stiffness matrix (accumulated), nullptr for the explicit update.
     * @param time Time data forwarded to materials.
     * @param dT Time increment.
     */
    void integrateInternalForce( const double* dQ, double* Pe, double* Ke, double time, double dT );

    /**
     * @brief Material update and integration at one quadrature point.
     * @tparam nComponents Strain components of the section: 1 (uniaxial), 3 (plane), 4 (axisymmetric), 6 (solid).
     * @param qp The quadrature point.
     * @param B Strain operator of the section.
     * @param dQ Incremental displacement.
     * @param Pe Internal force vector (accumulated).
     * @param Ke Tangent stiffness matrix (accumulated), nullptr for the explicit update.
     * @param energyGradients Whether the energy gradients are accumulated into the assigned outputs.
     * @param time Time data forwarded to materials.
     * @param dT Time increment.
     */
    template < ContinuumMechanics::Voigt::VoigtSize nComponents >
    void integrateQuadraturePoint( QuadraturePoint&                                            qp,
                                   const Eigen::Matrix< double, nComponents, sizeLoadVector >& B,
                                   const double*                                               dQ,
                                   double*                                                     Pe,
                                   double*                                                     Ke,
                                   bool                                                        energyGradients,
                                   double                                                      time,
                                   double                                                      dT );

    /**
     * @brief Material update through the 3D interface, shared by the solid, plane strain and axisymmetric sections,
     * implicit and explicit; stores the stress.
     * @param qp The quadrature point.
     * @param dE Strain increment in 3D Voigt notation.
     * @param C Algorithmic tangent in 3D Voigt notation; nullptr for the explicit update.
     * @param dElasticEnergyDensity_dStrain Derivative of the elastic energy density, nullptr unless requested
     * (implicit only).
     * @param dDissipation_dStrain Derivative of the dissipation, given together with
     * @p dElasticEnergyDensity_dStrain.
     * @param time Time data forwarded to materials.
     * @param dT Time increment.
     * @return The updated stress and energy densities.
     */
    Marmot::Materials::MaterialHypoElastic::state3D updateMaterial3D(
      QuadraturePoint&              qp,
      const Marmot::Math::Vector6d& dE,
      Marmot::Math::Matrix6d*       C,
      Marmot::Math::Vector6d*       dElasticEnergyDensity_dStrain,
      Marmot::Math::Vector6d*       dDissipation_dStrain,
      double                        time,
      double                        dT );

    /**
     * @brief Artificial bulk viscous stress at a quadrature point, see
     * Marmot::Discretizations::Elements::BulkViscosity.
     * @param qp The quadrature point, after its material update.
     * @param volumetricStrainIncrement Trace of the strain increment.
     * @param elasticEnergyDensity Updated elastic energy density.
     * @param dissipation Updated dissipation density.
     * @param dT Time increment.
     * @return The viscous stress added to the normal components.
     */
    double bulkViscousStress( QuadraturePoint& qp,
                              double           volumetricStrainIncrement,
                              double           elasticEnergyDensity,
                              double           dissipation,
                              double           dT );

    /**
     * @brief Components of a 3D Voigt vector carried by the section: (11, 22, 12) plane, (rr, zz, \f$\theta\theta\f$,
     * rz) axisymmetric.
     * @tparam nComponents Strain components of the section.
     * @param voigt Vector in 3D Voigt notation.
     * @return The section components.
     */
    template < ContinuumMechanics::Voigt::VoigtSize nComponents >
    static Eigen::Matrix< double, nComponents, 1 > sectionComponents( const Marmot::Math::Vector6d& voigt );

    /**
     * @brief Block of a 3D tangent carried by the section, see sectionComponents().
     * @tparam nComponents Strain components of the section.
     * @param C Tangent in 3D Voigt notation.
     * @return The section block.
     */
    template < ContinuumMechanics::Voigt::VoigtSize nComponents >
    static Eigen::Matrix< double, nComponents, nComponents > sectionTangent( const Marmot::Math::Matrix6d& C );

    /**
     * @brief Axisymmetric strain operator, rows (rr, zz, \f$\theta\theta\f$, rz): the planar operator of @p qp with
     * the hoop row \f$N_i / r\f$, see Spatial2D::Axisymmetric::BFromPlanarB().
     * @param qp The quadrature point.
     * @return The operator.
     */
    BSizedAxisymmetric axisymmetricB( const QuadraturePoint& qp ) const;
  };

  template < int nDim, int nNodes >
  DisplacementFiniteElement< nDim, nNodes >::DisplacementFiniteElement(
    int                                                   elementID,
    Discretizations::Elements::Quadrature::QuadratureRule quadratureRule,
    SectionType                                           sectionType )
    : ParentGeometryElement(),
      elementProperties( Eigen::Map< const Eigen::VectorXd >( nullptr, 0 ) ),
      elLabel( elementID ),
      sectionType( sectionType )
  {
    if constexpr ( nDim != 2 )
      if ( sectionType == SectionType::Axisymmetric )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__ << ": the axisymmetric section is two-dimensional." );

    for ( const auto& qpInfo :
          Discretizations::Elements::Quadrature::getQuadraturePoints( this->referenceDomain, quadratureRule ) ) {
      QuadraturePoint qp( qpInfo.xi, qpInfo.weight );
      qps.push_back( std::move( qp ) );
    }
  }

  template < int nDim, int nNodes >
  int DisplacementFiniteElement< nDim, nNodes >::getNumberOfRequiredStateVars()
  {
    return qps[0].getNumberOfRequiredStateVars() * static_cast< int >( qps.size() );
  }

  template < int nDim, int nNodes >
  std::vector< std::vector< std::string > > DisplacementFiniteElement< nDim, nNodes >::getNodeFields()
  {
    using namespace std;

    static const vector< vector< string > > nodeFields = [] {
      vector< vector< string > > fields;
      for ( int i = 0; i < nNodes; i++ ) {
        fields.push_back( vector< string >() );
        fields[i].push_back( "displacement" );
      }
      return fields;
    }();

    return nodeFields;
  }

  template < int nDim, int nNodes >
  std::vector< int > DisplacementFiniteElement< nDim, nNodes >::getDofIndicesPermutationPattern()
  {
    static const std::vector< int > permutationPattern = [] {
      std::vector< int > pattern;
      for ( int i = 0; i < nNodes * nDim; i++ )
        pattern.push_back( i );
      return pattern;
    }();

    return permutationPattern;
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::assignStateVars( double* stateVars, int nStateVars )
  {
    const int nQpStateVars = nStateVars / static_cast< int >( qps.size() );

    for ( size_t i = 0; i < qps.size(); i++ ) {
      auto&   qp          = qps[i];
      double* qpStateVars = stateVars + ( i * nQpStateVars );
      qp.assignStateVars( qpStateVars, nQpStateVars );
    }
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::assignProperty(
    const Marmot::Discretizations::ElementProperties& elementPropertiesInfo )
  {
    new ( &elementProperties ) Eigen::Map< const Eigen::VectorXd >( elementPropertiesInfo.elementProperties,
                                                                    elementPropertiesInfo.nElementProperties );
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::assignProperty(
    const Marmot::Discretizations::MaterialSection& section )
  {
    // the new materials replace the current ones only once all of them are accepted
    std::vector< std::unique_ptr< Marmot::Materials::MaterialHypoElastic > > materials;
    materials.reserve( qps.size() );
    for ( const auto& qp : qps ) {
      auto material = std::unique_ptr< Marmot::Materials::MaterialHypoElastic >(
        Marmot::Registry::MaterialHypoElasticFactory::createMaterial( section.materialName,
                                                                      section.materialProperties,
                                                                      section.nMaterialProperties,
                                                                      elLabel ) );

      if ( !material )
        throw std::invalid_argument(
          Marmot::Utilities::MakeString()
          << __PRETTY_FUNCTION__
          << ": invalid material assigned; cannot cast to Marmot::Materials::MaterialHypoElastic!" );

      /* Deliberately VOLUME based, and deliberately not the length
       * computeCriticalTimeStepForExplicitDynamics uses. The material's characteristic length is a
       * regularisation length -- it sets the width over which a softening law dissipates its
       * fracture energy -- and the volume-equivalent length is the right measure for that. The
       * stability estimate needs the opposite: the element's SMALLEST extent, because that is what
       * bounds the highest frequency. The two definitions are not interchangeable; do not unify
       * them.
       */
      material->setCharacteristicElementLength( ParentGeometryElement::volumeEquivalentLength( qp.detJ ) );

      if ( dElasticEnergy_dQOutput &&
           !dynamic_cast< const Marmot::Materials::HypoElasticEnergyTangentProvider* >( material.get() ) )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__ << ": element " << elLabel
                                     << ": the material provides no energy tangents "
                                        "(Marmot::Materials::HypoElasticEnergyTangentProvider), but energy gradient "
                                        "outputs are assigned; the element keeps its material." );

      materials.push_back( std::move( material ) );
    }

    for ( size_t i = 0; i < qps.size(); i++ ) {
      qps[i].material              = std::move( materials[i] );
      qps[i].energyTangentProvider = dynamic_cast< const Marmot::Materials::HypoElasticEnergyTangentProvider* >(
        qps[i].material.get() );
    }
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::assignProperty( const std::string& propertyName,
                                                                  const double*      properties,
                                                                  int                nProperties )
  {
    if ( propertyName == "bulk viscosity" ) {
      if ( nProperties != 2 )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__
                                     << ": the named property 'bulk viscosity' takes exactly "
                                        "2 values, the linear coefficient b1 and the quadratic coefficient b2, but "
                                     << nProperties << " were given." );

      /* Validated BEFORE anything is committed, so a rejected assignment leaves the element as it
       * was. Non-finite values are rejected alongside negative ones: every comparison against a
       * NaN is false, so a NaN would pass a `< 0.0` test and propagate into the viscous stress.
       */
      if ( !std::isfinite( properties[0] ) || !std::isfinite( properties[1] ) )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__
                                     << ": both bulk viscosity coefficients must be finite numbers." );

      if ( properties[0] < 0.0 || properties[1] < 0.0 )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__
                                     << ": both bulk viscosity coefficients must be non-negative, a negative one "
                                        "would feed energy into the solution rather than remove it." );

      /* Artificial bulk viscosity is a NUMERICAL device, so it is an element property and not a
       * material one: the same concrete integrated implicitly needs none of it, and two meshes of
       * the same material may want different amounts. Both coefficients are dimensionless; see
       * Marmot::Discretizations::Elements::BulkViscosity for what they multiply.
       */
      bulkViscosityCoefficients.linear    = properties[0];
      bulkViscosityCoefficients.quadratic = properties[1];
    }
    else if ( propertyName == "bulk viscosity damage degradation" ) {
      if ( nProperties != 1 )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__
                                     << ": the named property 'bulk viscosity damage degradation' takes exactly "
                                        "1 value, the exponent n of (c/c_0)^n, but "
                                     << nProperties << " were given." );

      /* Opt-in and separate from 'bulk viscosity' itself, so switching it on does not disturb the
       * coefficients: it costs a constitutive evaluation per quadrature point per increment. See
       * BulkViscosity::degradationFactor for the exponent, and for what it degrades with -- the
       * current tangent, not a damage variable.
       */
      if ( !std::isfinite( properties[0] ) )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__
                                     << ": the bulk viscosity damage degradation exponent must be a finite number; a "
                                        "NaN passes every ordering test and would silently disable the option." );

      if ( properties[0] < 0.0 )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__
                                     << ": the bulk viscosity damage degradation exponent must be non-negative, a "
                                        "negative one would AMPLIFY the viscous stress as the material fails." );

      bulkViscosityCoefficients.degradation = properties[0];
    }
    else {
      Marmot::Discretizations::Elements::Element::assignProperty( propertyName, properties, nProperties );
    }
  }

  template < int nDim, int nNodes >
  std::vector< std::string > DisplacementFiniteElement< nDim, nNodes >::getPropertyNames() const
  {
    return { "bulk viscosity", "bulk viscosity damage degradation" };
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::assignNodeCoordinates( const double* coordinates )
  {
    ParentGeometryElement::assignNodeCoordinates( coordinates );
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::initializeYourself()
  {
    const bool axisymmetric = sectionType == SectionType::Axisymmetric;

    // a mesh drawn in r < 0, or straddling the axis, would give J0xW <= 0 and an infinite hoop strain N / r; the
    // tolerance relative to the element size h (diagonal of its bounding box) accepts round-off of nodes on the axis
    if ( axisymmetric ) {
      const Eigen::Map< const Eigen::Matrix< double, nDim, nNodes > > nodes( this->coordinates.data() );
      const double h = ( nodes.rowwise().maxCoeff() - nodes.rowwise().minCoeff() ).norm();
      int          i = 0;
      const double r = nodes.row( 0 ).minCoeff( &i );
      if ( r < -1e-12 * h )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__ << ": axisymmetric element " << elLabel << ": node "
                                     << i + 1 << " (local) lies at r = x_1 = " << r
                                     << " < 0; axisymmetric elements are meshed in r >= 0." );
    }

    for ( QuadraturePoint& qp : qps ) {
      const dNdXiSized    dNdXi = this->dNdXi( qp.xi );
      const JacobianSized J     = this->Jacobian( dNdXi );
      const JacobianSized JInv  = J.inverse();
      const dNdXiSized    dNdX  = this->dNdX( dNdXi, JInv );
      qp.detJ                   = J.determinant();
      qp.B                      = this->B( dNdX );

      qp.characteristicLengthForCriticalTimeStep = this->characteristicLengthForCriticalTimeStepAt( qp.xi );
      qp.referenceWaveSpeed.reset();

      const double sectionFactor = sectionFactorAt( qp.xi );
      if ( axisymmetric && sectionFactor <= 0.0 )
        throw std::invalid_argument(
          Marmot::Utilities::MakeString()
          << __PRETTY_FUNCTION__ << ": axisymmetric element " << elLabel << ": a quadrature point lies at r = "
          << sectionFactor / ( 2.0 * Math::Constants::Pi ) << " <= 0; axisymmetric elements are meshed in r >= 0." );

      qp.J0xW = qp.weight * qp.detJ * sectionFactor;
    }
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::computeKernels( [[maybe_unused]] const double* QTotal,
                                                                  const double*                  dQ,
                                                                  double*                        Pe,
                                                                  double*                        Ke,
                                                                  double                         time,
                                                                  double                         dT )
  {
    integrateInternalForce( dQ, Pe, Ke, time, dT );
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::computeKernelsExplicit( [[maybe_unused]] const double* QTotal,
                                                                          const double*                  dQ,
                                                                          double*                        Pe,
                                                                          double                         time,
                                                                          double                         dT )
  {
    integrateInternalForce( dQ, Pe, nullptr, time, dT );
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::integrateInternalForce( const double* dQ,
                                                                          double*       Pe,
                                                                          double*       Ke,
                                                                          double        time,
                                                                          double        dT )
  {
    // energy gradients, if assigned, in the implicit update only
    const bool energyGradients = Ke != nullptr && dElasticEnergy_dQOutput != nullptr;
    if ( energyGradients ) {
      Eigen::Map< RhsSized >( dElasticEnergy_dQOutput ).setZero();
      Eigen::Map< RhsSized >( dDissipation_dQOutput ).setZero();
    }

    for ( QuadraturePoint& qp : qps ) {
      if constexpr ( nDim == 2 )
        if ( sectionType == SectionType::Axisymmetric ) {
          integrateQuadraturePoint< ContinuumMechanics::Voigt::Axial >( qp,
                                                                        axisymmetricB( qp ),
                                                                        dQ,
                                                                        Pe,
                                                                        Ke,
                                                                        energyGradients,
                                                                        time,
                                                                        dT );
          continue;
        }

      integrateQuadraturePoint< ParentGeometryElement::voigtSize >( qp, qp.B, dQ, Pe, Ke, energyGradients, time, dT );
    }
  }

  template < int nDim, int nNodes >
  template < ContinuumMechanics::Voigt::VoigtSize nComponents >
  void DisplacementFiniteElement< nDim, nNodes >::integrateQuadraturePoint(
    QuadraturePoint&                                            qp,
    const Eigen::Matrix< double, nComponents, sizeLoadVector >& B,
    const double*                                               dQ_,
    double*                                                     Pe_,
    double*                                                     Ke_,
    bool                                                        energyGradients,
    double                                                      time,
    double                                                      dT )
  {
    using namespace ContinuumMechanics::Voigt;
    using SectionVector = Eigen::Matrix< double, nComponents, 1 >;
    using SectionMatrix = Eigen::Matrix< double, nComponents, nComponents >;

    const Eigen::Map< const RhsSized > dQ( dQ_ );
    Eigen::Map< RhsSized >             Pe( Pe_ );

    const SectionVector dE = B * dQ;
    SectionVector       S  = SectionVector::Zero();
    SectionMatrix       C  = SectionMatrix::Zero();

    const Materials::MaterialHypoElastic::timeInfo timeInfo{ time, dT };
    double  elasticEnergyDensity = qp.J0xW > 0.0 ? qp.managedStateVars->elasticStrainEnergy / qp.J0xW : 0.0;
    double  dissipation          = qp.J0xW > 0.0 ? qp.managedStateVars->dissipation / qp.J0xW : 0.0;
    double* materialStateVars    = qp.managedStateVars->materialStateVars.data();

    bool planeStress = false;
    if constexpr ( nComponents == TwoD )
      planeStress = sectionType == SectionType::PlaneStress;

    if constexpr ( nComponents == OneD ) {
      if ( !Ke_ )
        throw std::runtime_error( "Explicit uniaxial stress not implemented yet" );

      Materials::MaterialHypoElastic::state1D state;
      state.stress               = qp.managedStateVars->stress( 0 );
      state.elasticEnergyDensity = elasticEnergyDensity;
      state.dissipation          = dissipation;
      state.stateVars            = materialStateVars;
      qp.material->computeUniaxialStress( state, C( 0, 0 ), dE( 0 ), timeInfo );

      S( 0 )                      = state.stress;
      qp.managedStateVars->stress = make3DVoigt< OneD >( S );
      elasticEnergyDensity        = state.elasticEnergyDensity;
      dissipation                 = state.dissipation;
    }
    else if ( planeStress ) {
      if constexpr ( nComponents == TwoD ) {
        Materials::MaterialHypoElastic::state2D state;
        state.stress               = voigtToPlaneVoigt( qp.managedStateVars->stress );
        state.elasticEnergyDensity = elasticEnergyDensity;
        state.dissipation          = dissipation;
        state.stateVars            = materialStateVars;
        qp.material->computePlaneStress( state, C, dE, timeInfo );

        S                           = state.stress;
        qp.managedStateVars->stress = make3DVoigt< TwoD >( S );
        elasticEnergyDensity        = state.elasticEnergyDensity;
        dissipation                 = state.dissipation;
      }
    }
    else {
      Math::Matrix6d C66;
      Math::Vector6d dElasticEnergyDensity_dStrain, dDissipation_dStrain;

      const auto state = updateMaterial3D( qp,
                                           make3DVoigt< nComponents >( dE ),
                                           Ke_ ? &C66 : nullptr,
                                           energyGradients ? &dElasticEnergyDensity_dStrain : nullptr,
                                           energyGradients ? &dDissipation_dStrain : nullptr,
                                           time,
                                           dT );

      S                    = sectionComponents< nComponents >( state.stress );
      elasticEnergyDensity = state.elasticEnergyDensity;
      dissipation          = state.dissipation;

      if ( Ke_ )
        C = sectionTangent< nComponents >( C66 );

      if ( energyGradients ) {
        Eigen::Map< RhsSized >(
          dElasticEnergy_dQOutput ) += B.transpose() *
                                       sectionComponents< nComponents >( dElasticEnergyDensity_dStrain ) * qp.J0xW;
        Eigen::Map< RhsSized >(
          dDissipation_dQOutput ) += B.transpose() * sectionComponents< nComponents >( dDissipation_dStrain ) * qp.J0xW;
      }
    }

    qp.managedStateVars->elasticStrainEnergy = elasticEnergyDensity * qp.J0xW;
    qp.managedStateVars->dissipation         = dissipation * qp.J0xW;
    qp.managedStateVars->totalStrainEnergy   = ( elasticEnergyDensity + dissipation ) * qp.J0xW;
    qp.managedStateVars->strain += make3DVoigt< nComponents >( dE );

    /* The artificial bulk viscosity is added to the stress that is INTEGRATED, never to the one
     * that is STORED. It is a numerical device: the constitutive law must not see it, it must
     * leave no trace in the state, and it must not appear in the reported stress. With inactive
     * coefficients nothing here is evaluated at all, so a run that does not ask for bulk
     * viscosity is bit-identical to one built before it existed.
     */
    if ( !Ke_ && bulkViscosityCoefficients.areActive() ) {

      // normal components: 11 (uniaxial); 11, 22 (plane); rr, zz, tt (axisymmetric); 11, 22, 33 (solid)
      constexpr int nNormalComponents = nComponents == OneD ? 1 : ( nComponents == TwoD ? 2 : 3 );

      /* In plane stress the out-of-plane strain is not carried by the element's kinematics, so
       * the trace is taken over the in-plane components only and the term is approximate there.
       * In 3D, in plane strain and axisymmetric, this IS the volumetric strain increment.
       */
      S.template head< nNormalComponents >()
        .array() += bulkViscousStress( qp,
                                       dE.template head< nNormalComponents >().sum(),
                                       elasticEnergyDensity,
                                       dissipation,
                                       dT );
    }

    if ( Ke_ )
      Eigen::Map< KeSizedMatrix >( Ke_ ) += B.transpose() * C * B * qp.J0xW;
    Pe += B.transpose() * S * qp.J0xW;
  }

  template < int nDim, int nNodes >
  Marmot::Materials::MaterialHypoElastic::state3D DisplacementFiniteElement< nDim, nNodes >::updateMaterial3D(
    QuadraturePoint&              qp,
    const Marmot::Math::Vector6d& dE,
    Marmot::Math::Matrix6d*       C,
    Marmot::Math::Vector6d*       dElasticEnergyDensity_dStrain,
    Marmot::Math::Vector6d*       dDissipation_dStrain,
    double                        time,
    double                        dT )
  {
    Materials::MaterialHypoElastic::state3D state( qp.managedStateVars->stress,
                                                   qp.J0xW > 0.0 ? qp.managedStateVars->elasticStrainEnergy / qp.J0xW
                                                                 : 0.0,
                                                   qp.J0xW > 0.0 ? qp.managedStateVars->dissipation / qp.J0xW : 0.0,
                                                   qp.managedStateVars->materialStateVars.data() );

    const Materials::MaterialHypoElastic::timeInfo timeInfo{ time, dT };

    if ( !C )
      qp.material->computeStressExplicit( state, dE, timeInfo );
    else if ( dElasticEnergyDensity_dStrain )
      qp.energyTangentProvider->computeStressWithEnergyTangents( state,
                                                                 *C,
                                                                 dE,
                                                                 timeInfo,
                                                                 *dElasticEnergyDensity_dStrain,
                                                                 *dDissipation_dStrain );
    else
      qp.material->computeStress( state, *C, dE, timeInfo );

    qp.managedStateVars->stress = state.stress;
    return state;
  }

  template < int nDim, int nNodes >
  double DisplacementFiniteElement< nDim, nNodes >::bulkViscousStress( QuadraturePoint& qp,
                                                                       double           volumetricStrainIncrement,
                                                                       double           elasticEnergyDensity,
                                                                       double           dissipation,
                                                                       double           dT )
  {
    const Marmot::Materials::MaterialHypoElastic::state3D stateForWaveSpeed( qp.managedStateVars->stress,
                                                                             elasticEnergyDensity,
                                                                             dissipation,
                                                                             qp.managedStateVars->materialStateVars
                                                                               .data() );

    if ( !qp.referenceWaveSpeed )
      qp.referenceWaveSpeed = qp.material->getMaximumWaveSpeed( stateForWaveSpeed );
    const double referenceWaveSpeed = *qp.referenceWaveSpeed;

    /* The optional degradation needs the material's CURRENT tangent, and asking for that costs
     * a full constitutive evaluation per quadrature point per increment unless the material
     * overrides getMaximumWaveSpeed() with a closed form. With the exponent at its default of
     * zero the material is never asked, so a deck that does not request the degradation
     * integrates exactly the stress it integrated before.
     */
    const double degradation = bulkViscosityCoefficients.isDegraded()
                                 ? Discretizations::Elements::BulkViscosity::
                                     degradationFactor( qp.material->getMaximumWaveSpeed( stateForWaveSpeed ),
                                                        referenceWaveSpeed,
                                                        bulkViscosityCoefficients.degradation )
                                 : 1.0;

    return degradation *
           Discretizations::Elements::BulkViscosity::
             viscousStressFromIncrement( volumetricStrainIncrement,
                                         dT,
                                         qp.material->getDensity( qp.managedStateVars->materialStateVars.data() ),
                                         referenceWaveSpeed,
                                         qp.characteristicLengthForCriticalTimeStep,
                                         bulkViscosityCoefficients );
  }

  template < int nDim, int nNodes >
  template < ContinuumMechanics::Voigt::VoigtSize nComponents >
  Eigen::Matrix< double, nComponents, 1 > DisplacementFiniteElement< nDim, nNodes >::sectionComponents(
    const Marmot::Math::Vector6d& voigt )
  {
    if constexpr ( nComponents == ContinuumMechanics::Voigt::Axial )
      return ContinuumMechanics::Voigt::voigtToAxisymmetricVoigt( voigt );
    else
      return ContinuumMechanics::Voigt::reduce3DVoigt< nComponents >( voigt );
  }

  template < int nDim, int nNodes >
  template < ContinuumMechanics::Voigt::VoigtSize nComponents >
  Eigen::Matrix< double, nComponents, nComponents > DisplacementFiniteElement< nDim, nNodes >::sectionTangent(
    const Marmot::Math::Matrix6d& C )
  {
    // gamma_rt = gamma_zt = 0 (axisymmetric) and eps_33 = 0 (plane strain) are kinematic, so the block is exact
    if constexpr ( nComponents == ContinuumMechanics::Voigt::TwoD )
      return ContinuumMechanics::PlaneStrain::getPlaneStrainTangent( C );
    else if constexpr ( nComponents == ContinuumMechanics::Voigt::Axial )
      return C.template topLeftCorner< 4, 4 >();
    else
      return C;
  }

  template < int nDim, int nNodes >
  typename DisplacementFiniteElement< nDim, nNodes >::BSizedAxisymmetric DisplacementFiniteElement< nDim, nNodes >::
    axisymmetricB( const QuadraturePoint& qp ) const
  {
    const NSized N = this->N( qp.xi );
    return Spatial2D::Axisymmetric::BFromPlanarB< nNodes >( qp.B, N, ( this->NB( N ) * this->coordinates )( 0 ) );
  }

  template < int nDim, int nNodes >
  double DisplacementFiniteElement< nDim, nNodes >::sectionFactorAt( const XiSized& xi ) const
  {
    if constexpr ( nDim == 3 )
      return 1.0;
    else if constexpr ( nDim == 2 ) {
      if ( sectionType == SectionType::Axisymmetric )
        return 2.0 * Marmot::Math::Constants::Pi * ( this->NB( this->N( xi ) ) * this->coordinates )( 0 );
      return elementProperties[0]; // thickness
    }
    else
      return elementProperties[0]; // cross section
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::setInitialConditions( StateTypes state, const double* values )
  {
    switch ( state ) {
    case Marmot::Discretizations::Elements::Element::MarmotMaterialInitialization: {
      for ( QuadraturePoint& qp : qps ) {
        qp.material->initializeYourself( qp.managedStateVars->materialStateVars.data(),
                                         static_cast< int >( qp.managedStateVars->materialStateVars.size() ) );
      }
      break;
    }
    case Marmot::Discretizations::Elements::Element::GeostaticStress: {
      if constexpr ( nDim >= 2 )
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
    case Marmot::Discretizations::Elements::Element::MarmotMaterialStateVars: {
      throw std::invalid_argument( "Please use initializeStateVars directly on material" );
    }
    default:
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__ << ": invalid initial condition" );
    }
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::computeDistributedLoad(
    Marmot::Discretizations::Elements::Element::DistributedLoadTypes loadType,
    double*                                                          P,
    [[maybe_unused]] double*                                         K,
    int                                                              elementFace,
    const double*                                                    load,
    [[maybe_unused]] const double*                                   QTotal,
    [[maybe_unused]] double                                          time,
    [[maybe_unused]] double                                          dT )
  {
    Eigen::Map< RhsSized > fU( P );

    // the boundary element of an axisymmetric element integrates over the ring; plane elements multiply by thickness
    const bool axisymmetric = sectionType == SectionType::Axisymmetric;
    const bool plane        = nDim == 2 && !axisymmetric;

    switch ( loadType ) {

    case Marmot::Discretizations::Elements::Element::Pressure: {
      const double p = load[0];

      Discretizations::Elements::BoundaryElement boundaryEl( this->shape,
                                                             elementFace,
                                                             nDim,
                                                             this->coordinates,
                                                             axisymmetric );

      Eigen::VectorXd Pk = -p * boundaryEl.computeSurfaceNormalVectorialLoadVector();

      if ( plane )
        Pk *= elementProperties[0]; // thickness

      boundaryEl.assembleIntoParentVectorial( Pk, fU );

      break;
    }
    case Marmot::Discretizations::Elements::Element::SurfaceTraction: {

      Discretizations::Elements::BoundaryElement boundaryEl( this->shape,
                                                             elementFace,
                                                             nDim,
                                                             this->coordinates,
                                                             axisymmetric );

      const Eigen::Map< const Eigen::VectorXd > tractionVector( load, nDim );

      auto Pk = boundaryEl.computeVectorialLoadVector( tractionVector );
      if ( plane )
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
  void DisplacementFiniteElement< nDim, nNodes >::computeBodyForce( double*                        P_,
                                                                    [[maybe_unused]] double*       K,
                                                                    const double*                  load,
                                                                    [[maybe_unused]] const double* QTotal,
                                                                    [[maybe_unused]] double        time,
                                                                    [[maybe_unused]] double        dT )
  {
    Eigen::Map< RhsSized >                                     Pe( P_ );
    const Eigen::Map< const Eigen::Matrix< double, nDim, 1 > > f( load );

    for ( const auto& qp : qps )
      Pe += this->NB( this->N( qp.xi ) ).transpose() * f * qp.J0xW;
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::computeConsistentInertia( double* M )
  {
    Eigen::Map< KeSizedMatrix > Me( M );
    Me.setZero();

    for ( const auto& point : Discretizations::Elements::ConsistentMass::integrationRule( this->shape ) ) {
      const XiSized xi = point.xi;
      const auto&   qp = qps[Discretizations::Elements::ConsistentMass::nearestQuadraturePoint( point.xi, qps )];
      // at the point of the mass rule, not at the nearest quadrature point: 2 pi r varies over the element
      const double sectionFactor = sectionFactorAt( xi );
      const double rho           = qp.material->getDensity( qp.managedStateVars->materialStateVars.data() );
      const double detJ          = this->Jacobian( this->dNdXi( xi ) ).determinant();
      const auto   N_            = this->NB( this->N( xi ) );
      Me += N_.transpose() * N_ * ( point.weight * detJ * sectionFactor ) * rho;
    }
  }
  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::computeLumpedInertia( double* M )
  {
    Eigen::Map< RhsSized > LMM( M );
    LMM.setZero();

    /* Row sums of the consistent mass matrix, from this element's own shape functions and from the
     * linear (corner-node) shape functions of the same element. Density-free: the blend weight and
     * the mass distribution are properties of the element geometry alone, and deriving them in one
     * place is what keeps this function and computeCriticalTimeStepForExplicitDynamics consistent
     * with each other -- a time step derived from a different mass distribution than the one
     * actually assembled is exactly the kind of inconsistency that surfaces as an unexplained
     * instability rather than as a clean failure.
     */
    constexpr int   nNodesLinear     = Marmot::Discretizations::Elements::MassLumping::nLinearNodes< nDim, nNodes >();
    auto            linGeometryEl    = Marmot::Discretizations::Elements::GeometryElement< nDim, nNodesLinear >();
    Eigen::VectorXd rowSumsHighOrder = Eigen::VectorXd::Zero( nNodes );
    Eigen::VectorXd rowSumsLinear    = Eigen::VectorXd::Zero( nNodesLinear );
    for ( const auto& qp : qps ) {
      rowSumsHighOrder += Eigen::VectorXd( this->N( qp.xi ) ) * qp.J0xW;
      rowSumsLinear += Eigen::VectorXd( linGeometryEl.N( qp.xi ) ) * qp.J0xW;
    }
    const double weight = Discretizations::Elements::MassLumping::manifoldBlendWeight( rowSumsHighOrder,
                                                                                       rowSumsLinear );

    for ( const auto& qp : qps ) {
      const auto N_    = this->N( qp.xi );
      const auto N_lin = linGeometryEl.N( qp.xi );

      Eigen::VectorXd N_weighted = weight * ( N_ );
      N_weighted.template head< nNodesLinear >() += ( 1.0 - weight ) * N_lin;

      const double    rho = qp.material->getDensity( qp.managedStateVars->materialStateVars.data() );
      Eigen::VectorXd m_  = N_weighted * qp.J0xW * rho;
      for ( int i = 0; i < nNodes; i++ ) {
        for ( int d = 0; d < nDim; d++ )
          LMM( i * nDim + d ) += m_( i );
      }
    }
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::computeCriticalTimeStepForExplicitDynamics(
    double&                        criticalTimeStep,
    [[maybe_unused]] const double* QTotal )
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
    constexpr int   nNodesLinear     = Marmot::Discretizations::Elements::MassLumping::nLinearNodes< nDim, nNodes >();
    auto            linGeometryEl    = Marmot::Discretizations::Elements::GeometryElement< nDim, nNodesLinear >();
    Eigen::VectorXd rowSumsHighOrder = Eigen::VectorXd::Zero( nNodes );
    Eigen::VectorXd rowSumsLinear    = Eigen::VectorXd::Zero( nNodesLinear );
    for ( const auto& qp : qps ) {
      rowSumsHighOrder += Eigen::VectorXd( this->N( qp.xi ) ) * qp.J0xW;
      rowSumsLinear += Eigen::VectorXd( linGeometryEl.N( qp.xi ) ) * qp.J0xW;
    }
    const double weight = Discretizations::Elements::MassLumping::manifoldBlendWeight( rowSumsHighOrder,
                                                                                       rowSumsLinear );
    const double lumpedMassTimeStepFactor = Discretizations::Elements::MassLumping::timeStepFactorFromMassDistribution(
      Discretizations::Elements::MassLumping::manifoldMassFractions( rowSumsHighOrder, rowSumsLinear, weight ) );

    criticalTimeStep = std::numeric_limits< double >::max();
    for ( const auto& qp : qps ) {
      /* The characteristic element length of GeometryElement: the smallest physical extent of a
       * quadrilateral or hexahedron, 1 / |grad N|_F of a linear simplex; see
       * GeometryElement::characteristicLengthForCriticalTimeStepAt.
       *
       * Read from the cache initializeYourself() fills rather than recomputed here: the artificial
       * bulk viscosity already needs this value on every explicit increment and caches it for that
       * reason, and the nodal coordinates it depends on do not change, so a second computation here
       * would only disagree with the cached one under floating-point noise, never under a real
       * update.
       */
      const double characteristicLengthForCriticalTimeStep = qp.characteristicLengthForCriticalTimeStep;

      Marmot::Materials::MaterialHypoElastic::state3D state( qp.managedStateVars->stress,
                                                             qp.managedStateVars->elasticStrainEnergy / qp.J0xW,
                                                             qp.managedStateVars->dissipation / qp.J0xW,
                                                             qp.managedStateVars->materialStateVars.data() );

      const double c = qp.material->getMaximumWaveSpeed( state );
      if ( c <= 0.0 )
        throw std::runtime_error( "Non-positive wave speed encountered in computeCriticalTimeStepForExplicitDynamics" );
      const double dt = lumpedMassTimeStepFactor * characteristicLengthForCriticalTimeStep / c;
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
  void DisplacementFiniteElement< nDim, nNodes >::assignEnergyGradientOutput( double* dElasticEnergy_dQ,
                                                                              double* dDissipation_dQ )
  {
    if ( ( dElasticEnergy_dQ == nullptr ) != ( dDissipation_dQ == nullptr ) )
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__ << ": assign both energy gradient outputs, or none." );

    if ( dElasticEnergy_dQ )
      checkEnergyGradientSupport();

    dElasticEnergy_dQOutput = dElasticEnergy_dQ;
    dDissipation_dQOutput   = dDissipation_dQ;
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::checkEnergyGradientSupport() const
  {
    bool sectionSupported = false;
    if constexpr ( nDim == 3 )
      sectionSupported = sectionType == SectionType::Solid;
    else if constexpr ( nDim == 2 )
      sectionSupported = sectionType == SectionType::PlaneStrain || sectionType == SectionType::Axisymmetric;

    if ( !sectionSupported )
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__ << ": element " << elLabel
                                   << ": energy gradients are available for solid, plane strain and axisymmetric "
                                      "sections only." );

    for ( const auto& qp : qps ) {
      if ( !qp.material )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__ << ": element " << elLabel
                                     << ": assign the material section before the energy gradient outputs." );
      if ( !qp.energyTangentProvider )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__ << ": element " << elLabel
                                     << ": the material provides no energy tangents "
                                        "(Marmot::Materials::HypoElasticEnergyTangentProvider)." );
    }
  }

  template < int nDim, int nNodes >
  std::vector< double > DisplacementFiniteElement< nDim, nNodes >::getCoordinatesAtCenter()
  {
    std::vector< double > coords( nDim );

    Eigen::Map< XiSized > coordsMap( &coords[0] );
    coordsMap = this->center();
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
    return static_cast< int >( qps.size() );
  }
} // namespace Marmot::Discretizations::Elements
