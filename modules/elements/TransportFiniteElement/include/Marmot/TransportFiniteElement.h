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
#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryElement.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotTransportMaterial.h"
#include "Marmot/MarmotTransportMaterialFactory.h"
#include "Marmot/MarmotTypedefs.h"
#include <array>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace Marmot::Discretizations::Elements {

  /**
   * @class Marmot::Discretizations::Elements::TransportFiniteElement
   * @tparam nDim Number of spatial dimensions (2 or 3).
   * @tparam nNodes Number of element nodes.
   * @tparam nFields Number of transported scalar fields per node (1 for heat conduction or moisture diffusion, 2 for
   * coupled heat and moisture transport).
   * @brief Continuum element for (coupled) transient scalar transport, e.g. heat conduction or moisture diffusion.
   *
   * @details The element integrates the balance \f$\dot m_a + \nabla\cdot\mathbf{q}_a = r_a + s_a\f$ of each field
   * \f$\phi_a\f$, with the content \f$m_a\f$, the flux \f$\mathbf{q}_a\f$ and the material source \f$r_a\f$ from a
   * Marmot::Materials::TransportMaterial, and the external body source \f$s_a\f$ from computeBodyForce(). The weak form
   * \f[
   *   \int_\Omega \delta\phi_a\,\dot m_a\,\mathrm{d}V - \int_\Omega \nabla\delta\phi_a\cdot\mathbf{q}_a\,\mathrm{d}V
   *   - \int_\Omega \delta\phi_a\,r_a\,\mathrm{d}V
   *   = \int_\Omega \delta\phi_a\,s_a\,\mathrm{d}V - \int_{\Gamma} \delta\phi_a\,q_{n,a}\,\mathrm{d}\Gamma
   * \f]
   * with the outward normal flux \f$q_{n,a}\f$ on the boundary is integrated in time with the backward Euler method
   * in conservative form, \f$\dot m_a \approx (m_a(\boldsymbol{\phi}_{n+1}) - m_a(\boldsymbol{\phi}_n))/\Delta t\f$,
   * with \f$\boldsymbol{\phi}_n\f$ = QTotal - dQ and \f$\boldsymbol{\phi}_{n+1}\f$ = QTotal.
   *
   * A zero time increment \f$\Delta t = 0\f$ (e.g. the zero increment a solver makes at the beginning of a step) is
   * the limit \f$\Delta t \to 0\f$ of the balance multiplied by \f$\Delta t\f$: no time passes, so nothing flows
   * and the content cannot change, \f$\mathbf{P}_e = \sum_{qp}\mathbf{N}_\phi^\mathsf{T}\Delta\mathbf{m}\,J_0 w\f$,
   * and the distributed and body loads vanish. Prescribed field values then spread by the capacity-weighted
   * projection only. A field without storage (Marmot::Materials::TransportMaterial::hasStorage(), e.g. a zero heat
   * capacity) instead obeys the steady state equations, whatever the time increment. The distinction is made field by
   * field: in a coupled element, the rows of a field without storage are steady, the others instantaneous.
   *
   * The dofs are ordered node by node, the fields of a node in the order of the field names given to the
   * constructor: dof \f$i\,n_f + a\f$ is field \f$a\f$ at node \f$i\f$. With the interpolation operator
   * \f$\mathbf{N}_\phi\f$ (\f$\boldsymbol{\phi} = \mathbf{N}_\phi\mathbf{Q}\f$) and the gradient operator
   * \f$\mathbf{B}_\phi\f$ (\f$\nabla\boldsymbol{\phi} = \mathbf{B}_\phi\mathbf{Q}\f$, three rows per field), the
   * residual (internal flux vector) and its tangent are
   * \f[
   *   \mathbf{P}_e = \sum_{qp} \left[ \mathbf{N}_\phi^\mathsf{T}\left(\frac{\Delta\mathbf{m}}{\Delta t}
   *                  - \mathbf{r}\right) - \mathbf{B}_\phi^\mathsf{T}\mathbf{q} \right] J_0 w,\qquad
   *   \mathbf{K}_e = \sum_{qp} \left[ \mathbf{N}_\phi^\mathsf{T}\left(\frac{1}{\Delta t}
   *                  \frac{\partial\Delta\mathbf{m}}{\partial\boldsymbol{\phi}}
   *                  - \frac{\partial\mathbf{r}}{\partial\boldsymbol{\phi}}\right)\mathbf{N}_\phi
   *                  - \mathbf{B}_\phi^\mathsf{T}\left(\frac{\partial\mathbf{q}}{\partial\boldsymbol{\phi}}
   *                  \mathbf{N}_\phi + \frac{\partial\mathbf{q}}{\partial\nabla\boldsymbol{\phi}}
   *                  \mathbf{B}_\phi\right) \right] J_0 w,
   * \f]
   * \f$J_0 w = w \det\mathbf{J}\, s\f$ with the section factor \f$s\f$: the thickness (plane), \f$2\pi r\f$
   * (axisymmetric, \f$r = x_1\f$), or 1 (solid).
   */
  template < int nDim, int nNodes, int nFields >
  class TransportFiniteElement : public Marmot::Discretizations::Elements::Element,
                                 public Marmot::Discretizations::Elements::GeometryElement< nDim, nNodes > {

    static_assert( nDim == 2 || nDim == 3, "Transport elements are two- or three-dimensional." );
    static_assert( nFields >= 1 && nFields <= Marmot::Materials::TransportMaterial::maxNumberOfFields,
                   "Unsupported number of fields." );

  public:
    /**
     * @brief The section assumption.
     * - Plane: 2D, integrated over the thickness (element property 0).
     * - Axisymmetric: 2D, \f$r = x_1 \geq 0\f$, \f$z = x_2\f$, integrated over the full ring.
     * - Solid: 3D.
     */
    enum SectionType {
      Plane,
      Axisymmetric,
      Solid,
    };

    static constexpr int sizeLoadVector = nNodes * nFields; ///< dofs per element
    static constexpr int nGradient      = 3 * nFields;      ///< gradient (and flux) components, three per field

    using ParentGeometryElement = Marmot::Discretizations::Elements::GeometryElement< nDim, nNodes >;
    using JacobianSized         = typename ParentGeometryElement::JacobianSized;
    using dNdXiSized            = typename ParentGeometryElement::dNdXiSized;
    using XiSized               = typename ParentGeometryElement::XiSized;
    using NSized                = typename ParentGeometryElement::NSized;
    using RhsSized              = Eigen::Matrix< double, sizeLoadVector, 1 >;
    using KeSizedMatrix         = Eigen::Matrix< double, sizeLoadVector, sizeLoadVector >;
    using NFieldSized           = Eigen::Matrix< double, nFields, sizeLoadVector >;   ///< \f$\mathbf{N}_\phi\f$
    using BFieldSized           = Eigen::Matrix< double, nGradient, sizeLoadVector >; ///< \f$\mathbf{B}_\phi\f$
    using FieldNames            = std::array< std::string, nFields >;

    /// Element properties: the thickness of a plane element
    Eigen::Map< const Eigen::VectorXd > elementProperties;
    /// Element label (ID)
    const int elLabel;
    /// Section assumption
    const SectionType sectionType;
    /// The names of the nodal fields, e.g. "temperature"
    const FieldNames fieldNames;
    /// Whether the material stores content of a field; false: steady state transport of that field, see the class
    /// description
    std::array< bool, nFields > withStorage;

    /**
     * @brief Whether a time increment lets no time pass for a field: then, with storage, nothing of it flows (see the
     * class description).
     * @param dT The time increment.
     * @param field The index of the field.
     * @return True for a zero time increment and a field with storage.
     */
    bool isInstantaneous( double dT, int field ) const { return withStorage[field] && !( dT > 0.0 ); }

    /// Data and state of a quadrature point
    struct QuadraturePoint {

      const XiSized xi;     ///< parametric coordinates
      const double  weight; ///< quadrature weight

      double      J0xW;     ///< \f$J_0 w\f$ including the section factor
      NFieldSized NPhi;     ///< interpolation operator
      BFieldSized BPhi;     ///< gradient operator

      /**
       * @brief The state variables of a quadrature point: [field (nFields), field gradient (3 nFields), flux
       * (3 nFields), material state (...)].
       */
      class QPStateVarManager : public Marmot::Utilities::StateVarVectorManager {

        /// \hideinitializer
        inline const static auto layout = makeLayout( {
          { .name = "field", .length = nFields },
          { .name = "field gradient", .length = 3 * nFields },
          { .name = "flux", .length = 3 * nFields },
          { .name = "begin of material state", .length = 0 },
        } );

      public:
        Eigen::Map< Eigen::Matrix< double, nFields, 1 > >   field;             ///< fields at the end of the increment
        Eigen::Map< Eigen::Matrix< double, nGradient, 1 > > gradient;          ///< their gradients
        Eigen::Map< Eigen::Matrix< double, nGradient, 1 > > flux;              ///< the fluxes
        Eigen::Map< Eigen::VectorXd >                       materialStateVars; ///< state of the material

        /// @return The number of state variables of the quadrature point without the material's
        static int getNumberOfRequiredStateVarsQuadraturePointOnly() { return layout.nRequiredStateVars; };

        /**
         * @brief Map the state variables of a quadrature point.
         * @param theStateVarVector The state variables.
         * @param nStateVars Their number.
         */
        QPStateVarManager( double* theStateVarVector, int nStateVars )
          : Marmot::Utilities::StateVarVectorManager( theStateVarVector, layout ),
            field( &find( "field" ) ),
            gradient( &find( "field gradient" ) ),
            flux( &find( "flux" ) ),
            materialStateVars( &find( "begin of material state" ),
                               nStateVars - getNumberOfRequiredStateVarsQuadraturePointOnly() ){};
      };

      std::unique_ptr< QPStateVarManager >                    managedStateVars; ///< the state variables
      std::unique_ptr< Marmot::Materials::TransportMaterial > material;         ///< the material

      /// @return The number of state variables of the quadrature point including the material's
      int getNumberOfRequiredStateVars()
      {
        return QPStateVarManager::getNumberOfRequiredStateVarsQuadraturePointOnly() +
               material->getNumberOfRequiredStateVars();
      };

      /**
       * @brief Map the state variables.
       * @param stateVars The state variables of the quadrature point.
       * @param nStateVars Their number.
       */
      void assignStateVars( double* stateVars, int nStateVars )
      {
        managedStateVars = std::make_unique< QPStateVarManager >( stateVars, nStateVars );
      }

      /**
       * @brief Construct the quadrature point.
       * @param xi Parametric coordinates.
       * @param weight Quadrature weight.
       */
      QuadraturePoint( XiSized xi, double weight )
        : xi( xi ), weight( weight ), J0xW( 0.0 ), NPhi( NFieldSized::Zero() ), BPhi( BFieldSized::Zero() ){};
    };

    /// Quadrature points
    std::vector< QuadraturePoint > qps;

    /**
     * @brief Construct the element.
     * @param elementID Element label.
     * @param quadratureRule Quadrature rule.
     * @param sectionType Section assumption.
     * @param fieldNames The names of the nodal fields, e.g. { "temperature" }.
     * @throws std::invalid_argument for a section not matching nDim.
     */
    TransportFiniteElement( int                                                   elementID,
                            Discretizations::Elements::Quadrature::QuadratureRule quadratureRule,
                            SectionType                                           sectionType,
                            const FieldNames&                                     fieldNames );

    /** @brief Total number of state variables of the element. */
    int getNumberOfRequiredStateVars() override
    {
      return qps[0].getNumberOfRequiredStateVars() * static_cast< int >( qps.size() );
    }

    /** @brief The field names, the same at each node. */
    std::vector< std::vector< std::string > > getNodeFields() override
    {
      return std::vector< std::vector< std::string > >( nNodes,
                                                        std::vector< std::string >( fieldNames.begin(),
                                                                                    fieldNames.end() ) );
    }

    /** @brief Identity: the element orders its dofs node by node already. */
    std::vector< int > getDofIndicesPermutationPattern() override
    {
      std::vector< int > pattern( sizeLoadVector );
      for ( int i = 0; i < sizeLoadVector; i++ )
        pattern[i] = i;
      return pattern;
    }

    /** @brief Number of nodes. */
    int getNNodes() override { return nNodes; }

    /** @brief Number of spatial dimensions. */
    int getNSpatialDimensions() override { return nDim; }

    /** @brief Number of dofs (nNodes * nFields). */
    int getNDofPerElement() override { return sizeLoadVector; }

    /** @brief Ensight shape of the element. */
    std::string getElementShape() override { return ParentGeometryElement::getElementShape(); }

    /** @brief Map the element state variables to the quadrature points. */
    void assignStateVars( double* stateVars, int nStateVars ) override;

    /** @brief Assign the element properties: the thickness of a plane element. */
    void assignProperty( const Marmot::Discretizations::ElementProperties& elementPropertiesInfo ) override;

    /**
     * @brief Create the transport materials at the quadrature points.
     * @throws std::invalid_argument if the material is no transport material or describes another number of fields.
     */
    void assignProperty( const Marmot::Discretizations::MaterialSection& section ) override;

    /** @brief Map the nodal coordinates. */
    void assignNodeCoordinates( const double* coordinates ) override
    {
      ParentGeometryElement::assignNodeCoordinates( coordinates );
    }

    /**
     * @brief Precompute the operators and \f$J_0 w\f$ at the quadrature points.
     * @throws std::invalid_argument for a plane element without thickness, an axisymmetric element with a node at
     * \f$r < 0\f$, or a non-positive \f$J_0 w\f$.
     */
    void initializeYourself() override;

    /**
     * @brief The section factor at @p xi.
     * @param xi Parametric coordinates.
     * @return The thickness (plane), \f$2\pi r\f$ (axisymmetric), or 1 (solid).
     */
    double sectionFactorAt( const XiSized& xi ) const;

    /**
     * @brief Initial conditions: only MarmotMaterialInitialization is supported.
     * @param state The state type.
     * @param values Not used.
     */
    void setInitialConditions( StateTypes state, const double* values ) override;

    /**
     * @brief Residual (internal flux vector) and its tangent, see the class description.
     * @param QTotal Nodal fields at the end of the increment.
     * @param dQ Their increment.
     * @param Pe Residual (accumulated).
     * @param Ke Tangent (accumulated).
     * @param time Time at the end of the increment (the current time).
     * @param dT Time increment, see the class description for a zero one.
     */
    void computeKernels( const double* QTotal, const double* dQ, double* Pe, double* Ke, double time, double dT )
      override;

    /**
     * @brief Surface flux loads: HeatFlux, HeatFilm and HeatRadiation, acting on the field "temperature", and
     * MoistureFlux and MoistureFilm, acting on the field "relative humidity".
     * @details The EXTERNAL load \f$\mathbf{P} = -\int_\Gamma \mathbf{N}^\mathsf{T} q_n\,\mathrm{d}\Gamma\f$ (inflow
     * positive) is added to @p P, the tangent of the residual \f$-\partial\mathbf{P}/\partial\mathbf{Q}\f$ to @p K:
     * - HeatFlux, load = [\f$\bar q\f$]: \f$q_n = -\bar q\f$ (\f$\bar q\f$ flows into the body),
     * - HeatFilm, load = [\f$h\f$, \f$T_\infty\f$]: \f$q_n = h\,(T - T_\infty)\f$,
     * - HeatRadiation, load = [\f$c\f$, \f$T_\infty\f$, \f$T_0\f$]:
     *   \f$q_n = c\,((T - T_0)^4 - (T_\infty - T_0)^4)\f$, \f$c = \varepsilon\sigma\f$, \f$T_0\f$ the absolute zero
     *   of the temperature scale used,
     * - MoistureFlux, load = [\f$\bar q\f$]: \f$q_n = -\bar q\f$ (\f$\bar q\f$ flows into the body),
     * - MoistureFilm, load = [\f$\beta\f$, \f$h_\infty\f$]: \f$q_n = \beta\,(h - h_\infty)\f$, the moisture
     *   exchange with an environment of relative humidity \f$h_\infty\f$ (surface emissivity / film coefficient
     *   \f$\beta\f$, in units of the moisture flux of the material per unit relative humidity).
     *
     * Plane elements multiply by the thickness, axisymmetric ones integrate over the ring. Nothing is added for a
     * zero time increment of an element with storage (see the class description).
     * @param loadType The load type.
     * @param P External load vector (accumulated).
     * @param K Tangent of the residual (accumulated).
     * @param elementFace Face (1-based).
     * @param load The load values, see above.
     * @param QTotal Nodal fields at the end of the increment.
     * @param time Not used.
     * @param dT Time increment.
     * @throws std::invalid_argument for another load type, or if the element has no field the load acts on.
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
     * @brief Volumetric source: \f$\mathbf{P} = \int_\Omega \mathbf{N}_\phi^\mathsf{T}\mathbf{s}\,\mathrm{d}V\f$.
     * @param P External load vector (accumulated).
     * @param K Not used.
     * @param load The source per unit volume, one value per field.
     * @param QTotal Not used.
     * @param time Not used.
     * @param dT Time increment: nothing is added to the fields with storage for a zero time increment.
     */
    void computeBodyForce( double* P, double* K, const double* load, const double* QTotal, double time, double dT )
      override;

    /**
     * @brief Access a state of a quadrature point: "field", "field gradient", "flux", or a state of the material.
     * @param stateName The name.
     * @param qpNumber The quadrature point.
     * @return The view.
     */
    Marmot::Utilities::StateView getStateView( const std::string& stateName, int qpNumber ) override
    {
      const auto& qp = qps[qpNumber];
      if ( qp.managedStateVars->contains( stateName ) )
        return qp.managedStateVars->getStateView( stateName );
      return qp.material->getStateView( stateName, qp.managedStateVars->materialStateVars.data() );
    }

    /** @brief Coordinates of the element center. */
    std::vector< double > getCoordinatesAtCenter() override;

    /** @brief Coordinates of the quadrature points. */
    std::vector< std::vector< double > > getCoordinatesAtQuadraturePoints() override;

    /** @brief Number of quadrature points. */
    int getNumberOfQuadraturePoints() override { return static_cast< int >( qps.size() ); }

    /**
     * @brief The index of the field a distributed load acts on.
     * @param loadType The load type.
     * @return The index into fieldNames.
     * @throws std::invalid_argument for a load type that is no transport load, or if the element has no such field.
     */
    int fieldIndexOfLoad( Marmot::Discretizations::Elements::Element::DistributedLoadTypes loadType ) const;
  };

  template < int nDim, int nNodes, int nFields >
  TransportFiniteElement< nDim, nNodes, nFields >::TransportFiniteElement(
    int                                                   elementID,
    Discretizations::Elements::Quadrature::QuadratureRule quadratureRule,
    SectionType                                           sectionType,
    const FieldNames&                                     fieldNames )
    : ParentGeometryElement(),
      elementProperties( Eigen::Map< const Eigen::VectorXd >( nullptr, 0 ) ),
      elLabel( elementID ),
      sectionType( sectionType ),
      fieldNames( fieldNames )
  {
    withStorage.fill( true );
    if ( ( nDim == 3 ) != ( sectionType == SectionType::Solid ) )
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__
                                   << ": plane and axisymmetric sections are two-dimensional, the solid section is "
                                      "three-dimensional." );

    for ( const auto& qpInfo :
          Discretizations::Elements::Quadrature::getQuadraturePoints( this->referenceDomain, quadratureRule ) )
      qps.push_back( QuadraturePoint( qpInfo.xi, qpInfo.weight ) );
  }

  template < int nDim, int nNodes, int nFields >
  void TransportFiniteElement< nDim, nNodes, nFields >::assignStateVars( double* stateVars, int nStateVars )
  {
    const int nQpStateVars = nStateVars / static_cast< int >( qps.size() );
    for ( size_t i = 0; i < qps.size(); i++ )
      qps[i].assignStateVars( stateVars + i * nQpStateVars, nQpStateVars );
  }

  template < int nDim, int nNodes, int nFields >
  void TransportFiniteElement< nDim, nNodes, nFields >::assignProperty(
    const Marmot::Discretizations::ElementProperties& elementPropertiesInfo )
  {
    new ( &elementProperties ) Eigen::Map< const Eigen::VectorXd >( elementPropertiesInfo.elementProperties,
                                                                    elementPropertiesInfo.nElementProperties );
  }

  template < int nDim, int nNodes, int nFields >
  void TransportFiniteElement< nDim, nNodes, nFields >::assignProperty(
    const Marmot::Discretizations::MaterialSection& section )
  {
    // the new materials replace the current ones only once all of them are accepted
    std::vector< std::unique_ptr< Marmot::Materials::TransportMaterial > > materials;
    materials.reserve( qps.size() );
    for ( size_t i = 0; i < qps.size(); i++ ) {
      auto material = std::unique_ptr< Marmot::Materials::TransportMaterial >(
        Marmot::Registry::TransportMaterialFactory::createMaterial( section.materialName,
                                                                    section.materialProperties,
                                                                    section.nMaterialProperties,
                                                                    elLabel ) );
      if ( material->getNumberOfFields() != nFields )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__ << ": element " << elLabel << " transports " << nFields
                                     << " field(s), but material " << section.materialName << " describes "
                                     << material->getNumberOfFields() << "." );
      materials.push_back( std::move( material ) );
    }

    withStorage.fill( false );
    for ( size_t i = 0; i < qps.size(); i++ ) {
      for ( int a = 0; a < nFields; a++ )
        withStorage[a] = withStorage[a] || materials[i]->hasStorage( a );
      qps[i].material = std::move( materials[i] );
    }
  }

  template < int nDim, int nNodes, int nFields >
  double TransportFiniteElement< nDim, nNodes, nFields >::sectionFactorAt( const XiSized& xi ) const
  {
    switch ( sectionType ) {
    case SectionType::Plane: {
      if ( elementProperties.size() < 1 )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__ << ": plane element " << elLabel
                                     << " has no thickness assigned (element property 0)." );
      return elementProperties[0];
    }
    case SectionType::Axisymmetric:
      return 2.0 * Marmot::Math::Constants::Pi * ( this->NB( this->N( xi ) ) * this->coordinates )( 0 );
    default: return 1.0;
    }
  }

  template < int nDim, int nNodes, int nFields >
  void TransportFiniteElement< nDim, nNodes, nFields >::initializeYourself()
  {
    if ( sectionType == SectionType::Axisymmetric ) {
      const Eigen::Map< const Eigen::Matrix< double, nDim, nNodes > > nodes( this->coordinates.data() );
      const double h = ( nodes.rowwise().maxCoeff() - nodes.rowwise().minCoeff() ).norm();
      if ( nodes.row( 0 ).minCoeff() < -1e-12 * h )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__ << ": axisymmetric element " << elLabel
                                     << " has a node at r = x_1 < 0; axisymmetric elements are meshed in r >= 0." );
    }

    for ( QuadraturePoint& qp : qps ) {
      const NSized        N     = this->N( qp.xi );
      const dNdXiSized    dNdXi = this->dNdXi( qp.xi );
      const JacobianSized J     = this->Jacobian( dNdXi );
      const dNdXiSized    dNdX  = this->dNdX( dNdXi, J.inverse() );

      qp.J0xW = qp.weight * J.determinant() * sectionFactorAt( qp.xi );
      if ( !( qp.J0xW > 0.0 ) )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__ << ": element " << elLabel
                                     << ": non-positive integration weight J0 w = " << qp.J0xW
                                     << " (distorted or inverted element, zero thickness, or r <= 0)." );

      qp.NPhi.setZero();
      qp.BPhi.setZero();
      for ( int j = 0; j < nNodes; j++ )
        for ( int a = 0; a < nFields; a++ ) {
          qp.NPhi( a, j * nFields + a ) = N( j );
          for ( int d = 0; d < nDim; d++ )
            qp.BPhi( 3 * a + d, j * nFields + a ) = dNdX( d, j );
        }
    }
  }

  template < int nDim, int nNodes, int nFields >
  void TransportFiniteElement< nDim, nNodes, nFields >::setInitialConditions( StateTypes                     state,
                                                                              [[maybe_unused]] const double* values )
  {
    switch ( state ) {
    case Marmot::Discretizations::Elements::Element::MarmotMaterialInitialization: {
      for ( QuadraturePoint& qp : qps )
        qp.material->initializeYourself( qp.managedStateVars->materialStateVars.data(),
                                         static_cast< int >( qp.managedStateVars->materialStateVars.size() ) );
      break;
    }
    default:
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__ << ": invalid initial condition" );
    }
  }

  template < int nDim, int nNodes, int nFields >
  void TransportFiniteElement< nDim, nNodes, nFields >::computeKernels( const double* QTotal_,
                                                                        const double* dQ_,
                                                                        double*       Pe_,
                                                                        double*       Ke_,
                                                                        double        time,
                                                                        double        dT )
  {
    using Material = Marmot::Materials::TransportMaterial;

    const Eigen::Map< const RhsSized > QTotal( QTotal_ );
    const Eigen::Map< const RhsSized > dQ( dQ_ );
    Eigen::Map< RhsSized >             Pe( Pe_ );
    Eigen::Map< KeSizedMatrix >        Ke( Ke_ );

    const RhsSized QOld = QTotal - dQ;

    // field by field: an instantaneous field (zero time increment, with storage) balances its content increment and
    // does not flow; the others balance rate, source and flux, with rate 1 / dT (with storage, dT > 0) or 0 (steady)
    Eigen::Matrix< double, nFields, 1 >   rateFactor;
    Eigen::Matrix< double, nFields, 1 >   flowing;
    Eigen::Matrix< double, nGradient, 1 > flowingGradient;
    for ( int a = 0; a < nFields; a++ ) {
      const bool instantaneous = isInstantaneous( dT, a );
      rateFactor( a )          = instantaneous ? 1.0 : ( withStorage[a] ? 1.0 / dT : 0.0 );
      flowing( a )             = instantaneous ? 0.0 : 1.0;
      flowingGradient.template segment< 3 >( 3 * a ).setConstant( flowing( a ) );
    }

    Material::FieldState                fieldState;
    Material::Response                  response;
    Material::AlgorithmicTangent        tangent;
    const Material::TimeIncrement       timeIncrement{ time, dT };
    Eigen::Matrix< double, nFields, 1 > rate;

    for ( QuadraturePoint& qp : qps ) {
      auto& state = *qp.managedStateVars;

      fieldState.field    = qp.NPhi * QTotal;
      fieldState.fieldOld = qp.NPhi * QOld;
      fieldState.gradient = qp.BPhi * QTotal;

      qp.material->initializeResponse( response, tangent );
      response.stateVars = state.materialStateVars.data();

      qp.material->computeTransport( response, tangent, fieldState, timeIncrement );

      state.field    = fieldState.field;
      state.gradient = fieldState.gradient;
      state.flux     = response.flux;

      // rows of instantaneous fields (dT = 0): the balance multiplied by dT, nothing flows, the content does not change
      const Eigen::Matrix< double, nFields, nFields > capacity = rateFactor.asDiagonal() *
                                                                   tangent.dContentIncrement_dField -
                                                                 flowing.asDiagonal() * tangent.dSource_dField;
      const Eigen::Matrix< double, nGradient, nFields >   dQ_dPhi = flowingGradient.asDiagonal() * tangent.dFlux_dField;
      const Eigen::Matrix< double, nGradient, nGradient > dQ_dGradPhi = flowingGradient.asDiagonal() *
                                                                        tangent.dFlux_dGradient;
      const Eigen::Matrix< double, nGradient, 1 > flux = flowingGradient.cwiseProduct( response.flux );

      rate = rateFactor.cwiseProduct( response.contentIncrement ) - flowing.cwiseProduct( response.source );

      Pe += ( qp.NPhi.transpose() * rate - qp.BPhi.transpose() * flux ) * qp.J0xW;
      Ke += ( qp.NPhi.transpose() * capacity * qp.NPhi -
              qp.BPhi.transpose() * ( dQ_dPhi * qp.NPhi + dQ_dGradPhi * qp.BPhi ) ) *
            qp.J0xW;
    }
  }

  template < int nDim, int nNodes, int nFields >
  int TransportFiniteElement< nDim, nNodes, nFields >::fieldIndexOfLoad(
    Marmot::Discretizations::Elements::Element::DistributedLoadTypes loadType ) const
  {
    // each transport load acts on the field of its physics; extend for the loads of further fields here
    std::string fieldName;
    switch ( loadType ) {
    case Marmot::Discretizations::Elements::Element::HeatFlux:
    case Marmot::Discretizations::Elements::Element::HeatFilm:
    case Marmot::Discretizations::Elements::Element::HeatRadiation: fieldName = "temperature"; break;
    case Marmot::Discretizations::Elements::Element::MoistureFlux:
    case Marmot::Discretizations::Elements::Element::MoistureFilm: fieldName = "relative humidity"; break;
    default:
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__ << ": element " << elLabel
                                   << ": load type not supported by transport elements." );
    }

    for ( int a = 0; a < nFields; a++ )
      if ( fieldNames[a] == fieldName )
        return a;

    throw std::invalid_argument( Marmot::Utilities::MakeString()
                                 << __PRETTY_FUNCTION__ << ": element " << elLabel << " has no field '" << fieldName
                                 << "' the load could act on." );
  }

  template < int nDim, int nNodes, int nFields >
  void TransportFiniteElement< nDim, nNodes, nFields >::computeDistributedLoad(
    Marmot::Discretizations::Elements::Element::DistributedLoadTypes loadType,
    double*                                                          P_,
    double*                                                          K_,
    int                                                              elementFace,
    const double*                                                    load,
    const double*                                                    QTotal_,
    [[maybe_unused]] double                                          time,
    double                                                           dT )
  {
    const int field = fieldIndexOfLoad( loadType );

    if ( isInstantaneous( dT, field ) )
      return; // no time passes: nothing flows in

    Eigen::Map< RhsSized >             P( P_ );
    Eigen::Map< KeSizedMatrix >        K( K_ );
    const Eigen::Map< const RhsSized > QTotal( QTotal_ );

    const bool   axisymmetric = sectionType == SectionType::Axisymmetric;
    const double factor       = sectionType == SectionType::Plane ? sectionFactorAt( XiSized::Zero() ) : 1.0;

    Discretizations::Elements::BoundaryElement boundaryEl( this->shape,
                                                           elementFace,
                                                           nDim,
                                                           this->coordinates,
                                                           axisymmetric );

    const Eigen::VectorXi& nodes      = boundaryEl.getBoundaryNodeIndicesInParent();
    const int              nFaceNodes = static_cast< int >( nodes.size() );

    Eigen::VectorXd Pb;
    Eigen::MatrixXd dPb_dQ = Eigen::MatrixXd::Zero( nFaceNodes, nFaceNodes );

    switch ( loadType ) {
    case Marmot::Discretizations::Elements::Element::HeatFlux:
    case Marmot::Discretizations::Elements::Element::MoistureFlux: {
      Pb = load[0] * boundaryEl.computeScalarLoadVector();
      break;
    }
    case Marmot::Discretizations::Elements::Element::HeatFilm:
    case Marmot::Discretizations::Elements::Element::MoistureFilm: {
      // q_n = h (T - T_inf), or q_n = beta (h - h_inf) for moisture
      const double h = load[0], TInf = load[1];
      // the external load density is the inflow, the negative of the outward flux
      const auto inflow = [h, TInf]( double T, double& f, double& df_dT ) {
        f     = -h * ( T - TInf );
        df_dT = -h;
      };
      Eigen::VectorXd Tb( nFaceNodes );
      for ( int k = 0; k < nFaceNodes; k++ )
        Tb( k ) = QTotal( nodes( k ) * nFields + field );
      boundaryEl.computeFieldDependentScalarLoadVector( Tb, inflow, Pb, dPb_dQ );
      break;
    }
    case Marmot::Discretizations::Elements::Element::HeatRadiation: {
      const double c = load[0], TInf = load[1], T0 = load[2];
      const double TInfAbs4 = std::pow( TInf - T0, 4 );
      const auto   inflow   = [c, T0, TInfAbs4]( double T, double& f, double& df_dT ) {
        const double TAbs = T - T0;
        f                 = -c * ( TAbs * TAbs * TAbs * TAbs - TInfAbs4 );
        df_dT             = -4.0 * c * TAbs * TAbs * TAbs;
      };
      Eigen::VectorXd Tb( nFaceNodes );
      for ( int k = 0; k < nFaceNodes; k++ )
        Tb( k ) = QTotal( nodes( k ) * nFields + field );
      boundaryEl.computeFieldDependentScalarLoadVector( Tb, inflow, Pb, dPb_dQ );
      break;
    }
    default: throw std::invalid_argument( "Invalid load type for a transport element" );
    }

    for ( int k = 0; k < nFaceNodes; k++ ) {
      const int row = nodes( k ) * nFields + field;
      P( row ) += factor * Pb( k );
      for ( int l = 0; l < nFaceNodes; l++ )
        K( row, nodes( l ) * nFields + field ) -= factor * dPb_dQ( k, l );
    }
  }

  template < int nDim, int nNodes, int nFields >
  void TransportFiniteElement< nDim, nNodes, nFields >::computeBodyForce( double*                        P_,
                                                                          [[maybe_unused]] double*       K,
                                                                          const double*                  load,
                                                                          [[maybe_unused]] const double* QTotal,
                                                                          [[maybe_unused]] double        time,
                                                                          double                         dT )
  {
    Eigen::Map< RhsSized >              P( P_ );
    Eigen::Matrix< double, nFields, 1 > s;
    for ( int a = 0; a < nFields; a++ )
      s( a ) = isInstantaneous( dT, a ) ? 0.0 : load[a]; // no time passes: the source produces nothing

    for ( const auto& qp : qps )
      P += qp.NPhi.transpose() * s * qp.J0xW;
  }

  template < int nDim, int nNodes, int nFields >
  std::vector< double > TransportFiniteElement< nDim, nNodes, nFields >::getCoordinatesAtCenter()
  {
    std::vector< double > coords( nDim );
    Eigen::Map< XiSized > coordsMap( coords.data() );
    coordsMap = this->center();
    return coords;
  }

  template < int nDim, int nNodes, int nFields >
  std::vector< std::vector< double > > TransportFiniteElement< nDim, nNodes, nFields >::
    getCoordinatesAtQuadraturePoints()
  {
    std::vector< std::vector< double > > listedCoords;
    std::vector< double >                coords( nDim );
    Eigen::Map< XiSized >                coordsMap( coords.data() );
    for ( const auto& qp : qps ) {
      coordsMap = this->NB( this->N( qp.xi ) ) * this->coordinates;
      listedCoords.push_back( coords );
    }
    return listedCoords;
  }

} // namespace Marmot::Discretizations::Elements
