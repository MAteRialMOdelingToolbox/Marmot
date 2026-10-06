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
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryElement.h"
#include "Marmot/MarmotGeostaticStress.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrain.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrainFactory.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotTypedefs.h"
#include <memory>
#include <vector>

namespace Marmot::Elements {

  /**
   * @class Marmot::Elements::GradientEnhancedFiniteStrainDisplacementElement
   * @brief Gradient-enhanced (implicit-gradient, nonlocal) finite-strain displacement element.
   *
   * The displacement field @f$ \boldsymbol{u} @f$ (nDim dofs per node) is coupled to a scalar nonlocal field
   * @f$ \bar{N} @f$ (one dof per node), which regularizes softening materials by the implicit-gradient balance in the
   * reference configuration,
   * @f[
   *   \int_{\Omega_0} \boldsymbol{\tau} : \nabla_x \delta\boldsymbol{u} \, dV = \text{external work}, \qquad
   *   \int_{\Omega_0} \left( \bar{N}\,\delta\bar{N} + c\,\nabla_X\bar{N}\cdot\nabla_X\delta\bar{N}
   *   - L\,\delta\bar{N} \right) dV = 0,
   * @f]
   * with the Kirchhoff stress @f$ \boldsymbol{\tau} @f$, the local driving force @f$ L @f$ and @f$ c = l^2 @f$
   * from the material. The kinematics are those of DisplacementFiniteStrainULElement; the element is the
   * non-micropolar sibling of GradientEnhancedMicropolarULFiniteElement. It consumes a
   * MarmotMaterialGradientEnhancedFiniteStrain.
   *
   * The residual is the internal force vector (the same sign convention as DisplacementFiniteStrainULElement), the
   * dofs are ordered node by node, @f$ \{u_1, \dots, u_{nDim}, \bar{N}\} @f$, and internally field by field (see
   * getDofIndicesPermutationPattern()). Plane stress is not supported.
   *
   * @tparam nDim Number of spatial dimensions (2: plane strain, 3: solid).
   * @tparam nNodes Number of nodes.
   */
  template < int nDim, int nNodes >
  class GradientEnhancedFiniteStrainDisplacementElement : public MarmotElement,
                                                          public MarmotGeometryElement< nDim, nNodes > {

  public:
    /// @brief Section type of the element.
    enum SectionType {
      PlaneStress, ///< plane stress (not supported)
      PlaneStrain, ///< plane strain (nDim = 2)
      Solid,       ///< three-dimensional solid (nDim = 3)
    };

    static constexpr int nDofPerNodeU = nDim;          ///< dofs per node of the displacement field
    static constexpr int nDofPerNodeN = 1;             ///< dofs per node of the nonlocal field

    static constexpr int nCoordinates = nNodes * nDim; ///< number of nodal coordinates

    static constexpr int bsU = nNodes * nDofPerNodeU;  ///< size of the displacement block
    static constexpr int bsN = nNodes * nDofPerNodeN;  ///< size of the nonlocal block

    static constexpr int sizeLoadVector = bsU + bsN;   ///< number of dofs of the element

    static constexpr int idxU = 0;                     ///< first index of the displacement block
    static constexpr int idxN = idxU + bsU;            ///< first index of the nonlocal block

    using ParentGeometryElement = MarmotGeometryElement< nDim, nNodes >;
    using Material              = MarmotMaterialGradientEnhancedFiniteStrain;

    using JacobianSized = typename ParentGeometryElement::JacobianSized;
    using NSized        = typename ParentGeometryElement::NSized;
    using dNdXiSized    = typename ParentGeometryElement::dNdXiSized;
    using XiSized       = typename ParentGeometryElement::XiSized;
    using RhsSized      = Eigen::Matrix< double, sizeLoadVector, 1 >;
    using KSizedMatrix  = Eigen::Matrix< double, sizeLoadVector, sizeLoadVector >;
    using USizedVector  = Eigen::Matrix< double, bsU, 1 >;

    Eigen::Map< const Eigen::VectorXd > elementProperties;   ///< element properties: thickness (2D)
    const int                           elLabel;             ///< element label
    const SectionType                   sectionType;         ///< section type
    bool                                hasEigenDeformation; ///< whether a geostatic eigen deformation is applied

    /**
     * @struct QuadraturePoint
     * @brief A quadrature point: its location, the cached shape functions and derivatives, state and material.
     */
    struct QuadraturePoint {

      const XiSized xi;     ///< parametric coordinates
      const double  weight; ///< quadrature weight

      NSized     N;         ///< shape functions at @c xi
      dNdXiSized dNdX;      ///< shape function derivatives in the reference configuration
      double     detJ;      ///< determinant of the reference Jacobian
      double     J0xW;      ///< integration weight in the reference configuration (times the thickness in 2D)

      /**
       * @class QPStateVarManager
       * @brief State of a quadrature point: stress, energies, eigen deformation, followed by the material state.
       */
      class QPStateVarManager : public MarmotStateVarVectorManager {

        inline const static auto layout = makeLayout( {
          { .name = "stress", .length = 9 },
          { .name = "elastic energy density", .length = 1 },
          { .name = "dissipation density", .length = 1 },
          { .name = "F0 XX", .length = 1 },
          { .name = "F0 YY", .length = 1 },
          { .name = "F0 ZZ", .length = 1 },
          { .name = "begin of material state", .length = 0 },
        } );

      public:
        Eigen::Map< Marmot::Vector9d > stress; ///< Kirchhoff stress (3D, also in plane strain)
        double& elasticEnergyDensity;          ///< elastic energy per undeformed volume, as returned by the material
        double& dissipationDensity;            ///< dissipation per undeformed volume, accumulated by the material
        double& F0_XX;                         ///< eigen deformation, XX component
        double& F0_YY;                         ///< eigen deformation, YY component
        double& F0_ZZ;                         ///< eigen deformation, ZZ component
        Eigen::Map< Eigen::VectorXd > materialStateVars; ///< state variables of the material

        /// @brief Number of state variables of the quadrature point without those of the material.
        static int getNumberOfRequiredStateVarsQuadraturePointOnly() { return layout.nRequiredStateVars; };

        /**
         * @brief Map the state of a quadrature point.
         * @param[in] theStateVarVector State variables of the quadrature point.
         * @param[in] nStateVars Their number, including those of the material.
         */
        QPStateVarManager( double* theStateVarVector, int nStateVars )
          : MarmotStateVarVectorManager( theStateVarVector, layout ),
            stress( &find( "stress" ) ),
            elasticEnergyDensity( find( "elastic energy density" ) ),
            dissipationDensity( find( "dissipation density" ) ),
            F0_XX( find( "F0 XX" ) ),
            F0_YY( find( "F0 YY" ) ),
            F0_ZZ( find( "F0 ZZ" ) ),
            materialStateVars( &find( "begin of material state" ),
                               nStateVars - getNumberOfRequiredStateVarsQuadraturePointOnly() ){};
      };

      std::unique_ptr< QPStateVarManager > managedStateVars; ///< state of the quadrature point
      std::unique_ptr< Material >          material;         ///< material of the quadrature point

      /// @brief Number of state variables of the quadrature point without those of the material.
      int getNumberOfRequiredStateVarsQuadraturePointOnly()
      {
        return QPStateVarManager::getNumberOfRequiredStateVarsQuadraturePointOnly();
      };

      /// @brief Number of state variables of the quadrature point, including those of the material.
      int getNumberOfRequiredStateVars()
      {
        return getNumberOfRequiredStateVarsQuadraturePointOnly() + material->getNumberOfRequiredStateVars();
      };

      /**
       * @brief Assign the state variables of the quadrature point.
       * @param[in] stateVars State variables of the quadrature point.
       * @param[in] nStateVars Their number.
       */
      void assignStateVars( double* stateVars, int nStateVars )
      {
        managedStateVars = std::make_unique< QPStateVarManager >( stateVars, nStateVars );
      }

      /**
       * @brief Construct a quadrature point.
       * @param[in] xi Parametric coordinates.
       * @param[in] weight Quadrature weight.
       * @param[in] N Shape functions at @p xi.
       */
      QuadraturePoint( XiSized xi, double weight, const NSized& N )
        : xi( xi ), weight( weight ), N( N ), dNdX( dNdXiSized::Zero() ), detJ( 0.0 ), J0xW( 0.0 ){};
    };

    std::vector< QuadraturePoint > qps; ///< the quadrature points

    /**
     * @brief Construct the element.
     * @param[in] elementID Element label.
     * @param[in] integrationType Full or reduced integration.
     * @param[in] sectionType Section type: PlaneStrain for nDim = 2, Solid for nDim = 3.
     * @throws std::invalid_argument for a section type that does not match nDim, e.g. plane stress.
     */
    GradientEnhancedFiniteStrainDisplacementElement(
      int                                                 elementID,
      Marmot::FiniteElement::Quadrature::IntegrationTypes integrationType,
      SectionType                                         sectionType );

    /// @brief Number of state variables of the element (of all quadrature points).
    int getNumberOfRequiredStateVars();

    /// @brief Fields per node: "displacement" and "nonlocal damage".
    std::vector< std::vector< std::string > > getNodeFields();

    /// @brief Permutation from the node-by-node dof order of the host to the field-by-field order of the element.
    std::vector< int > getDofIndicesPermutationPattern();

    /// @brief Number of nodes.
    int getNNodes() { return nNodes; }

    /// @brief Number of spatial dimensions.
    int getNSpatialDimensions() { return nDim; }

    /// @brief Number of dofs of the element.
    int getNDofPerElement() { return sizeLoadVector; }

    /// @brief Shape of the element.
    std::string getElementShape() { return ParentGeometryElement::getElementShape(); }

    /**
     * @brief Assign the state variables, split evenly among the quadrature points.
     * @param[in] managedStateVars State variables of the element.
     * @param[in] nStateVars Their number.
     */
    void assignStateVars( double* managedStateVars, int nStateVars );

    /**
     * @brief Assign the element properties.
     * @param[in] MarmotElementProperty Element properties: the thickness in 2D.
     */
    void assignProperty( const ElementProperties& MarmotElementProperty );

    /**
     * @brief Create the material of every quadrature point.
     * @param[in] MarmotElementProperty Material name and properties.
     */
    void assignProperty( const MarmotMaterialSection& MarmotElementProperty );

    /**
     * @brief Assign the nodal coordinates in the reference configuration.
     * @param[in] coordinates Nodal coordinates, node by node.
     */
    void assignNodeCoordinates( const double* coordinates );

    /// @brief Compute the reference shape function derivatives and integration weights of the quadrature points.
    void initializeYourself();

    /**
     * @brief Set initial conditions.
     * @param[in] state MarmotMaterialInitialization (initializes the material state) or GeostaticStress (finds the
     * eigen deformation of a linearly distributed geostatic stress).
     * @param[in] values Definition of the geostatic stress distribution.
     * @throws std::invalid_argument for another initial condition.
     */
    void setInitialConditions( StateTypes state, const double* values );

    /**
     * @brief Compute a distributed load on an element face.
     * @param[in] loadType Pressure (follower load, with its stiffness @f$ K = -\partial P/\partial Q @f$) or
     * SurfaceTraction (in the reference configuration).
     * @param[in,out] P Load vector, added to.
     * @param[in,out] K Stiffness matrix, added to.
     * @param[in] elementFace Face of the element.
     * @param[in] load Pressure, or traction vector.
     * @param[in] QTotal Total dofs.
     * @param[in] time Time.
     * @param[in] dT Time increment.
     * @throws std::invalid_argument for another load type.
     */
    void computeDistributedLoad( MarmotElement::DistributedLoadTypes loadType,
                                 double*                             P,
                                 double*                             K,
                                 int                                 elementFace,
                                 const double*                       load,
                                 const double*                       QTotal,
                                 double                              time,
                                 double                              dT );

    /**
     * @brief Compute a body force per reference volume.
     * @param[in,out] P Load vector, added to.
     * @param[in,out] K Stiffness matrix (not changed, the load does not depend on the dofs).
     * @param[in] load Body force vector.
     * @param[in] QTotal Total dofs.
     * @param[in] time Time.
     * @param[in] dT Time increment.
     */
    void computeBodyForce( double* P, double* K, const double* load, const double* QTotal, double time, double dT );

    /**
     * @brief Compute the residual (internal forces) and the stiffness matrix, and update the state.
     * @param[in] QTotal Total dofs at the end of the increment.
     * @param[in] dQ Dof increment.
     * @param[in,out] Pe Residual, added to.
     * @param[in,out] Ke Stiffness matrix, added to.
     * @param[in] time Time.
     * @param[in] dT Time increment.
     */
    void computeKernels( const double* QTotal, const double* dQ, double* Pe, double* Ke, double time, double dT );

    /**
     * @brief Compute the residual (internal forces) only, and update the state, for explicit time integration.
     * @details The material is updated as in computeKernels(); its tangents are discarded.
     * @param[in] QTotal Total dofs at the end of the increment.
     * @param[in] dQ Dof increment.
     * @param[in,out] Pe Residual, added to.
     * @param[in] time Time.
     * @param[in] dT Time increment.
     */
    void computeKernelsExplicit( const double* QTotal, const double* dQ, double* Pe, double time, double dT );

    /**
     * @brief A view of a state of a quadrature point, of the element or of its material.
     * @param[in] stateName Name of the state.
     * @param[in] qpNumber Number of the quadrature point.
     * @return The view.
     */
    StateView getStateView( const std::string& stateName, int qpNumber );

    /// @brief Coordinates of the element center in the reference configuration.
    std::vector< double > getCoordinatesAtCenter();

    /// @brief Coordinates of the quadrature points in the reference configuration.
    std::vector< std::vector< double > > getCoordinatesAtQuadraturePoints();

    /// @brief Number of quadrature points.
    int getNumberOfQuadraturePoints();

  private:
    /**
     * @brief Update the material of a quadrature point, in plane strain through its 3D response.
     * @param[in,out] qp Quadrature point; its stress and energies are updated.
     * @param[in] F Deformation gradient at the end of the increment.
     * @param[in] nonlocalField Nonlocal field at the quadrature point.
     * @param[in] time Time.
     * @param[in] dT Time increment.
     * @param[out] response Kirchhoff stress, local driving force, nonlocal radius and energies.
     * @param[out] tangents Algorithmic tangents.
     */
    void computeMaterialResponse( QuadraturePoint&                                 qp,
                                  const Fastor::Tensor< double, nDim, nDim >&      F,
                                  double                                           nonlocalField,
                                  double                                           time,
                                  double                                           dT,
                                  typename Material::ConstitutiveResponse< nDim >& response,
                                  typename Material::AlgorithmicModuli< nDim >&    tangents );
  };

  template < int nDim, int nNodes >
  StateView GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::getStateView( const std::string& stateName,
                                                                                           int                qpNumber )
  {
    const auto& qp = qps[qpNumber];
    if ( qp.managedStateVars->contains( stateName ) ) {
      return qp.managedStateVars->getStateView( stateName );
    }
    else {
      return qp.material->getStateView( stateName, qp.managedStateVars->materialStateVars.data() );
    }
  }

  template < int nDim, int nNodes >
  GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::GradientEnhancedFiniteStrainDisplacementElement(
    int                                                 elementID,
    Marmot::FiniteElement::Quadrature::IntegrationTypes integrationType,
    SectionType                                         sectionType )
    : ParentGeometryElement(),
      elementProperties( Eigen::Map< const Eigen::VectorXd >( nullptr, 0 ) ),
      elLabel( elementID ),
      sectionType( sectionType ),
      hasEigenDeformation( false )
  {
    if ( ( nDim == 2 && sectionType != SectionType::PlaneStrain ) ||
         ( nDim == 3 && sectionType != SectionType::Solid ) )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__
                                                << ": only plane strain (2D) and solid (3D) sections are supported" );

    for ( const auto& qpInfo : Marmot::FiniteElement::Quadrature::getGaussPointInfo( this->shape, integrationType ) ) {
      QuadraturePoint qp( qpInfo.xi, qpInfo.weight, this->N( qpInfo.xi ) );
      qps.push_back( std::move( qp ) );
    }
  }

  template < int nDim, int nNodes >
  int GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::getNumberOfRequiredStateVars()
  {
    return qps[0].getNumberOfRequiredStateVars() * static_cast< int >( qps.size() );
  }

  template < int nDim, int nNodes >
  std::vector< std::vector< std::string > > GradientEnhancedFiniteStrainDisplacementElement< nDim,
                                                                                             nNodes >::getNodeFields()
  {
    using namespace std;

    static vector< vector< string > > nodeFields;
    if ( nodeFields.empty() )
      for ( int i = 0; i < nNodes; i++ ) {
        nodeFields.push_back( vector< string >() );
        nodeFields[i].push_back( "displacement" );
        nodeFields[i].push_back( "nonlocal damage" );
      }

    return nodeFields;
  }

  template < int nDim, int nNodes >
  std::vector< int > GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::getDofIndicesPermutationPattern()
  {
    static std::vector< int > permutationPattern;
    if ( permutationPattern.empty() ) {
      for ( int i = 0; i < nNodes; i++ )
        for ( int j = 0; j < nDim; j++ )
          permutationPattern.push_back( i * ( nDim + 1 ) + j );
      for ( int i = 0; i < nNodes; i++ )
        permutationPattern.push_back( i * ( nDim + 1 ) + nDim );
    }

    return permutationPattern;
  }

  template < int nDim, int nNodes >
  void GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::assignStateVars( double* managedStateVars,
                                                                                         int     nStateVars )
  {
    const int nQpStateVars = nStateVars / static_cast< int >( qps.size() );

    for ( size_t i = 0; i < qps.size(); i++ ) {
      auto&   qp          = qps[i];
      double* qpStateVars = &managedStateVars[i * nQpStateVars];
      qp.assignStateVars( qpStateVars, nQpStateVars );
    }
  }

  template < int nDim, int nNodes >
  void GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::assignProperty(
    const ElementProperties& elementPropertiesInfo )
  {
    new ( &elementProperties ) Eigen::Map< const Eigen::VectorXd >( elementPropertiesInfo.elementProperties,
                                                                    elementPropertiesInfo.nElementProperties );
  }

  template < int nDim, int nNodes >
  void GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::assignProperty(
    const MarmotMaterialSection& section )
  {
    for ( auto& qp : qps ) {
      qp.material = std::unique_ptr< Material >(
        MarmotLibrary::MarmotMaterialGradientEnhancedFiniteStrainFactory::createMaterial( section.materialName,
                                                                                          section.materialProperties,
                                                                                          section.nMaterialProperties,
                                                                                          elLabel ) );
    }
  }

  template < int nDim, int nNodes >
  void GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::assignNodeCoordinates(
    const double* coordinates )
  {
    ParentGeometryElement::assignNodeCoordinates( coordinates );
  }

  template < int nDim, int nNodes >
  void GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::initializeYourself()
  {
    for ( QuadraturePoint& qp : qps ) {

      const dNdXiSized dNdXi_ = this->dNdXi( qp.xi );

      const JacobianSized J    = this->Jacobian( dNdXi_ );
      const JacobianSized JInv = J.inverse();
      const double        detJ = J.determinant();
      qp.detJ                  = detJ;

      qp.dNdX = this->dNdX( dNdXi_, JInv );

      if constexpr ( nDim == 3 ) {
        qp.J0xW = qp.weight * detJ;
      }
      if constexpr ( nDim == 2 ) {
        const double& thickness = elementProperties[0];
        qp.J0xW                 = qp.weight * detJ * thickness;
      }
      if constexpr ( nDim == 1 ) {
        const double& crossSection = elementProperties[0];
        qp.J0xW                    = qp.weight * detJ * crossSection;
      }
    }
  }

  template < int nDim, int nNodes >
  void GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::computeMaterialResponse(
    QuadraturePoint&                                 qp,
    const Fastor::Tensor< double, nDim, nDim >&      F,
    double                                           nonlocalField,
    double                                           time,
    double                                           dT,
    typename Material::ConstitutiveResponse< nDim >& response,
    typename Material::AlgorithmicModuli< nDim >&    tangents )
  {
    using namespace Fastor;
    using namespace Marmot::FastorIndices;

    const typename Material::Deformation< nDim > deformation = { F, nonlocalField };

    const typename Material::TimeIncrement timeIncrement{ time, dT };

    if constexpr ( nDim == 2 ) {
      // plane strain (the constructor admits no other 2D section): the 3D response with F_33 = 1

      typename Material::ConstitutiveResponse< 3 >
        response3D( FastorStandardTensors::Tensor33d( qp.managedStateVars->stress.data(), Fastor::ColumnMajor ),
                    0.0,
                    0.0,
                    qp.managedStateVars->elasticEnergyDensity,
                    qp.managedStateVars->dissipationDensity,
                    qp.managedStateVars->materialStateVars.data() );

      typename Material::AlgorithmicModuli< 3 > algorithmicModuli3D;

      typename Material::Deformation< 3 > deformation3D{ expandTo3D( deformation.F ), deformation.N };
      deformation3D.F( 2, 2 ) = 1.0;

      if ( hasEigenDeformation )
        qp.material->computePlaneStrain( response3D,
                                         algorithmicModuli3D,
                                         deformation3D,
                                         timeIncrement,
                                         { qp.managedStateVars->F0_XX,
                                           qp.managedStateVars->F0_YY,
                                           qp.managedStateVars->F0_ZZ } );
      else
        qp.material->computePlaneStrain( response3D, algorithmicModuli3D, deformation3D, timeIncrement );

      response.tau                  = reduceTo2D< U, U >( response3D.tau );
      response.L                    = response3D.L;
      response.nonLocalRadius       = response3D.nonLocalRadius;
      response.elasticEnergyDensity = response3D.elasticEnergyDensity;
      response.dissipation          = response3D.dissipation;

      tangents.dTau_dF = reduceTo2D< U, U, U, U >( algorithmicModuli3D.dTau_dF );
      tangents.dTau_dN = reduceTo2D< U, U >( algorithmicModuli3D.dTau_dN );
      tangents.dL_dF   = reduceTo2D< U, U >( algorithmicModuli3D.dL_dF );
      tangents.dL_dN   = algorithmicModuli3D.dL_dN;

      qp.managedStateVars->stress = Marmot::mapEigenToFastor( response3D.tau ).reshaped();
    }
    else {
      response = typename Material::ConstitutiveResponse< nDim >( Tensor< double, nDim, nDim >( qp.managedStateVars
                                                                                                  ->stress.data(),
                                                                                                ColumnMajor ),
                                                                  0.0,
                                                                  0.0,
                                                                  qp.managedStateVars->elasticEnergyDensity,
                                                                  qp.managedStateVars->dissipationDensity,
                                                                  qp.managedStateVars->materialStateVars.data() );

      // as for plane strain above: a geostatic initial state lives in the eigen deformation
      if ( hasEigenDeformation )
        qp.material->computeStress( response,
                                    tangents,
                                    deformation,
                                    timeIncrement,
                                    { qp.managedStateVars->F0_XX,
                                      qp.managedStateVars->F0_YY,
                                      qp.managedStateVars->F0_ZZ } );
      else
        qp.material->computeStress( response, tangents, deformation, timeIncrement );
      qp.managedStateVars->stress = Marmot::mapEigenToFastor( response.tau ).reshaped();
    }

    // the materials accumulate the dissipation onto the incoming value: keep both with the state
    qp.managedStateVars->elasticEnergyDensity = response.elasticEnergyDensity;
    qp.managedStateVars->dissipationDensity   = response.dissipation;
  }

  template < int nDim, int nNodes >
  void GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::computeKernels( const double* qTotal,
                                                                                        const double* dQ,
                                                                                        double*       rightHandSide,
                                                                                        double*       stiffnessMatrix,
                                                                                        double        time,
                                                                                        double        dT )
  {
    using namespace Fastor;

    static const Tensor< double, nDim, nDim > I(
      ( Eigen::Matrix< double, nDim, nDim >() << Eigen::Matrix< double, nDim, nDim >::Identity() ).finished().data() );

    // in  ...
    const auto qU_np = TensorMap< const double, nNodes, nDim >( qTotal );
    const auto qN_np = TensorMap< const double, nNodes >( qTotal + idxN );

    // ... and out: residuals
    TensorMap< double, nNodes, nDim > r_U( rightHandSide );
    TensorMap< double, nNodes >       r_N( rightHandSide + idxN );

    // temporary stiffness sub-blocks
    Tensor< double, nDim, nNodes, nDim, nNodes > k_UU( 0.0 );
    Tensor< double, nDim, nNodes, nNodes >       k_UN( 0.0 );
    Tensor< double, nNodes, nDim, nNodes >       k_NU( 0.0 );
    Tensor< double, nNodes, nNodes >             k_NN( 0.0 );

    for ( auto& qp : qps ) {

      using namespace Marmot::FastorIndices;

      const auto N    = Tensor< double, nNodes >( qp.N.data() );
      const auto dNdX = Tensor< double, nDim, nNodes >( qp.dNdX.data(), ColumnMajor );

      const auto F_np = evaluate( einsum< Ai, jA >( qU_np, dNdX ) + I );

      const double nonlocalField = inner( N, qN_np );

      typename Material::ConstitutiveResponse< nDim > response;
      typename Material::AlgorithmicModuli< nDim >    tangents;
      computeMaterialResponse( qp, F_np, nonlocalField, time, dT, response, tangents );

      const auto dNdx = evaluate( einsum< ji, jA >( inv( F_np ), dNdX ) );

      const double& J0xW = qp.J0xW;

      const auto&  S          = response.tau;
      const double localField = response.L;
      const double c          = response.nonLocalRadius * response.nonLocalRadius;

      const auto& t = tangents;

      // aux stiffness tensors
      const auto dS_dqU = evaluate( +einsum< ijkl, lB >( t.dTau_dF, dNdX ) );
      const auto dS_dqN = evaluate( +einsum< ij, B >( t.dTau_dN, N ) );
      const auto dL_dqU = evaluate( +einsum< kl, lB >( t.dL_dF, dNdX ) );

      // residuals (Pe = +internal force, matching DisplacementFiniteStrainULElement)
      r_U += ( +einsum< iA, ij >( dNdx, S ) ) * J0xW;
      r_N += ( N * nonlocalField + c * einsum< iA, iB, B >( dNdX, dNdX, qN_np ) - N * localField ) * J0xW;

      // stiffness
      k_UU += ( +einsum< iA, ijkB, to_jAkB >( dNdx, dS_dqU ) ) * J0xW;
      k_UU += ( -einsum< kA, ij, iB, to_jAkB >( dNdx, S, dNdx ) ) * J0xW; // geometric
      k_UN += ( +einsum< iA, ijB, to_jAB >( dNdx, dS_dqN ) ) * J0xW;
      k_NU += ( -einsum< A, kB >( N, dL_dqU ) ) * J0xW;
      k_NN += ( +einsum< A, B >( N, N ) + einsum< iA, iB >( dNdX, dNdX ) * c ) * J0xW;
      k_NN += ( -einsum< A, B >( N, N ) * t.dL_dN ) * J0xW;
    }

    using namespace Eigen;

    Map< KSizedMatrix > K( stiffnessMatrix );

    K.template block< bsU, bsU >( idxU, idxU ) += Map< Matrix< double, bsU, bsU > >( torowmajor( k_UU ).data() );
    K.template block< bsU, bsN >( idxU, idxN ) += Map< Matrix< double, bsU, bsN > >( torowmajor( k_UN ).data() );
    K.template block< bsN, bsU >( idxN, idxU ) += Map< Matrix< double, bsN, bsU > >( torowmajor( k_NU ).data() );
    K.template block< bsN, bsN >( idxN, idxN ) += Map< Matrix< double, bsN, bsN > >( torowmajor( k_NN ).data() );
  }

  template < int nDim, int nNodes >
  void GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::computeKernelsExplicit( const double* qTotal,
                                                                                                const double* dQ,
                                                                                                double* rightHandSide,
                                                                                                double  time,
                                                                                                double  dT )
  {
    using namespace Fastor;

    static const Tensor< double, nDim, nDim > I(
      ( Eigen::Matrix< double, nDim, nDim >() << Eigen::Matrix< double, nDim, nDim >::Identity() ).finished().data() );

    const auto qU_np = TensorMap< const double, nNodes, nDim >( qTotal );
    const auto qN_np = TensorMap< const double, nNodes >( qTotal + idxN );

    TensorMap< double, nNodes, nDim > r_U( rightHandSide );
    TensorMap< double, nNodes >       r_N( rightHandSide + idxN );

    for ( auto& qp : qps ) {

      using namespace Marmot::FastorIndices;

      const auto N    = Tensor< double, nNodes >( qp.N.data() );
      const auto dNdX = Tensor< double, nDim, nNodes >( qp.dNdX.data(), ColumnMajor );

      const auto F_np = evaluate( einsum< Ai, jA >( qU_np, dNdX ) + I );

      const double nonlocalField = inner( N, qN_np );

      // the implicit material update; its tangents are not needed
      typename Material::ConstitutiveResponse< nDim > response;
      typename Material::AlgorithmicModuli< nDim >    tangents;
      computeMaterialResponse( qp, F_np, nonlocalField, time, dT, response, tangents );

      const auto dNdx = evaluate( einsum< ji, jA >( inv( F_np ), dNdX ) );

      const double c = response.nonLocalRadius * response.nonLocalRadius;

      r_U += ( +einsum< iA, ij >( dNdx, response.tau ) ) * qp.J0xW;
      r_N += ( N * nonlocalField + c * einsum< iA, iB, B >( dNdX, dNdX, qN_np ) - N * response.L ) * qp.J0xW;
    }
  }

  template < int nDim, int nNodes >
  void GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::computeDistributedLoad(
    MarmotElement::DistributedLoadTypes loadType,
    double*                             rightHandSide,
    double*                             stiffnessMatrix,
    int                                 elementFace,
    const double*                       load,
    const double*                       QTotal_,
    double                              time,
    double                              dT )
  {
    Eigen::Map< USizedVector > r_U( rightHandSide );

    switch ( loadType ) {

    case MarmotElement::Pressure: {
      const double                                                          p = load[0];
      const Eigen::Map< const RhsSized >                                    QTotal( QTotal_ );
      const Eigen::Ref< const USizedVector >                                qU( QTotal.head( bsU ) );
      Eigen::Map< Eigen::Matrix< double, sizeLoadVector, sizeLoadVector > > K( stiffnessMatrix );

      const USizedVector             coordinates_np = this->coordinates + qU;
      FiniteElement::BoundaryElement boundaryEl( this->shape, elementFace, nDim, coordinates_np );

      Eigen::VectorXd Pb = -p * boundaryEl.computeSurfaceNormalVectorialLoadVector();
      Eigen::MatrixXd Kb = -p * boundaryEl.computeDSurfaceNormalVectorialLoadVector_dCoordinates();

      if constexpr ( nDim == 2 ) {
        Pb *= elementProperties[0]; // thickness
        Kb *= elementProperties[0];
      }

      boundaryEl.assembleIntoParentVectorial( Pb, r_U );
      boundaryEl.assembleIntoParentStiffnessVectorial( Kb, K );

      break;
    }
    case MarmotElement::SurfaceTraction: {

      FiniteElement::BoundaryElement boundaryEl( this->shape, elementFace, nDim, this->coordinates );

      const XiSized tractionVector( load );

      auto Pk = boundaryEl.computeVectorialLoadVector( tractionVector );
      if constexpr ( nDim == 2 )
        Pk *= elementProperties[0]; // thickness
      boundaryEl.assembleIntoParentVectorial( Pk, r_U );

      break;
    }
    default: {
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid load type" );
    }
    }
  }

  template < int nDim, int nNodes >
  void GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::setInitialConditions(
    StateTypes    state,
    const double* initialConditionDefinition )
  {
    if constexpr ( nDim > 1 ) {
      switch ( state ) {

      case MarmotElement::MarmotMaterialInitialization: {
        for ( QuadraturePoint& qp : qps ) {

          qp.managedStateVars->F0_XX = 1.0;
          qp.managedStateVars->F0_YY = 1.0;
          qp.managedStateVars->F0_ZZ = 1.0;

          qp.material->initializeYourself( qp.managedStateVars->materialStateVars.data(),
                                           static_cast< int >( qp.managedStateVars->materialStateVars.size() ) );
        }
        break;
      }

      case MarmotElement::GeostaticStress: {

        for ( QuadraturePoint& qp : qps ) {

          XiSized coordAtGauss = this->NB( qp.N ) * this->coordinates;

          const auto geostaticNormalStressComponents = Marmot::GeostaticStress::
            getGeostaticStressFromLinearDistribution( initialConditionDefinition, coordAtGauss[1] );

          const auto [F0_XX,
                      F0_YY,
                      F0_ZZ] = qp.material
                                 ->findEigenDeformationForEigenStress( { qp.managedStateVars->F0_XX,
                                                                         qp.managedStateVars->F0_YY,
                                                                         qp.managedStateVars->F0_ZZ },
                                                                       geostaticNormalStressComponents,
                                                                       qp.managedStateVars->materialStateVars.data() );

          qp.managedStateVars->F0_XX = F0_XX;
          qp.managedStateVars->F0_YY = F0_YY;
          qp.managedStateVars->F0_ZZ = F0_ZZ;

          hasEigenDeformation = true;
        }
        break;
      }
      default: throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid initial condition" );
      }
    }
  }

  template < int nDim, int nNodes >
  void GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::computeBodyForce( double*       rightHandSide,
                                                                                          double*       stiffnessMatrix,
                                                                                          const double* load,
                                                                                          const double* qTotal,
                                                                                          double        time,
                                                                                          double        dT )
  {
    Eigen::Map< RhsSized >                                     r( rightHandSide );
    Eigen::Ref< USizedVector >                                 r_U( r.head( bsU ) );
    const Eigen::Map< const Eigen::Matrix< double, nDim, 1 > > f( load );

    for ( const auto& qp : qps )
      r_U += this->NB( qp.N ).transpose() * f * qp.J0xW;
  }

  template < int nDim, int nNodes >
  std::vector< double > GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::getCoordinatesAtCenter()
  {
    std::vector< double > coords( nDim );

    Eigen::Map< XiSized > coordsMap( &coords[0] );
    const auto            centerXi = XiSized::Zero();
    coordsMap                      = this->NB( this->N( centerXi ) ) * this->coordinates;
    return coords;
  }

  template < int nDim, int nNodes >
  std::vector< std::vector< double > > GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::
    getCoordinatesAtQuadraturePoints()
  {
    std::vector< std::vector< double > > listedCoords;

    std::vector< double > coords( nDim );
    Eigen::Map< XiSized > coordsMap( &coords[0] );

    for ( const auto& qp : qps ) {
      coordsMap = this->NB( qp.N ) * this->coordinates;
      listedCoords.push_back( coords );
    }

    return listedCoords;
  }

  template < int nDim, int nNodes >
  int GradientEnhancedFiniteStrainDisplacementElement< nDim, nNodes >::getNumberOfQuadraturePoints()
  {
    return static_cast< int >( qps.size() );
  }

} // namespace Marmot::Elements
