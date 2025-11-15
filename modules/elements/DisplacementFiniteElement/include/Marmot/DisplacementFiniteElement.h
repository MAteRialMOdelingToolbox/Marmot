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
 * Matthias Neuner matthias.neuner@uibk.ac.at
 * Magdalena Schreter magdalena.schreter@uibk.ac.at
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
#include "Marmot/Marmot.h"
#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryElement.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotLowerDimensionalStress.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include <iostream>
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

        inline const static auto layout = makeLayout( {
          { .name = "stress", .length = 6 },
          { .name = "strain", .length = 6 },
          { .name = "begin of material state", .length = 0 },
        } );

      public:
        mVector6d                     stress;
        mVector6d                     strain;
        Eigen::Map< Eigen::VectorXd > materialStateVars;

        static int getNumberOfRequiredStateVarsQuadraturePointOnly() { return layout.nRequiredStateVars; };

        QPStateVarManager( double* theStateVarVector, int nStateVars )
          : MarmotStateVarVectorManager( theStateVarVector, layout ),
            stress( &find( "stress" ) ),
            strain( &find( "strain" ) ),
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
        // material->assignStateVars( managedStateVars->materialStateVars.data(),
        //                            managedStateVars->materialStateVars.size() );
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
                                 const double*                       time,
                                 double                              dT );

    /**
     * @brief Assemble body force contribution.
     * @details Integrates \f$\mathbf{P}_e^{(b)} = \int_{\Omega_e} \mathbf{N}^\mathsf{T} \mathbf{f}\,
     * \mathrm{d}\Omega\f$.
     */
    void computeBodyForce( double*       P,
                           double*       K,
                           const double* load,
                           const double* QTotal,
                           const double* time,
                           double        dT );

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
     * @param pNewdT Suggested scaling of dT by the material; if reduced (<1), the routine returns early.
     */
    void computeYourself( const double* QTotal,
                          const double* dQ,
                          double*       Pe,
                          double*       Ke,
                          const double* time,
                          double        dT,
                          double&       pNewdT );

    /**
     * @brief Compute consistent mass matrix using material density.
     * @details \f$\mathbf{M}_e = \sum_{qp} \rho\, \mathbf{N}^\mathsf{T}\mathbf{N}\, J_0 w\f$.
     */
    void computeConsistentInertia( double* M );

    /**
     * @brief Compute lumped mass vector via row-sum of consistent mass.
     * @details \f$\mathbf{m}_e = \mathrm{rowsum}(\mathbf{M}_e)\f$.
     */
    void computeLumpedInertia( double* M );

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

    static vector< vector< string > > nodeFields;
    if ( nodeFields.empty() )
      for ( int i = 0; i < nNodes; i++ ) {
        nodeFields.push_back( vector< string >() );
        nodeFields[i].push_back( "displacement" );
      }

    return nodeFields;
  }

  template < int nDim, int nNodes >
  std::vector< int > DisplacementFiniteElement< nDim, nNodes >::getDofIndicesPermutationPattern()
  {
    static std::vector< int > permutationPattern;
    if ( permutationPattern.empty() )
      for ( int i = 0; i < nNodes * nDim; i++ )
        permutationPattern.push_back( i );

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
      qp.material = std::unique_ptr< MarmotMaterialHypoElastic >( dynamic_cast< MarmotMaterialHypoElastic* >(
        MarmotLibrary::MarmotMaterialFactory::createMaterial( section.materialCode,
                                                              section.materialProperties,
                                                              section.nMaterialProperties,
                                                              elLabel ) ) );

      if ( !qp.material )
        throw std::invalid_argument( MakeString()
                                     << __PRETTY_FUNCTION__
                                     << ": invalid material assigned; cannot cast to MarmotMaterialHypoElastic!" );

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
  void DisplacementFiniteElement< nDim, nNodes >::computeYourself( const double* QTotal_,
                                                                   const double* dQ_,
                                                                   double*       Pe_,
                                                                   double*       Ke_,
                                                                   const double* time,
                                                                   double        dT,
                                                                   double&       pNewDT )
  {
    using namespace Marmot;
    using namespace ContinuumMechanics::VoigtNotation;

    Map< const RhsSized > QTotal( QTotal_ );
    Map< const RhsSized > dQ( dQ_ );
    Map< KeSizedMatrix >  Ke( Ke_ );
    Map< RhsSized >       Pe( Pe_ );

    Voigt  S, dE;
    CSized C;

    for ( QuadraturePoint& qp : qps ) {

      const BSized& B = qp.B;
      dE              = B * dQ;

      if constexpr ( nDim == 1 ) {

        MarmotMaterialHypoElastic::state1D  state;
        MarmotMaterialHypoElastic::timeInfo timeInfo;

        // set state info
        state.stress       = reduce3DVoigt< ParentGeometryElement::voigtSize >( qp.managedStateVars->stress )( 0 );
        state.strainEnergy = 0.0;
        state.stateVars    = qp.managedStateVars->materialStateVars.data();

        // set time info
        timeInfo.time = time[1];
        timeInfo.dT   = dT;
        try {
          qp.material->computeUniaxialStress( state, C.data(), dE.data(), timeInfo );
        }
        catch ( const std::runtime_error& e ) {
          pNewDT = 0.5;
          return;
        }
        Eigen::VectorXd stress1D( 1 );
        stress1D( 0 )               = state.stress;
        qp.managedStateVars->stress = make3DVoigt< ParentGeometryElement::voigtSize >( stress1D );
      }

      else if constexpr ( nDim == 2 ) {

        if ( sectionType == SectionType::PlaneStress ) {

          MarmotMaterialHypoElastic::state2D  state;
          MarmotMaterialHypoElastic::timeInfo timeInfo;

          // set state info
          state.stress       = reduce3DVoigt< ParentGeometryElement::voigtSize >( qp.managedStateVars->stress );
          state.strainEnergy = 0.0;
          state.stateVars    = qp.managedStateVars->materialStateVars.data();

          // set time info
          timeInfo.time = time[1];
          timeInfo.dT   = dT;
          try {
            qp.material->computePlaneStress( state, C.data(), dE.data(), timeInfo );
          }
          catch ( const std::runtime_error& e ) {
            pNewDT = 0.5;
            return;
          }
          qp.managedStateVars->stress = make3DVoigt< ParentGeometryElement::voigtSize >( state.stress );
          S                           = state.stress;
        }

        else if ( sectionType == SectionType::PlaneStrain ) {

          Vector6d dE6 = planeVoigtToVoigt( dE );
          Matrix6d C66;

          // Vector6d S6 = qp.managedStateVars->stress;
          MarmotMaterialHypoElastic::state3D  state;
          MarmotMaterialHypoElastic::timeInfo timeInfo;

          // set state info
          state.stress       = qp.managedStateVars->stress;
          state.strainEnergy = 0.0;
          state.stateVars    = qp.managedStateVars->materialStateVars.data();

          // set time info
          timeInfo.time = time[1];
          timeInfo.dT   = dT;
          try {
            qp.material->computeStress( state, C66.data(), dE6.data(), timeInfo );
          }
          catch ( const std::runtime_error& e ) {
            pNewDT = 0.5;
            return;
          }
          qp.managedStateVars->stress = state.stress;

          S = reduce3DVoigt< ParentGeometryElement::voigtSize >( state.stress );
          C = ContinuumMechanics::PlaneStrain::getPlaneStrainTangent( C66 );
        }
      }

      else if constexpr ( nDim == 3 ) {
        if ( sectionType == SectionType::Solid ) {

          MarmotMaterialHypoElastic::state3D  state;
          MarmotMaterialHypoElastic::timeInfo timeInfo;

          // set state info
          state.stress       = qp.managedStateVars->stress;
          state.strainEnergy = 0.0;
          state.stateVars    = qp.managedStateVars->materialStateVars.data();

          // set time info
          timeInfo.time = time[1];
          timeInfo.dT   = dT;
          try {
            qp.material->computeStress( state, C.data(), dE.data(), timeInfo );
          }
          catch ( const std::runtime_error& e ) {
            pNewDT = 0.5;
            return;
          }
          qp.managedStateVars->stress = state.stress;
          S                           = state.stress;
        }
      }

      qp.managedStateVars->strain += make3DVoigt< ParentGeometryElement::voigtSize >( dE );

      Ke += B.transpose() * C * B * qp.J0xW;
      Pe -= B.transpose() * S * qp.J0xW;
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
                                                                          const double* time,
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
                                                                    const double* time,
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
      const double rho = qp.material->getDensity();
      Me += N_.transpose() * N_ * qp.detJ * qp.weight * rho;
    }
  }
  template < int nDim, int nNodes >
  void DisplacementFiniteElement< nDim, nNodes >::computeLumpedInertia( double* M )
  {
    Map< RhsSized > Me( M );
    Me.setZero();

    KeSizedMatrix CMM;
    CMM.setZero();
    computeConsistentInertia( CMM.data() );

    Me = CMM.rowwise().sum();
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
