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
 * Alexander Dummer alexander.dummer@uibk.ac.at
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
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryElement.h"
#include "Marmot/MarmotGeostaticStress.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotTypedefs.h"
#include <iostream>
#include <memory>
#include <type_traits>
#include <vector>

namespace Marmot::Elements {
  /** @brief Implementation of a displacement-based finite element
   * for geometrically nonlinear analysis.
   *
   * @tparam nDim   Number of spatial dimensions (1, 2, 3)
   * @tparam nNodes Number of nodes of the element
   */
  template < int nDim, int nNodes >
  class DisplacementFiniteStrainULElement : public MarmotElement, public MarmotGeometryElement< nDim, nNodes > {

  public:
    /**
     * #SectionType is an enum class which involves the following cases:*/

    enum SectionType {
      PlaneStress, /**< Plane stress section for 2D elements. */
      PlaneStrain, /**< Plane strain section for 2D elements. */
      Solid,       /**< Solid (homogeneous) section for 1D, 2D and 3D elements. */
    };

    /// @brief Displacement degrees of freedom per node
    static constexpr int nDofPerNodeU = nDim; // Displacement   field U
    
    /// @brief Total number of coordinates of the element
    static constexpr int nCoordinates = nNodes * nDim;

    /// @brief Block size of element stiffness matrix and load vector for displacement field U
    static constexpr int bsU = nNodes * nDofPerNodeU;

    /// @brief Size of element stiffness matrix and load vector
    static constexpr int sizeLoadVector = bsU;

    /// @brief Starting index of displacement field U in element stiffness matrix and load vector
    static constexpr int idxU = 0;

    /// @brief Parent element class for geometry related operations as, e.g., shape functions
    using ParentGeometryElement = MarmotGeometryElement< nDim, nNodes >;
    /// @brief Material class for finite strain material formulations
    using Material              = MarmotMaterialFiniteStrain;
    /// @brief Sized matrix type used for Jacobian matrix (inherited from ParentGeometryElement)
    using JacobianSized = typename ParentGeometryElement::JacobianSized;
    /// @brief Sized matrix type used for shape function matrix (inherited from ParentGeometryElement)
    using NSized        = typename ParentGeometryElement::NSized;
    /// @brief Sized matrix typed used for shape function derivatives (inherited from ParentGeometryElement)
    using dNdXiSized    = typename ParentGeometryElement::dNdXiSized;
    /// @brief Sized vector type used for local coordinates (inherited from ParentGeometryElement)
    using XiSized       = typename ParentGeometryElement::XiSized;
    /// @brief Sized vector type used for the right hand side of the global equation system (negative element residual vector)
    using RhsSized      = Eigen::Matrix< double, sizeLoadVector, 1 >;
    /// @brief Sized matrix type used for the element stiffness matrix
    using KSizedMatrix  = Eigen::Matrix< double, sizeLoadVector, sizeLoadVector >;
    /// @brief Sized vector type used for the element displacement vector
    using USizedVector  = Eigen::Matrix< double, bsU, 1 >;
    /// @brief Sized vector type used for force vectors
    using ForceSized = Eigen::Matrix< double, nDim, 1 >;

    /// @brief Element properties as provided in the input file
    Eigen::Map< const Eigen::VectorXd > elementProperties;

    /// @brief Element label (ID)
    const int                           elLabel;
    /// @brief Section type of the element
    const SectionType                   sectionType;
    /// @brief Boolean for indicating whether initial deformation is considered or not
    bool hasEigenDeformation;

    /// @struct QuadraturePoint
    /// @brief Structure for storing quadrature point related information
    struct QuadraturePoint {

      const XiSized xi; /**< Local coordinates of the quadrature point */
      const double  weight; /**< Weight of the quadrature point */

      dNdXiSized dNdX; /**< Shape function derivatives w.r.t. material (undeformed) coordinates evaluated at the quadrature point */
      double     J0xW; /**< Determinant of the undeformed Jacobian times quadrature weight */

      /// @class QPStateVarManager
      /// @brief Manager class for handling state variables at the quadrature point
      class QPStateVarManager : public MarmotStateVarVectorManager {

        /// @brief Layout of the state variable vector at the quadrature point
        inline const static auto layout = makeLayout( {
          { .name = "stress", .length = 9 },
          { .name = "F0 XX", .length = 1 }, 
          { .name = "F0 YY", .length = 1 }, 
          { .name = "F0 ZZ", .length = 1 }, 
          { .name = "begin of material state", .length = 0 },
        } );

      public:
        Eigen::Map< Marmot::Vector9d > stress; /**< Stress tensor at the quadrature point */
        double&                        F0_XX; /**< Deformation gradient component XX for prescribing an initial deformation state*/
        double&                        F0_YY; /**< Deformation gradient component YY for prescribing an initial deformation state*/
        double&                        F0_ZZ; /**< Deformation gradient component ZZ for prescribing an initial deformation state*/
        Eigen::Map< Eigen::VectorXd >  materialStateVars; /**< Material state variables at the quadrature point */

        /// @brief Get number of required state variables at the quadrature point only (without material state variables)
        static int getNumberOfRequiredStateVarsQuadraturePointOnly() { return layout.nRequiredStateVars; };

        /** @brief Constructor of the state variable manager at the quadrature point
         * @param theStateVarVector[in] Pointer to the state variable vector at the quadrature point
         * @param nStateVars[in] Number of state variables at the quadrature point
         */
        QPStateVarManager( double* theStateVarVector, int nStateVars )
          : MarmotStateVarVectorManager( theStateVarVector, layout ),
            stress( &find( "stress" ) ), 
            F0_XX( find( "F0 XX" ) ),    
            F0_YY( find( "F0 YY" ) ),    
            F0_ZZ( find( "F0 ZZ" ) ),    
            materialStateVars( &find( "begin of material state" ),  
                               nStateVars - getNumberOfRequiredStateVarsQuadraturePointOnly() ){};
      };

      /// @brief Managed state variables at the quadrature point
      std::unique_ptr< QPStateVarManager > managedStateVars;

      /// @brief Material at the quadrature point
      std::unique_ptr< Material > material;

      /** @brief Get number of required state variables at the quadrature point only
       * @return Number of required state variables at the quadrature point only (without material state variables)
       */
      int getNumberOfRequiredStateVarsQuadraturePointOnly()
      {
        return QPStateVarManager::getNumberOfRequiredStateVarsQuadraturePointOnly();
      };

      /** @brief Get total number of required state variables at the quadrature point
       * @return Total number of required state variables at the quadrature point (including material state variables)
       */
      int getNumberOfRequiredStateVars()
      {
        return getNumberOfRequiredStateVarsQuadraturePointOnly() + material->getNumberOfRequiredStateVars();
      };

      /** @brief Assign the state variable vector at the quadrature point
       * @param stateVars[in] Pointer to the state variable vector at the quadrature point
       * @param nStateVars[in] Number of state variables at the quadrature point
       */
      void assignStateVars( double* stateVars, int nStateVars )
      {
        managedStateVars = std::make_unique< QPStateVarManager >( stateVars, nStateVars );
        material->assignStateVars( managedStateVars->materialStateVars.data(),
                                   managedStateVars->materialStateVars.size() );
      }

      /** @brief Constructor of the quadrature point
       * @param xi[in] Local coordinates of the quadrature point
       * @param weight[in] Weight of the quadrature point
       * @note The shape function derivatives w.r.t. material (undeformed) coordinates and the determinant of the undeformed Jacobian times quadrature weight are initialized with zero values.
       */
      QuadraturePoint( XiSized xi, double weight )
        : xi( xi ), weight( weight ), dNdX( dNdXiSized::Zero() ), J0xW( 0.0 ){};
    };

    /// @brief List of quadrature points of the element
    std::vector< QuadraturePoint > qps;

    /** @brief Constructor of the displacement-based finite strain element
     * @param elementID[in] Element ID (label) of the element
     * @param integrationType[in] Integration type of the element
     * @param sectionType[in] Section type of the element
     */
    DisplacementFiniteStrainULElement( int                                                 elementID,
                                       Marmot::FiniteElement::Quadrature::IntegrationTypes integrationType,
                                       SectionType                                         sectionType );

    /// @brief Get the total number of required state variables of the element
    int getNumberOfRequiredStateVars();

    /// @brief Get the nodal fields of the element
    std::vector< std::vector< std::string > > getNodeFields();

    /// @brief Get the permutation pattern of the element degrees of freedom
    std::vector< int > getDofIndicesPermutationPattern();

    /** @brief Get the number of nodes of the element
     * @return Number of nodes of the element
     */
    int getNNodes() { return nNodes; }

    /** @brief Get the number of spatial dimensions of the element
     * @return Number of spatial dimensions of the element
     */
    int getNSpatialDimensions() { return nDim; }

    /** @brief Get the number of degrees of freedom per node of the element
     * @return Number of degrees of freedom per node of the element
     */
    int getNDofPerElement() { return sizeLoadVector; }

    /** @brief Get the shape of the element
     * @return Shape of the element
     */
    std::string getElementShape() { return ParentGeometryElement::getElementShape(); }

    /** @brief Assign the state variable vector of the element
     * @param managedStateVars[in] Pointer to the state variable vector of the element
     * @param nStateVars[in] Number of state variables of the element
     */
    void assignStateVars( double* managedStateVars, int nStateVars );

    /** @brief Assign the element properties of the element
     * @param MarmotElementProperty[in] Element properties
     */
    void assignProperty( const ElementProperties& MarmotElementProperty );

    /** @brief Assign the material section of the element
     * @param MarmotElementProperty[in] Material section
     */
    void assignProperty( const MarmotMaterialSection& MarmotElementProperty );

    /** @brief Assign the nodal coordinates of the element
     * @param coordinates[in] Pointer to the nodal coordinates of the element
     */
    void assignNodeCoordinates( const double* coordinates );

    /// @brief Initialize the element
    void initializeYourself();

    /** @brief Set the initial conditions of the element
     * @param state[in] Type of the initial state
     * @param values[in] Pointer to the values defining the initial state
     */
    void setInitialConditions( StateTypes state, const double* values );

    /** @brief Compute the contributions of distributed loads to the element residual vector and stiffness matrix
     *
     *  For a given distributed load vector \f$\bar{\boldsymbol{t}}^{(n+1)}\f$ at the current time step \f$t^{(n+1)} = t^{(n)} + \Delta\,t\f$, compute the distributed load contribution to the negative element residual vector (right hand side of global newton) \f$\int_\bar{A}\,\mathbf{N}_A\,\bar{t}_j\,d\bar{A}\f$ and the element stiffness matrix \f$-\int_\bar{A}\,\mathbf{N}_{A}\,\bar{t}_i\left(\delta_{ij}\delta_{lk} - \delta_{ik}\delta_{lj}\right)\,\mathbf{N}_{B,l}\,d\bar{A}\f$.
     *
     * @param loadType[in] Type of the distributed load, e.g., pressure or surface traction
     * @param P[in,out] Pointer to the element residual vector (right hand side of the global equation system)
     * @param K[in,out] Pointer to the element stiffness matrix
     * @param elementFace[in] Local face number of the element where the distributed load is applied
     * @param load[in] Pointer to the distributed load vector
     * @param QTotal[in] Pointer to the total element displacement vector at the current time step
     * @param time[in] Pointer to the time at the beginning of the current time step
     * @param dT[in] Length of the current time step
     */
    void computeDistributedLoad( MarmotElement::DistributedLoadTypes loadType,
                                 double*                             P,
                                 double*                             K,
                                 const int                           elementFace,
                                 const double*                       load,
                                 const double*                       QTotal,
                                 const double*                       time,
                                 double                              dT );

    /** @brief Compute the contributions of body forces to the element residual vector and stiffness matrix
     *
     *  For a given body force vector \f$\mathbf{f}^{(n+1)}\f$ at the current time step \f$t^{(n+1)} = t^{(n)} + \Delta\,t\f$, compute the body force contribution for the negative element residual vector (right hand side of global newton) \f$\int_{V_0}\,\mathbf{N}_A\,f_j\,dV_0\f$. The stiffness matrix contribution is zero and thus not computed.
     *
     * @param P[in,out] Pointer to the element residual vector (right hand side of the global equation system)
     * @param K[in,out] Pointer to the element stiffness matrix
     * @param load[in] Pointer to the body force vector
     * @param QTotal[in] Pointer to the total element displacement vector at the current time step
     * @param time[in] Pointer to the time at the beginning of the current time step
     * @param dT[in] Length of the current time step
     */
    void computeBodyForce( double*       P,
                           double*       K,
                           const double* load,
                           const double* QTotal,
                           const double* time,
                           double        dT );

    /** @brief Compute the negative element residual vector (right hand side of global newton) and stiffness matrix
     *
     * For a given displacement \f$\mathbf{q}^{(n+1)}\f$ at the current time step \f$t^{(n+1)} = t^{(n)} + \Delta\,t\f$, compute the internal work contribution for the negative element residual vector (right hand side of global newton) \f$-\int_{V_0}\,\mathbf{N}_{A,i}\,\tau_{ij}\,dV_0\f$ and the element stiffness matrix \f$\int_{V_0}\,\mathbf{N}_{A,i}\,\frac{\partial \tau_{ij}}{\partial F_{kK}}\,\mathbf{N}_{B,K}\,-\,\mathbf{N}_{A,k\,}\mathbf{N}_{B,i}\,\tau_{ij}\,dV_0\f$.
     *
     * @param QTotal[in] Pointer to the total element displacement vector at the current time step
     * @param dQ[in] Pointer to the increment of the element displacement vector at the current time step
     * @param Pe[in,out] Pointer to the negative element residual vector (right hand side of global newton)
     * @param Ke[in,out] Pointer to the element stiffness matrix
     * @param time[in] Pointer to the time at the beginning of the current time step
     * @param dT[in] Length of the current time step
     * @param pNewdT[in,out] Suggested length of the next time step
     */
    void computeYourself( const double* QTotal,
                          const double* dQ,
                          double*       Pe,
                          double*       Ke,
                          const double* time,
                          double        dT,
                          double&       pNewdT );

    /** @brief Get a view to a state variable at a specific quadrature point of the element
     * @param stateName[in] Name of the state variable
     * @param qpNumber[in] Number of the quadrature point where the state variable is stored
     * @return View to the requested state variable at the specified quadrature point
     */
    StateView getStateView( const std::string& stateName, int qpNumber );

    /** @brief Get the coordinates of the element center
     * @return Coordinates of the element center
     */
    std::vector< double > getCoordinatesAtCenter();

    /** @brief Get the coordinates of all quadrature points of the element
     * @return Coordinates of all quadrature points of the element
     */
    std::vector< std::vector< double > > getCoordinatesAtQuadraturePoints();

    /** @brief Get the number of quadrature points of the element
     * @return Number of quadrature points of the element
     */
    int getNumberOfQuadraturePoints();
  };

  template < int nDim, int nNodes >
  StateView DisplacementFiniteStrainULElement< nDim, nNodes >::getStateView( const std::string& stateName,
                                                                             int                qpNumber )
  {
    const auto& qp = qps[qpNumber];
    if ( qp.managedStateVars->contains( stateName ) ) {
      return qp.managedStateVars->getStateView( stateName );
    }
    else {
      return qp.material->getStateView( stateName );
    }
  }

  template < int nDim, int nNodes >
  DisplacementFiniteStrainULElement< nDim, nNodes >::DisplacementFiniteStrainULElement(
    int                                                 elementID,
    Marmot::FiniteElement::Quadrature::IntegrationTypes integrationType,
    SectionType                                         sectionType )
    : ParentGeometryElement(),
      elementProperties( Eigen::Map< const Eigen::VectorXd >( nullptr, 0 ) ),
      elLabel( elementID ),
      sectionType( sectionType ),
      hasEigenDeformation( false )
  {
    for ( const auto& qpInfo : Marmot::FiniteElement::Quadrature::getGaussPointInfo( this->shape, integrationType ) ) {
      QuadraturePoint qp( qpInfo.xi, qpInfo.weight );
      qps.push_back( std::move( qp ) );
    }
  }

  template < int nDim, int nNodes >
  int DisplacementFiniteStrainULElement< nDim, nNodes >::getNumberOfRequiredStateVars()
  {
    return qps[0].getNumberOfRequiredStateVars() * qps.size();
  }

  template < int nDim, int nNodes >
  std::vector< std::vector< std::string > > DisplacementFiniteStrainULElement< nDim, nNodes >::getNodeFields()
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
  std::vector< int > DisplacementFiniteStrainULElement< nDim, nNodes >::getDofIndicesPermutationPattern()
  {
    static std::vector< int > permutationPattern;
    if ( permutationPattern.empty() ) {
      for ( int i = 0; i < nNodes; i++ )
        for ( int j = 0; j < nDim; j++ )
          permutationPattern.push_back( i * nDim + j );
    }

    return permutationPattern;
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteStrainULElement< nDim, nNodes >::assignStateVars( double* managedStateVars, int nStateVars )
  {
    const int nQpStateVars = nStateVars / qps.size();

    for ( size_t i = 0; i < qps.size(); i++ ) {
      auto&   qp          = qps[i];
      double* qpStateVars = &managedStateVars[i * nQpStateVars];
      qp.assignStateVars( qpStateVars, nQpStateVars );
    }
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteStrainULElement< nDim, nNodes >::assignProperty(
    const ElementProperties& elementPropertiesInfo )
  {
    new ( &elementProperties ) Eigen::Map< const Eigen::VectorXd >( elementPropertiesInfo.elementProperties,
                                                                    elementPropertiesInfo.nElementProperties );
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteStrainULElement< nDim, nNodes >::assignProperty( const MarmotMaterialSection& section )
  {

    for ( auto& qp : qps ) {
      qp.material = std::unique_ptr< Material >(
        dynamic_cast< Material* >( MarmotLibrary::MarmotMaterialFactory::createMaterial( section.materialCode,
                                                                                         section.materialProperties,
                                                                                         section.nMaterialProperties,
                                                                                         elLabel ) ) );
    }
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteStrainULElement< nDim, nNodes >::assignNodeCoordinates( const double* coordinates )
  {
    ParentGeometryElement::assignNodeCoordinates( coordinates );
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteStrainULElement< nDim, nNodes >::initializeYourself()
  {
    for ( QuadraturePoint& qp : qps ) {

      const dNdXiSized dNdXi_ = this->dNdXi( qp.xi );

      const JacobianSized J    = this->Jacobian( dNdXi_ );
      const JacobianSized JInv = J.inverse();
      const double        detJ = J.determinant();

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
  void DisplacementFiniteStrainULElement< nDim, nNodes >::computeYourself( const double* qTotal,
                                                                           const double* dQ,
                                                                           double*       rightHandSide,
                                                                           double*       stiffnessMatrix,
                                                                           const double* time,
                                                                           double        dT,
                                                                           double&       pNewDT )
  {
    using namespace Fastor;

    const static Tensor< double, nDim, nDim > I(
      ( Eigen::Matrix< double, nDim, nDim >() << Eigen::Matrix< double, nDim, nDim >::Identity() ).finished().data() );

    // in  ...
    const auto qU_np = TensorMap< const double, nNodes, nDim >( qTotal );

    // ... and out: residuals and stiffness
    TensorMap< double, nNodes, nDim > r_U( rightHandSide );

    // temporary stiffness matrices, which are assembled into the large one at the end of the method
    Tensor< double, nDim, nNodes, nDim, nNodes > k_UU( 0.0 );

    Eigen::Map< Eigen::VectorXd > rhs( rightHandSide, sizeLoadVector );

    for ( auto& qp : qps ) {

      using namespace Marmot::FastorIndices;

      const auto& dNdX_ = qp.dNdX;

      const auto dNdX = Tensor< double, nDim, nNodes >( dNdX_.data(), ColumnMajor );

      const auto F_np = evaluate( einsum< Ai, jA >( qU_np, dNdX ) + I );

      const Material::Deformation< nDim > deformation = { F_np };

      const Material::TimeIncrement timeIncrement{ time[0], dT };

      Material::ConstitutiveResponse< nDim > response;
      Material::AlgorithmicModuli< nDim >    tangents;
      try {
        if constexpr ( nDim == 2 ) {

          if ( sectionType == SectionType::PlaneStrain ) {

            using namespace Marmot;

            Material::ConstitutiveResponse< 3 >
              response3D{ FastorStandardTensors::Tensor33d( qp.managedStateVars->stress.data(), Fastor::ColumnMajor ),
                          -1.0,
                          -1.0 };

            Material::AlgorithmicModuli< 3 > algorithmicModuli3D;

            Material::Deformation< 3 > deformation3D{
              expandTo3D( deformation.F ),
            };

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

            response = { reduceTo2D< U, U >( response3D.tau ), response3D.rho, response3D.elasticEnergyDensity };

            tangents = {
              reduceTo2D< U, U, U, U >( algorithmicModuli3D.dTau_dF ),
            };

            qp.managedStateVars->stress = Marmot::mapEigenToFastor( response3D.tau ).reshaped();
          }
        }
        else {
          response = { Marmot::FastorStandardTensors::Tensor33d( qp.managedStateVars->stress.data(), ColumnMajor ),
                       -1.0,
                       -1.0 };

          qp.material->computeStress( response, tangents, deformation, timeIncrement );

          // implicit conversion to col major
          qp.managedStateVars->stress = Marmot::mapEigenToFastor( response.tau ).reshaped();
        }
      }
      catch ( const std::runtime_error& ) {
        pNewDT = 0.25;
        return;
      }
      const auto dNdx = evaluate( einsum< ji, jA >( inv( F_np ), dNdX ) );

      const double& J0xW = qp.J0xW;

      const auto& tau = response.tau;

      const auto& t = tangents;

      // aux stiffness tensors
      const auto dTau_dqU = evaluate( +einsum< ijkl, lB >( t.dTau_dF, dNdX ) );

      // r[ node, dim ] (swap to abuse directly colmajor layout)
      // directly operate via TensorMap
      r_U -= ( +einsum< iA, ij >( dNdx, tau ) ) * J0xW;

      // K [dim, node, dim, node ]
      k_UU += ( +einsum< iA, ijkB, to_jAkB >( dNdx, dTau_dqU ) ) * J0xW;

      // geometric contribution
      k_UU += ( -einsum< kA, ij, iB, to_jAkB >( dNdx, tau, dNdx ) ) * J0xW;
    }
    // copy back to the subblocks using mighty Eigen block access,
    // note the layout swap rowmajor -> colmajor
    // Fastor offers similar functionality, but performancy is (slightly) inferior

    using namespace Eigen;

    Map< KSizedMatrix > K( stiffnessMatrix );

    K.template block< bsU, bsU >( idxU, idxU ) += Map< Matrix< double, bsU, bsU > >( torowmajor( k_UU ).data() );
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteStrainULElement< nDim, nNodes >::computeDistributedLoad(
    MarmotElement::DistributedLoadTypes loadType,
    double*                             rightHandSide,
    double*                             stiffnessMatrix,
    const int                           elementFace,
    const double*                       load,
    const double*                       QTotal_,
    const double*                       time,
    double                              dT )
  {

    Eigen::Map< USizedVector > r_U( rightHandSide );

    switch ( loadType ) {

    case MarmotElement::Pressure: {
      const double                                                          p = load[0];
      const Eigen::Map< const RhsSized >                                    QTotal( QTotal_ );
      const Eigen::Ref< const USizedVector >                                qU( QTotal.head( bsU ) );
      Eigen::Map< Eigen::Matrix< double, sizeLoadVector, sizeLoadVector > > K( stiffnessMatrix );
      Eigen::Ref< Eigen::Matrix< double, bsU, bsU > >                       kUU( K.topLeftCorner( bsU, bsU ) );

      const USizedVector             coordinates_np = this->coordinates + qU;
      FiniteElement::BoundaryElement boundaryEl( this->shape, elementFace, nDim, coordinates_np );

      Eigen::VectorXd Pb = -p * boundaryEl.computeSurfaceNormalVectorialLoadVector();
      Eigen::MatrixXd Kb = -p * boundaryEl.computeDSurfaceNormalVectorialLoadVector_dCoordinates();

      if ( nDim == 2 ) {
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
      if ( nDim == 2 )
        Pk *= elementProperties[0]; // thickness
      boundaryEl.assembleIntoParentVectorial( Pk, r_U );

      break;
    }
    default: {
      throw std::invalid_argument( "Invalid Load Type specified" );
    }
    }
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteStrainULElement< nDim, nNodes >::setInitialConditions(
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

          qp.material->initializeYourself();
        }
        break;
      }

      case MarmotElement::GeostaticStress: {

        for ( QuadraturePoint& qp : qps ) {

          XiSized coordAtGauss = this->NB( this->N( qp.xi ) ) * this->coordinates;

          const auto geostaticNormalStressComponents = Marmot::GeostaticStress::
            getGeostaticStressFromLinearDistribution( initialConditionDefinition, coordAtGauss[1] );

          const auto [F0_XX,
                      F0_YY,
                      F0_ZZ] = qp.material->findEigenDeformationForEigenStress( { qp.managedStateVars->F0_XX,
                                                                                  qp.managedStateVars->F0_YY,
                                                                                  qp.managedStateVars->F0_ZZ },
                                                                                geostaticNormalStressComponents );

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
  void DisplacementFiniteStrainULElement< nDim, nNodes >::computeBodyForce( double*       rightHandSide,
                                                                            double*       stiffnessMatrix,
                                                                            const double* load,

                                                                            const double* qTotal,
                                                                            const double* time,
                                                                            double        dT )
  {
    Eigen::Map< RhsSized >                                     r( rightHandSide );
    Eigen::Ref< USizedVector >                                 r_U( r.head( bsU ) );
    const Eigen::Map< const Eigen::Matrix< double, nDim, 1 > > f( load );

    for ( const auto& qp : qps )
      r_U += this->NB( this->N( qp.xi ) ).transpose() * f * qp.J0xW;
  }

  template < int nDim, int nNodes >
  std::vector< double > DisplacementFiniteStrainULElement< nDim, nNodes >::getCoordinatesAtCenter()
  {
    std::vector< double > coords( nDim );

    Eigen::Map< XiSized > coordsMap( &coords[0] );
    const auto            centerXi = XiSized::Zero();
    coordsMap                      = this->NB( this->N( centerXi ) ) * this->coordinates;
    return coords;
  }

  template < int nDim, int nNodes >
  std::vector< std::vector< double > > DisplacementFiniteStrainULElement< nDim,
                                                                          nNodes >::getCoordinatesAtQuadraturePoints()
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
  int DisplacementFiniteStrainULElement< nDim, nNodes >::getNumberOfQuadraturePoints()
  {
    return qps.size();
  }

  template < int nNodes >
  class AxiSymmetricDisplacementFiniteStrainULElement : public DisplacementFiniteStrainULElement< 2, nNodes > {

    using DisplacementFiniteStrainULElement< 2, nNodes >::DisplacementFiniteStrainULElement;

    void computeYourself( const double* QTotal,
                          const double* dQ,
                          double*       Pe,
                          double*       Ke,
                          const double* time,
                          double        dT,
                          double&       pNewdT );
  };

  template < int nNodes >
  void AxiSymmetricDisplacementFiniteStrainULElement< nNodes >::computeYourself( const double* qTotal,
                                                                                 const double* dQ,
                                                                                 double*       rightHandSide,
                                                                                 double*       stiffnessMatrix,
                                                                                 const double* time,
                                                                                 double        dT,
                                                                                 double&       pNewDT )
  {
    constexpr int nDim = 2;

    using Parent              = DisplacementFiniteStrainULElement< 2, nNodes >;
    const auto idxU           = Parent::idxU;
    const auto sizeLoadVector = DisplacementFiniteStrainULElement< 2, nNodes >::sizeLoadVector;
    using Material            = MarmotMaterialFiniteStrain;

    using namespace Fastor;

    const static Tensor< double, nDim, nDim > I(
      ( Eigen::Matrix< double, nDim, nDim >() << Eigen::Matrix< double, nDim, nDim >::Identity() ).finished().data() );

    // in  ...
    const auto qU_np = TensorMap< const double, nNodes, nDim >( qTotal );

    // ... and out: residuals and stiffness
    TensorMap< double, nNodes, nDim > r_U( rightHandSide );

    // temporary stiffness matrices, which are assembled into the large one at the end of the method
    Tensor< double, nDim, nNodes, nDim, nNodes > k_UU( 0.0 );

    Eigen::Map< Eigen::VectorXd > rhs( rightHandSide, sizeLoadVector );

    for ( auto& qp : Parent::qps ) {

      using namespace Marmot::FastorIndices;

      auto        N_    = this->N( qp.xi );
      const auto& dNdX_ = qp.dNdX;

      Eigen::Vector2d coords = this->NB( N_ ) * this->coordinates;
      const double    r      = coords[0];

      const auto N    = Tensor< double, nNodes >( N_.data() );
      const auto dNdX = Tensor< double, nDim, nNodes >( dNdX_.data(), ColumnMajor );

      const auto u_np = evaluate( einsum< A, Ai >( N, qU_np ) );

      const auto F_np = evaluate( einsum< Ai, jA >( qU_np, dNdX ) + I );

      const Material::Deformation< nDim > deformation = {
        F_np,
      };

      const Material ::TimeIncrement timeIncrement{ time[0], dT };

      Material::ConstitutiveResponse< nDim > response;
      Material::AlgorithmicModuli< nDim >    tangents;

      using namespace Marmot;
      Material::ConstitutiveResponse< 3 >
        response3D{ FastorStandardTensors::Tensor33d( qp.managedStateVars->stress.data(), Fastor::ColumnMajor ),
                    -1.0,
                    -1.0 };

      Material::AlgorithmicModuli< 3 > algorithmicModuli3D;

      Material::Deformation< 3 > deformation3D{ expandTo3D( deformation.F ) };

      deformation3D.F( 2, 2 ) = 1 + u_np[0] / r;

      try {
        qp.material->computePlaneStrain( response3D, algorithmicModuli3D, deformation3D, timeIncrement );
      }
      catch ( const std::runtime_error& ) {
        pNewDT = 0.25;
        return;
      }
      response = { reduceTo2D< U, U >( response3D.tau ), response3D.rho, response3D.elasticEnergyDensity };

      tangents = {
        reduceTo2D< U, U, U, U >( algorithmicModuli3D.dTau_dF ),
      };

      qp.managedStateVars->stress = Marmot::mapEigenToFastor( response3D.tau ).reshaped();

      const auto dNdx = evaluate( einsum< ji, jA >( inv( F_np ), dNdX ) );

      const double J0xWxRx2Pi = qp.J0xW * 2 * Constants::Pi * r;

      const auto& tau = response.tau;

      const auto& t = tangents;

      // aux stiffness tensors
      auto dTau_dqU = evaluate( +einsum< ijkl, lB >( t.dTau_dF, dNdX ) );

      Fastor::Tensor< double, 2, 2 > dTau33_dF_2D( 0.0 );
      for ( int i = 0; i < 2; i++ )
        for ( int j = 0; j < 2; j++ ) {
          dTau33_dF_2D( i, j ) = algorithmicModuli3D.dTau_dF( 2, 2, i, j );
        }
      auto dTau33_dqU = evaluate( +einsum< kl, lB >( dTau33_dF_2D, dNdX ) );

      // additional terms due to axisymmetry:

      for ( int B = 0; B < nNodes; B++ ) {
        for ( int i = 0; i < 2; i++ )
          for ( int j = 0; j < 2; j++ ) {
            dTau_dqU( i, j, 0, B ) += algorithmicModuli3D.dTau_dF( i, j, 2, 2 ) * N( B ) / r;
          }
        dTau33_dqU( 0, B ) += algorithmicModuli3D.dTau_dF( 2, 2, 2, 2 ) * N( B ) / r;
      }

      // r[ node, dim ] (swap to abuse directly colmajor layout)
      // directly operate via TensorMap
      r_U -= ( +einsum< iA, ij >( dNdx, tau ) ) * J0xWxRx2Pi;

      const double F33          = 1 + u_np[0] / r;
      const double invF33       = 1. / F33;
      const double dInvF33_dF33 = -1. / ( F33 * F33 );

      for ( int A = 0; A < nNodes; A++ ) {
        r_U( A, 0 ) -= ( +N( A ) * invF33 * response3D.tau( 2, 2 ) / r ) * J0xWxRx2Pi;

        for ( int B = 0; B < nNodes; B++ ) {
          for ( int j = 0; j < 2; j++ ) {
            k_UU( 0, A, j, B ) += ( +N( A ) * invF33 * dTau33_dqU( j, B ) / r ) * J0xWxRx2Pi;
          }
          const double dF33_dN_qU_0 = ( N( B ) / r );
          k_UU( 0, A, 0, B ) += ( +N( A ) * dInvF33_dF33 * dF33_dN_qU_0 * response3D.tau( 2, 2 ) / r ) * J0xWxRx2Pi;
        }
      }

      // K [dim, node, dim, node ]
      k_UU += ( +einsum< iA, ijkB, to_jAkB >( dNdx, dTau_dqU ) ) * J0xWxRx2Pi;
    }
    // copy back to the subblocks using mighty Eigen block access,
    // note the layout swap rowmajor -> colmajor
    // Fastor offers similar functionality, but performancy is (slightly) inferior

    using namespace Eigen;

    Map< typename Parent::KSizedMatrix > K( stiffnessMatrix );

    K.template block< Parent::bsU, Parent::bsU >( idxU, idxU ) += Map< Matrix< double, Parent::bsU, Parent::bsU > >(
      torowmajor( k_UU ).data() );
  }

} // namespace Marmot::Elements
