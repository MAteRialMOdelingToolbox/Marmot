/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck
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
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotGeostaticStress.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotLowerDimensionalStress.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrain.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrainFactory.h"
#include "Marmot/MarmotMaterialPoint.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotOutOfPlaneTangents.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

#include <array>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace Marmot::Discretizations::MaterialPoints {

  /**
   * @brief Names of the nonlocal fields of the gradient-enhanced finite-strain discretizations.
   * @tparam nNonlocal Number of nonlocal variables.
   *
   * The node field follows the gradient-enhanced finite elements
   * (Marmot::Discretizations::Elements::GeneralGradientEnhancedDisplacementFiniteElement): six nonlocal variables are a
   * nonlocal strain tensor, the field `strain symmetric`; any other number is a set of scalar fields, `nonlocal
   * damage`, `nonlocal damage 2`, ... The state variables of a material point are named alike: `nonlocal strain` and
   * `local strain` for six, `nonlocal damage` and `local damage` otherwise.
   */
  template < int nNonlocal >
  struct GradientEnhancedFiniteStrainNames {

    /**
     * @brief The node fields of the nonlocal variables.
     * @return `strain symmetric` for six variables, else `nonlocal damage`, `nonlocal damage 2`, ...
     */
    static std::vector< std::string > nonlocalNodeFields()
    {
      if constexpr ( nNonlocal == 6 )
        return { "strain symmetric" };
      else {
        std::vector< std::string > fields = { "nonlocal damage" };
        for ( int K = 1; K < nNonlocal; K++ )
          fields.push_back( "nonlocal damage " + std::to_string( K + 1 ) );
        return fields;
      }
    }

    /**
     * @brief The state variable holding the nonlocal fields @f$ \bar{N}_K @f$.
     * @return `nonlocal strain` for six variables, else `nonlocal damage`.
     */
    static std::string nonlocalState() { return nNonlocal == 6 ? "nonlocal strain" : "nonlocal damage"; }

    /**
     * @brief The state variable holding the local driving forces @f$ L_K @f$.
     * @return `local strain` for six variables, else `local damage`.
     */
    static std::string localState() { return nNonlocal == 6 ? "local strain" : "local damage"; }
  };

  /**
   * @class Marmot::Discretizations::MaterialPoints::GradientEnhancedFiniteStrainMaterialPoint
   * @brief Material point for gradient-enhanced (implicit-gradient) finite-strain materials without a micropolar
   *        continuum.
   *
   * The material point carries the displacement field @f$ \boldsymbol{u} @f$ and @p nNonlocal nonlocal fields
   * @f$ \bar{N}_K @f$ and drives a Marmot::Materials::MaterialGradientEnhancedFiniteStrain, which returns the Kirchhoff
   * stress
   * @f$ \boldsymbol{\tau} @f$, the local driving forces @f$ L_K @f$ (the sources of @f$ \bar{N}_K - c_K\,\nabla^2
   * \bar{N}_K = L_K @f$, @f$ c_K = R_K^2 @f$ with the nonlocal radius @f$ R_K @f$) and the four algorithmic tangents.
   * It is the material point both of the MPM cell GradientEnhancedFiniteStrainCell and of the RKPM particle
   * GradientEnhancedFiniteStrainParticle (and its SQCNI / NSNI variants). There is no micro-rotation field and hence
   * no couple stress.
   *
   * **Kinematics (semi-Lagrangian).** The state variable `deformation gradient` (@c dY_dX) holds the deformation
   * gradient @f$ \boldsymbol{F}_n = \partial\boldsymbol{Y}/\partial\boldsymbol{X} @f$ of the last accepted state,
   * i.e. of the intermediate reference configuration @f$ \boldsymbol{Y} @f$ with respect to the undeformed
   * configuration @f$ \boldsymbol{X} @f$, and `delta deformation gradient` (@c dx_dY) the increment
   * @f$ \Delta\boldsymbol{F} = \partial\boldsymbol{x}/\partial\boldsymbol{Y} @f$ accumulated since. The material is
   * evaluated at
   * @f[
   *   \boldsymbol{F} = \Delta\boldsymbol{F}\,\boldsymbol{F}_n, \qquad \bar{N}_K = \bar{N}_{K,n} + \Delta\bar{N}_K,
   * @f]
   * where @f$ \Delta F_{ij} = \delta_{ij} + \partial\Delta u_i/\partial Y_j @f$ and @f$ \Delta\bar{N}_K @f$ are
   * accumulated by incrementDeformation() from the increments interpolated by the cell or particle. Because the
   * unknowns of an increment are
   * @f$ \Delta\boldsymbol{F} @f$ and @f$ \Delta\bar{N}_K @f$, computeYourself() converts the stress and driving-force
   * tangents with respect to @f$ \boldsymbol{F} @f$ to tangents with respect to @f$ \Delta\boldsymbol{F} @f$ by the
   * chain rule
   * @f[
   *   \frac{\partial F_{iI}}{\partial \Delta F_{jJ}} = \delta_{ij}\,F_{n,JI},
   * @f]
   * and stores them in #tangents together with the response in #response.
   *
   * **Totals for the explicit scheme.** The material point also keeps the running sum of
   * @f$ \nabla_X\Delta\bar{N}_K @f$ (state `<nonlocal> gradient`, i.e. @f$ \nabla_X\bar{N}_K @f$) and the rate
   * @f$ \dot{\bar{N}}_K @f$ of the last explicit increment (state `<nonlocal> rate`): the explicit residual of the
   * nonlocal balance is written in total form, and the momentum of a nonlocal field with a micro-inertia needs its
   * rate. The gradient is accumulated by incrementDeformation() in the implicit and the explicit scheme alike, so
   * that an explicit step may follow an implicit one.
   *
   * In 2D (@c nDim = 2) all tensors are stored in 3D, the material is evaluated through
   * Marmot::Materials::MaterialGradientEnhancedFiniteStrain::computePlaneStrain, and the in-plane parts are handed out.
   * The out-of-plane stretch @f$ \Delta F_{33} @f$ of the increment is 1 in plane strain; an axisymmetric particle sets
   * it to the hoop stretch increment (addToOutOfPlaneStretchIncrement()), so that @f$ F_{33} @f$ is the hoop stretch,
   * as for the axisymmetric elements. In 2D, the response and the tangents also hold the out-of-plane parts the hoop
   * terms of an axisymmetric particle need (Response::tau33, Tangents::outOfPlane, Tangents::dTau33_dN and
   * Tangents::dL_dDeltaF33).
   *
   * Registered as `GradientEnhancedFiniteStrain/PlaneStrain` and `GradientEnhancedFiniteStrain/3D` (one nonlocal
   * variable) and as `GradientEnhanced6FiniteStrain/3D` (six).
   *
   * @tparam nDim      Spatial dimension (2: plane strain, 3: 3D).
   * @tparam nNonlocal Number of nonlocal variables (1: a scalar nonlocal damage driver, 6: a nonlocal strain).
   */
  template < int nDim, int nNonlocal >
  class GradientEnhancedFiniteStrainMaterialPoint : public Marmot::Discretizations::MaterialPoints::MaterialPoint {

  public:
    using Material = Marmot::Materials::MaterialGradientEnhancedFiniteStrain< nNonlocal >; ///< the material interface
    using NonlocalVector = typename Material::NonlocalVector;              ///< a value per nonlocal variable
    using NonlocalMatrix = typename Material::NonlocalMatrix;              ///< a value per pair of nonlocal variables
    using Names          = GradientEnhancedFiniteStrainNames< nNonlocal >; ///< names of the nonlocal fields

  protected:
    constexpr static int _nVertices = 1;                     ///< a material point is represented by a single vertex

    using TensorD    = Fastor::Tensor< double, nDim >;       ///< vector of size nDim
    using TensorDD   = Fastor::Tensor< double, nDim, nDim >; ///< second-order tensor of size nDim
    using TensorDDDD = Fastor::Tensor< double, nDim, nDim, nDim, nDim >; ///< fourth-order tensor of size nDim

    int _mpNumber;                        ///< label of the material point, also passed to the material as its number

    TensorD _x0;                          ///< coordinates @f$ \boldsymbol{X} @f$ in the undeformed configuration

    double _vol0;                         ///< volume @f$ V_0 @f$ in the undeformed configuration
    double _density;                      ///< mass density in the reference configuration, from initializeYourself()

    std::unique_ptr< Material > material; ///< the assigned material

    /**
     * @class Marmot::Discretizations::MaterialPoints::GradientEnhancedFiniteStrainMaterialPoint::MPStateVarManager
     * @brief State variable layout of the material point, followed by the state of the material.
     *
     * All kinematic quantities are stored with 3 (vectors) or 9 (tensors) entries, also in plane strain. The entries
     * are packed one after the other, without padding.
     */
    class MPStateVarManager : public Marmot::Utilities::StateVarVectorManager {

      /// @cond INTERNAL (designated initializers are not parseable by the Sphinx C++ domain)

      // names and lengths of the state variables; the material state starts at `begin of material state`
      inline const static auto layout = makeLayout( {
        { .name = "displacement", .length = 3 },
        { .name = "velocity", .length = 3 },
        { .name = "acceleration", .length = 3 },
        { .name = "delta displacement", .length = 3 },
        { .name = "delta deformation gradient", .length = 9 },
        { .name = "deformation gradient", .length = 9 },
        { .name = Names::nonlocalState(), .length = nNonlocal },
        { .name = Names::localState(), .length = nNonlocal },
        { .name = "stress", .length = 9 },
        { .name = "F0 XX", .length = 1 },
        { .name = "F0 YY", .length = 1 },
        { .name = "F0 ZZ", .length = 1 },
        { .name = Names::nonlocalState() + " gradient", .length = 3 * nNonlocal },
        { .name = Names::nonlocalState() + " rate", .length = nNonlocal },
        { .name = "begin of material state", .length = 0 },
      } );

      /// @endcond

    public:
      Marmot::Math::Tensors::FastorTensors::TensorMap3d u;  ///< displacement of the last accepted state
      Marmot::Math::Tensors::FastorTensors::TensorMap3d v;  ///< velocity (set by the particle)
      Marmot::Math::Tensors::FastorTensors::TensorMap3d a;  ///< acceleration (set by the particle)
      Marmot::Math::Tensors::FastorTensors::TensorMap3d du; ///< displacement increment of the current step
      Marmot::Math::Tensors::FastorTensors::TensorMap33d
        dx_dY;                                    ///< deformation gradient increment @f$ \Delta\boldsymbol{F} @f$
      Marmot::Math::Tensors::FastorTensors::TensorMap33d
                                   dY_dX;         ///< deformation gradient @f$ \boldsymbol{F}_n @f$, last accepted
      Eigen::Map< NonlocalVector > nonlocalField; ///< nonlocal fields @f$ \bar{N}_K @f$ (totals; incrementDeformation()
                                                  ///< adds the increment of the step to the committed values, see the
                                                  ///< state contract of
                                                  ///< Marmot::Discretizations::MaterialPoints::MaterialPoint)
      Eigen::Map< NonlocalVector > localField;    ///< local driving forces @f$ L_K @f$ of the last material evaluation
      Marmot::Math::Tensors::FastorTensors::TensorMap33d stress; ///< Kirchhoff stress @f$ \boldsymbol{\tau} @f$ (3x3)
      double&                                            F0_XX;  ///< eigen deformation (geostatic stress), XX component
      double&                                            F0_YY;  ///< eigen deformation (geostatic stress), YY component
      double&                                            F0_ZZ;  ///< eigen deformation (geostatic stress), ZZ component
      /// @f$ \nabla_X\bar{N}_K @f$ in column K (3 rows, also in plane strain): the running sum of
      /// @f$ \nabla_X\Delta\bar{N}_K @f$ added by incrementDeformation()
      Eigen::Map< Eigen::Matrix< double, 3, nNonlocal > > nonlocalFieldGradient;
      Eigen::Map< NonlocalVector >  nonlocalFieldRate; ///< @f$ \dot{\bar{N}}_K @f$ of the last explicit increment
      Eigen::Map< Eigen::VectorXd > materialState;     ///< state variables of the material

      /**
       * @brief Number of state variables of the material point itself (without the material).
       * @return Size of the fixed part of the layout.
       */
      static int getNumberOfRequiredStateVars() { return layout.nRequiredStateVars; };

      /**
       * @brief Map the layout onto a state variable vector.
       * @param[in,out] theStateVarVector State variable vector of the material point.
       * @param[in]     nStateVars        Total size of the vector (material point + material).
       */
      MPStateVarManager( double* theStateVarVector, int nStateVars )
        : Marmot::Utilities::StateVarVectorManager( theStateVarVector, layout ),
          u( &find( "displacement" ) ),
          v( &find( "velocity" ) ),
          a( &find( "acceleration" ) ),
          du( &find( "delta displacement" ) ),
          dx_dY( &find( "delta deformation gradient" ) ),
          dY_dX( &find( "deformation gradient" ) ),
          nonlocalField( &find( Names::nonlocalState() ) ),
          localField( &find( Names::localState() ) ),
          stress( &find( "stress" ) ),
          F0_XX( find( "F0 XX" ) ),
          F0_YY( find( "F0 YY" ) ),
          F0_ZZ( find( "F0 ZZ" ) ),
          nonlocalFieldGradient( &find( Names::nonlocalState() + " gradient" ) ),
          nonlocalFieldRate( &find( Names::nonlocalState() + " rate" ) ),
          materialState( &find( "begin of material state" ), nStateVars - getNumberOfRequiredStateVars() ){};
    };

    std::unique_ptr< MPStateVarManager > state; ///< view on the assigned state variable vector

  public:
    bool hasEigenDeformation; ///< true once a geostatic initial condition has set an eigen deformation

    /**
     * @brief Construct a material point.
     * @param[in] mpNumber           Label of the material point.
     * @param[in] vertexCoordinates  Undeformed coordinates @f$ \boldsymbol{X} @f$ (nDim values).
     * @param[in] nVertexCoordinates Number of coordinates (not checked).
     * @param[in] volume             Undeformed volume @f$ V_0 @f$.
     */
    GradientEnhancedFiniteStrainMaterialPoint( int           mpNumber,
                                               const double* vertexCoordinates,
                                               int           nVertexCoordinates,
                                               double        volume )
      : _mpNumber( mpNumber ), hasEigenDeformation( false )
    {
      assignVertexCoordinates( vertexCoordinates );
      assignVolume( volume );
    };

    /**
     * @brief Assign the state variable vector (material point layout followed by the material state).
     * @param[in,out] stateVars  State variable vector, owned by the host.
     * @param[in]     nStateVars Its size, getNumberOfRequiredStateVars().
     */
    void assignStateVars( double* stateVars, int nStateVars ) override;

    /**
     * @brief Access a state variable of the material point or, if it has none of that name, of the material.
     * @param[in] stateName Name of the state variable (see MPStateVarManager).
     * @return View on the state variable.
     */
    Marmot::Utilities::StateView getStateView( const std::string& stateName ) const override;

    /**
     * @brief Shape of the material point.
     * @return Always "point".
     */
    std::string getMaterialPointShape() const override { return "point"; };

    /**
     * @brief Create the material from the Marmot::Registry::MaterialGradientEnhancedFiniteStrainFactory of
     * @p nNonlocal nonlocal variables.
     * @param[in] property Material section (name and properties).
     * @throws std::invalid_argument if no such material is registered.
     */
    void assignMaterial( const Marmot::Discretizations::MaterialSection& property ) override;

    /**
     * @brief Initialize the state: @f$ \boldsymbol{F}_n = \boldsymbol{I} @f$, unit eigen deformation, a zero
     * increment, the material state, and the density.
     *
     * The density is queried from the material here already, so that the inertia can be assembled before the first
     * computation. Requires assignMaterial() and assignStateVars() to have been called.
     */
    void initializeYourself() override;

    /**
     * @brief Label of the material point.
     * @return The label passed to the constructor.
     */
    int getMaterialPointNumber() const override { return _mpNumber; }

    /**
     * @brief Spatial dimension.
     * @return nDim.
     */
    int getDimension() const override { return nDim; }

    /**
     * @brief Number of vertices.
     * @return 1.
     */
    int getNumberOfVertices() const override { return _nVertices; };

    /**
     * @brief Number of state variables: material point layout plus material.
     * @return Required size of the state variable vector; requires an assigned material.
     */
    int getNumberOfRequiredStateVars() const override
    {
      return MPStateVarManager::getNumberOfRequiredStateVars() + material->getNumberOfRequiredStateVars();
    };

    /**
     * @brief Set the undeformed volume.
     * @param[in] volume Undeformed volume @f$ V_0 @f$.
     */
    void assignVolume( double volume ) { _vol0 = volume; };

    /**
     * @brief Undeformed volume.
     * @return @f$ V_0 @f$.
     */
    double getVolumeUndeformed() const override { return _vol0; }

    /**
     * @brief Set the undeformed coordinates.
     * @param[in] coordinates Undeformed coordinates @f$ \boldsymbol{X} @f$ (nDim values).
     */
    void assignVertexCoordinates( const double* coordinates ) { _x0 = TensorD( coordinates ); };

    /**
     * @brief Coordinates of the single vertex, identical to getCoordinatesAtCenter().
     * @param[out] coordinates Coordinates @f$ \boldsymbol{X} + \boldsymbol{u} @f$ (nDim values).
     */
    void getVertexCoordinates( double* coordinates ) const override { return getCoordinatesAtCenter( coordinates ); };

    /**
     * @brief Coordinates of the last accepted state, @f$ \boldsymbol{X} + \boldsymbol{u} @f$.
     *
     * The displacement increment of the current step is not included, so this is the position in the intermediate
     * reference configuration @f$ \boldsymbol{Y} @f$.
     *
     * @param[out] coordinates Coordinates (nDim values).
     */
    void getCoordinatesAtCenter( double* coordinates ) const override
    {
      Eigen::Map< const Eigen::Matrix< double, nDim, 1 > > x0( _x0.data() );
      Eigen::Map< Eigen::Matrix< double, nDim, 1 > >       newCoords( coordinates );
      Eigen::Map< Eigen::Matrix< double, nDim, 1 > >       u( state->u.data() );

      newCoords = x0 + u;
    };

    /**
     * @brief Displacement of the last accepted state.
     * @param[out] displacement Displacement @f$ \boldsymbol{u} @f$ (nDim values).
     */
    void getCenterDisplacement( double* displacement ) const override
    {
      for ( int i = 0; i < nDim; i++ )
        displacement[i] = state->u( i );
    };

    /**
     * @brief Mass density in the reference configuration.
     * @return Density reported by the material at initializeYourself(), the mass of the material point over its
     * reference volume: an evaluation does not change it.
     */
    double getDensityUndeformed() const override { return _density; };

    /**
     * @brief Undeformed coordinates.
     * @return @f$ \boldsymbol{X} @f$.
     */
    const TensorD& coordinates() const { return _x0; };

    /**
     * @brief Reset the increment of the current step: @f$ \Delta\boldsymbol{u} = \boldsymbol{0} @f$,
     * @f$ \Delta\boldsymbol{F} = \boldsymbol{I} @f$.
     *
     * The nonlocal fields are not reset; incrementDeformation() adds to the stored totals @f$ \bar{N}_K @f$.
     *
     * @param[in] timeNew Time at the end of the increment (unused).
     * @param[in] dT      Time increment (unused).
     */
    void prepareYourself( double timeNew, double dT ) override;

    /**
     * @brief Evaluate the material at @f$ \boldsymbol{F} = \Delta\boldsymbol{F}\,\boldsymbol{F}_n @f$ and
     * @f$ \bar{N}_K @f$ and fill #response and #tangents.
     *
     * In 3D through Marmot::Materials::MaterialGradientEnhancedFiniteStrain::computeStress; in plane strain through
     * Marmot::Materials::MaterialGradientEnhancedFiniteStrain::computePlaneStrain with the 3D deformation gradient
     * whose out-of-plane entries are those of the identity, and the response and tangents are reduced to their in-plane
     * components. With the eigen deformation of a geostatic initial condition, if one is set.
     *
     * @param[in] timeNew Time at the end of the increment.
     * @param[in] dT      Time increment.
     */
    void computeYourself( double timeNew, double dT ) override;

    /**
     * @brief Accept the current increment: @f$ \boldsymbol{u} \leftarrow \boldsymbol{u} + \Delta\boldsymbol{u} @f$,
     * @f$ \boldsymbol{F}_n \leftarrow \Delta\boldsymbol{F}\,\boldsymbol{F}_n @f$.
     *
     * The increment itself is not reset here; the next prepareYourself() does that. A second accept before it would
     * apply the increment twice: a particle, whose accept must be idempotent, calls prepareYourself() right after
     * (GradientEnhancedFiniteStrainParticle::acceptStateAndPosition()).
     */
    void acceptStateAndPosition() override
    {
      const auto&                                          u_n  = state->u;
      const Marmot::Math::Tensors::FastorTensors::Tensor3d u_np = u_n + state->du;

#if defined( __GNUC__ ) && !defined( __clang__ )
// GCC's -Wmaybe-uninitialized misfires on Fastor's AVX copies of the mapped state tensors: the
// mapped memory is initialized by the host, which GCC cannot see (observed with GCC 14 at x86-64-v3).
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
      const Marmot::Math::Tensors::FastorTensors::Tensor33d dx_dX_n  = state->dY_dX;
      const Marmot::Math::Tensors::FastorTensors::Tensor33d dx_dX_np = state->dx_dY % dx_dX_n;
#if defined( __GNUC__ ) && !defined( __clang__ )
#  pragma GCC diagnostic pop
#endif

      Marmot::Math::Tensors::mapEigenToFastor( state->u )     = Marmot::Math::Tensors::mapEigenToFastor( u_np );
      Marmot::Math::Tensors::mapEigenToFastor( state->dY_dX ) = Marmot::Math::Tensors::mapEigenToFastor( dx_dX_np );
    };

    /**
     * @brief Add an interpolated increment to the state.
     *
     * @f$ \Delta\boldsymbol{u} \mathrel{+}= \delta\boldsymbol{u} @f$,
     * @f$ \Delta\boldsymbol{F} \mathrel{+}= \partial\delta\boldsymbol{u}/\partial\boldsymbol{Y} @f$,
     * @f$ \bar{N}_K \mathrel{+}= \delta\bar{N}_K @f$ and
     * @f$ \nabla_X\bar{N}_K \mathrel{+}= \boldsymbol{F}_n^T\,\partial\delta\bar{N}_K/\partial\boldsymbol{Y} @f$; in
     * plane strain the inputs are expanded to 3D.
     *
     * @param[in] displacementIncrement         Displacement increment @f$ \delta\boldsymbol{u} @f$.
     * @param[in] displacementGradientIncrement Its gradient with respect to the intermediate reference
     *                                          configuration, @f$ \partial\delta u_i/\partial Y_j @f$.
     * @param[in] nonlocalFieldIncrement         Increments @f$ \delta\bar{N}_K @f$ of the nonlocal fields.
     * @param[in] nonlocalFieldGradientIncrement Their gradients with respect to @f$ \boldsymbol{Y} @f$, one per
     *                                           nonlocal variable.
     */
    void incrementDeformation( const TensorD&                          displacementIncrement,
                               const TensorDD&                         displacementGradientIncrement,
                               const NonlocalVector&                   nonlocalFieldIncrement,
                               const std::array< TensorD, nNonlocal >& nonlocalFieldGradientIncrement );

    /**
     * @brief 2D: add to the out-of-plane stretch increment, @f$ \Delta F_{33} \mathrel{+}= \delta @f$ (@f$ \Delta
     * F_{33} @f$ is 1 after prepareYourself(), and stays 1 in plane strain): the hoop stretch increment of an
     * axisymmetric particle.
     * @param[in] delta @f$ \delta @f$, e.g. @f$ \Delta u_r / \bar r @f$.
     * @throws std::logic_error for a 3D material point.
     */
    void addToOutOfPlaneStretchIncrement( [[maybe_unused]] double delta )
    {
      if constexpr ( nDim == 2 )
        state->dx_dY( 2, 2 ) += delta;
      else
        throw std::logic_error( "GradientEnhancedFiniteStrainMaterialPoint: the out-of-plane stretch is 2D only" );
    };

    /**
     * @brief 2D: the out-of-plane stretch increment, the 33 component of @f$ \Delta\boldsymbol{F} @f$ (not of
     * @f$ \boldsymbol{F} @f$).
     * @return @f$ \Delta F_{33} @f$ (1 in plane strain, the hoop stretch increment of an axisymmetric particle).
     * @throws std::logic_error for a 3D material point.
     */
    double outOfPlaneStretchIncrement() const
    {
      if constexpr ( nDim == 2 )
        return state->dx_dY( 2, 2 );
      else
        throw std::logic_error( "GradientEnhancedFiniteStrainMaterialPoint: the out-of-plane stretch is 2D only" );
    };

    /// @brief Response of the last computeYourself(), in nDim components.
    struct Response {
      Fastor::Tensor< double, nDim, nDim > S;  ///< Kirchhoff stress @f$ \boldsymbol{\tau} @f$
      NonlocalVector                       dL; ///< @f$ \Delta L_K @f$: the local driving forces of this evaluation
                                               ///< minus the values stored in the local state variable (which is
                                               ///< then overwritten by the new values)
      NonlocalVector nonLocalRadius;           ///< nonlocal radii @f$ R_K @f$ of the material, @f$ c_K = R_K^2 @f$
      double         tau33 = 0.0;              ///< 2D: the out-of-plane (hoop) Kirchhoff stress @f$ \tau_{33} @f$
    } response;

    /// @brief Algorithmic tangents of the last computeYourself(), with respect to @f$ \Delta\boldsymbol{F} @f$ and
    /// @f$ \bar{N}_K @f$, in nDim components.
    struct Tangents {
      TensorDDDD                        dS_dDeltaF;         ///< @f$ \partial\tau_{ij}/\partial\Delta F_{kL} @f$
      std::array< TensorDD, nNonlocal > dS_dN;              ///< @f$ \partial\tau_{ij}/\partial\bar{N}_K @f$, one per K
      std::array< TensorDD, nNonlocal > dL_dDeltaF;         ///< @f$ \partial L_K/\partial\Delta F_{kL} @f$, one per K
      NonlocalMatrix                    dL_dN;              ///< @f$ \partial L_K/\partial\bar{N}_M @f$
      OutOfPlaneTangents                outOfPlane;         ///< 2D: the out-of-plane parts with respect to
                                                            ///< @f$ \Delta\boldsymbol{F} @f$
      NonlocalVector dTau33_dN    = NonlocalVector::Zero(); ///< 2D: @f$ \partial\tau_{33}/\partial\bar{N}_K @f$
      NonlocalVector dL_dDeltaF33 = NonlocalVector::Zero(); ///< 2D: @f$ \partial L_K/\partial\Delta F_{33} @f$
    } tangents;

    /**
     * @brief Deformation gradient increment.
     * @return @f$ \Delta\boldsymbol{F} = \partial\boldsymbol{x}/\partial\boldsymbol{Y} @f$ (nDim x nDim).
     */
    TensorDD dx_dY() const { return state->dx_dY( Fastor::seq( 0, nDim ), Fastor::seq( 0, nDim ) ); };

    /**
     * @brief Deformation gradient of the last accepted state.
     * @return @f$ \boldsymbol{F}_n = \partial\boldsymbol{Y}/\partial\boldsymbol{X} @f$ (nDim x nDim).
     */
    TensorDD dY_dX() const { return state->dY_dX( Fastor::seq( 0, nDim ), Fastor::seq( 0, nDim ) ); };

    /**
     * @brief Velocity.
     * @return Stored velocity (nDim components).
     */
    TensorD getVelocity() const { return state->v( Fastor::seq( 0, nDim ) ); };

    /**
     * @brief Acceleration.
     * @return Stored acceleration (nDim components).
     */
    TensorD getAcceleration() const { return state->a( Fastor::seq( 0, nDim ) ); };

    /**
     * @brief Store the velocity.
     * @param[in] velocity Velocity (nDim components).
     */
    void setVelocity( const TensorD& velocity )
    {
      for ( int i = 0; i < nDim; i++ )
        state->v( i ) = velocity( i );
    };

    /**
     * @brief Store the acceleration.
     * @param[in] acceleration Acceleration (nDim components).
     */
    void setAcceleration( const TensorD& acceleration )
    {
      for ( int i = 0; i < nDim; i++ )
        state->a( i ) = acceleration( i );
    };

    /**
     * @brief The nonlocal fields, totals.
     * @return @f$ \bar{N}_K @f$, including the increment of the current step.
     */
    NonlocalVector getNonlocalField() const { return state->nonlocalField; };

    /**
     * @brief The local driving forces of the last computeYourself(), totals.
     * @return @f$ L_K @f$.
     */
    NonlocalVector getLocalField() const { return state->localField; };

    /**
     * @brief Gradient of the nonlocal fields with respect to the undeformed configuration, totals.
     * @return @f$ \nabla_X\bar{N}_K @f$, one vector of nDim components per nonlocal variable.
     */
    std::array< TensorD, nNonlocal > getNonlocalFieldGradient() const
    {
      std::array< TensorD, nNonlocal > gradient;
      for ( int K = 0; K < nNonlocal; K++ )
        for ( int i = 0; i < nDim; i++ )
          gradient[K]( i ) = state->nonlocalFieldGradient( i, K );
      return gradient;
    };

    /**
     * @brief Rate of the nonlocal fields.
     * @return @f$ \dot{\bar{N}}_K @f$ as stored by setNonlocalFieldRate().
     */
    NonlocalVector getNonlocalFieldRate() const { return state->nonlocalFieldRate; };

    /**
     * @brief Store the rate of the nonlocal fields.
     * @param[in] rate @f$ \dot{\bar{N}}_K @f$.
     */
    void setNonlocalFieldRate( const NonlocalVector& rate ) { state->nonlocalFieldRate = rate; };

    /**
     * @brief Nonlocal viscosities of the assigned material at the current state.
     * @return @f$ \eta_K @f$.
     */
    NonlocalVector getNonlocalViscosity() const
    {
      return material->getNonlocalViscosity( state->materialState.data() );
    };

    /**
     * @brief Nonlocal micro-inertias of the assigned material at the current state.
     * @return @f$ m_K @f$; zero for a nonlocal field that is first order in time.
     */
    NonlocalVector getNonlocalMicroInertia() const
    {
      return material->getNonlocalMicroInertia( state->materialState.data() );
    };

    /**
     * @brief Apply an initial condition.
     *
     * Only `geostaticstress` is supported: the eigen deformation that produces the hydrostatic stress
     * @f$ \tau_{XX} = \tau_{YY} = \tau_{ZZ} = \mathrm{value}[0] @f$ is found by
     * Marmot::Materials::MaterialGradientEnhancedFiniteStrain::findEigenDeformationForEigenStress and applied in every
     * subsequent material evaluation.
     *
     * @param[in] conditionName Name of the condition.
     * @param[in] value         Its value(s); `value[0]` is the geostatic normal stress.
     * @throws std::invalid_argument for any other condition.
     */
    void setInitialCondition( const std::string& conditionName, const double* value ) override
    {
      if ( conditionName == "geostaticstress" ) {
        std::tuple< double, double, double > geostaticNormalStressComponents = { value[0], value[0], value[0] };
        const auto [F0_XX,
                    F0_YY,
                    F0_ZZ] = material->findEigenDeformationForEigenStress( { state->F0_XX, state->F0_YY, state->F0_ZZ },
                                                                           geostaticNormalStressComponents,
                                                                           state->materialState.data() );

        state->F0_XX = F0_XX;
        state->F0_YY = F0_YY;
        state->F0_ZZ = F0_ZZ;

        hasEigenDeformation = true;
      }
      else {
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__ << ": invalid initial condition" );
      }
    }

  protected:
    /**
     * @brief The total deformation gradient @f$ \boldsymbol{F} = \Delta\boldsymbol{F}\,\boldsymbol{F}_n @f$ (3x3).
     * @return The deformation gradient the material is evaluated at.
     */
    Marmot::Math::Tensors::FastorTensors::Tensor33d deformationGradient() const;
  };

  extern template class GradientEnhancedFiniteStrainMaterialPoint< 2, 1 >;
  extern template class GradientEnhancedFiniteStrainMaterialPoint< 3, 1 >;
  extern template class GradientEnhancedFiniteStrainMaterialPoint< 3, 6 >;

} // namespace Marmot::Discretizations::MaterialPoints
