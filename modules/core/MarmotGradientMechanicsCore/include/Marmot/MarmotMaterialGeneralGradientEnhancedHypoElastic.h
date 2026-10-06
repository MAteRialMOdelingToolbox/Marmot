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
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotStateHelpers.h"
#include "Marmot/MarmotTypedefs.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

/**
 * @brief Base class for general gradient-enhanced hypoelastic material models.
 * @details This class defines the interface for general gradient-enhanced hypoelastic material models,
 * including methods for computing stress, plane stress response, and managing state variables.
 * The template parameter `nNonlocalVariables` specifies the number of nonlocal variables used in the material model.
 *
 * In addition to the standard balance of linear momentum each nonlocal variable introduces an additional balance
 * equation, which is solved simultaneously with the balance of linear momentum. This balance equation is defined as:
 * \f[ \knl_i - \nabla ( c( \knl_i )\, \nabla \knl_i ) = s_i (\boldsymbol \varepsilon,\, \knl_i ) \f]
 * where \f$ \knl_i \f$ is the nonlocal variable, \f$ c( \knl_i ) \f$ is the nonlocal interaction parameter, and
 * \f$ s_i \f$ is the local driving variable for the nonlocal variable \f$ \knl_i \f$. The interaction
 * parameter
 * \f$ c( \knl_i ) \f$ defines the influence of the nonlocal variable on its own gradient and can be used to model
 * phenomena such as damage-dependent interactions. Also phase-field models can be implemented in this framework by
 * defining the nonlocal variable as the phase-field variable and the local driving variable as a function of the strain
 * tensor.
 */
template < int nNonlocalVariables >
class MarmotMaterialGeneralGradientEnhancedHypoElastic {

protected:
  /// @brief Pointer to the array of material properties.
  const double* materialProperties;
  /// @brief Number of material properties.
  const int nMaterialProperties;

public:
  /// @brief Material number (identifier for the material).
  const int materialNumber;

  /**
   * @brief Constructor for the general gradient-enhanced hypoelastic material model.
   * @param matProperties_ Pointer to the array of material properties.
   * @param nMaterialProperties_ Number of material properties.
   * @param materialNumber_ Material number (identifier for the material).
   */
  MarmotMaterialGeneralGradientEnhancedHypoElastic( const double* matProperties_,
                                                    int           nMaterialProperties_,
                                                    int           materialNumber_ )
    : materialProperties( matProperties_ ),
      nMaterialProperties( nMaterialProperties_ ),
      materialNumber( materialNumber_ )
  {
  }

  /**
   * @brief Virtual destructor.
   *
   * This class is abstract and instances are owned through base-class pointers
   * (e.g. the std::unique_ptr held by the gradient-enhanced elements), so the
   * destructor must be virtual — destroying a derived object through a base
   * pointer with a non-virtual destructor is undefined behaviour. libstdc++
   * silently invokes the wrong destructor, while libc++/clang diagnoses it
   * (-Wdelete-abstract-non-virtual-dtor) and emits a trap instruction, which
   * aborted every simulation using such a material on macOS/arm64.
   */
  virtual ~MarmotMaterialGeneralGradientEnhancedHypoElastic() = default;

  /// @brief Struct to hold the increment information.
  struct increment {
    Marmot::Vector6d                            dStrain; ///< Increment of the strain tensor in Voigt notation
    Eigen::Vector< double, nNonlocalVariables > K;       ///< Nonlocal variables at the current increment
    Eigen::Vector< double, nNonlocalVariables > dK;      ///< Increment of the nonlocal variables
    double                                      time;    ///< Current time
    double                                      dT;      ///< Time increment
  };

  /// @brief Struct to hold the material response information.
  struct response {
    Marmot::Vector6d                            stress;    ///< Stress tensor in Voigt notation
    Eigen::Vector< double, nNonlocalVariables > KLocal;    ///< Local driving variables at the current increment
    Eigen::Vector< double, nNonlocalVariables > c;         ///< Nonlocal interaction parameters at the current increment
    double*                                     stateVars; ///< Pointer to the array of state variables
    double                                      elasticEnergyDensity; ///< Elastic strain energy density
    double                                      dissipation;          ///< Dissipation (if applicable)
  };

  /// @brief Struct to hold the algorithmic tangent matrices for the material model.
  struct tangents {
    /// @brief Algorithmic tangent matrix relating the increment of stress to the increment of strain.
    Marmot::Matrix6d dStressddStrain = Marmot::Matrix6d::Zero();
    /// @brief Algorithmic tangent matrix relating the increment of stress to the increment of nonlocal variables.
    Eigen::Matrix< double, 6, nNonlocalVariables > dStressddK = Eigen::Matrix< double, 6, nNonlocalVariables >::Zero();
    /// @brief Algorithmic tangent matrix relating the increment of local driving variables to the increment of strain.
    Eigen::Matrix< double, nNonlocalVariables, 6 >
      dKLocalddStrain = Eigen::Matrix< double, nNonlocalVariables, 6 >::Zero();
    /// @brief Algorithmic tangent matrix relating the increment of local driving variables to the increment of nonlocal
    /// variables.
    Eigen::Matrix< double, nNonlocalVariables, nNonlocalVariables >
      dKLocalddK = Eigen::Matrix< double, nNonlocalVariables, nNonlocalVariables >::Zero();
    /// @brief First derivative of nonlocal interaction parameters with respect to nonlocal variables.
    Eigen::Matrix< double, nNonlocalVariables, nNonlocalVariables >
      dcddK = Eigen::Matrix< double, nNonlocalVariables, nNonlocalVariables >::Zero();
    /// @brief Second derivative of nonlocal interaction parameters with respect to nonlocal variables.
    Eigen::Matrix< double, nNonlocalVariables, nNonlocalVariables >
      d2cddK2 = Eigen::Matrix< double, nNonlocalVariables, nNonlocalVariables >::Zero();
  };

  /**
   * @brief Layout of the state variables for the material model.
   * @note Must be defined in derived classes to specify the structure and organization of the state variables used in
   * the material model.
   */
  MarmotStateLayoutDynamic stateLayout;

  /**
   * @brief Compute the stress response of the material model including the algorithmic tangents.
   * @param[in,out] res Reference to the response struct to be filled with the computed stress and other response
   * variables.
   * @param[in,out] tan Reference to the tangents struct to be filled with the computed algorithmic tangent matrices.
   * @param[in] inc Reference to the increment struct containing the strain increment, nonlocal variable increments, and
   * time information.
   *
   * This method must be implemented in derived classes to compute the stress response based on the specific material
   * model.
   */
  virtual void computeStress( response& res, tangents& tan, const increment& inc ) const = 0;

  /**
   * @brief Compute the stress response of the material model.
   * @param[in,out] res Reference to the response struct to be filled with the computed stress and other response
   * variables.
   * @param[in] inc Reference to the increment struct containing the strain increment, nonlocal variable increments, and
   * time information.
   *
   * This method can be used when only the stress response is needed without the tangents, e.g., for explicit dynamics.
   * The default implementation calls the `computeStress` method and ignores the tangents, which can be computationally
   * expensive to compute.
   */
  virtual void computeStressExplicit( response& res, const increment& inc ) const
  {
    tangents tan;
    computeStress( res, tan, inc );
  }

  /**
   * @brief Compute the plane stress response of the material model.
   * @param[in,out] res Reference to the response struct to be filled with the computed plane stress and other response
   * variables.
   * @param[in,out] tan Reference to the tangents struct to be filled with the computed algorithmic tangent matrices for
   * plane stress.
   * @param[in] inc Reference to the increment struct containing the strain increment, nonlocal variable increments, and
   * time information.
   *
   * This method provides a default implementation for computing the plane stress response using an iterative approach.
   * It assumes an initial guess of isochoric deformation and iteratively corrects the out-of-plane strain until
   * convergence is achieved. Derived classes can override this method if a different approach for plane stress
   * computation is desired.
   */
  virtual void computePlaneStress( response& res, tangents& tan, const increment& inc ) const
  {
    using namespace Marmot;
    using namespace Eigen;

    Map< VectorXd > stateVars( res.stateVars, stateLayout.totalSize() );

    VectorXd  stateVarsOld = stateVars;
    response  resTemp = { res.stress, res.KLocal, res.c, res.stateVars, res.elasticEnergyDensity, res.dissipation };
    increment incTemp = { inc.dStrain, inc.K, inc.dK, inc.time, inc.dT };

    double residual          = 1;
    double tangentCompliance = 1.;
    // assumption of isochoric deformation for initial guess
    double strainIncrement = ( -incTemp.dStrain( 0 ) - incTemp.dStrain( 1 ) );
    incTemp.dStrain( 2 )   = strainIncrement;

    int planeStressCount = 1;
    while ( true ) {

      // set old response
      resTemp = { .stress               = res.stress,
                  .KLocal               = res.KLocal,
                  .c                    = res.c,
                  .stateVars            = res.stateVars,
                  .elasticEnergyDensity = res.elasticEnergyDensity,
                  .dissipation          = res.dissipation };
      // set old state variables
      stateVars = stateVarsOld;
      // compute stress
      computeStress( resTemp, tan, incTemp );

      // evauate residual
      residual = std::abs( resTemp.stress.array().abs()[2] / std::max( resTemp.stress.array().abs().maxCoeff(), 1. ) );
      if ( ( residual < 1e-10 && std::abs( strainIncrement ) < 1e-8 ) || ( planeStressCount > 7 && residual < 1e-5 ) ) {
        break;
      }

      // correct strain increment
      tangentCompliance = 1. / tan.dStressddStrain( 2, 2 );
      if ( Math::isNaN( tangentCompliance ) || std::abs( tangentCompliance ) > 1e10 ) {
        tangentCompliance = 1e10;
      }

      strainIncrement = -resTemp.stress( 2 ) * tangentCompliance;
      incTemp.dStrain( 2 ) += strainIncrement;

      planeStressCount += 1;
      if ( planeStressCount > 13 ) {
        throw Marmot::StressUpdateFailed( "PlaneStressWrapper requires cutback" );
      }
    }

    res = resTemp;
  }

  /**
   * @brief Compute the plane stress response of the material model without computing tangents.
   * @param[in,out] res Reference to the response struct to be filled with the computed plane stress and other response
   * variables.
   * @param[in] inc Reference to the increment struct containing the strain increment, nonlocal variable increments, and
   * time information.
   *
   * This method can be used when only the plane stress response is needed without the tangents, e.g., for explicit
   * dynamics. The default implementation calls the `computePlaneStress` method ignoring the tangents.
   * This method maybe overridden in derived classes if a different approach for plane stress computation is desired.
   */
  virtual void computePlaneStressExplicit( response& res, const increment& inc ) const
  {
    tangents tan;
    computePlaneStress( res, tan, inc );
  }

  /**
   * @brief Get a view to the state variables.
   * @param stateName Name of the state variable
   * @param stateVars Pointer to the state variable array
   * @return StatView to access the state variable
   */
  StateView getStateView( const std::string& stateName, double* stateVars ) const
  {
    return stateLayout.getStateView( stateVars, stateName );
  }

  /**
   * @brief Get the total number of required state variables.
   * @return Total number of required state variables
   */
  int getNumberOfRequiredStateVars() const { return stateLayout.totalSize(); }

  /**
   * @brief Initialize the state variables at a material point.
   * @param stateVars Pointer to the state variable array
   * @param nStateVars Number of state variables
   *
   * @note The default implementation initializes all state variables to zero.
   */
  virtual void initializeYourself( double* stateVars, int nStateVars )
  {
    for ( int i = 0; i < nStateVars; ++i ) {
      stateVars[i] = 0.0;
    }
  }
  /**
   * @brief Get the density of the material.
   * @param stateVars Pointer to the array of state variables
   * @return Density of the material
   *
   * This method must be implemented in derived classes to return the density of the material, which is required for
   * dynamic analyses.
   */
  virtual double getDensity( const double* stateVars ) const = 0;

  /**
   * @brief Get the maximum wave speed for the current response state.
   * @param currentResponse Current response state
   * @param K The non-local field values the tangent is to be evaluated at; zero (the default)
   *          asks for the response of the material with its non-local field undamaged.
   * @return Maximum wave speed
   * @details The default implementation computes the 3D algorithmic tangent and returns
   *          `sqrt(max(C_ii) / rho)` with `C_ii` from the Voigt tangent diagonal entries.
   *
   * @note The field is passed in because it is an INPUT to the constitutive law, not a state the
   *       response carries. Leaving it at zero therefore means "zero", not "whatever it currently
   *       is" -- and where damage is driven by that field alone (`m = 1` in GCDP) that is the
   *       virgin tangent however damaged the point is. Pass the current field for the current wave
   *       speed; take the default for an undamaged reference or a conservative critical time step.
   */
  virtual double getMaximumWaveSpeed(
    const response&                                    currentResponse,
    const Eigen::Vector< double, nNonlocalVariables >& K = Eigen::Vector< double, nNonlocalVariables >::Zero() ) const
  {
    const int nStateVars = getNumberOfRequiredStateVars();

    std::vector< double > stateVarsCopy( nStateVars, 0.0 );
    if ( currentResponse.stateVars != nullptr && nStateVars > 0 ) {
      std::copy_n( currentResponse.stateVars, nStateVars, stateVarsCopy.begin() );
    }

    response responseCopy  = currentResponse;
    responseCopy.stateVars = stateVarsCopy.data();

    tangents  tan;
    increment inc;
    inc.dStrain = Marmot::Vector6d::Zero();
    inc.K       = K;
    inc.dK      = Eigen::Vector< double, nNonlocalVariables >::Zero();
    inc.time    = 0.0;
    inc.dT      = 1.0;

    computeStress( responseCopy, tan, inc );

    const double maxStiffnessDiagonal = std::max( { tan.dStressddStrain( 0, 0 ),
                                                    tan.dStressddStrain( 1, 1 ),
                                                    tan.dStressddStrain( 2, 2 ),
                                                    tan.dStressddStrain( 3, 3 ),
                                                    tan.dStressddStrain( 4, 4 ),
                                                    tan.dStressddStrain( 5, 5 ) } );
    const double density              = getDensity( responseCopy.stateVars );

    return density > 0.0 ? std::sqrt( std::max( 0.0, maxStiffnessDiagonal ) / density ) : 0.0;
  }

  /**
   * @brief Get the nonlocal viscosity of the material.
   * @param stateVars Pointer to the array of state variables
   * @return Vector containing the nonlocal viscosity values for each nonlocal variable
   *
   * This method must be implemented in derived classes to return the nonlocal viscosity values, which are required for
   * dynamic analyses involving nonlocal variables. The returned vector should have a length equal to
   * `nNonlocalVariables`.
   */
  virtual std::vector< double > getNonlocalViscosity( const double* stateVars ) const = 0;

  /**
   * @brief Get the micro-inertia \f$m_k\f$ of each nonlocal variable.
   * @param stateVars Pointer to the array of state variables
   * @return Vector containing the micro-inertia of each nonlocal variable, in seconds squared
   *
   * @details It multiplies the second time derivative of the nonlocal variable, turning its balance from the
   * parabolic (viscous) equation
   * \f$ \eta\,\dot{\bar\varepsilon} + \bar\varepsilon - c\,\nabla^2\bar\varepsilon = \tilde\varepsilon \f$
   * into the damped hyperbolic one
   * \f$ m_k\,\ddot{\bar\varepsilon} + \eta\,\dot{\bar\varepsilon} + \bar\varepsilon -
   * c\,\nabla^2\bar\varepsilon = \tilde\varepsilon \f$, whose stable increment falls off with \f$h\f$ rather than
   * with \f$h^2\f$. The nonlocal viscosity keeps its meaning exactly: what was the coefficient of the highest time
   * derivative becomes the damping.
   *
   * It sits here, next to getNonlocalViscosity(), because the two are **one parameter and not two**: the zeroth-order
   * reaction mode does not ring only while \f$m_k \le \eta^2/4\f$, and the largest admissible value is the best one
   * because the stable increment grows with \f$\sqrt{m_k}\f$. Nothing could check that while the two lived at
   * different levels of the stack, and the deck author computed \f$\eta^2/4\f$ by hand; see
   * validatedNonlocalMicroInertia(), which derived classes should return through.
   *
   * @note The default is zero for every nonlocal variable, which is the quasi-static model and the behaviour of every
   * material that predates this interface: the field is then first order in time and is integrated by its viscosity
   * alone. A derived class opts in by overriding this, and a run opts in by giving the material the property the
   * override reads -- so this is deliberately not pure virtual, unlike getNonlocalViscosity(), which a
   * gradient-enhanced material must answer for any explicit path to exist at all.
   */
  virtual std::vector< double > getNonlocalMicroInertia( const double* stateVars ) const
  {
    return std::vector< double >( nNonlocalVariables, 0.0 );
  }

protected:
  /**
   * @brief Validate a micro-inertia against this material's own nonlocal viscosity.
   * @param microInertia The micro-inertia of each nonlocal variable, as read from the material properties
   * @param stateVars Pointer to the array of state variables, to evaluate the nonlocal viscosity at
   * @return @p microInertia unchanged, if every entry is admissible
   * @throws std::invalid_argument if an entry is non-finite, negative, or exceeds \f$\eta^2/4\f$
   *
   * @details Three failures, each of which is silent without this check:
   *
   * - a **non-finite** micro-inertia passes every ordering test, so a NaN would propagate into the lumped inertia
   *   and from there into every increment;
   * - a **negative** one would make the nonlocal field integrate backwards in time;
   * - one **above \f$\eta^2/4\f$** leaves the zeroth-order reaction mode underdamped, so the regularisation rings.
   *   That does not present as a failure but as noise on the nonlocal field, which is indistinguishable by eye from
   *   the mesh-scale oscillation the gradient enhancement exists to remove.
   *
   * A viscosity of zero admits no micro-inertia at all: an undamped second-order field rings forever.
   */
  std::vector< double > validatedNonlocalMicroInertia( std::vector< double > microInertia,
                                                       const double*         stateVars ) const
  {
    const std::vector< double > eta = getNonlocalViscosity( stateVars );

    for ( size_t n = 0; n < microInertia.size(); n++ ) {

      if ( !std::isfinite( microInertia[n] ) )
        throw std::invalid_argument( MakeString()
                                     << __PRETTY_FUNCTION__ << ": the micro-inertia of nonlocal variable " << n
                                     << " is not a finite number; a NaN passes every ordering test and would "
                                        "propagate into the lumped inertia." );

      if ( microInertia[n] < 0.0 )
        throw std::invalid_argument( MakeString()
                                     << __PRETTY_FUNCTION__ << ": the micro-inertia of nonlocal variable " << n
                                     << " is negative, which would make the nonlocal field integrate backwards in "
                                        "time." );

      /* The relative slack is what makes the RECOMMENDED value reachable. m_k = eta^2/4 is the best
       * choice, and a deck states both numbers in decimal -- eta = 2e-7 with m_k = 1e-14, say --
       * which rounds to 9.999999999999998e-15 here and would miss a bare `>` by one ulp. Rejecting
       * the one pairing the documentation asks for is not a defensible reading of a bound whose
       * violation is a matter of degree.
       */
      const double admissible = 0.25 * eta[n] * eta[n];
      if ( microInertia[n] > admissible * ( 1.0 + 1e-9 ) )
        throw std::invalid_argument(
          MakeString() << __PRETTY_FUNCTION__ << ": the micro-inertia of nonlocal variable " << n << " is "
                       << microInertia[n] << ", above the largest admissible eta^2/4 = " << admissible
                       << " for a nonlocal viscosity of " << eta[n]
                       << ". Above it the zeroth-order reaction mode is underdamped and the regularisation rings, "
                          "which looks like the mesh-scale oscillation the gradient enhancement exists to remove. "
                          "Take exactly eta^2/4: it is the largest admissible value and hence the best one, since "
                          "the stable increment grows with sqrt(m_k)." );
    }

    return microInertia;
  }
};
