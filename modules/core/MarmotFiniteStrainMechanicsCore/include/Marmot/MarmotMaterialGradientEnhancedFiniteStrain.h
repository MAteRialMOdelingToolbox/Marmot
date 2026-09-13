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
#include "Marmot/MarmotStateHelpers.h"
#include "Marmot/MarmotUtils.h"
#include <Fastor/tensor/Tensor.h>
#include <string>
#include <tuple>

/**
 * @class MarmotMaterialGradientEnhancedFiniteStrain
 * @brief Abstract base class for gradient-enhanced (implicit-gradient) materials in the finite strain regime.
 *
 * The non-micropolar sibling of MarmotMaterialGradientEnhancedMicropolar, and the finite-strain
 * counterpart of MarmotMaterialGeneralGradientEnhancedHypoElastic: the displacement field is coupled to a
 * single scalar nonlocal field @f$ \bar{N} @f$ governed by the additional balance equation
 * @f[
 *   \bar{N} - c\,\nabla^2\bar{N} = L(\boldsymbol{F},\,\bar{N}),
 *   \qquad c = R^2,
 * @f]
 * where @f$ L @f$ is the local driving force reported in ConstitutiveResponse::L and @f$ R @f$ the
 * nonlocal radius reported in ConstitutiveResponse::nonLocalRadius.
 *
 * @note **The nonlocal balance is formulated in the material (reference) configuration.** Consumers take
 * the gradients with respect to @f$ \boldsymbol{X} @f$ and integrate over the reference volume, so
 * @f$ c @f$ is a material constant and is neither pushed forward nor weighted by @f$ J @f$. The stress,
 * in contrast, is the Kirchhoff stress @f$ \boldsymbol{\tau} = J\,\boldsymbol{\sigma} @f$, which is what
 * integrating against the reference volume requires.
 *
 * @note The interface carries a **single** scalar nonlocal field, unlike the small-strain
 * MarmotMaterialGeneralGradientEnhancedHypoElastic, which is templated on the number of nonlocal
 * variables. It also has **no slot for @f$ \partial c/\partial\bar{N} @f$**: a material with
 * damage-dependent interactions (see MarmotDecreasingInteractions) can report the current @f$ c @f$, but
 * not its derivative, so its consistent tangent would be incomplete.
 *
 * @warning **This header is a reconstruction.** The original was written by the author of
 * GradientEnhancedFiniteStrainMaterialPoint / ...Particle / ...DisplacementElement but was never
 * published to any repository. Everything below that those consumers observe -- member names, struct
 * shapes, call signatures, the six-argument ConstitutiveResponse constructor -- is fixed by them. What
 * they do not observe is inferred from the sibling interfaces: the *order* of the four doubles in
 * ConstitutiveResponse (the consumers only ever pass zeros there and read the members by name), the
 * names and semantics of elasticEnergyDensity and dissipation, and the defaulted methods. Should the
 * original resurface, diff it against this file before adopting either.
 */
class MarmotMaterialGradientEnhancedFiniteStrain {

protected:
  const double* materialProperties;  ///< Pointer to the array of material property values.
  const int     nMaterialProperties; ///< Number of material property values.

public:
  const int materialNumber; ///< Unique identifier for this material instance.

  /**
   * @brief Construct a MarmotMaterialGradientEnhancedFiniteStrain.
   * @param[in] matProperties_       Pointer to the array of material property values.
   * @param[in] nMaterialProperties_ Number of material property values.
   * @param[in] materialNumber_      Unique identifier for this material instance.
   */
  MarmotMaterialGradientEnhancedFiniteStrain( const double* matProperties_,
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
   * Instances are owned through base-class pointers (the std::unique_ptr held by the gradient-enhanced
   * elements and material points), so destroying a derived object through a base pointer with a
   * non-virtual destructor would be undefined behaviour.
   */
  virtual ~MarmotMaterialGradientEnhancedFiniteStrain() = default;

  /// Layout of the state variables.
  MarmotStateLayoutDynamic stateLayout;

  /**
   * @struct ConstitutiveResponse
   * @brief Constitutive response of a gradient-enhanced material at a given state.
   * @tparam nDim Number of spatial dimensions (2 or 3).
   */
  template < int nDim >
  struct ConstitutiveResponse {
    Fastor::Tensor< double, nDim, nDim > tau; ///< Kirchhoff stress
    double                               L;   ///< local driving force, the source of the nonlocal
                                              ///< balance; a total, not an increment
    double  nonLocalRadius;                   ///< nonlocal radius @f$ R @f$; @f$ c = R^2 @f$
    double  elasticEnergyDensity;             ///< elastic energy per unit reference volume
    double  dissipation;                      ///< dissipation per unit reference volume
    double* stateVars;                        ///< pointer to the state variables

    /// @brief Default constructor; zeroes everything and leaves stateVars null.
    ConstitutiveResponse()
      : tau( Fastor::Tensor< double, nDim, nDim >( 0.0 ) ),
        L( 0.0 ),
        nonLocalRadius( 0.0 ),
        elasticEnergyDensity( 0.0 ),
        dissipation( 0.0 ),
        stateVars( nullptr )
    {
    }

    /**
     * @brief Constructor for initializing the constitutive response.
     * @param[in] tau_                  Kirchhoff stress.
     * @param[in] L_                    Local driving force.
     * @param[in] nonLocalRadius_       Nonlocal radius.
     * @param[in] elasticEnergyDensity_ Elastic energy density.
     * @param[in] dissipation_          Dissipation per unit reference volume.
     * @param[in] stateVars_            Pointer to the state variables.
     */
    ConstitutiveResponse( const Fastor::Tensor< double, nDim, nDim >& tau_,
                          double                                      L_,
                          double                                      nonLocalRadius_,
                          double                                      elasticEnergyDensity_,
                          double                                      dissipation_,
                          double*                                     stateVars_ )
      : tau( tau_ ),
        L( L_ ),
        nonLocalRadius( nonLocalRadius_ ),
        elasticEnergyDensity( elasticEnergyDensity_ ),
        dissipation( dissipation_ ),
        stateVars( stateVars_ )
    {
    }
  };

  /**
   * @struct AlgorithmicModuli
   * @brief Algorithmic tangent moduli of a gradient-enhanced material.
   * @tparam nDim Number of spatial dimensions (2 or 3).
   *
   * The four blocks of the linearized two-field problem. They are zero-initialized, so a material that
   * legitimately has no coupling in one block may simply leave it alone.
   */
  template < int nDim >
  struct AlgorithmicModuli {
    /// @f$ \partial\boldsymbol{\tau}/\partial\boldsymbol{F} @f$
    Fastor::Tensor< double, nDim, nDim, nDim, nDim > dTau_dF = 0.0;
    /// @f$ \partial\boldsymbol{\tau}/\partial\bar{N} @f$
    Fastor::Tensor< double, nDim, nDim > dTau_dN = 0.0;
    /// @f$ \partial L/\partial\boldsymbol{F} @f$
    Fastor::Tensor< double, nDim, nDim > dL_dF = 0.0;
    /// @f$ \partial L/\partial\bar{N} @f$
    double dL_dN = 0.0;
  };

  /**
   * @struct Deformation
   * @brief Deformation state handed to a gradient-enhanced material.
   * @tparam nDim Number of spatial dimensions (2 or 3).
   */
  template < int nDim >
  struct Deformation {
    Fastor::Tensor< double, nDim, nDim > F; ///< deformation gradient
    double                               N; ///< nonlocal field @f$ \bar{N} @f$; a total, not an increment
  };

  /**
   * @struct TimeIncrement
   * @brief Time and time increment of the current step.
   */
  struct TimeIncrement {
    const double time; ///< time at the beginning of the increment
    const double dT;   ///< size of the time increment
  };

  /**
   * @brief Update the material state.
   * @param[in,out] response      ConstitutiveResponse instance; carries the state variables in and the
   *                              Kirchhoff stress, local driving force and nonlocal radius out.
   * @param[out]    tangents      AlgorithmicModuli instance.
   * @param[in]     deformation   Deformation gradient and nonlocal field.
   * @param[in]     timeIncrement Current time and time increment.
   */
  virtual void computeStress( ConstitutiveResponse< 3 >& response,
                              AlgorithmicModuli< 3 >&    tangents,
                              const Deformation< 3 >&    deformation,
                              const TimeIncrement&       timeIncrement ) const = 0;

  /**
   * @brief Explicit version of computeStress, for explicit time integration.
   * @param[in,out] response      ConstitutiveResponse instance.
   * @param[in]     deformation   Deformation gradient and nonlocal field.
   * @param[in]     timeIncrement Current time and time increment.
   *
   * @note The default implementation calls computeStress() and discards the algorithmic tangent.
   * Derived classes may override it for efficiency.
   */
  virtual void computeStressExplicit( ConstitutiveResponse< 3 >& response,
                                      const Deformation< 3 >&    deformation,
                                      const TimeIncrement&       timeIncrement ) const
  {
    AlgorithmicModuli< 3 > tangents;
    computeStress( response, tangents, deformation, timeIncrement );
  }

  /**
   * @brief Update the material state, accounting for an eigen deformation (e.g. a geostatic stress state).
   * @param[in,out] response         ConstitutiveResponse instance.
   * @param[out]    tangents         AlgorithmicModuli instance; scaled to the eigen deformation.
   * @param[in]     deformation      Deformation gradient and nonlocal field.
   * @param[in]     timeIncrement    Current time and time increment.
   * @param[in]     eigenDeformation Eigen deformation in each spatial direction.
   */
  virtual void computeStress( ConstitutiveResponse< 3 >&                  response,
                              AlgorithmicModuli< 3 >&                     tangents,
                              const Deformation< 3 >&                     deformation,
                              const TimeIncrement&                        timeIncrement,
                              const std::tuple< double, double, double >& eigenDeformation ) const;

  /**
   * @brief Compute the response under plane strain conditions.
   * @param[in,out] response          ConstitutiveResponse instance.
   * @param[out]    algorithmicModuli AlgorithmicModuli instance.
   * @param[in]     deformation       Deformation gradient and nonlocal field.
   * @param[in]     timeIncrement     Current time and time increment.
   *
   * @note Plane strain is the 3D response of a deformation gradient whose out-of-plane entries are
   * those of the identity, so the default implementation forwards to computeStress(). Materials that
   * need to do more may override it.
   */
  virtual void computePlaneStrain( ConstitutiveResponse< 3 >& response,
                                   AlgorithmicModuli< 3 >&    algorithmicModuli,
                                   const Deformation< 3 >&    deformation,
                                   const TimeIncrement&       timeIncrement ) const;

  /**
   * @brief Compute the response under plane strain conditions, accounting for an eigen deformation.
   * @param[in,out] response          ConstitutiveResponse instance.
   * @param[out]    algorithmicModuli AlgorithmicModuli instance.
   * @param[in]     deformation       Deformation gradient and nonlocal field.
   * @param[in]     timeIncrement     Current time and time increment.
   * @param[in]     eigenDeformation  Eigen deformation in each spatial direction.
   */
  virtual void computePlaneStrain( ConstitutiveResponse< 3 >&                  response,
                                   AlgorithmicModuli< 3 >&                     algorithmicModuli,
                                   const Deformation< 3 >&                     deformation,
                                   const TimeIncrement&                        timeIncrement,
                                   const std::tuple< double, double, double >& eigenDeformation ) const;

  /**
   * @brief Explicit version of computePlaneStrain, for explicit time integration.
   * @param[in,out] response      ConstitutiveResponse instance.
   * @param[in]     deformation   Deformation gradient and nonlocal field.
   * @param[in]     timeIncrement Current time and time increment.
   *
   * @note The default implementation calls computePlaneStrain() and discards the algorithmic tangent.
   */
  virtual void computePlaneStrainExplicit( ConstitutiveResponse< 3 >& response,
                                           const Deformation< 3 >&    deformation,
                                           const TimeIncrement&       timeIncrement ) const
  {
    AlgorithmicModuli< 3 > algorithmicModuli;
    computePlaneStrain( response, algorithmicModuli, deformation, timeIncrement );
  }

  /**
   * @brief Find the eigen deformation that corresponds to a given eigen stress.
   * @param[in]     initialGuess Initial guess for the eigen deformation.
   * @param[in]     eigenStress  Target eigen normal stress components.
   * @param[in,out] stateVars    State variable array used during the iteration; restored afterwards.
   * @return Eigen deformation corresponding to the given eigen stress.
   *
   * Used for geostatic stress initialization. The nonlocal field is held at zero throughout, which is
   * what an undamaged initial state means.
   */
  std::tuple< double, double, double > findEigenDeformationForEigenStress(
    const std::tuple< double, double, double >& initialGuess,
    const std::tuple< double, double, double >& eigenStress,
    double*                                     stateVars ) const;

  /**
   * @brief Get a view of a state variable.
   * @param[in] stateName Name of the state variable.
   * @param[in] stateVars Pointer to the state variable array.
   * @return StateView giving access to the requested state variable.
   */
  virtual StateView getStateView( const std::string& stateName, double* stateVars ) const
  {
    return stateLayout.getStateView( stateVars, stateName );
  }

  /**
   * @brief Number of state variables this material requires.
   * @return Total size of the state layout.
   */
  virtual int getNumberOfRequiredStateVars() const { return stateLayout.totalSize(); }

  /**
   * @brief Initialize the state variables.
   * @param[in,out] stateVars  Pointer to the state variable array.
   * @param[in]     nStateVars Number of state variables.
   *
   * @note The default implementation zeroes them. A material whose pristine state is not all-zero
   * (a stored deformation gradient, say) must override this.
   */
  virtual void initializeYourself( double* stateVars, int nStateVars )
  {
    for ( int i = 0; i < nStateVars; ++i )
      stateVars[i] = 0.0;
  }

  /**
   * @brief Mass density in the reference configuration.
   * @param[in] stateVars Pointer to the state variable array.
   * @return Mass density.
   */
  virtual double getDensity( const double* stateVars ) const = 0;

  /**
   * @brief Assign the characteristic element length.
   * @param[in] length Characteristic element length at the considered evaluation point.
   *
   * @note A gradient-enhanced material regularises through its own nonlocal radius and normally has no
   * use for a mesh-derived length, so this does nothing by default. It exists for materials that wrap a
   * small-strain model which reads one.
   */
  virtual void setCharacteristicElementLength( double length ) {}

  /**
   * @brief Nonlocal viscosity @f$ \eta @f$ of the nonlocal balance.
   * @param[in] stateVars Pointer to the state variable array.
   * @return Nonlocal viscosity, the damping coefficient of the damped hyperbolic (or, with zero
   * micro-inertia, parabolic) balance
   * @f[
   *   m_k\,\ddot{\bar{N}} + \eta\,\dot{\bar{N}} + \bar{N} - c\,\nabla^2\bar{N} = L(\boldsymbol{F},\,\bar{N}) .
   * @f]
   *
   * @note Zero by default, which is the quasi-static model. The small-strain
   * MarmotMaterialGeneralGradientEnhancedHypoElastic::getNonlocalViscosity is pure virtual instead, on
   * the grounds that a gradient-enhanced material must answer it for any explicit path to exist -- but
   * this header is a reconstruction (see the class note), and a reconstruction must not impose a
   * breaking requirement on implementors it cannot see. Appended at the very end of the class, after
   * every other virtual, so that a stale-library mismatch is confined to callers of this method and of
   * getNonlocalMicroInertia() below -- the reasoning Marmot #84 applied to MarmotElement.
   */
  virtual double getNonlocalViscosity( const double* stateVars ) const { return 0.0; }

  /**
   * @brief Nonlocal micro-inertia @f$ m_k @f$ of the nonlocal balance.
   * @param[in] stateVars Pointer to the state variable array.
   * @return Nonlocal micro-inertia, in units of [time]^2.
   *
   * @note The default is zero, the quasi-static model in which the nonlocal field is first order in
   * time and integrated by its viscosity alone -- the behaviour of every material that predates this
   * method. A derived class opts in to the damped hyperbolic balance by overriding this, mirroring
   * MarmotMaterialGeneralGradientEnhancedHypoElastic::getNonlocalMicroInertia.
   */
  virtual double getNonlocalMicroInertia( const double* stateVars ) const { return 0.0; }
};
