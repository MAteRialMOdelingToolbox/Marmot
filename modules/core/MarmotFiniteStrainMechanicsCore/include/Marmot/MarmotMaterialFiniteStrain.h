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
#include <Fastor/tensor/Tensor.h>
#include <algorithm>
#include <cmath>
#include <vector>

/**
 * @class MarmotMaterialFiniteStrain
 * @brief Abstract base class for mechanical materials in the finite strain regime.
 *
 * Derived classes implement computeStress() to provide the constitutive response
 * (Kirchhoff stress, density, elastic energy density) and the algorithmic tangent.
 */
class MarmotMaterialFiniteStrain {

protected:
  const double* materialProperties;  ///< Pointer to the array of material property values.
  const int     nMaterialProperties; ///< Number of material property values.

public:
  const int materialNumber; ///< Unique identifier for this material instance.

  /**
   * @brief Construct a MarmotMaterialFiniteStrain.
   * @param[in] matProperties_       Pointer to the array of material property values.
   * @param[in] nMaterialProperties_ Number of material property values.
   * @param[in] materialNumber_      Unique identifier for this material instance.
   */
  MarmotMaterialFiniteStrain( const double* matProperties_, int nMaterialProperties_, int materialNumber_ )
    : materialProperties( matProperties_ ),
      nMaterialProperties( nMaterialProperties_ ),
      materialNumber( materialNumber_ )
  {
  }

  /// Default destructor
  virtual ~MarmotMaterialFiniteStrain() = default;

  /// Layout of the state variables
  MarmotStateLayoutDynamic stateLayout;

  /**
   * @struct ConstitutiveResponse
   * @brief Constitutive response of a material at given state.
   * @tparam nDim Number of spatial dimensions (2 or 3).
   *
   *  Contains stress, density and elastic energy density.
   */
  template < int nDim >
  struct ConstitutiveResponse {
    Fastor::Tensor< double, nDim, nDim > tau;                  ///< Kirchhoff stress
    double                               elasticEnergyDensity; ///< elastic energy per unit volume
    double                               dissipation;          ///< dissipation per unit volume
    double*                              stateVars;            ///< pointer to state variables

    /**
     * @brief Default constructor.
     * Initializes stress to zero, energy and dissipation to zero, and stateVars to nullptr.
     */
    ConstitutiveResponse()
      : tau( Fastor::Tensor< double, nDim, nDim >( 0.0 ) ),
        elasticEnergyDensity( 0.0 ),
        dissipation( 0.0 ),
        stateVars( nullptr )
    {
    }
    /**
     * @brief Constructor for initializing the constitutive response.
     * @param tau_ Kirchhoff stress tensor.
     * @param elasticEnergyDensity_ Elastic energy density.
     * @param dissipation_ Dissipation per unit volume.
     * @param stateVars_ Pointer to state variables.
     */
    ConstitutiveResponse( const Fastor::Tensor< double, nDim, nDim >& tau_,
                          double                                      elasticEnergyDensity_,
                          double                                      dissipation_,
                          double*                                     stateVars_ )
      : tau( tau_ ), elasticEnergyDensity( elasticEnergyDensity_ ), dissipation( dissipation_ ), stateVars( stateVars_ )
    {
    }
  };

  /**
   * @struct AlgorithmicModuli
   * @brief Algorithmic tangent moduli of a material.
   * @tparam nDim Number of spatial dimensions (2 or 3).
   *
   * Contains the algorithmic tangent moduli \f$\frac{\partial \boldsymbol{\tau}}{\partial \boldsymbol{F}}\f$
   * with respect to the deformation gradient \f$\boldsymbol{F}\f$.
   * */
  template < int nDim >
  struct AlgorithmicModuli {
    Fastor::Tensor< double, nDim, nDim, nDim, nDim > dTau_dF; ///< tangent operator w.r.t. deformation gradient
  };

  /**
   * @struct Deformation
   * @brief Represents the deformation state of a material.
   *
   * This struct holds the deformation gradient \f$\boldsymbol{F}\f$ which describes the local deformation of a
   * material.
   */
  template < int nDim >
  struct Deformation {
    Fastor::Tensor< double, nDim, nDim > F; ///< deformation gradient
  };

  /**
   * @struct TimeIncrement
   * @brief Represents a time increment in a simulation.
   *
   * This struct holds information about the current time
   * and the time step size (dT).
   */
  struct TimeIncrement {
    const double time; ///< time at the beginning of the increment
    const double dT;   ///< size of the time increment
  };

  /**
   * @brief Updates the material state.
   * @param[inout] response ConstitutiveResponse instance
   * @param[out] tangents AlgorithmicModuli instance
   * @param[in] deformation Deformation instance
   * @param[in] timeIncrement TimeIncrement instance
   *
   * Computes the Kirchhoff \f$\boldsymbol{\tau}\f$ stress
   * from the deformation gradient \f$\boldsymbol{F}\f$ and the
   * time increment \f$\Delta t\f$.
   * It further updates the mass density \f$\rho\f$ and the elastic energy density.
   * Additionally, computes the algorithmic tangent moduli
   * \f$\frac{\partial \boldsymbol{\tau}}{\partial \boldsymbol{F}}\f$.
   *
   * */
  virtual void computeStress( ConstitutiveResponse< 3 >& response,
                              AlgorithmicModuli< 3 >&    tangents,
                              const Deformation< 3 >&    deformation,
                              const TimeIncrement&       timeIncrement ) const = 0;

  /**
   * @brief Explicit version of computeStress for use in explicit time integration schemes.
   * @param[inout] response ConstitutiveResponse instance
   * @param[in] deformation Deformation instance
   * @param[in] timeIncrement TimeIncrement instance
   *
   * @note The default implementation calls computeStress and ignores the algorithmic tangent.
   * @note Derived classes may override this method for efficiency reasons.
   */
  virtual void computeStressExplicit( ConstitutiveResponse< 3 >& response,
                                      const Deformation< 3 >&    deformation,
                                      const TimeIncrement&       timeIncrement ) const
  {
    AlgorithmicModuli< 3 > tangents;
    computeStress( response, tangents, deformation, timeIncrement );
  }

  /**
   * @brief Computes the Kirchhoff stress given the deformation, time increment, and eigen deformation.
   * @param[inout] response ConstitutiveResponse instance
   * @param[out] tangents AlgorithmicModuli instance
   * @param[in] deformation Deformation instance
   * @param[in] timeIncrement TimeIncrement instance
   * @param[in] eigenDeformation Tuple representing eigen deformation in each spatial direction.
   *
   */
  virtual void computeStress( ConstitutiveResponse< 3 >&                  response,
                              AlgorithmicModuli< 3 >&                     tangents,
                              const Deformation< 3 >&                     deformation,
                              const TimeIncrement&                        timeIncrement,
                              const std::tuple< double, double, double >& eigenDeformation ) const;

  /**
   * @brief Compute stress under plane strain conditions.
   * @param[inout] response ConstitutiveResponse instance
   * @param[out] algorithmicModuli AlgorithmicModuli instance
   * @param[in] deformation Deformation instance
   * @param[in] timeIncrement TimeIncrement instance
   *
   * It uses the general 3D computeStress function for a plane strain Deformation.
   * The algorithmic tangent is modified according to plane strain conditions.
   */
  virtual void computePlaneStrain( ConstitutiveResponse< 3 >& response,
                                   AlgorithmicModuli< 3 >&    algorithmicModuli,
                                   const Deformation< 3 >&    deformation,
                                   const TimeIncrement&       timeIncrement ) const;

  /**
   * @brief Explicit version of computePlaneStrain for use in explicit time integration schemes.
   * @note The default implementation calls computePlaneStrain and ignores the algorithmic tangent.
   */
  virtual void computePlaneStrainExplicit( ConstitutiveResponse< 3 >& response,
                                           const Deformation< 3 >&    deformation,
                                           const TimeIncrement&       timeIncrement ) const
  {
    AlgorithmicModuli< 3 > algorithmicModuli;
    computePlaneStrain( response, algorithmicModuli, deformation, timeIncrement );
  }

  /**
   * @brief Compute stress under plane strain conditions with eigen deformation.
   * @param[inout] response ConstitutiveResponse instance
   * @param[out] algorithmicModuli AlgorithmicModuli instance
   * @param[in] deformation Deformation instance
   * @param[in] timeIncrement TimeIncrement instance
   * @param[in] eigenDeformation Tuple representing eigen deformation in each spatial direction.
   *
   * It uses the general 3D computeStress function for a plane strain Deformation.
   * The algorithmic tangent is modified according to plane strain conditions.
   */
  virtual void computePlaneStrain( ConstitutiveResponse< 3 >&                  response,
                                   AlgorithmicModuli< 3 >&                     algorithmicModuli,
                                   const Deformation< 3 >&                     deformation,
                                   const TimeIncrement&                        timeIncrement,
                                   const std::tuple< double, double, double >& eigenDeformation ) const;

  /**
   * @brief Explicit version of computePlaneStrain with eigen deformation for use in explicit time integration schemes.
   * @note The default implementation calls computePlaneStrain and ignores the algorithmic tangent.
   */
  virtual void computePlaneStrainExplicit( ConstitutiveResponse< 3 >&                  response,
                                           const Deformation< 3 >&                     deformation,
                                           const TimeIncrement&                        timeIncrement,
                                           const std::tuple< double, double, double >& eigenDeformation ) const
  {
    AlgorithmicModuli< 3 > algorithmicModuli;
    computePlaneStrain( response, algorithmicModuli, deformation, timeIncrement, eigenDeformation );
  }

  /**
   * @brief Compute stress under plane stress conditions.
   * @param[inout] response ConstitutiveResponse instance
   * @param[out] algorithmicModuli AlgorithmicModuli instance
   * @param[in] deformation Deformation instance
   * @param[in] timeIncrement TimeIncrement instance
   *
   * It uses the general 3D computeStress function and iteratively finds the out-of-plane deformation.
   * The algorithmic tangent is modified according to plane stress conditions.
   */
  virtual void computePlaneStress( ConstitutiveResponse< 2 >& response,
                                   AlgorithmicModuli< 2 >&    algorithmicModuli,
                                   const Deformation< 2 >&    deformation,
                                   const TimeIncrement&       timeIncrement ) const;

  /**
   * @brief Explicit version of computePlaneStress for use in explicit time integration schemes.
   * @note The default implementation calls computePlaneStress and ignores the algorithmic tangent.
   */
  virtual void computePlaneStressExplicit( ConstitutiveResponse< 2 >& response,
                                           const Deformation< 2 >&    deformation,
                                           const TimeIncrement&       timeIncrement ) const
  {
    AlgorithmicModuli< 2 > algorithmicModuli;
    computePlaneStress( response, algorithmicModuli, deformation, timeIncrement );
  }

  /**
   * @brief Find the eigen deformation that corresponds to a given eigen stress.
   * @param initialGuess Initial guess for the eigen deformation.
   * @param eigenStress Target eigen stress.
   * @param stateVars Pointer to the state variable array used during iteration.
   * @return Eigen deformation that corresponds to the given eigen stress.
   *
   * This function iteratively finds the eigen deformation that corresponds to a given eigen stress.
   * This is used e.g. for geostatic stress initialization.
   */
  std::tuple< double, double, double > findEigenDeformationForEigenStress(
    const std::tuple< double, double, double >& initialGuess,
    const std::tuple< double, double, double >& eigenStress,
    double*                                     stateVars ) const;

  /**
   * @brief Get a view to the state variables.
   * @param stateName Name of the state variable
   * @param stateVars Pointer to the state variable array
   * @return StatView to access the state variable
   */
  virtual StateView getStateView( const std::string& stateName, double* stateVars ) const
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
   * @brief Get the maximum wave speed of the material for a given deformation gradient.
   * @param[in] stateVars Pointer to the state variable array
   * @param[in] deformationGradient Current deformation gradient used as perturbation reference
   * @return Maximum wave speed
   * @details The default implementation computes diagonal Voigt tangent entries by centered
   * finite differences and returns `sqrt(max(C_ii) / rho)`.
   */
  virtual double getMaximumWaveSpeed( const double*                         stateVars,
                                      const Fastor::Tensor< double, 3, 3 >& deformationGradient ) const
  {
    constexpr double dEps = 1e-6;

    const int nStateVars = getNumberOfRequiredStateVars();

    double     maxStiffnessDiagonal = 0.0;
    const auto timeIncrement        = TimeIncrement{ 0.0, 1.0 };

    for ( int i = 0; i < 3; ++i ) {
      std::vector< double > stateVarsPlus( nStateVars, 0.0 );
      std::vector< double > stateVarsMinus( nStateVars, 0.0 );
      if ( stateVars != nullptr && nStateVars > 0 ) {
        std::copy_n( stateVars, nStateVars, stateVarsPlus.begin() );
        std::copy_n( stateVars, nStateVars, stateVarsMinus.begin() );
      }

      auto FPlus  = deformationGradient;
      auto FMinus = deformationGradient;
      FPlus( i, i ) *= ( 1.0 + dEps );
      FMinus( i, i ) *= ( 1.0 - dEps );

      Deformation< 3 > deformationPlus{ FPlus };
      Deformation< 3 > deformationMinus{ FMinus };

      ConstitutiveResponse< 3 > responsePlus{ Fastor::Tensor< double, 3, 3 >( 0.0 ), 0.0, 0.0, stateVarsPlus.data() };
      ConstitutiveResponse< 3 > responseMinus{ Fastor::Tensor< double, 3, 3 >( 0.0 ), 0.0, 0.0, stateVarsMinus.data() };

      computeStressExplicit( responsePlus, deformationPlus, timeIncrement );
      computeStressExplicit( responseMinus, deformationMinus, timeIncrement );

      const double dLogStrain = std::log( 1.0 + dEps ) - std::log( 1.0 - dEps );
      const double Cii        = ( responsePlus.tau( i, i ) - responseMinus.tau( i, i ) ) / dLogStrain;
      maxStiffnessDiagonal    = std::max( maxStiffnessDiagonal, Cii );
    }

    const double density = getDensity( stateVars );
    return density > 0.0 ? std::sqrt( std::max( 0.0, maxStiffnessDiagonal ) / density ) : 0.0;
  }

  /**
   * @brief Set the characteristic element length at the considered evaluation point.
   *
   * Needed by materials whose softening behaviour is regularised via a mesh-adjusted softening modulus.
   * The default implementation does nothing, so materials that do not depend on a length are unaffected.
   *
   * @param[in] length Characteristic element length.
   */
  virtual void setCharacteristicElementLength( double length ) {}

  /**
   * @brief Get the mass density of the material.
   * @param[in] stateVars Pointer to the state variable array
   * @return Mass density
   */
  virtual double getDensity( const double* stateVars ) const = 0;
};
