/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck.
 *
 * Thomas Mader thomas.mader@boku.ac.at
 *
 * This file is part of the MAteRialMOdellingToolbox (marmot).
 * LGPL v2.1+, see LICENSE.md at the top level directory of marmot.
 * ---------------------------------------------------------------------
 */
#pragma once
#include "Marmot/MarmotStateHelpers.h"
#include <Fastor/tensor/Tensor.h>
#include <algorithm>
#include <cmath>
#include <tuple>
#include <vector>

/**
 * @class MarmotMaterialGradientEnhancedFiniteStrain
 * @brief Abstract base class for gradient-enhanced (implicit-gradient / nonlocal)
 *        mechanical materials in the finite strain regime, WITHOUT a micropolar
 *        continuum.
 *
 * This is the non-micropolar sibling of @ref MarmotMaterialGradientEnhancedMicropolar:
 * it couples the standard displacement field (deformation gradient @f$\boldsymbol{F}@f$)
 * to a single scalar nonlocal field @f$\bar{N}@f$ used to regularise softening damage.
 * The material returns the Kirchhoff stress @f$\boldsymbol{\tau}@f$, a scalar local
 * driving force @f$L@f$ (the source term of the Helmholtz-type nonlocal averaging
 * equation @f$\bar{N} - c\,\nabla^2\bar{N} = L@f$, with @f$c = @f$ nonLocalRadius@f$^2@f$),
 * and the algorithmic tangents
 * @f$\partial\boldsymbol{\tau}/\partial\boldsymbol{F}@f$,
 * @f$\partial\boldsymbol{\tau}/\partial\bar{N}@f$,
 * @f$\partial L/\partial\boldsymbol{F}@f$ and
 * @f$\partial L/\partial\bar{N}@f$.
 */
class MarmotMaterialGradientEnhancedFiniteStrain {

protected:
  const double* materialProperties;  ///< Pointer to the array of material property values.
  const int     nMaterialProperties; ///< Number of material property values.

public:
  const int materialNumber; ///< Unique identifier for this material instance.

  MarmotMaterialGradientEnhancedFiniteStrain( const double* matProperties_,
                                              int           nMaterialProperties_,
                                              int           materialNumber_ )
    : materialProperties( matProperties_ ),
      nMaterialProperties( nMaterialProperties_ ),
      materialNumber( materialNumber_ )
  {
  }

  virtual ~MarmotMaterialGradientEnhancedFiniteStrain() = default;

  /// Layout of the state variables
  MarmotStateLayoutDynamic stateLayout;

  /**
   * @struct ConstitutiveResponse
   * @brief Constitutive response of a gradient-enhanced finite-strain material.
   */
  template < int nDim >
  struct ConstitutiveResponse {
    Fastor::Tensor< double, nDim, nDim > tau;                  ///< Kirchhoff stress
    double                               L;                    ///< local damage driving force (nonlocal source term)
    double                               nonLocalRadius;       ///< nonlocal length (c = nonLocalRadius^2)
    double                               elasticEnergyDensity; ///< elastic energy per unit volume
    double                               dissipation;          ///< dissipation per unit volume
    double*                              stateVars;            ///< pointer to state variables

    ConstitutiveResponse()
      : tau( Fastor::Tensor< double, nDim, nDim >( 0.0 ) ),
        L( 0.0 ),
        nonLocalRadius( 0.0 ),
        elasticEnergyDensity( 0.0 ),
        dissipation( 0.0 ),
        stateVars( nullptr )
    {
    }

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
   * @brief Algorithmic tangent moduli of a gradient-enhanced finite-strain material.
   */
  template < int nDim >
  struct AlgorithmicModuli {
    Fastor::Tensor< double, nDim, nDim, nDim, nDim > dTau_dF; ///< d tau / d F
    Fastor::Tensor< double, nDim, nDim >             dTau_dN; ///< d tau / d Nbar
    Fastor::Tensor< double, nDim, nDim >             dL_dF;   ///< d L   / d F
    double                                           dL_dN;   ///< d L   / d Nbar
  };

  /**
   * @struct Deformation
   * @brief Deformation state fed to a gradient-enhanced finite-strain material:
   *        the deformation gradient plus the nonlocal field value at the point.
   */
  template < int nDim >
  struct Deformation {
    Fastor::Tensor< double, nDim, nDim > F; ///< deformation gradient
    double                               N; ///< nonlocal field value at the point
  };

  struct TimeIncrement {
    const double time; ///< time at the beginning of the increment
    const double dT;   ///< size of the time increment
  };

  /**
   * @brief Update the material state (3D).
   *
   * Computes the Kirchhoff stress, the local driving force L, the nonlocal
   * radius and the algorithmic tangents from the deformation gradient F and
   * the nonlocal field value N. A failed local iteration must be signalled by
   * throwing Marmot::StressUpdateFailed.
   */
  virtual void computeStress( ConstitutiveResponse< 3 >& response,
                              AlgorithmicModuli< 3 >&    tangents,
                              const Deformation< 3 >&    deformation,
                              const TimeIncrement&       timeIncrement ) const = 0;

  /**
   * @brief computeStress accounting for an eigen deformation (e.g. geostatic states).
   */
  virtual void computeStress( ConstitutiveResponse< 3 >&                  response,
                              AlgorithmicModuli< 3 >&                     tangents,
                              const Deformation< 3 >&                     deformation,
                              const TimeIncrement&                        timeIncrement,
                              const std::tuple< double, double, double >& eigenDeformation ) const;

  /**
   * @brief Compute stress under plane strain conditions.
   *
   * Default implementation forwards to the 3D computeStress (the caller is
   * expected to have set F(2,2)=1).
   */
  virtual void computePlaneStrain( ConstitutiveResponse< 3 >& response,
                                   AlgorithmicModuli< 3 >&    algorithmicModuli,
                                   const Deformation< 3 >&    deformation,
                                   const TimeIncrement&       timeIncrement ) const;

  virtual void computePlaneStrain( ConstitutiveResponse< 3 >&                  response,
                                   AlgorithmicModuli< 3 >&                     algorithmicModuli,
                                   const Deformation< 3 >&                     deformation,
                                   const TimeIncrement&                        timeIncrement,
                                   const std::tuple< double, double, double >& eigenDeformation ) const;

  /**
   * @brief Find the eigen deformation that corresponds to a given eigen stress
   *        (used for geostatic stress initialization). The nonlocal field is
   *        held at zero during the search.
   */
  std::tuple< double, double, double > findEigenDeformationForEigenStress(
    const std::tuple< double, double, double >& initialGuess,
    const std::tuple< double, double, double >& eigenStress,
    double*                                     stateVars ) const;

  StateView getStateView( const std::string& stateName, double* stateVars ) const
  {
    return stateLayout.getStateView( stateVars, stateName );
  }

  int getNumberOfRequiredStateVars() const { return stateLayout.totalSize(); }

  virtual void initializeYourself( double* stateVars, int nStateVars )
  {
    for ( int i = 0; i < nStateVars; ++i )
      stateVars[i] = 0.0;
  }

  virtual double getDensity( const double* stateVars ) const { return 0.0; }
};
