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

#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialGeneralGradientEnhancedHypoElastic.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrain.h"
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

namespace Marmot::Materials {

  /**
   * @brief Number of nonlocal variables of a small-strain gradient-enhanced material type.
   * @tparam nNonlocalVariables Deduced from the base class Marmot::Materials::MaterialGradientEnhancedHypoElastic.
   * @return @p nNonlocalVariables.
   * @details Called unevaluated, on a null pointer of the material type, to read the template argument of its base.
   */
  template < int nNonlocalVariables >
  constexpr int nNonlocalVariablesOf( const MaterialGradientEnhancedHypoElastic< nNonlocalVariables >* )
  {
    return nNonlocalVariables;
  }

  /**
   * @class GradientEnhancedHypoElasticWrapper
   * @brief Common base of the decorators presenting a small-strain gradient-enhanced material as a finite-strain one.
   *
   * @tparam BaseMaterialType A concrete class derived from Marmot::Materials::MaterialGradientEnhancedHypoElastic<n>;
   *                          the wrapper has the same number n of nonlocal variables.
   *
   * What the decorators share: they own the wrapped material, keep its state variables in the slot `materialstate`
   * of their own layout, forward its density, nonlocal viscosity and micro-inertia, resolve state names they do not
   * know themselves in the wrapped material, and hand it the nonlocal fields together with their increment. The
   * finite-strain interface passes the TOTAL nonlocal fields, while the wrapped material wants the fields **and
   * their increment**; the increment is formed against the value of the last call, which each decorator keeps in a
   * slot of its own layout (baseMaterialIncrement() below). The kinematics -- which strain increment the wrapped
   * material is driven by, and which stress it returns -- is the derived class' own:
   * GradientEnhancedHughesWingetWrapper (an objective rate, for a material incremental in stress) and
   * GradientEnhancedCorotationalWrapper (the logarithmic strain in the frame of the polar rotation, for any material).
   *
   * Each derived class builds its layout in its constructor: its own slots, a slot for the nonlocal fields of
   * the last call, and baseMaterialSlot, of size baseMaterial->getNumberOfRequiredStateVars().
   */
  template < typename BaseMaterialType >
  class GradientEnhancedHypoElasticWrapper : public MaterialGradientEnhancedFiniteStrain< nNonlocalVariablesOf(
                                               static_cast< BaseMaterialType* >( nullptr ) ) > {

  public:
    /// The number of nonlocal variables, that of the wrapped material.
    static constexpr int nNonlocal = nNonlocalVariablesOf( static_cast< BaseMaterialType* >( nullptr ) );

    /// The finite-strain interface presented.
    using Interface = MaterialGradientEnhancedFiniteStrain< nNonlocal >;

    using typename Interface::NonlocalVector;

  protected:
    /// The small-strain interface of the wrapped material.
    using BaseMaterial = MaterialGradientEnhancedHypoElastic< nNonlocal >;

    /// The wrapped small-strain material, held as the base interface.
    std::unique_ptr< BaseMaterial > baseMaterial;

    /// Name of the state layout slot holding the state variables of the wrapped material.
    static constexpr const char* baseMaterialSlot = "materialstate";

  public:
    /**
     * @brief Construct the wrapper and the wrapped material.
     * @param[in] matProperties_       Material property values; forwarded **unchanged** to the wrapped material.
     * @param[in] nMaterialProperties_ Number of material property values.
     * @param[in] materialNumber_      Unique identifier for this material instance.
     */
    GradientEnhancedHypoElasticWrapper( const double* matProperties_, int nMaterialProperties_, int materialNumber_ )
      : Interface( matProperties_, nMaterialProperties_, materialNumber_ ),
        baseMaterial( std::make_unique< BaseMaterialType >( matProperties_, nMaterialProperties_, materialNumber_ ) )
    {
    }

    /**
     * @brief Resolve a state variable by name: a slot of the wrapper, or else a state of the wrapped material.
     * @param[in] stateName Name of the state variable.
     * @param[in] stateVars Pointer to the state variable array of this wrapper.
     * @return A view of the requested state variable.
     */
    Marmot::Utilities::StateView getStateView( const std::string& stateName, double* stateVars ) const override
    {
      if ( this->stateLayout.contains( stateName ) )
        return this->stateLayout.getStateView( stateVars, stateName );

      return baseMaterial->getStateView( stateName, baseMaterialState( stateVars ) );
    }

    /**
     * @brief Mass density of the wrapped material.
     * @param[in] stateVars Pointer to the state variable array of this wrapper.
     * @return Mass density.
     */
    double getDensity( const double* stateVars ) const override
    {
      return baseMaterial->getDensity( baseMaterialState( const_cast< double* >( stateVars ) ) );
    }

    /**
     * @brief Nonlocal viscosities of the wrapped small-strain material.
     * @param[in] stateVars Pointer to the state variable array of this wrapper.
     * @return The wrapped material's own getNonlocalViscosity(), one entry per nonlocal variable.
     */
    NonlocalVector getNonlocalViscosity( const double* stateVars ) const override
    {
      return toNonlocalVector(
        baseMaterial->getNonlocalViscosity( baseMaterialState( const_cast< double* >( stateVars ) ) ) );
    }

    /**
     * @brief Nonlocal micro-inertias of the wrapped small-strain material.
     * @param[in] stateVars Pointer to the state variable array of this wrapper.
     * @return The wrapped material's own getNonlocalMicroInertia(), one entry per nonlocal variable.
     */
    NonlocalVector getNonlocalMicroInertia( const double* stateVars ) const override
    {
      return toNonlocalVector(
        baseMaterial->getNonlocalMicroInertia( baseMaterialState( const_cast< double* >( stateVars ) ) ) );
    }

  protected:
    /**
     * @brief Convert an Eigen 3x3 matrix into a Fastor tensor.
     * @param[in] m The matrix.
     * @return The same values as a Fastor tensor.
     */
    static Marmot::Math::Tensors::FastorTensors::Tensor33d fromEigen( const Eigen::Matrix3d& m )
    {
      Marmot::Math::Tensors::FastorTensors::Tensor33d t;
      Marmot::Math::Tensors::mapEigenToFastor( t ) = m;
      return t;
    }

    /**
     * @brief Convert a Fastor 3x3 tensor into an Eigen matrix.
     * @param[in] t The tensor.
     * @return The same values as an Eigen matrix.
     */
    static Eigen::Matrix3d toEigen( const Marmot::Math::Tensors::FastorTensors::Tensor33d& t )
    {
      return Eigen::Matrix3d( Marmot::Math::Tensors::mapEigenToFastor( t ) );
    }

    /**
     * @brief The state variables of the wrapped material inside those of this wrapper.
     * @param[in] stateVars Pointer to the state variable array of this wrapper.
     * @return Pointer to the slot baseMaterialSlot.
     */
    double* baseMaterialState( double* stateVars ) const
    {
      return this->stateLayout.getPtr( stateVars, baseMaterialSlot );
    }

    /**
     * @brief Initialise the state variables of the wrapped material.
     * @param[in,out] stateVars Pointer to the state variable array of this wrapper.
     */
    void initializeBaseMaterial( double* stateVars )
    {
      baseMaterial->initializeYourself( baseMaterialState( stateVars ), baseMaterial->getNumberOfRequiredStateVars() );
    }

    /**
     * @brief The increment the wrapped material is handed for a total nonlocal field, and the field it is formed
     * against.
     * @param[in] N  Total nonlocal fields of this call.
     * @param[in] Nn Nonlocal fields of the last call, as kept in the wrapper's own slot.
     * @param[in] dStrain Strain increment, in Voigt notation with engineering shear.
     * @param[in] timeIncrement Time and time increment of the call.
     * @return The increment for the wrapped material: @p dStrain, @f$ K = \bar{N} @f$,
     *         @f$ \Delta K = \bar{N} - \bar{N}_n @f$ and the time.
     */
    static typename BaseMaterial::increment baseMaterialIncrement(
      const NonlocalVector&                    N,
      const NonlocalVector&                    Nn,
      const Marmot::Math::Vector6d&            dStrain,
      const typename Interface::TimeIncrement& timeIncrement )
    {
      typename BaseMaterial::increment inc;
      inc.dStrain = dStrain;
      inc.K       = N;
      inc.dK      = N - Nn;
      inc.time    = timeIncrement.time;
      inc.dT      = timeIncrement.dT;
      return inc;
    }

    /**
     * @brief Copy a per-nonlocal-variable std::vector, as the small-strain interface reports one, into a vector of
     * the finite-strain interface.
     * @param[in] values One value per nonlocal variable.
     * @return The same values.
     */
    static NonlocalVector toNonlocalVector( const std::vector< double >& values )
    {
      NonlocalVector v;
      for ( int K = 0; K < nNonlocal; K++ )
        v( K ) = values[K];
      return v;
    }
  };

} // namespace Marmot::Materials
