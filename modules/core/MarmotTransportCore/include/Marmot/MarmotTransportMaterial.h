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
#include <Eigen/Core>
#include <string>

namespace Marmot::Materials {

  /**
   * @class Marmot::Materials::TransportMaterial
   * @brief Abstract base class for the constitutive laws of (coupled) scalar transport problems: heat conduction,
   * moisture diffusion, coupled heat and moisture transport.
   *
   * A transport material describes \f$n_f\f$ scalar fields \f$\phi_a\f$ (e.g. the temperature, or the temperature
   * and the pore relative humidity) through three constitutive quantities per field:
   *
   * - the **content** \f$m_a(\boldsymbol{\phi})\f$, the conserved quantity per unit volume whose rate enters the
   *   balance (the enthalpy for heat, the water content for moisture),
   * - the **flux** \f$\mathbf{q}_a(\boldsymbol{\phi}, \nabla\boldsymbol{\phi})\f$, the physical flux (e.g.
   *   \f$\mathbf{q} = -\mathbf{k}\,\nabla T\f$),
   * - the **source** \f$r_a(\boldsymbol{\phi})\f$, a volumetric production rate of the material itself (e.g. the
   *   heat of hydration, or a self-desiccation sink); loads applied from outside are body loads of the element.
   *
   * The balance of each field reads \f$\dot m_a + \nabla\cdot\mathbf{q}_a = r_a\f$. Elements integrate it in time
   * with the backward Euler method in its conservative form: the material returns the content INCREMENT
   * \f$\Delta m_a = m_a(\boldsymbol{\phi}_{n+1}) - m_a(\boldsymbol{\phi}_n)\f$ over the time increment, which keeps
   * the balance exact for a field dependent capacity (the "enthalpy method"), and leaves path dependent contents
   * (e.g. hysteretic sorption) to the material and its state variables.
   *
   * Gradients and fluxes are always three-dimensional, also in two-dimensional elements, which pad the gradient with
   * zeros and use the in-plane components of the flux only. Field \f$a\f$ occupies the components \f$3a, 3a+1,
   * 3a+2\f$ of the gradient and flux vectors.
   *
   * The interface has its own small structs instead of the state3D/timeInfo of the mechanical materials, because
   * none of its quantities is a stress. The sizes are dynamic with a fixed maximum (maxNumberOfFields), so that
   * they live on the stack.
   */
  class TransportMaterial {

  protected:
    const double* materialProperties;  ///< Pointer to the array of material properties
    const int     nMaterialProperties; ///< Number of material properties

  public:
    const int materialNumber; ///< Integer identifier for this material instance

    /// Maximum number of transported fields a material may describe
    static constexpr int maxNumberOfFields = 4;

    /// One value per field
    using FieldVector = Eigen::Matrix< double, Eigen::Dynamic, 1, Eigen::ColMajor, maxNumberOfFields, 1 >;
    /// Derivative of one value per field with respect to the fields
    using FieldMatrix = Eigen::
      Matrix< double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, maxNumberOfFields, maxNumberOfFields >;
    /// Three components per field (gradients and fluxes), field after field
    using GradientVector = Eigen::Matrix< double, Eigen::Dynamic, 1, Eigen::ColMajor, 3 * maxNumberOfFields, 1 >;
    /// Derivative of the fluxes with respect to the fields
    using FluxFieldMatrix = Eigen::
      Matrix< double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, 3 * maxNumberOfFields, maxNumberOfFields >;
    /// Derivative of the fluxes with respect to the gradients
    using FluxGradientMatrix = Eigen::
      Matrix< double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, 3 * maxNumberOfFields, 3 * maxNumberOfFields >;

    /**
     * @brief Construct the material with a given set of properties and an identifier.
     * @param[in] materialProperties_  Pointer to the array of material properties.
     * @param[in] nMaterialProperties_ Number of entries in @p materialProperties_.
     * @param[in] materialNumber_      Integer identifying this material instance.
     */
    TransportMaterial( const double* materialProperties_, int nMaterialProperties_, int materialNumber_ )
      : materialProperties( materialProperties_ ),
        nMaterialProperties( nMaterialProperties_ ),
        materialNumber( materialNumber_ )
    {
    }

    /// Default destructor
    virtual ~TransportMaterial() = default;

    /// Layout of the state variables
    Marmot::Utilities::StateLayoutDynamic stateLayout;

    /// The fields at a material point
    struct FieldState {
      FieldVector    field;    ///< \f$\boldsymbol{\phi}_{n+1}\f$, the fields at the end of the increment
      FieldVector    fieldOld; ///< \f$\boldsymbol{\phi}_n\f$, the fields at the beginning of the increment
      GradientVector gradient; ///< \f$\nabla\boldsymbol{\phi}_{n+1}\f$, three components per field
    };

    /// The response of the material at the end of the increment
    struct Response {
      FieldVector    contentIncrement; ///< \f$\Delta m_a = m_a(\boldsymbol{\phi}_{n+1}) - m_a(\boldsymbol{\phi}_n)\f$
      GradientVector flux;             ///< \f$\mathbf{q}_a\f$ at the end of the increment, three components per field
      FieldVector    source;           ///< \f$r_a\f$, volumetric production rate at the end of the increment
      double*        stateVars;        ///< state variables: values at the beginning of the increment on entry
    };

    /// The algorithmic tangent of the response with respect to the fields and gradients at the end of the increment
    struct AlgorithmicTangent {
      FieldMatrix        dContentIncrement_dField; ///< \f$\partial\Delta m_a / \partial\phi_b\f$
      FluxFieldMatrix    dFlux_dField;             ///< \f$\partial\mathbf{q}_a / \partial\phi_b\f$
      FluxGradientMatrix dFlux_dGradient;          ///< \f$\partial\mathbf{q}_a / \partial\nabla\phi_b\f$
      FieldMatrix        dSource_dField;           ///< \f$\partial r_a / \partial\phi_b\f$
    };

    /// Time at the beginning of the increment and size of the time increment
    struct TimeIncrement {
      double time; ///< time at the beginning of the increment
      double dT;   ///< size of the time increment
    };

    /**
     * @brief The number of transported fields \f$n_f\f$ the material describes.
     * @return The number of fields.
     */
    virtual int getNumberOfFields() const = 0;

    /**
     * @brief Whether the content depends on the fields at all.
     * @return False for a material without storage (steady state transport), whose content increment is always zero.
     * @details Elements use it to tell the two meanings of a zero time increment apart: with storage, no time passes
     * and nothing may flow, without storage the steady state equations hold at every time.
     */
    virtual bool hasStorage() const { return true; }

    /**
     * @brief Compute the content increment, the flux and the source at the end of the increment, and their
     * algorithmic tangent.
     * @param[inout] response Sized to the number of fields by the caller; the state variables are updated in place.
     * @param[out] tangent Sized to the number of fields by the caller.
     * @param[in] fieldState The fields at the beginning and the end of the increment, the gradients at its end.
     * @param[in] timeIncrement Time information.
     */
    virtual void computeTransport( Response&            response,
                                   AlgorithmicTangent&  tangent,
                                   const FieldState&    fieldState,
                                   const TimeIncrement& timeIncrement ) const = 0;

    /**
     * @brief Size the response and the tangent to the number of fields of the material, all entries zero.
     * @param[out] response The response to size.
     * @param[out] tangent The tangent to size.
     */
    void initializeResponse( Response& response, AlgorithmicTangent& tangent ) const
    {
      const int nF = getNumberOfFields();
      response.contentIncrement.setZero( nF );
      response.flux.setZero( 3 * nF );
      response.source.setZero( nF );
      tangent.dContentIncrement_dField.setZero( nF, nF );
      tangent.dFlux_dField.setZero( 3 * nF, nF );
      tangent.dFlux_dGradient.setZero( 3 * nF, 3 * nF );
      tangent.dSource_dField.setZero( nF, nF );
    }

    /**
     * @brief Get a view to the state variables.
     * @param stateName Name of the state variable
     * @param stateVars Pointer to the state variable array
     * @return Marmot::Utilities::StateView to access the state variable
     */
    Marmot::Utilities::StateView getStateView( const std::string& stateName, double* stateVars ) const
    {
      return stateLayout.getStateView( stateVars, stateName );
    }

    /**
     * @brief Get the total number of required state variables.
     * @return Total number of required state variables
     */
    int getNumberOfRequiredStateVars() const { return stateLayout.totalSize(); }

    /**
     * @brief Initialize the state variables.
     * @param stateVars Pointer to the state variable array
     * @param nStateVars Number of state variables
     *
     * @note The default implementation initializes all state variables to zero.
     */
    virtual void initializeYourself( double* stateVars, int nStateVars )
    {
      for ( int i = 0; i < nStateVars; ++i )
        stateVars[i] = 0.0;
    }
  };

} // namespace Marmot::Materials
