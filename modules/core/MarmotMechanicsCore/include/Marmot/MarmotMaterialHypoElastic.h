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
#include "Marmot/MarmotPortability.h"
#include "Marmot/MarmotStateHelpers.h"
#include "Marmot/MarmotTypedefs.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace Marmot::Materials {

  /**
   *
   * Derived abstract base class for elastic materials expressed purely in rate form.
   *
   * In general, the nominal stress rate tensor \f$ \sigRate \f$ can be written as a function of the nominal stress
   * tensor \f$ \sig \f$, the stretching rate tensor \f$ \epsRate \f$ and the time \f$ t \f$.
   *
   * \f[  \displaystyle \sigRate = f( \sig, \epsRate, t, ...) \f]
   *
   * In course of numerical time integration, this relation will be formulated incrementally as
   *
   * \f[  \displaystyle \Delta \sig = f ( \sig_n, \Delta\eps, \Delta t, t_n, ...) \f]
   *
   * with
   *
   * \f[  \displaystyle \Delta\eps =  \epsRate\, \Delta t \f]
   *
   * and the algorithmic tangent
   *
   * \f[ \displaystyle \frac{d \sig }{d \eps } =  \frac{d \Delta \sig }{d \Delta \eps } \f]
   *
   * This formulation is compatible with an Abaqus interface.
   */
  class MaterialHypoElastic {

  protected:
    const double* materialProperties;  ///< Pointer to the array of material properties
    const int     nMaterialProperties; ///< Number of material properties

  public:
    const int materialNumber; ///< Integer identifier for this material instance
    /**
     * @brief Constructs the material with a given set of material properties and an identifier.
     * @param[in] matProperties_       Pointer to the array of material properties.
     * @param[in] nMaterialProperties_ Number of entries in @p matProperties_.
     * @param[in] materialNumber_      Integer identifying this material instance.
     */
    MaterialHypoElastic( const double* matProperties_, int nMaterialProperties_, int materialNumber_ )
      : materialProperties( matProperties_ ),
        nMaterialProperties( nMaterialProperties_ ),
        materialNumber( materialNumber_ )
    {
    }

    /// Default destructor
    virtual ~MaterialHypoElastic() = default;

    /// Layout of the state variables
    Marmot::Utilities::StateLayoutDynamic stateLayout;

    /// Structure to hold the material state at a material point in 3D
    struct state3D {
      Marmot::Math::Vector6d stress;               ///< Cauchy stress tensor in Voigt notation
      double                 elasticEnergyDensity; ///< Elastic strain energy density
      double                 dissipation;          ///< Dissipation
      double*                stateVars;            ///< Pointer to array of state variables
      /**
       * @brief Default constructor for state3D
       * Initializes stress to zero, energy and dissipation to zero, and stateVars to nullptr.
       */
      state3D()
        : stress( Marmot::Math::Vector6d::Zero() ),
          elasticEnergyDensity( 0.0 ),
          dissipation( 0.0 ),
          stateVars( nullptr )
      {
      }
      /**
       * @brief Constructor for initializing state3D
       * @param stress_ Cauchy stress tensor in Voigt notation
       * @param elasticEnergyDensity_ Elastic strain energy density
       * @param dissipation_ Dissipation
       * @param stateVars_ Pointer to array of state variables
       */
      state3D( Marmot::Math::Vector6d stress_, double elasticEnergyDensity_, double dissipation_, double* stateVars_ )
        : stress( stress_ ),
          elasticEnergyDensity( elasticEnergyDensity_ ),
          dissipation( dissipation_ ),
          stateVars( stateVars_ )
      {
      }
    };

    // Structure to hold the material state at a material point for 2D plane stress
    /// @brief Structure holding the material state at a material point for 2D plane stress.
    struct state2D {
      Marmot::Math::Vector3d stress;               ///< 2D Cauchy stress tensor in Voigt notation
      double                 elasticEnergyDensity; ///< Elastic strain energy density
      double                 dissipation;          ///< Dissipation
      double*                stateVars;            ///< Pointer to array of state variables

      state2D()
        : stress( Marmot::Math::Vector3d::Zero() ),
          elasticEnergyDensity( 0.0 ),
          dissipation( 0.0 ),
          stateVars( nullptr )
      {
      }
    };

    // Structure to hold the material state at a material point for 1D uniaxial stress
    /// @brief Structure holding the material state at a material point for 1D uniaxial stress.
    struct state1D {
      double  stress;               ///< 1D Cauchy stress
      double  elasticEnergyDensity; ///< Elastic strain energy density
      double  dissipation;          ///< Dissipation
      double* stateVars;            ///< Pointer to array of state variables

      state1D() : stress( 0.0 ), elasticEnergyDensity( 0.0 ), dissipation( 0.0 ), stateVars( nullptr ) {}
    };

    /// @brief Structure holding the material state at a material point (fiber) of a beam: the axial and the
    /// transverse shear stresses, \f$[ \sigma_{11}, \sigma_{12}, \sigma_{13} ]\f$, all other stresses vanish.
    struct stateBeam {
      Marmot::Math::Vector3d stress;               ///< [sigma11, sigma12, sigma13]
      double                 elasticEnergyDensity; ///< Elastic strain energy density
      double                 dissipation;          ///< Dissipation
      double*                stateVars;            ///< Pointer to array of state variables

      stateBeam()
        : stress( Marmot::Math::Vector3d::Zero() ),
          elasticEnergyDensity( 0.0 ),
          dissipation( 0.0 ),
          stateVars( nullptr )
      {
      }
    };

    /// @brief Structure carrying (pseudo-)time information passed to the material routines.
    struct timeInfo {
      double time; ///< Current (pseudo-)time
      double dT;   ///< (Pseudo-)time increment from the old (pseudo-)time to the current (pseudo-)time
    };

    /// Characteristic element length
    double characteristicElementLength;
    /**
     * Set the characteristic element length at the considered quadrature point.
     * It is needed for the regularization of materials with softening behavior based on the mesh-adjusted softening
     * modulus.
     *
     * @param[in] length characteristic length; will be assigned to @ref characteristicElementLength
     */
    MARMOT_API void setCharacteristicElementLength( double length );

    /**
     * For a given linearized strain increment \f$\Delta\boldsymbol{\varepsilon}\f$ at the old and the current time,
     * compute the Cauchy stress and the algorithmic tangent
     * \f$\frac{\partial\boldsymbol{\sigma}^{(n+1)}}{\partial\boldsymbol{\varepsilon}^{(n+1)}}\f$.
     *
     * @param[in,out]	state  A state3D instance carrying stress, strain energy, and state variables
     * @param[in,out]	dStress_dStrain	Algorithmic tangent representing the derivative of the Cauchy stress tensor with
     * respect to the linearized strain
     * @param[in]	dStrain linearized strain increment
     * @param[in]	timeInfo Structure carrying the current (pseudo-)time and the (pseudo-)time increment
     */
    virtual void computeStress( state3D&                      state,
                                Marmot::Math::Matrix6d&       dStress_dStrain,
                                const Marmot::Math::Vector6d& dStrain,
                                const timeInfo&               timeInfo ) const = 0;

    /**
     * Explicit version of @ref computeStress for use in explicit time integration schemes.
     * The algorithmic tangent is not needed in explicit schemes and will therefore not be computed.
     * @param[in,out] state  A state3D instance carrying stress, strain energy, and state variables
     * @param[in]   dStrain linearized strain increment
     * @param[in]   timeInfo Structure carrying time information
     *
     * @note The default implementation calls @ref computeStress and ignores the algorithmic tangent.
     * @note Derived classes may override this method for efficiency reasons.
     */
    virtual void computeStressExplicit( state3D&                      state,
                                        const Marmot::Math::Vector6d& dStrain,
                                        const timeInfo&               timeInfo ) const
    {
      Marmot::Math::Matrix6d dStress_dStrain = Marmot::Math::Matrix6d::Zero();
      computeStress( state, dStress_dStrain, dStrain, timeInfo );
    }

    /**
     * Plane stress implementation of @ref computeStress.
     */
    virtual void computePlaneStress( state2D&                      stress2D,
                                     Marmot::Math::Matrix3d&       dStress_dStrain2D,
                                     const Marmot::Math::Vector3d& dStrain2D,
                                     const timeInfo&               timeInfo ) const;

    /**
     * Uniaxial stress implementation of @ref computeStress.
     */
    virtual void computeUniaxialStress( state1D&        stress1D,
                                        double&         dStress_dStrain1D,
                                        const double    dStrain,
                                        const timeInfo& timeInfo ) const;

    /**
     * Beam stress implementation of @ref computeStress. The stress state of a fiber of a beam, with the prescribed
     * axial strain and transverse shear strains \f$[ \Delta\varepsilon_{11}, \Delta\gamma_{12}, \Delta\gamma_{13} ]\f$
     * and \f$\sigma_{22} = \sigma_{33} = \sigma_{23} = 0\f$ (computeStressWithFreeComponents()). The tangent is the
     * condensed 3x3 tangent \f$\partial [ \sigma_{11}, \sigma_{12}, \sigma_{13} ] / \partial [ \varepsilon_{11},
     * \gamma_{12}, \gamma_{13} ]\f$. Throws StressUpdateFailed if the iteration does not converge.
     *
     * @param[in,out] stateBeam           stresses, energy, dissipation and state variables of the fiber
     * @param[out]    dStress_dStrainBeam the condensed tangent
     * @param[in]     dStrainBeam         the increment \f$[ \Delta\varepsilon_{11}, \Delta\gamma_{12},
     * \Delta\gamma_{13} ]\f$
     * @param[in]     timeInfo            (pseudo-)time and its increment
     */
    virtual void computeBeamStress( stateBeam&                    stateBeam,
                                    Marmot::Math::Matrix3d&       dStress_dStrainBeam,
                                    const Marmot::Math::Vector3d& dStrainBeam,
                                    const timeInfo&               timeInfo ) const;

    /**
     * @brief A reduced stress state by static condensation of the free components: the strain components
     * \f$\boldsymbol{\varepsilon}_p\f$ are prescribed, the stress components \f$\boldsymbol{\sigma}_f = \mathbf{0}\f$,
     * all other strain increments are zero.
     *
     * @details Indices in Marmot's Voigt order \f$[11, 22, 33, 12, 13, 23]\f$ (shear strains as engineering strains),
     * e.g. plane stress \f$p = \{0, 1, 3\}, f = \{2\}\f$, uniaxial stress \f$p = \{0\}, f = \{1, 2\}\f$, beam stress
     * \f$p = \{0, 3, 4\}, f = \{1, 2, 5\}\f$. Newton's method on the free strains: the 3D material (computeStress())
     * is evaluated from the start-of-increment state every iteration (the state variables are reset), the free strains
     * are corrected by \f$\Delta\boldsymbol{\varepsilon}_f \leftarrow \Delta\boldsymbol{\varepsilon}_f -
     * \mathbf{C}_{ff}^{-1} \boldsymbol{\sigma}_f\f$, until \f$\max|\sigma_f| \leq 10^{-12} s\f$ (or \f$10^{-9} s\f$
     * after 8 iterations), with the stress scale \f$s = \max(1, \max|\sigma|)\f$ of the prescribed components: the
     * round-off of \f$\boldsymbol{\sigma}_f\f$ is relative to the size of the stresses, so an absolute tolerance cannot
     * be met for large stresses. The tangent is the condensed one, \f$\mathbf{C}_{pp} - \mathbf{C}_{pf}
     * \mathbf{C}_{ff}^{-1} \mathbf{C}_{fp}\f$.
     *
     * @tparam nP number of prescribed components
     * @tparam nF number of free components
     * @param[in,out] state            in: the start-of-increment stress (its free components are ignored), energy,
     * dissipation and state variables; out: the updated ones
     * @param[out]    dStressP_dStrainP the condensed tangent of the prescribed components
     * @param[in,out] dStrain          in: the prescribed strain increments at p and an initial guess at f (others
     * zero); out: the converged increment
     * @param[in]     p                the prescribed components
     * @param[in]     f                the free components
     * @param[in]     timeInfo         (pseudo-)time and its increment
     */
    template < int nP, int nF >
    void computeStressWithFreeComponents( state3D&                         state,
                                          Eigen::Matrix< double, nP, nP >& dStressP_dStrainP,
                                          Marmot::Math::Vector6d&          dStrain,
                                          const std::array< int, nP >&     p,
                                          const std::array< int, nF >&     f,
                                          const timeInfo&                  timeInfo ) const;

    /**
     * @brief Get a view to the state variables.
     * @param stateName Name of the state variable
     * @param stateVars Pointer to the state variable array
     * @return StatView to access the state variable
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
     * @brief Get the maximum wave speed of the material.
     * @param[in] state Current material state
     * @return Maximum wave speed
     * @details The default implementation computes the 3D algorithmic tangent and returns
     *          `sqrt(max(C_ii) / rho)` with `C_ii` from the Voigt tangent diagonal entries.
     */
    virtual double getMaximumWaveSpeed( const state3D& state ) const;

    /**
     * @brief Get the mass density of the material.
     * @param stateVars Pointer to the state variable array
     * @return Mass density of the material
     */
    virtual double getDensity( const double* stateVars ) const = 0;
  };

} // namespace Marmot::Materials
