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
#include "Marmot/MarmotKelvinChain.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotSolidification.h"
#include <string>

namespace Marmot::Materials {

  /**
   * \brief Implementation of a linear elastic material
   * according to the B4 model by Bazant et al. (2015)
   * generalized for 3D stress states.
   *
   * For further information see \ref b4.
   */
  class B4 : public MarmotMaterialHypoElastic {

    /// \brief Poisson's ratio
    const double& nu;
    /**< #nu represents Poisson's ratio for isotropic linear elasticity.
     * It is a reference variable to #materialProperties[0]. */

    /// \brief asymptotic elastic compliance parameter
    const double& q1;
    /**< #q1 represents the first compliance parameter for %Solidification Theory.
     * It is a reference variable to #materialProperties[1]. */

    /// \brief viscoelastic compliance parameter
    const double& q2;
    /**< #q2 represents the second compliance parameter for %Solidification Theory.
     * It is a reference variable to #materialProperties[2]. */

    /// \brief viscoelastic compliance parameter
    const double& q3;
    /**< #q3 represents the third compliance parameter for %Solidification Theory.
     * It is a reference variable to #materialProperties[3]. */

    /// \brief flow compliance parameter
    const double& q4;
    /**< #q4 represents the fourth compliance parameter for %Solidification Theory.
     * It is a reference variable to #materialProperties[4]. */

    /// \brief log-power law exponent
    const double& n;
    /**< #n represents the exponent of the log-power law used in the compliance function of the hardened constituent.
     * It is a reference variable to #materialProperties[5]. */

    /// \brief solidified volume exponent
    const double& m;
    /**< #m represents the exponent in the function describing the solidification of the material.
     * It is a reference variable to #materialProperties[6]. */

    /// \brief number of Kelvin units to approximate the viscoelastic compliance
    const size_t nKelvinBasic;
    /**< #nKelvinBasic the number of Kelvin units used to approximate the viscoelastic compliance.
     * It is a reference variable to #materialProperties[7]. */

    /// \brief minimal retardation time used in the viscoelastic Kelvin chain
    const double& minTauBasic;
    /**< #minTauBasic represents the minimal retardation time in days used in the Kelvin chain to approximate the
     * viscoelastic compliance. It is a reference variable to #materialProperties[8]. */

    /// \brief ultimate autogenous shrinkage strain
    const double& ultimateAutogenousShrinkageStrain;
    /**< #ultimateAutogenousShrinkageStrain represents the ultimate autogenous shrinkage strain.
     * It is a reference variable to #materialProperties[9]. */

    /// \brief autogenous shrinkage half time
    const double& autogenousShrinkageHalfTime;
    /**< #autogenousShrinkageHalfTime represents the autogenous shrinkage half time.
     * It is a reference variable to #materialProperties[10]. */

    /// \brief autogenous shrinkage material parameter
    const double& alpha;
    /**< #alpha represents a cement type and concrete composition dependent material parameter for autogenous shrinkage.
     * It is a reference variable to #materialProperties[11]. */

    /// \brief autogenous shrinkage material parameter
    const double& rt;
    /**< #rt represents a cement type dependent material parameter for autogenous shrinkage, whose value is #rt = -4.50
     * for all cements.
     * It is a reference variable to #materialProperties[12]. */

    /// \brief ultimate drying shrinkage strain
    const double& ultimateDryingShrinkageStrain;
    /**< #ultimateDryingShrinkageStrain represents the ultimate drying shrinkage strain at zero relative ambient
     * humidity. It is a reference variable to #materialProperties[13]. */

    /// \brief drying shrinkage half time
    const double& dryingShrinkageHalfTime;
    /**< #dryingShrinkageHalfTime represents the drying shrinkage half time.
     * It is a reference variable to #materialProperties[14]. */

    /// \brief drying start time
    const double& dryingStart;
    /**< #dryingStart represents start of drying in days.
     * It is a reference variable to #materialProperties[15]. */

    /// \brief relative ambient humidity
    const double& hEnv;
    /**< #hEnv represents relative ambient humidity.
     * It is a reference variable to #materialProperties[16]. */

    /// \brief drying creep compliance parameter
    const double& q5;
    /**< #q5 represents drying creep compliance parameter.
     * It is a reference variable to #materialProperties[17]. */

    /// \brief number of Kelvin units to approximate the drying creep compliance
    const size_t nKelvinDrying;
    /**< #q5 represents the number of Kelvin units to approximate the drying creep compliance function.
     * It is a reference variable to #materialProperties[18]. */

    /// \brief minimal retardation time used in the drying creep Kelvin chain
    const double& minTauDrying;
    /**< #minTauDrying represents the minimal retardation time used in the Kelvin chain to approximate the drying creep
     * compliance function. It is a reference variable to #materialProperties[19]. */

    /// \brief time to start hydration
    const double& castTime;
    /**< #castTime represents the time at which hydration starts with respect to simulation time.
     * It is a reference variable to #materialProperties[20]. */

    /// \brief ratio of simulation time to days
    const double& timeToDays;
    /**< #timeToDays represents the ratio of simulation time to days.
     * It is a reference variable to #materialProperties[21]. */

  public:
    using MarmotMaterialHypoElastic::MarmotMaterialHypoElastic;

    B4( const double* materialProperties, int nMaterialProperties, int materialLabel );

    void computeStress( state3D& state, double* dStressDDStrain, const double* dStrain, const timeInfo& timeInfo );

    int getNumberOfRequiredStateVars();

    StateView getStateView( const std::string& stateName, double* stateVars );

  private:
    /// \brief material parameters for Solidification Theory
    SolidificationTheory::Parameters solidificationParameters;
    /**
     * \brief properties of the Kelvin chain for approximating
     * the viscoelastic compliance of the %Solidification Theory */
    SolidificationTheory::KelvinChainProperties solidificationKelvinProperties;

    /// \brief Young's modulus of the Kelvin units representing basic creep
    KelvinChain::Properties basicCreepElasticModuli;
    /// \brief retardation times of the Kelvin units representing basic creep
    KelvinChain::Properties basicCreepRetardationTimes;

    /// \brief approximation order of the Post-Widder formula for drying creep
    static constexpr int dryingCreepComplianceApproximationOrder = 5;
    /// \brief approximation order of the Post-Widder formula for basic creep
    static constexpr int basicCreepComplianceApproximationOrder = 2;

    /// \brief drying creep compliance function
    template < typename T_ >
    T_ phi( T_ xi, double b, double xiZero )
    {
      T_ val = sqrt( exp( tanh( sqrt( xi - xiZero ) ) * b ) - exp( tanh( sqrt( -xiZero ) ) * b ) );
      return val;
    }
  };
} // namespace Marmot::Materials
