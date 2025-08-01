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
#include "Marmot/MarmotStateVarVectorManager.h"
#include <iostream>
#include <string>
#include <vector>

namespace Marmot::Materials {

  /**
   * \brief Implementation of an orthotropic linear viscoelastic material model
   * following a Power Law compliance function and assuming constant Poisson's ratios
   * generalized for 3D stress states.
   *
   */
  class LinearViscoelasticOrthotropicPowerLaw : public MarmotMaterialHypoElastic {

    /// \brief Young's modulus in x1 direction
    const double& E1;

    /// \brief Young's modulus in x2 direction
    const double& E2;

    /// \brief Young's modulus in x3 direction
    const double& E3;

    /// \brief Poisson's ratio
    const double& nu12;

    /// \brief Poisson's ratio
    const double& nu23;

    /// \brief Poisson's ratio
    const double& nu13;

    /// \brief Shear modulus in x1-x2 plane
    const double& G12;

    /// \brief Shear modulus in x2-x3 plane
    const double& G23;

    /// \brief Shear modulus in x1-x3 plane
    const double& G13;

    /// \brief power law compliance parameter
    const double& m;

    /// \brief power law exponent
    const double& n;

    /// \brief approximation order for the retardation spectrum
    const int powerLawApproximationOrder;

    /// \brief number of Kelvin units to approximate the viscoelastic compliance
    const size_t nKelvin;

    /// \brief minimal retardation time used in the viscoelastic Kelvin chain
    const double& minTau;

    /// \brief log spacing between the retardation times of the Kelvin Chain
    const double& spacing;

    /// \brief ratio of simulation time to days
    const double& timeToDays;

    /// \brief direction x1 w.r.t the global coordinate system
    const Vector3d direction1;

    /// \brief direction x2 w.r.t the global coordinate system
    const Vector3d direction2;

    class LinearViscoelasticOrthotropicPowerLawStateVarManager : public MarmotStateVarVectorManager {

    public:
      inline const static auto layout = makeLayout( {
        { .name = "kelvinStateVars", .length = 0 },
      } );

      KelvinChain::mapStateVarMatrix kelvinStateVars;

      LinearViscoelasticOrthotropicPowerLawStateVarManager( double* theStateVarVector, int nKelvinUnits )
        : MarmotStateVarVectorManager( theStateVarVector, layout ),
          kelvinStateVars( &find( "kelvinStateVars" ), 6, nKelvinUnits ){};
    };
    std::unique_ptr< LinearViscoelasticOrthotropicPowerLawStateVarManager > stateVarManager;

  public:
    using MarmotMaterialHypoElastic::MarmotMaterialHypoElastic;

    LinearViscoelasticOrthotropicPowerLaw( const double* materialProperties,
                                           int           nMaterialProperties,
                                           int           materialLabel );

    void computeStress( double* stress,
                        double* dStressDDStrain,

                        const double* dStrain,
                        const double* timeOld,
                        const double  dT,
                        double&       pNewDT );

    int getNumberOfRequiredStateVars();

    void assignStateVars( double* stateVars_, int nStateVars );

    StateView getStateView( const std::string& stateName );

  private:
    KelvinChain::Properties elasticModuli;
    KelvinChain::Properties retardationTimes;
    double                  zerothKelvinChainCompliance;

    Matrix6d CInv, Cel, CelUnit, CelUnitInv, CelUnitGlobal;

    Matrix3d localCoordinateSystem;
  };
} // namespace Marmot::Materials
