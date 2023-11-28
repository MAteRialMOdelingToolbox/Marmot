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
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include <string>

namespace Marmot::Materials {

  class CompressibleNeoHooke : public MarmotMaterialFiniteStrain {
  public:
    using MarmotMaterialFiniteStrain::MarmotMaterialFiniteStrain;

    CompressibleNeoHooke( const double* materialProperties, int nMaterialProperties, int materialLabel );

    static constexpr int nStateVarsRequired = 0;

    void computeStress( ConstitutiveResponse< 3 >&,
                        AlgorithmicModuli< 3 >&,
                        const Deformation< 3 >&,
                        const TimeIncrement& );

    int getNumberOfRequiredStateVars() { return this->nStateVarsRequired; }

    void assignStateVars( double* stateVars, int nStateVars )
    {
      this->stateVars  = stateVars;
      this->nStateVars = nStateVars;
    };

    StateView getStateView( const std::string& result );
  };

  template < typename T >
  T psi( const Fastor::Tensor< T, 3, 3 >& C, const double K, const double G )
  {

    const T J  = sqrt( Fastor::determinant( C ) );
    const T I1 = Fastor::trace( C );

    T res = K / 8. * pow( J - 1. / J, 2. ) + G / 2. * ( I1 * pow( J, -2. / 3 ) - 3. );

    return res;
  }
} // namespace Marmot::Materials
