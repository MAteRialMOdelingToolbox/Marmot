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
 * Matthias Neuner matthias.neuner@uibk.ac.at
 * Magdalena Schreter magdalena.schreter@uibk.ac.at
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
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTypedefs.h"

#define VOIGTFROMDIM( x ) ( ( ( x * x ) + x ) >> 1 )

namespace Marmot {
  namespace ContinuumMechanics::VoigtNotation {

    constexpr int VoigtSize = 6;

    extern const Marmot::Vector6d P;
    extern const Marmot::Vector6d PInv;

    extern const Marmot::Vector6d I;
    extern const Marmot::Vector6d IHyd;
    extern const Matrix6d         IDev;

    // Plane Stress handling
    Eigen::Vector3d  voigtToPlaneVoigt( const Marmot::Vector6d& voigt );
    Marmot::Vector6d planeVoigtToVoigt( const Eigen::Vector3d& voigtPlane );

    template < int voigtSize >
    Eigen::Matrix< double, voigtSize, 1 > reduce3DVoigt( const Marmot::Vector6d& Voigt3D )
    {
      if constexpr ( voigtSize == 1 )
        return ( Eigen::Matrix< double, 1, 1 >() << Voigt3D( 0 ) ).finished();
      else if constexpr ( voigtSize == 3 )
        return voigtToPlaneVoigt( Voigt3D );
      else if constexpr ( voigtSize == 6 )
        return Voigt3D;
      else
        throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid dimension specified" );
    }

    template < int voigtSize >
    Marmot::Vector6d make3DVoigt( const Eigen::Matrix< double, voigtSize, 1 >& Voigt )
    {
      if constexpr ( voigtSize == 1 )
        return ( Marmot::Vector6d() << Voigt( 0 ), 0, 0, 0, 0, 0 ).finished();
      else if constexpr ( voigtSize == 3 )
        return planeVoigtToVoigt( Voigt );
      else if constexpr ( voigtSize == 6 )
        return Voigt;
      else
        throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid dimension specified" );
    }

    // function prototypes for  Marmot::Vector6d handling
    Eigen::Matrix3d  voigtToStrain( const Marmot::Vector6d& strainVector );
    Eigen::Matrix3d  voigtToStress( const Marmot::Vector6d& stressVector );
    Marmot::Vector6d strainToVoigt( const Eigen::Matrix3d& strainTensor );
    Marmot::Vector6d stressToVoigt( const Eigen::Matrix3d& stressTensor );

    template < int nDim >
    Eigen::Matrix< double, nDim, nDim > stressMatrixFromVoigt(
      const Eigen::Matrix< double, VOIGTFROMDIM( nDim ), 1 >& Voigt )
    {
      if constexpr ( nDim == 1 )
        return ( Eigen::Matrix< double, nDim, nDim >() << Voigt( 0 ) ).finished();
      else if constexpr ( nDim == 2 )
        return ( Eigen::Matrix< double, nDim, nDim >() << Voigt( 0 ), Voigt( 2 ), Voigt( 2 ), Voigt( 1 ) ).finished();
      else if constexpr ( nDim == 3 )
        return voigtToStress( Voigt );
      else
        throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid dimension specified" );
    }

    template < int nDim >
    Eigen::Matrix< double, VOIGTFROMDIM( nDim ), 1 > voigtFromStrainMatrix(
      const Eigen::Matrix< double, nDim, nDim >& strain )
    {
      if constexpr ( nDim == 1 )
        return ( Eigen::Matrix< double, VOIGTFROMDIM( nDim ), 1 >() << strain( 0, 0 ) ).finished();
      else if constexpr ( nDim == 2 )
        return ( Eigen::Matrix< double, VOIGTFROMDIM( nDim ), 1 >() << strain( 0, 0 ),
                 strain( 1, 1 ),
                 2 * strain( 0, 1 ) )
          .finished();
      else if constexpr ( nDim == 3 )
        return strainToVoigt( strain );
      else
        throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid dimension specified" );
    }

    namespace Invariants {

      // principal strains calculated by solving eigenvalue problem ( !NOT sorted! )
      Eigen::Vector3d principalStrains( const Marmot::Vector6d& strain );
      // principal stresses calculated by solving eigenvalue problem ( !NOT sorted! )
      Eigen::Vector3d principalStresses( const Marmot::Vector6d& stress );
      // principal strains calculated from haigh westergaard strains ( sorted --> e1 > e2 > e3 )
      Eigen::Vector3d sortedPrincipalStrains( const Marmot::Vector6d& strain );
      // principal stressDirections calculated by solving eigenvalue problem ( !NOT sorted! )
      Eigen::Matrix3d principalStressesDirections( const Marmot::Vector6d& stress );

      // equivalent von Mises stress
      double vonMisesEquivalentStress( const Marmot::Vector6d& stress );
      // equivalent von Mises strain
      double vonMisesEquivalentStrain( const Marmot::Vector6d& strain );
      // Euclidian norm of strain
      double normStrain( const Marmot::Vector6d& strain );
      // Euclidian norm of stress
      double normStress( const Marmot::Vector6d& stress );
      // Trace of compressive strains
      double StrainVolumetricNegative( const Marmot::Vector6d& strain );

      // Invariants - keep atention: different for stress/strain tensor
      double I1( const Marmot::Vector6d& stress );
      double I2( const Marmot::Vector6d& stress );
      double I2Strain( const Marmot::Vector6d& strain );
      double I3( const Marmot::Vector6d& stress );
      double I3Strain( const Marmot::Vector6d& strain );

      // Invariants of the deviatoric part of the stress/strain tensor
      double J2( const Marmot::Vector6d& stress );
      double J2Strain( const Marmot::Vector6d& strain );
      double J3( const Marmot::Vector6d& stress );
      double J3Strain( const Marmot::Vector6d& strain );

      // principal values in voigt
      std::pair< Eigen::Vector3d, Eigen::Matrix< double, 3, 6 > > principalValuesAndDerivatives(
        const Eigen::Matrix< double, 6, 1 >& S );

    } // namespace Invariants

    namespace Derivatives {

      // derivatives of Haigh Westergaard stresses with respect to cauchy stress in eng. notation
      Marmot::Vector6d dStressMean_dStress();
      Marmot::Vector6d dRho_dStress( double rho, const Marmot::Vector6d& stress );
      Marmot::Vector6d dTheta_dStress( double theta, const Marmot::Vector6d& stress );
      // derivatives of Haigh Westergaard stresses with respect to deviatoric invariants
      double dTheta_dJ2( const Marmot::Vector6d& stress );
      double dTheta_dJ3( const Marmot::Vector6d& stress );
      // derivatives of Haigh Westergaard strains with respect to deviatoric invariants
      double dThetaStrain_dJ2Strain( const Marmot::Vector6d& strain );
      double dThetaStrain_dJ3Strain( const Marmot::Vector6d& strain );
      // derivatives of deviatoric invariants with respect to eng. stresses
      Marmot::Vector6d dJ2_dStress( const Marmot::Vector6d& stress );
      Marmot::Vector6d dJ3_dStress( const Marmot::Vector6d& stress );
      // derivatives of deviatoric invariants with respect to eng. strains
      Marmot::Vector6d dJ2Strain_dStrain( const Marmot::Vector6d& strain );
      Marmot::Vector6d dJ3Strain_dStrain( const Marmot::Vector6d& strain );
      // derivatives of Haigh Westergaard strains with respect to eng. strains
      Marmot::Vector6d dThetaStrain_dStrain( const Marmot::Vector6d& strain );

      // derivatives of principalStess with respect to stress
      Marmot::Matrix36 dStressPrincipals_dStress( const Marmot::Vector6d& stress );

      // derivatives of plastic strains with respect to strains
      Eigen::Vector3d  dStrainVolumetricNegative_dStrainPrincipal( const Marmot::Vector6d& strain );
      Matrix6d         dEp_dE( const Matrix6d& CelInv, const Matrix6d& Cep );
      RowVector6d      dDeltaEpv_dE( const Matrix6d& CelInv, const Matrix6d& Cep );
      Marmot::Matrix36 dSortedStrainPrincipal_dStrain( const Marmot::Vector6d& dEp );
      RowVector6d      dDeltaEpvneg_dE( const Marmot::Vector6d& dEp, const Matrix6d& CelInv, const Matrix6d& Cep );

    } // namespace Derivatives

    namespace Transformations {

      Matrix6d  transformationMatrixStrainVoigt( const Matrix3d& transformedCoordinateSystem );
      Matrix6d  transformationMatrixStressVoigt( const Matrix3d& transformedCoordinateSystem );
      Matrix36d projectVoigtStressToPlane( const Vector3d& normalVector );
      Matrix36d projectVoigtStrainToPlane( const Vector3d& normalVector );
      /// rotate a 2nd order tensor T in voigt notation by : T' = Q * T * Q^T
      Marmot::Vector6d rotateVoigtStress( const Eigen::Matrix3d& Q, const Marmot::Vector6d& stress );

    } // namespace Transformations

  } // namespace ContinuumMechanics::VoigtNotation
} // namespace Marmot
