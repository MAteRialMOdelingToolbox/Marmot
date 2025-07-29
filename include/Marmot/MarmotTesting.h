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
 * Alexander Dummer alexander.dummer@uibk.ac.at
 * Paul Hofer Paul.Hofer@uibk.ac.at
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

#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "autodiff/forward/dual/dual.hpp"
#include "unsupported/Eigen/CXX11/Tensor"
#include <Eigen/Core>
#include <Fastor/Fastor.h>
#include <iostream>

namespace Marmot::Testing {

  bool checkIfEqual( const double a, const double b, const double tol = 1e-15 );

  bool checkIfEqual( const autodiff::dual a, const autodiff::dual b, const double tol = 1e-15 );

  bool checkIfEqual( const std::complex< double > a, const std::complex< double > b, const double tol = 1e-15 );

  std::string getString( const double a );
  std::string getString( const autodiff::dual a );

  template < typename T >
  bool checkIfEqual( const Eigen::Matrix< T, -1, -1 >& a,
                     const Eigen::Matrix< T, -1, -1 >& b,
                     const double                      tol = 1e-15 )
  {
    if ( a.rows() != b.rows() || a.cols() != b.cols() ) {
      return false;
    }
    for ( int i = 0; i < a.rows(); i++ ) {
      for ( int j = 0; j < a.cols(); j++ ) {
        auto cond = checkIfEqual( a( i, j ), b( i, j ), tol );
        if ( !cond ) {
          std::cout << "  -> HINT:  a(" << i << "," << j << ") = " << getString( a( i, j ) ) << " !=  b(" << i << ","
                    << j << ") =" << getString( b( i, j ) ) << std::endl;
          return false;
        }
      }
    }
    return true;
  }

  template < typename T, long int... Rest >
  bool checkIfEqual( const Eigen::TensorFixedSize< T, Eigen::Sizes< Rest... > >& a,
                     const Eigen::TensorFixedSize< T, Eigen::Sizes< Rest... > >& b,
                     const double                                                tol = 1e-15 )
  {
    const T* a_data = a.data();
    const T* b_data = b.data();

    for ( int i = 0; i < a.size(); i++ ) {
      auto cond = checkIfEqual( a_data[i], b_data[i], tol );
      if ( !cond ) {
        std::cout << "  -> HINT:  a(" << i << ") = " << getString( a_data[i] ) << " !=  b(" << i
                  << ") =" << getString( b_data[i] ) << std::endl;
        return false;
      }
    }
    return true;
  }

  template < typename T, size_t... Rest >
  bool checkIfEqual( const Fastor::Tensor< T, Rest... >& a,
                     const Fastor::Tensor< T, Rest... >& b,
                     const double                        tol = 1e-15 )
  {
    const T* a_data = a.data();
    const T* b_data = b.data();

    for ( size_t i = 0; i < a.size(); i++ ) {
      auto cond = checkIfEqual( a_data[i], b_data[i], tol );
      if ( !cond ) {
        std::cout << "  -> HINT:  a(" << i << ") = " << getString( a_data[i] ) << " !=  b(" << i
                  << ") =" << getString( b_data[i] ) << std::endl;
        return false;
      }
    }
    return true;
  }
  void throwExceptionOnFailure( const bool condition, const std::string& message = "" );

  void executeTestsAndCollectExceptions( const std::vector< std::function< void() > >& testFunctions );

  template < int N >
  Eigen::Matrix< double, N, 2 > fibonacciLatticeHemisphere()
  {
    double                        phi;
    double                        theta;
    Eigen::Matrix< double, 1, 2 > pt;
    Eigen::Matrix< double, N, 2 > pts;

    // i1 = 1, 2, ..., N
    for ( int i1 = 1; i1 <= N; i1++ ) {
      theta = acos( -( i1 - 1. ) / N );
      phi   = ( i1 - 1. ) * Constants::GoldenAngle;

      pts.row( i1 - 1 ) = Eigen::Vector2d( phi, theta );
    }

    return pts;
  };

  bool spinTurbokreisel( const std::unique_ptr< MarmotMaterialHypoElastic >& material,
                         double*                                             stress,
                         double*                                             dStress_dStrain,
                         const double*                                       dStrain,
                         const double*                                       timeOld,
                         const double                                        dT,
                         double&                                             pNewDT,
                         const double                                        stressTol    = 1e-15,
                         const double                                        stiffnessTol = 1e-15 );

} // namespace Marmot::Testing
