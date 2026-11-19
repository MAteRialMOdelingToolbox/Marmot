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

/// Test utilities shared by the beam element tests: a single element with its own state (Beam), a chain of elements
/// solved by Newton iterations (Chain), the tangent check against central differences, section point generators and
/// the elastica of a cantilever as reference solution. A new beam test starts at the physics, after including this.

#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numbers>
#include <string>
#include <utility>
#include <vector>

namespace Marmot::Testing::BeamTests {

  /// A single beam element with its own state
  struct Beam {
    std::unique_ptr< Marmot::Discretizations::Elements::Element > element;
    std::vector< double >                                         coordinates;
    std::vector< double >                                         properties;
    std::vector< double >                                         materialProperties;
    std::vector< double >                                         stateVars;
    int                                                           nDof;

    Beam( const std::string&    type,
          std::vector< double > coordinates_,
          std::vector< double > properties_,
          const std::string&    material,
          std::vector< double > materialProperties_ )
      : element( Marmot::Registry::ElementFactory::createElement( type, 1 ) ),
        coordinates( std::move( coordinates_ ) ),
        properties( std::move( properties_ ) ),
        materialProperties( std::move( materialProperties_ ) )
    {
      nDof = element->getNDofPerElement();
      element->assignNodeCoordinates( coordinates.data() );
      element->assignProperty(
        Marmot::Discretizations::ElementProperties( properties.data(), static_cast< int >( properties.size() ) ) );
      element->initializeYourself();
      element->assignProperty(
        Marmot::Discretizations::MaterialSection( material,
                                                  materialProperties.data(),
                                                  static_cast< int >( materialProperties.size() ) ) );
      stateVars.assign( element->getNumberOfRequiredStateVars(), 0.0 );
      element->assignStateVars( stateVars.data(), static_cast< int >( stateVars.size() ) );
      element->setInitialConditions( Marmot::Discretizations::Elements::Element::MarmotMaterialInitialization,
                                     nullptr );
    }

    Beam( Beam&& )            = default;
    Beam& operator=( Beam&& ) = delete;

    /// residual and tangent for the total dofs U from the committed state at UOld; commits if requested. The
    /// committed state is restored if the evaluation is not committed, also if computeKernels throws
    std::pair< Eigen::VectorXd, Eigen::MatrixXd > kernels( const Eigen::VectorXd& U,
                                                           const Eigen::VectorXd& UOld,
                                                           bool                   commit )
    {
      const auto            backup = stateVars;
      Eigen::VectorXd       P      = Eigen::VectorXd::Zero( nDof );
      Eigen::MatrixXd       K      = Eigen::MatrixXd::Zero( nDof, nDof );
      const Eigen::VectorXd dU     = U - UOld;
      try {
        element->computeKernels( U.data(), dU.data(), P.data(), K.data(), 0.0, 1.0 );
      }
      catch ( ... ) {
        stateVars = backup;
        throw;
      }
      if ( !commit )
        stateVars = backup;
      return { P, K };
    }

    /// central difference tangent from the committed state
    Eigen::MatrixXd numericalTangent( const Eigen::VectorXd& U, const Eigen::VectorXd& UOld, double h )
    {
      Eigen::MatrixXd K( nDof, nDof );
      for ( int j = 0; j < nDof; j++ ) {
        Eigen::VectorXd Up = U, Um = U;
        Up( j ) += h;
        Um( j ) -= h;
        K.col( j ) = ( kernels( Up, UOld, false ).first - kernels( Um, UOld, false ).first ) / ( 2 * h );
      }
      return K;
    }

    /// a quadrature point state of the element
    Eigen::VectorXd state( const std::string& name, int qp = 0 )
    {
      const auto view = element->getStateView( name, qp );
      return Eigen::Map< const Eigen::VectorXd >( view.stateLocation, view.stateSize );
    }
  };

  /// |value - expected| <= tol (absolute)
  inline void check( double value, double expected, double tol, const std::string& what )
  {
    throwExceptionOnFailure( std::abs( value - expected ) <= tol,
                             Marmot::Utilities::MakeString() << std::setprecision( 12 ) << what << ": " << value
                                                             << " != " << expected << " (tol " << tol << ")" );
  }

  /// |value - expected| <= relTol |expected| (relative)
  inline void checkRelative( double value, double expected, double relTol, const std::string& what )
  {
    check( value, expected, relTol * std::abs( expected ), what );
  }

  /// the tangent against central differences: max |K - K_num| <= relTol max |K_num|; returns the relative error
  inline double checkTangent( Beam&                  beam,
                              const Eigen::VectorXd& U,
                              const Eigen::VectorXd& UOld,
                              double                 h,
                              double                 relTol,
                              const std::string&     what )
  {
    const auto [P, K]  = beam.kernels( U, UOld, false );
    const auto   numK  = beam.numericalTangent( U, UOld, h );
    const double scale = numK.cwiseAbs().maxCoeff();
    const double err   = ( K - numK ).cwiseAbs().maxCoeff() / scale;
    std::cout << what << ": tangent max rel. error " << err << "\n";
    throwExceptionOnFailure( err <= relTol,
                             Marmot::Utilities::MakeString()
                               << what << ": tangent does not match the numerical one\nK=\n"
                               << K << "\nnumK=\n"
                               << numK );
    return err;
  }

  /// 2D section points [n, (y, A)...] of a rectangle of height h, width b: n layers with 2 Gauss points each (exact
  /// A, I)
  inline std::vector< double > rectangle2D( double h, double b, int nLayers )
  {
    const double          g = 1 / std::sqrt( 3. ), dy = h / nLayers;
    std::vector< double > p = { double( 2 * nLayers ) };
    for ( int i = 0; i < nLayers; i++ )
      for ( double a : { -g, g } )
        p.insert( p.end(), { -h / 2 + ( i + 0.5 + 0.5 * a ) * dy, 0.5 * b * dy } );
    return p;
  }

  /// 3D properties [v, J, n, (y, z, A)...] of a rectangle of height h (along y), width b (along z): ny x nz cells
  /// with 2 x 2 Gauss points each (exact A, Iy, Iz); the torsion constant J is given
  inline std::vector< double > rectangle3D( double          h,
                                            double          b,
                                            int             ny,
                                            int             nz,
                                            double          J,
                                            Eigen::Vector3d v = { 0, 1, 0 } )
  {
    const double          g = 1 / std::sqrt( 3. ), dy = h / ny, dz = b / nz;
    std::vector< double > p = { v( 0 ), v( 1 ), v( 2 ), J, double( 4 * ny * nz ) };
    for ( int i = 0; i < ny; i++ )
      for ( int k = 0; k < nz; k++ )
        for ( double a : { -g, g } )
          for ( double c : { -g, g } )
            p.insert( p.end(),
                      { -h / 2 + ( i + 0.5 + 0.5 * a ) * dy, -b / 2 + ( k + 0.5 + 0.5 * c ) * dz, 0.25 * dy * dz } );
    return p;
  }

  /// A chain of 2-node beam elements through given nodes (2D or 3D), clamped at the first node, solved by Newton
  /// iterations on the assembled system under load control with dead nodal loads at the last node
  struct Chain {
    std::vector< Beam > beams;
    int                 nEl, nDof, nd;
    Eigen::VectorXd     U, UOld;
    Eigen::VectorXd     applied; ///< the load at the last node of the last converged state

    Chain( const std::string&                    type,
           const std::vector< Eigen::VectorXd >& nodes,
           const std::vector< double >&          section,
           const std::string&                    material,
           const std::vector< double >&          matProps )
      : nEl( static_cast< int >( nodes.size() ) - 1 )
    {
      const int dim = static_cast< int >( nodes[0].size() );
      nd            = dim == 2 ? 3 : 6;
      nDof          = nd * ( nEl + 1 );
      for ( int e = 0; e < nEl; e++ ) {
        std::vector< double > X( nodes[e].data(), nodes[e].data() + dim );
        X.insert( X.end(), nodes[e + 1].data(), nodes[e + 1].data() + dim );
        beams.emplace_back( type, X, section, material, matProps );
      }
      U       = Eigen::VectorXd::Zero( nDof );
      UOld    = U;
      applied = Eigen::VectorXd::Zero( nd );
    }

    std::pair< Eigen::VectorXd, Eigen::MatrixXd > assemble( bool commit )
    {
      Eigen::VectorXd P = Eigen::VectorXd::Zero( nDof );
      Eigen::MatrixXd K = Eigen::MatrixXd::Zero( nDof, nDof );
      for ( int e = 0; e < nEl; e++ ) {
        const auto [Pe, Ke] = beams[e].kernels( U.segment( nd * e, 2 * nd ), UOld.segment( nd * e, 2 * nd ), commit );
        P.segment( nd * e, 2 * nd ) += Pe;
        K.block( nd * e, nd * e, 2 * nd, 2 * nd ) += Ke;
      }
      return { P, K };
    }

    /// increase the load at the last node (dead) from the current one to tipLoad in nSteps equal increments; returns
    /// the number of Newton iterations
    int solve( const Eigen::VectorXd& tipLoad, int nSteps, double tol = 1e-10 )
    {
      int                   total = 0;
      Eigen::VectorXd       F     = Eigen::VectorXd::Zero( nDof );
      const Eigen::VectorXd start = applied;
      for ( int step = 1; step <= nSteps; step++ ) {
        F.tail( nd )   = start + ( tipLoad - start ) * double( step ) / nSteps;
        bool converged = false;
        for ( int it = 0; it < 50; it++, total++ ) {
          const auto [P, K]       = assemble( false );
          const Eigen::VectorXd r = ( P - F ).tail( nDof - nd );
          if ( r.norm() < tol * std::max( 1.0, F.norm() ) ) {
            converged = true;
            break;
          }
          U.tail( nDof - nd ) += K.bottomRightCorner( nDof - nd, nDof - nd ).lu().solve( -r );
        }
        throwExceptionOnFailure( converged,
                                 Marmot::Utilities::MakeString() << "chain: Newton did not converge in step " << step );
        assemble( true );
        UOld    = U;
        applied = F.tail( nd );
      }
      return total;
    }

    /// the dofs of the last node
    Eigen::VectorXd tip() const { return U.tail( nd ); }
  };

  /// nEl + 1 equidistant nodes from the origin along x, in dim dimensions
  inline std::vector< Eigen::VectorXd > straightNodes( int nEl, double L, int dim )
  {
    std::vector< Eigen::VectorXd > n;
    for ( int e = 0; e <= nEl; e++ ) {
      Eigen::VectorXd x = Eigen::VectorXd::Zero( dim );
      x( 0 )            = e * L / nEl;
      n.push_back( x );
    }
    return n;
  }

  /// the elastica of a cantilever under a dead tip force P normal to the undeformed axis: tip deflection v/L, axial
  /// shortening (L - x)/L and tip rotation for lambda = P L^2 / EI. theta'' = -lambda cos(theta) on s in [0, 1],
  /// theta(0) = 0, theta'(1) = 0; first integral theta'(0)^2 = 2 lambda sin(theta_tip): bisection on theta_tip, RK4
  inline std::array< double, 3 > elastica( double lambda )
  {
    const int n     = 20000;
    auto      shoot = [&]( double thetaTip, double& x, double& y, double& th ) {
      double t = 0, k = std::sqrt( 2 * lambda * std::sin( thetaTip ) );
      x = y          = 0;
      const double h = 1.0 / n;
      auto         f = [&]( double t_, double k_ ) {
        return std::array< double, 4 >{ k_, -lambda * std::cos( t_ ), std::cos( t_ ), std::sin( t_ ) };
      };
      for ( int i = 0; i < n; i++ ) {
        const auto a = f( t, k );
        const auto b = f( t + 0.5 * h * a[0], k + 0.5 * h * a[1] );
        const auto c = f( t + 0.5 * h * b[0], k + 0.5 * h * b[1] );
        const auto d = f( t + h * c[0], k + h * c[1] );
        t += h / 6 * ( a[0] + 2 * b[0] + 2 * c[0] + d[0] );
        k += h / 6 * ( a[1] + 2 * b[1] + 2 * c[1] + d[1] );
        x += h / 6 * ( a[2] + 2 * b[2] + 2 * c[2] + d[2] );
        y += h / 6 * ( a[3] + 2 * b[3] + 2 * c[3] + d[3] );
      }
      th = t;
      return k; // theta'(1): < 0 if theta_tip is too small, > 0 if too large
    };
    double x, y, th, lo = 1e-9, hi = std::numbers::pi / 2 - 1e-12;
    for ( int it = 0; it < 60; it++ ) {
      const double mid                         = 0.5 * ( lo + hi );
      ( shoot( mid, x, y, th ) < 0 ? lo : hi ) = mid;
    }
    shoot( 0.5 * ( lo + hi ), x, y, th );
    return { y, 1 - x, th };
  }

} // namespace Marmot::Testing::BeamTests
