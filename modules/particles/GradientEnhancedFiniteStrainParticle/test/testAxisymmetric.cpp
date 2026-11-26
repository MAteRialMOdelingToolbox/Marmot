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

// The axisymmetric particles (the page "Axisymmetric Particles" of the documentation): the displacement SQCNI particle
// and the gradient-enhanced SQCNIxNSNI particle in the meridian half plane (r, z), on patches of particles whose
// nodes are a grid of B-spline kernels.

#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotMeshfreeKernelFunctionBSpline2ndOrderBoxed.h"
#include "Marmot/MarmotMeshfreeReproducingKernelApproximation.h"
#include "Marmot/MarmotParticleLibrary.h"
#include "Marmot/MarmotParticleStiffnessVisitors.h"
#include "Marmot/MarmotParticleTesting.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

using namespace Marmot;
using namespace Marmot::Discretizations::Meshfree;
using namespace Marmot::Testing;

namespace {

  constexpr double Pi = Marmot::Math::Constants::Pi;

  // COMPRESSIBLENEOHOOKE: K, G, rho
  const std::vector< double > neoHookeProperties = { 3500., 1500., 2.0 };
  // GRADIENTENHANCEDCOMPRESSIBLENEOHOOKEDAMAGE: K, G, kappa0, kappaF, l, rho
  const std::vector< double > damageProperties = { 3500., 1500., 1e-3, 1e-2, 0.3, 2.0 };

  /// a particle type and its material
  struct ParticleType {
    std::string           name;
    std::string           material;
    std::vector< double > properties;
    int                   nBlock; ///< dofs per node: the displacement, then the nonlocal field
    double                density;
  };

  const ParticleType displacement = { "DisplacementSQCNI/Axisymmetric/Quad",
                                      "COMPRESSIBLENEOHOOKE",
                                      neoHookeProperties,
                                      2,
                                      neoHookeProperties[2] };

  const ParticleType gradientEnhanced = { "GradientEnhancedFiniteStrainSQCNIxNSNI/Axisymmetric/Quad",
                                          "GRADIENTENHANCEDCOMPRESSIBLENEOHOOKEDAMAGE",
                                          damageProperties,
                                          3,
                                          damageProperties[5] };

  /**
   * A patch of nr x nz square particles of the edge h on [r0, r0 + nr h] x [z0, z0 + nz h], one node per particle,
   * a 2nd-order B-spline kernel at its center (support radius 2.1 h, as for nodes that are the particles), the
   * approximation reproducing linear fields.
   * Each particle gets the kernels whose support covers one of its evaluation points, as the host's particle manager
   * does it.
   */
  struct Patch {
    static constexpr double                supportRadius = 2.1; ///< of the kernels, relative to h
    const ParticleType&                    type;
    double                                 r0, z0, h;
    int                                    nr, nz;
    std::vector< std::array< double, 2 > > nodes;
    std::vector< std::unique_ptr< KernelFunctionBSpline2ndOrderBoxed > > kernels;
    ReproducingKernelApproximation                                       approximation{ 2, 1 };
    std::vector< std::unique_ptr< Particle > >                           particles;
    std::vector< std::vector< double > >                                 stateVars;
    std::vector< std::vector< int > >                                    nodesOfParticle;

    Patch( const ParticleType& type_, double r0_, double z0_, double h_, int nr_, int nz_ )
      : type( type_ ), r0( r0_ ), z0( z0_ ), h( h_ ), nr( nr_ ), nz( nz_ )
    {
      for ( int j = 0; j < nz; j++ )
        for ( int i = 0; i < nr; i++ )
          nodes.push_back( { r0 + ( i + 0.5 ) * h, z0 + ( j + 0.5 ) * h } );
      for ( auto& node : nodes )
        kernels.emplace_back(
          std::make_unique< KernelFunctionBSpline2ndOrderBoxed >( node.data(), 2, supportRadius * h ) );

      for ( int j = 0; j < nz; j++ )
        for ( int i = 0; i < nr; i++ ) {
          const double                r = r0 + i * h, z = z0 + j * h;
          const std::vector< double > vertices = { r, z, r + h, z, r + h, z + h, r, z + h };
          particles.emplace_back(
            Marmot::Registry::ParticleFactory::createParticle( type.name,
                                                               static_cast< int >( particles.size() + 1 ),
                                                               vertices.data(),
                                                               8,
                                                               0.0,
                                                               type.material,
                                                               type.properties.data(),
                                                               static_cast< int >( type.properties.size() ),
                                                               approximation ) );
          stateVars.emplace_back( particles.back()->getNumberOfRequiredStateVars(), 0.0 );
          particles.back()->assignStateVars( stateVars.back().data(), static_cast< int >( stateVars.back().size() ) );
          particles.back()->initializeYourself();
        }
      assignKernels();
    }

    void assignKernels()
    {
      nodesOfParticle.clear();
      for ( auto& particle : particles ) {
        const int             nEval = particle->getNumberOfEvaluationPoints();
        std::vector< double > eval( 2 * nEval );
        particle->getEvaluationCoordinates( eval.data() );
        std::vector< const KernelFunction* > assigned;
        std::vector< int >                   indices;
        for ( size_t A = 0; A < kernels.size(); A++ )
          for ( int e = 0; e < nEval; e++ )
            if ( kernels[A]->isInSupport( &eval[2 * e] ) ) {
              assigned.push_back( kernels[A].get() );
              indices.push_back( static_cast< int >( A ) );
              break;
            }
        particle->assignMeshfreeKernelFunctions( assigned );
        nodesOfParticle.push_back( indices );
      }
    }

    int nDof() const { return static_cast< int >( nodes.size() ) * type.nBlock; }

    /// the nodal values of the particle p, gathered from the global vector
    Eigen::VectorXd gather( int p, const Eigen::VectorXd& Q ) const
    {
      const auto&     indices = nodesOfParticle[p];
      Eigen::VectorXd q( indices.size() * type.nBlock );
      for ( size_t a = 0; a < indices.size(); a++ )
        q.segment( a * type.nBlock, type.nBlock ) = Q.segment( indices[a] * type.nBlock, type.nBlock );
      return q;
    }

    /// add the values of the particle p to the global vector, or the global matrix
    void scatter( int p, const Eigen::VectorXd& f, Eigen::VectorXd& F ) const
    {
      const auto& indices = nodesOfParticle[p];
      for ( size_t a = 0; a < indices.size(); a++ )
        F.segment( indices[a] * type.nBlock, type.nBlock ) += f.segment( a * type.nBlock, type.nBlock );
    }

    void scatter( int p, const Eigen::MatrixXd& k, Eigen::MatrixXd& K ) const
    {
      const auto& indices = nodesOfParticle[p];
      for ( size_t a = 0; a < indices.size(); a++ )
        for ( size_t b = 0; b < indices.size(); b++ )
          K.block( indices[a] * type.nBlock,
                   indices[b] * type.nBlock,
                   type.nBlock,
                   type.nBlock ) += k.block( a * type.nBlock, b * type.nBlock, type.nBlock, type.nBlock );
    }

    /// the internal force and its tangent of the increment dQ; the state is put back unless accepted
    std::pair< Eigen::VectorXd, Eigen::MatrixXd > forceAndTangent( const Eigen::VectorXd& dQ, bool accept = false )
    {
      Eigen::VectorXd F = Eigen::VectorXd::Zero( nDof() );
      Eigen::MatrixXd K = Eigen::MatrixXd::Zero( nDof(), nDof() );
      for ( size_t p = 0; p < particles.size(); p++ ) {
        const auto      backup = stateVars[p];
        const auto      q      = gather( static_cast< int >( p ), dQ );
        Eigen::VectorXd f      = Eigen::VectorXd::Zero( q.size() );
        Eigen::MatrixXd k      = Eigen::MatrixXd::Zero( q.size(), q.size() );
        particles[p]->computePhysicsKernels( q.data(), f.data(), k.data(), 1.0, 1.0 );
        scatter( static_cast< int >( p ), f, F );
        scatter( static_cast< int >( p ), k, K );
        if ( accept )
          particles[p]->acceptStateAndPosition();
        else
          stateVars[p] = backup;
      }
      if ( accept )
        assignKernels();
      return { F, K };
    }

    /// a displacement field u(r, z) at the nodes (the nonlocal dofs zero)
    Eigen::VectorXd field( const std::function< std::array< double, 2 >( double, double ) >& u ) const
    {
      Eigen::VectorXd Q = Eigen::VectorXd::Zero( nDof() );
      for ( size_t A = 0; A < nodes.size(); A++ ) {
        const auto uA            = u( nodes[A][0], nodes[A][1] );
        Q( A * type.nBlock )     = uA[0];
        Q( A * type.nBlock + 1 ) = uA[1];
      }
      return Q;
    }

    /// the nodes whose support lies in the interior: no evaluation point of the boundary of the patch in it; the axis
    /// (r0 = 0) is no boundary of the body of revolution
    std::vector< int > interiorNodes() const
    {
      std::vector< int > interior;
      const double       R = supportRadius * h;
      for ( size_t A = 0; A < nodes.size(); A++ )
        if ( ( r0 == 0.0 || nodes[A][0] - R > r0 ) && nodes[A][0] + R < r0 + nr * h && nodes[A][1] - R > z0 &&
             nodes[A][1] + R < z0 + nz * h )
          interior.push_back( static_cast< int >( A ) );
      return interior;
    }

    /// a 3 x 3 state tensor of particle p (row-major, as the material points store it)
    Eigen::Matrix3d stateTensor( int p, const std::string& name ) const
    {
      const auto view = particles[p]->getStateView( name, 0 );
      return Eigen::Map< const Eigen::Matrix< double, 3, 3, Eigen::RowMajor > >( view.stateLocation );
    }
  };

  /// the largest displacement force of the interior nodes, relative to the largest of all nodes
  double interiorForce( const Patch& patch, const Eigen::VectorXd& F )
  {
    double largest = 0.0, interior = 0.0;
    for ( size_t A = 0; A < patch.nodes.size(); A++ )
      for ( int i = 0; i < 2; i++ )
        largest = std::max( largest, std::abs( F( A * patch.type.nBlock + i ) ) );
    for ( int A : patch.interiorNodes() )
      for ( int i = 0; i < 2; i++ )
        interior = std::max( interior, std::abs( F( A * patch.type.nBlock + i ) ) );
    return interior / largest;
  }

  // The patch tests of a uniform state, u_r = a r, u_z = b z (sigma_rr = sigma_tt, sigma_rz = 0, in equilibrium
  // without a body force): every particle sees the same deformation and stress, and the interior nodes have no force,
  // from the undeformed state and again from the deformed one after an accepted increment. On a patch away from the
  // axis (r0 > 0) and on one at the axis (r0 = 0), whose particles touch it with a face: a face on the axis has the
  // ring weight 0, and the nodes near the axis are interior nodes.
  void testUniformStatesArePatchTests( const ParticleType& type, double a, double b, double r0 )
  {
    Patch patch( type, r0, 0.0, 0.25, 8, 8 );
    throwExceptionOnFailure( patch.interiorNodes().size() == ( r0 == 0.0 ? 24u : 16u ),
                             "the patch has 16 interior nodes (24 at the axis)" );

    for ( auto [ar, bz] : std::vector< std::pair< double, double > >{ { a, 0.0 }, { 0.0, b }, { a, b } } ) {
      const Eigen::VectorXd Q = patch.field( [&]( double r, double z ) -> std::array< double, 2 > {
        return { ar * r, bz * z };
      } );
      const auto [F, K]       = patch.forceAndTangent( Q );
      const double rel        = interiorForce( patch, F );
      std::cout << type.name << ": patch test at r0 = " << r0 << ", u_r = " << ar << " r, u_z = " << bz
                << " z: interior force / largest " << rel << std::endl;
      throwExceptionOnFailure( rel < 1e-11,
                               Marmot::Utilities::MakeString()
                                 << type.name << ": u_r = " << ar << " r, u_z = " << bz << " z: the interior force is "
                                 << rel << " of the largest" );

      // every particle the same increment of the deformation gradient, with the hoop stretch 1 + a
      for ( size_t p = 0; p < patch.particles.size(); p++ ) {
        const auto      backup = patch.stateVars[p];
        const auto      q      = patch.gather( static_cast< int >( p ), Q );
        Eigen::VectorXd f      = Eigen::VectorXd::Zero( q.size() );
        Eigen::MatrixXd k      = Eigen::MatrixXd::Zero( q.size(), q.size() );
        patch.particles[p]->computePhysicsKernels( q.data(), f.data(), k.data(), 1.0, 1.0 );
        const Eigen::Matrix3d deltaF   = patch.stateTensor( static_cast< int >( p ), "delta deformation gradient" );
        patch.stateVars[p]             = backup;
        const Eigen::Matrix3d expected = Eigen::Vector3d( 1 + ar, 1 + bz, 1 + ar ).asDiagonal();
        throwExceptionOnFailure( ( deltaF - expected ).cwiseAbs().maxCoeff() < 1e-12,
                                 Marmot::Utilities::MakeString()
                                   << type.name << ": particle " << p << " sees the deformation gradient increment\n"
                                   << deltaF );
      }
    }

    // and again from the deformed state
    patch.forceAndTangent( patch.field( [&]( double r, double z ) -> std::array< double, 2 > {
      return { a * r, b * z };
    } ),
                           true );
    const auto [F, K] = patch.forceAndTangent( patch.field( [&]( double r, double z ) -> std::array< double, 2 > {
      return { 0.5 * a * r, -0.5 * b * z };
    } ) );
    const double rel  = interiorForce( patch, F );
    std::cout << type.name << ": patch test at r0 = " << r0
              << ", second increment from the deformed state: interior force / largest " << rel << std::endl;
    throwExceptionOnFailure( rel < 1e-11,
                             Marmot::Utilities::MakeString()
                               << type.name << ": the second increment: the interior force is " << rel
                               << " of the largest" );
  }

  // The fields with a shear, u_z = c r and u_r = c z, whose stress needs a body force: the smoothed gradient and the
  // hoop strain of every particle reproduce them exactly (the exact ring averages of linear fields), also at the axis.
  void testLinearFieldsAreReproduced( const ParticleType& type, double c, double r0 )
  {
    Patch patch( type, r0, 0.0, 0.25, 4, 4 );
    for ( int which = 0; which < 2; which++ ) {
      const Eigen::VectorXd Q = patch.field( [&]( double r, double z ) -> std::array< double, 2 > {
        return which == 0 ? std::array< double, 2 >{ 0.0, c * r } : std::array< double, 2 >{ c * z, 0.0 };
      } );
      for ( size_t p = 0; p < patch.particles.size(); p++ ) {
        std::array< double, 2 > center;
        patch.particles[p]->getCenterCoordinates( center.data() );
        const auto      backup = patch.stateVars[p];
        const auto      q      = patch.gather( static_cast< int >( p ), Q );
        Eigen::VectorXd f      = Eigen::VectorXd::Zero( q.size() );
        Eigen::MatrixXd k      = Eigen::MatrixXd::Zero( q.size(), q.size() );
        patch.particles[p]->computePhysicsKernels( q.data(), f.data(), k.data(), 1.0, 1.0 );
        const Eigen::Matrix3d deltaF = patch.stateTensor( static_cast< int >( p ), "delta deformation gradient" );
        patch.stateVars[p]           = backup;

        Eigen::Matrix3d expected = Eigen::Matrix3d::Identity();
        if ( which == 0 )
          expected( 1, 0 ) = c;
        else {
          expected( 0, 1 ) = c;
          expected( 2, 2 ) = 1 + c * center[1] / center[0];
        }
        throwExceptionOnFailure( ( deltaF - expected ).cwiseAbs().maxCoeff() < 1e-12,
                                 Marmot::Utilities::MakeString() << type.name << ": particle " << p << ", field "
                                                                 << which << ": the deformation gradient increment is\n"
                                                                 << deltaF << "\ninstead of\n"
                                                                 << expected );
      }
    }
  }

  // The patch test with a load on the boundary: u_r = a r, u_z = b z, and on the outer faces of the patch the traction
  // sigma n of its uniform stress as a PRESSURE, p = -sigma_nn. The load balances the internal force of every node,
  // the boundary nodes included, to round-off: the pressure is integrated with the face rule of the smoothing, the
  // rule of the integration constraint. On a patch at the axis (r0 = 0) the faces on the axis carry no load: the axis
  // is no boundary. Then the tangent of the pressure against central differences, after an accepted increment, for an
  // increment with all components.
  void testTractionPatchTest( const ParticleType& type, double a, double b, double r0 )
  {
    Patch      patch( type, r0, 0.0, 0.25, 4, 4 );
    const int  nBlock = type.nBlock;
    const auto Q = patch.field( [&]( double r, double z ) -> std::array< double, 2 > { return { a * r, b * z }; } );

    // the faces of particle (i, j) on the boundary of the patch (Abaqus: 1 z = z0, 2 r = r1, 3 z = z1, 4 r = r0)
    const auto loadedFaces = [&]( int i, int j ) {
      std::vector< int > faces;
      if ( j == 0 )
        faces.push_back( 1 );
      if ( i == patch.nr - 1 )
        faces.push_back( 2 );
      if ( j == patch.nz - 1 )
        faces.push_back( 3 );
      if ( i == 0 && r0 > 0.0 )
        faces.push_back( 4 );
      return faces;
    };

    Eigen::VectorXd FInt = Eigen::VectorXd::Zero( patch.nDof() ), FLoad = Eigen::VectorXd::Zero( patch.nDof() );
    for ( int j = 0; j < patch.nz; j++ )
      for ( int i = 0; i < patch.nr; i++ ) {
        const int       p      = j * patch.nr + i;
        const auto      backup = patch.stateVars[p];
        const auto      q      = patch.gather( p, Q );
        const int       n      = static_cast< int >( q.size() );
        Eigen::VectorXd f      = Eigen::VectorXd::Zero( n );
        Eigen::MatrixXd k      = Eigen::MatrixXd::Zero( n, n );
        patch.particles[p]->computePhysicsKernels( q.data(), f.data(), k.data(), 1.0, 1.0 );
        patch.scatter( p, f, FInt );
        // the Cauchy stress of the first increment, the Kirchhoff stress over det Delta F (3D)
        const Eigen::Matrix3d sigma = patch.stateTensor( p, "stress" ) /
                                      patch.stateTensor( p, "delta deformation gradient" ).determinant();
        for ( int face : loadedFaces( i, j ) ) {
          const double    pressure = face % 2 == 1 ? -sigma( 1, 1 ) : -sigma( 0, 0 );
          Eigen::VectorXd fExt     = Eigen::VectorXd::Zero( n );
          Eigen::MatrixXd kExt     = Eigen::MatrixXd::Zero( n, n );
          patch.particles[p]->computeDistributedLoad( 0, face, &pressure, fExt.data(), kExt.data(), 1.0, 1.0 );
          patch.scatter( p, fExt, FLoad );
        }
        patch.stateVars[p] = backup;
      }
    double largest = 0.0, residual = 0.0;
    for ( size_t A = 0; A < patch.nodes.size(); A++ )
      for ( int i = 0; i < 2; i++ ) {
        largest  = std::max( largest, std::abs( FInt( A * nBlock + i ) ) );
        residual = std::max( residual, std::abs( FInt( A * nBlock + i ) + FLoad( A * nBlock + i ) ) );
      }
    std::cout << type.name << ": traction patch test at r0 = " << r0 << ": largest nodal residual / largest force "
              << residual / largest << std::endl;
    throwExceptionOnFailure( residual < 1e-11 * largest,
                             Marmot::Utilities::MakeString()
                               << type.name << ": the traction patch test at r0 = " << r0 << ": the residual is "
                               << residual / largest << " of the largest force" );

    // the tangent of the pressure on the corner particle (r1, z1), faces 2 and 3
    patch.forceAndTangent( patch.field( [&]( double r, double z ) -> std::array< double, 2 > {
      return { 0.01 * r + 0.004 * z, -0.006 * z + 0.002 * r };
    } ),
                           true );
    const int  p        = patch.nz * patch.nr - 1;
    const int  n        = static_cast< int >( patch.nodesOfParticle[p].size() ) * nBlock;
    auto&      particle = *patch.particles[p];
    auto&      state    = patch.stateVars[p];
    const auto load     = [&]( const Eigen::VectorXd& q, Eigen::MatrixXd* K ) {
      const auto      backup = state;
      Eigen::VectorXd f = Eigen::VectorXd::Zero( n ), fExt = Eigen::VectorXd::Zero( n );
      Eigen::MatrixXd k = Eigen::MatrixXd::Zero( n, n ), kExt = Eigen::MatrixXd::Zero( n, n );
      particle.computePhysicsKernels( q.data(), f.data(), k.data(), 1.0, 1.0 );
      for ( int face : { 2, 3 } ) {
        const double pressure = face == 2 ? 0.37 : -0.21;
        particle.computeDistributedLoad( 0, face, &pressure, fExt.data(), kExt.data(), 1.0, 1.0 );
      }
      state = backup;
      if ( K )
        *K = kExt;
      return fExt;
    };
    Eigen::VectorXd dQ = Eigen::VectorXd::Zero( n );
    for ( int dof = 0; dof < n; dof++ )
      if ( dof % nBlock < 2 )
        dQ( dof ) = 0.01 * std::sin( 1.7 * ( dof + 1 ) );
    Eigen::MatrixXd K;
    load( dQ, &K );
    Eigen::MatrixXd KNumerical( n, n );
    for ( int c = 0; c < n; c++ ) {
      const double    step = 1e-7;
      Eigen::VectorXd e    = Eigen::VectorXd::Zero( n );
      e( c )               = step;
      KNumerical.col( c )  = ( load( dQ + e, nullptr ) - load( dQ - e, nullptr ) ) / ( 2 * step );
    }
    const double error = ( K - KNumerical ).norm() / KNumerical.norm();
    std::cout << type.name << ": tangent of the pressure at r0 = " << r0
              << " against central differences, relative error " << error << std::endl;
    throwExceptionOnFailure( error < 1e-7,
                             Marmot::Utilities::MakeString()
                               << type.name << ": the tangent of the pressure is inconsistent, relative error "
                               << error );
  }

  // The ring volume of every particle, 2 pi A r_c, and the lumped mass, rho 2 pi int R dA, node by node summing to it,
  // over the patch rho pi (r1^2 - r0^2) (z1 - z0).
  void testRingVolumeAndLumpedMass( const ParticleType& type )
  {
    Patch        patch( type, 2.0, 0.0, 0.25, 4, 4 );
    const double h         = patch.h;
    double       totalMass = 0.0;
    for ( size_t p = 0; p < patch.particles.size(); p++ ) {
      std::array< double, 2 > center;
      patch.particles[p]->getCenterCoordinates( center.data() );
      const double V0 = 2 * Pi * h * h * center[0];
      throwExceptionOnFailure( std::abs( patch.particles[p]->getVolumeUndeformed() - V0 ) < 1e-14 * V0,
                               type.name + ": the ring volume of a particle" );
      const int       n = static_cast< int >( patch.nodesOfParticle[p].size() ) * type.nBlock;
      Eigen::VectorXd M = Eigen::VectorXd::Zero( n );
      patch.particles[p]->computeLumpedInertia( M.data() );
      double massOfParticle = 0.0;
      for ( int a = 0; a < n; a += type.nBlock ) {
        massOfParticle += M( a );
        throwExceptionOnFailure( M( a + 1 ) == M( a ), type.name + ": the same lumped mass of u_r and u_z" );
      }
      throwExceptionOnFailure( std::abs( massOfParticle - type.density * V0 ) < 1e-13 * type.density * V0,
                               type.name + ": the lumped mass of a particle is rho 2 pi A r_c" );
      totalMass += massOfParticle;
    }
    const double r1 = patch.r0 + patch.nr * h, z1 = patch.z0 + patch.nz * h;
    const double expected = type.density * Pi * ( r1 * r1 - patch.r0 * patch.r0 ) * ( z1 - patch.z0 );
    throwExceptionOnFailure( std::abs( totalMass - expected ) < 1e-12 * expected,
                             type.name + ": the lumped mass of the patch is rho pi (r1^2 - r0^2) h" );
  }

  // The tangent of one particle against central differences, after an accepted increment (F_n != I, the hoop stretch
  // != 1), for an increment with all components (and a nonlocal field that damages the gradient-enhanced material).
  // The displacement particle's tangent is exact. The NSNI stabilization of the gradient-enhanced particle omits the
  // derivative of dtau/dDeltaF itself (its dependence on the deformation and on the damage), as in plane strain: its
  // nonlocal rows are exact, its displacement rows within the bound of the plane-strain test. For a displacement
  // increment that is a linear field, whose second gradient and hence stabilization vanish, its tangent is exact as
  // well, including the couplings of the hoop stress with the nonlocal field.
  void testTangent( const ParticleType& type )
  {
    Patch      patch( type, 1.0, 0.0, 0.25, 3, 3 );
    const auto stretch = patch.field( [&]( double r, double z ) -> std::array< double, 2 > {
      return { 0.01 * r + 0.004 * z, -0.006 * z + 0.002 * r };
    } );
    patch.forceAndTangent( stretch, true );

    const int   p        = 4; // the middle one
    const int   nBlock   = type.nBlock;
    const auto& nodes    = patch.nodesOfParticle[p];
    const int   n        = static_cast< int >( nodes.size() ) * nBlock;
    auto&       particle = *patch.particles[p];
    auto&       state    = patch.stateVars[p];

    const auto force = [&]( const Eigen::VectorXd& q, Eigen::MatrixXd* K ) {
      const auto      backup = state;
      Eigen::VectorXd f      = Eigen::VectorXd::Zero( n );
      Eigen::MatrixXd k      = Eigen::MatrixXd::Zero( n, n );
      particle.computePhysicsKernels( q.data(), f.data(), k.data(), 1.0, 1.0 );
      state = backup;
      if ( K )
        *K = k;
      return f;
    };

    // the relative errors of the tangent at dQ, of the blocks U-U, U-N, N-U, N-N, and of all of it
    const auto tangentErrors = [&]( const Eigen::VectorXd& dQ ) {
      Eigen::MatrixXd K;
      force( dQ, &K );
      Eigen::MatrixXd KNumerical( n, n );
      for ( int b = 0; b < n; b++ ) {
        const double    step = b % nBlock < 2 ? 1e-7 : 1e-9;
        Eigen::VectorXd e    = Eigen::VectorXd::Zero( n );
        e( b )               = step;
        KNumerical.col( b )  = ( force( dQ + e, nullptr ) - force( dQ - e, nullptr ) ) / ( 2 * step );
      }
      const auto blockError = [&]( bool uRows, bool uColumns ) {
        double num = 0, den = 0;
        for ( int a = 0; a < n; a++ )
          for ( int b = 0; b < n; b++ )
            if ( ( a % nBlock < 2 ) == uRows && ( b % nBlock < 2 ) == uColumns ) {
              num += std::pow( K( a, b ) - KNumerical( a, b ), 2 );
              den += std::pow( KNumerical( a, b ), 2 );
            }
        return den > 0 ? std::sqrt( num / den ) : std::sqrt( num );
      };
      return std::array< double, 5 >{ blockError( true, true ),
                                      blockError( true, false ),
                                      blockError( false, true ),
                                      blockError( false, false ),
                                      ( K - KNumerical ).norm() / KNumerical.norm() };
    };

    Eigen::VectorXd dQ( n ), dQLinear( n );
    for ( int a = 0; a < n; a++ ) {
      const auto& Y = patch.nodes[nodes[a / nBlock]];
      const auto  i = a % nBlock;
      dQ( a )       = i < 2 ? 0.01 * std::sin( 1.7 * ( a + 1 ) ) : 3e-3 * ( 1. + 0.3 * std::cos( 0.9 * ( a + 1 ) ) );
      dQLinear( a ) = i == 0 ? 0.01 * Y[0] + 0.003 * Y[1] : i == 1 ? -0.004 * Y[0] + 0.006 * Y[1] : dQ( a );
    }

    const auto e = tangentErrors( dQ );
    if ( nBlock == 2 ) {
      std::cout << type.name << ": tangent against central differences, relative error " << e[4] << std::endl;
      throwExceptionOnFailure( e[4] < 1e-7,
                               Marmot::Utilities::MakeString()
                                 << type.name << ": the tangent is inconsistent, relative error " << e[4] );
    }
    else {
      const auto eLinear = tangentErrors( dQLinear );
      std::cout << type.name << ": tangent against central differences, relative errors U-U " << e[0] << ", U-N "
                << e[1] << " (both NSNI), N-U " << e[2] << ", N-N " << e[3] << "; for a linear displacement increment "
                << eLinear[4] << std::endl;
      throwExceptionOnFailure( e[2] < 5e-5 && e[3] < 1e-7,
                               Marmot::Utilities::MakeString()
                                 << type.name << ": the nonlocal rows are inconsistent, " << e[2] << ", " << e[3] );
      throwExceptionOnFailure( e[0] < 5e-3,
                               Marmot::Utilities::MakeString()
                                 << type.name << ": the U-U block is beyond the NSNI approximation, " << e[0] );
      throwExceptionOnFailure( eLinear[4] < 1e-6,
                               Marmot::Utilities::MakeString()
                                 << type.name << ": the tangent of a linear displacement increment is inconsistent, "
                                 << eLinear[4] );
    }
  }

  // The bound of |K| W the particle describes (Particle::describeStiffness()) against its stiffness, and that the
  // description leaves the particle as it was; and the bound of the critical time step of the assembled patch, the
  // Collatz-Wielandt bound max_i (|K| 1)_i / m_i of the displacement block with the ring-weighted lumped mass, never
  // below the largest eigenvalue of M^-1 K.
  void testStiffnessBound( const ParticleType& type, bool explicitPath )
  {
    Patch patch( type, 1.0, 0.0, 0.25, 4, 4 );
    patch.forceAndTangent( patch.field( [&]( double r, double z ) -> std::array< double, 2 > {
      return { 0.01 * r + 0.004 * z, -0.006 * z };
    } ),
                           true );
    const int nBlock = type.nBlock;

    // one particle against its stiffness; its state is put back for the patch
    {
      const int             p         = 5;
      const auto            committed = patch.stateVars[p];
      const int             n         = static_cast< int >( patch.nodesOfParticle[p].size() ) * nBlock;
      const auto            backup    = patch.stateVars[p];
      Eigen::VectorXd       f         = Eigen::VectorXd::Zero( n );
      Eigen::MatrixXd       K         = Eigen::MatrixXd::Zero( n, n );
      const Eigen::VectorXd zero      = Eigen::VectorXd::Zero( n );
      patch.particles[p]->computePhysicsKernels( zero.data(), f.data(), K.data(), 1.0, 1.0 );
      patch.stateVars[p] = backup;
      Eigen::MatrixXd W  = Eigen::MatrixXd::Zero( n, nBlock == 2 ? 1 : 2 );
      for ( int a = 0; a < n; a++ )
        W( a, a % nBlock < 2 ? 0 : 1 ) = 1.0 + 0.5 * std::sin( 0.7 * a );
      checkAbsoluteStiffnessBoundAgainstStiffness( type.name, *patch.particles[p], K, W, 1.5 );
      Eigen::VectorXd dQ( n );
      for ( int a = 0; a < n; a++ )
        dQ( a ) = a % nBlock < 2 ? 0.003 * std::sin( 1.3 * a ) : 1e-4;
      checkDescribingTheStiffnessLeavesTheParticleUnchanged( type.name,
                                                             *patch.particles[p],
                                                             patch.stateVars[p],
                                                             dQ,
                                                             explicitPath );
      patch.stateVars[p] = committed;
    }

    // the assembled patch
    const int       nDof  = patch.nDof();
    Eigen::MatrixXd K     = Eigen::MatrixXd::Zero( nDof, nDof );
    Eigen::VectorXd M     = Eigen::VectorXd::Zero( nDof );
    Eigen::VectorXd bound = Eigen::VectorXd::Zero( nDof );
    for ( size_t p = 0; p < patch.particles.size(); p++ ) {
      const int       n = static_cast< int >( patch.nodesOfParticle[p].size() ) * nBlock;
      Eigen::MatrixXd k = Eigen::MatrixXd::Zero( n, n );
      Eigen::VectorXd m = Eigen::VectorXd::Zero( n ), w = Eigen::VectorXd::Zero( n ), b = Eigen::VectorXd::Zero( n );
      for ( int a = 0; a < n; a++ )
        w( a ) = a % nBlock < 2 ? 1.0 : 0.0;
      assembleStiffness( *patch.particles[p], k.data(), 1.0, 1.0 );
      patch.particles[p]->computeLumpedInertia( m.data() );
      addBoundOfAbsoluteStiffnessTimes( *patch.particles[p], w.data(), 1, b.data(), 1.0, 1.0 );
      patch.scatter( static_cast< int >( p ), k, K );
      patch.scatter( static_cast< int >( p ), m, M );
      patch.scatter( static_cast< int >( p ), b, bound );
    }
    std::vector< int > u;
    for ( int a = 0; a < nDof; a++ )
      if ( a % nBlock < 2 )
        u.push_back( a );
    const int       nU = static_cast< int >( u.size() );
    Eigen::MatrixXd MinvK( nU, nU );
    double          collatzWielandt = 0.0;
    for ( int i = 0; i < nU; i++ ) {
      for ( int j = 0; j < nU; j++ )
        MinvK( i, j ) = K( u[i], u[j] ) / M( u[i] );
      collatzWielandt = std::max( collatzWielandt, bound( u[i] ) / M( u[i] ) );
    }
    const auto   eigenvalues = Eigen::EigenSolver< Eigen::MatrixXd >( MinvK, false ).eigenvalues();
    const double largest     = eigenvalues.cwiseAbs().maxCoeff();
    throwExceptionOnFailure( collatzWielandt >= largest && collatzWielandt < 10 * largest,
                             Marmot::Utilities::MakeString()
                               << type.name << ": the bound of the critical time step, " << collatzWielandt
                               << ", against the largest eigenvalue of M^-1 K, " << largest );
    std::cout << type.name << ": Collatz-Wielandt bound / largest eigenvalue of M^-1 K = " << collatzWielandt / largest
              << std::endl;
  }

  // Lame's thick-walled cylinder in plane strain (u_z = 0), inner radius a, outer b, inner pressure p:
  // u_r = (1 + nu) p a^2 / (E (b^2 - a^2)) ((1 - 2 nu) r + b^2 / r).
  double lame( double r, double a, double b, double p )
  {
    const double K = neoHookeProperties[0], G = neoHookeProperties[1];
    const double E = 9 * K * G / ( 3 * K + G ), nu = ( 3 * K - 2 * G ) / ( 2 * ( 3 * K + G ) );
    return ( 1 + nu ) * p * a * a / ( E * ( b * b - a * a ) ) * ( ( 1 - 2 * nu ) * r + b * b / r );
  }

  /// solve for the radial dofs of a system with u_z = 0 at all nodes (the z dofs are the odd ones of a node block)
  Eigen::VectorXd solveRadial(
    const std::function< std::pair< Eigen::VectorXd, Eigen::MatrixXd >( const Eigen::VectorXd& ) >& residualAndTangent,
    int                                                                                             nDof,
    int                                                                                             nBlock )
  {
    std::vector< int > radial;
    for ( int a = 0; a < nDof; a += nBlock )
      radial.push_back( a );
    Eigen::VectorXd Q = Eigen::VectorXd::Zero( nDof );
    for ( int iteration = 0; iteration < 4; iteration++ ) {
      const auto [R, K] = residualAndTangent( Q );
      Eigen::VectorXd r( radial.size() );
      Eigen::MatrixXd k( radial.size(), radial.size() );
      for ( size_t i = 0; i < radial.size(); i++ ) {
        r( i ) = R( radial[i] );
        for ( size_t j = 0; j < radial.size(); j++ )
          k( i, j ) = K( radial[i], radial[j] );
      }
      const Eigen::VectorXd dq = k.fullPivLu().solve( -r );
      for ( size_t i = 0; i < radial.size(); i++ )
        Q( radial[i] ) += dq( i );
    }
    return Q;
  }

  // The thick-walled cylinder a = 1, b = 2 under the inner pressure p, in plane strain (u_z = 0 everywhere, a slice
  // of height 0.5): the displacement particles (PRESSURE on their inner faces) against Lame's solution and against a
  // mesh of CX8RUL elements (the ring pressure as consistent nodal forces), in the small strain limit.
  void testThickWalledCylinder()
  {
    const double a = 1.0, b = 2.0, height = 0.5, p = 1e-2;

    // the particles
    Patch      patch( displacement, a, 0.0, 0.125, 8, 4 );
    const auto particleSystem = [&]( const Eigen::VectorXd& Q ) {
      auto [R, K] = patch.forceAndTangent( Q );
      for ( int j = 0; j < patch.nz; j++ ) {
        const int       index  = j * patch.nr; // the particle on the inner surface of row j
        const int       n      = static_cast< int >( patch.nodesOfParticle[index].size() ) * 2;
        const auto      q      = patch.gather( index, Q );
        const auto      backup = patch.stateVars[index];
        Eigen::VectorXd f = Eigen::VectorXd::Zero( n ), pExt = Eigen::VectorXd::Zero( n );
        Eigen::MatrixXd k = Eigen::MatrixXd::Zero( n, n ), kExt = Eigen::MatrixXd::Zero( n, n );
        patch.particles[index]->computePhysicsKernels( q.data(), f.data(), k.data(), 1.0, 1.0 );
        patch.particles[index]->computeDistributedLoad( 0, 4, &p, pExt.data(), kExt.data(), 1.0, 1.0 );
        patch.stateVars[index] = backup;
        patch.scatter( index, pExt, R );
        patch.scatter( index, kExt, K );
      }
      return std::make_pair( R, K );
    };
    const Eigen::VectorXd QParticles = solveRadial( particleSystem, patch.nDof(), 2 );

    // the particles' displacement at a radius, in the middle of the slice
    const auto particleDisplacement = [&]( double r ) {
      const std::array< double, 2 >        point = { r, 0.5 * height };
      std::vector< const KernelFunction* > covering;
      std::vector< int >                   indices;
      for ( size_t A = 0; A < patch.kernels.size(); A++ )
        if ( patch.kernels[A]->isInSupport( point.data() ) ) {
          covering.push_back( patch.kernels[A].get() );
          indices.push_back( static_cast< int >( A ) );
        }
      Eigen::VectorXd N( covering.size() );
      patch.approximation.computeShapeFunctions( point.data(), covering, N.data() );
      double u = 0.0;
      for ( size_t i = 0; i < indices.size(); i++ )
        u += N( i ) * QParticles( 2 * indices[i] );
      return u;
    };

    // the CX8RUL mesh: 8 elements across, one in z; nodes of a row of corners, a row of midsides, ...
    const int                              nElements = 8;
    const double                           he        = ( b - a ) / nElements;
    std::vector< std::array< double, 2 > > feNodes;
    std::map< std::pair< int, int >, int > nodeAt; // (2 i, 2 j) on the grid of half the element size
    for ( int j = 0; j <= 2; j++ )
      for ( int i = 0; i <= 2 * nElements; i++ ) {
        if ( j == 1 && i % 2 == 1 )
          continue; // the center of an 8-node element
        nodeAt[{ i, j }] = static_cast< int >( feNodes.size() );
        feNodes.push_back( { a + 0.5 * he * i, 0.5 * height * j } );
      }
    std::vector< std::unique_ptr< Marmot::Discretizations::Elements::Element > > elements;
    std::vector< std::vector< double > >                                         elementStateVars;
    // an element keeps pointers to its coordinates and its section: they live as long as it
    std::vector< std::vector< double > >           elementCoordinates( nElements );
    const Marmot::Discretizations::MaterialSection materialSection( "COMPRESSIBLENEOHOOKE",
                                                                    neoHookeProperties.data(),
                                                                    static_cast< int >( neoHookeProperties.size() ) );
    std::vector< std::array< int, 8 > >            connectivity;
    for ( int e = 0; e < nElements; e++ ) {
      const int                  i           = 2 * e;
      const std::array< int, 8 > c           = { nodeAt[{ i, 0 }],
                                                 nodeAt[{ i + 2, 0 }],
                                                 nodeAt[{ i + 2, 2 }],
                                                 nodeAt[{ i, 2 }],
                                                 nodeAt[{ i + 1, 0 }],
                                                 nodeAt[{ i + 2, 1 }],
                                                 nodeAt[{ i + 1, 2 }],
                                                 nodeAt[{ i, 1 }] };
      auto&                      coordinates = elementCoordinates[e];
      for ( int node : c ) {
        coordinates.push_back( feNodes[node][0] );
        coordinates.push_back( feNodes[node][1] );
      }
      elements.emplace_back( Marmot::Registry::ElementFactory::createElement( "CX8RUL", e + 1 ) );
      elements.back()->assignNodeCoordinates( coordinates.data() );
      elements.back()->assignProperty( materialSection );
      elementStateVars.emplace_back( elements.back()->getNumberOfRequiredStateVars(), 0.0 );
      elements.back()->assignStateVars( elementStateVars.back().data(),
                                        static_cast< int >( elementStateVars.back().size() ) );
      elements.back()->initializeYourself();
      connectivity.push_back( c );
    }
    const int       nFeDof = 2 * static_cast< int >( feNodes.size() );
    Eigen::VectorXd feLoad = Eigen::VectorXd::Zero( nFeDof );
    // the ring pressure on the inner face r = a, a quadratic edge of the length height: 1/6, 4/6, 1/6
    feLoad( 2 * nodeAt[{ 0, 0 }] ) += 2 * Pi * a * p * height / 6;
    feLoad( 2 * nodeAt[{ 0, 1 }] ) += 2 * Pi * a * p * height * 4 / 6;
    feLoad( 2 * nodeAt[{ 0, 2 }] ) += 2 * Pi * a * p * height / 6;
    const auto feSystem = [&]( const Eigen::VectorXd& Q ) {
      Eigen::VectorXd R = -feLoad;
      Eigen::MatrixXd K = Eigen::MatrixXd::Zero( nFeDof, nFeDof );
      for ( size_t e = 0; e < elements.size(); e++ ) {
        Eigen::VectorXd q( 16 ), dq = Eigen::VectorXd::Zero( 16 ), Pe = Eigen::VectorXd::Zero( 16 );
        Eigen::MatrixXd Ke = Eigen::MatrixXd::Zero( 16, 16 );
        for ( int A = 0; A < 8; A++ )
          q.segment< 2 >( 2 * A ) = Q.segment< 2 >( 2 * connectivity[e][A] );
        elements[e]->computeKernels( q.data(), dq.data(), Pe.data(), Ke.data(), 1.0, 1.0 );
        for ( int A = 0; A < 8; A++ ) {
          R.segment< 2 >( 2 * connectivity[e][A] ) += Pe.segment< 2 >( 2 * A );
          for ( int B = 0; B < 8; B++ )
            K.block< 2, 2 >( 2 * connectivity[e][A], 2 * connectivity[e][B] ) += Ke.block< 2, 2 >( 2 * A, 2 * B );
        }
      }
      return std::make_pair( R, K );
    };
    const Eigen::VectorXd QElements = solveRadial( feSystem, nFeDof, 2 );

    for ( double r : { a, 1.5, b } ) {
      const double exact = lame( r, a, b, p );
      const double uFE = QElements( 2 * nodeAt[{ static_cast< int >( std::lround( ( r - a ) / ( 0.5 * he ) ) ), 1 }] );
      const double uRKPM   = particleDisplacement( r );
      const double errorFE = std::abs( uFE - exact ) / exact;
      const double errorRK = std::abs( uRKPM - exact ) / exact;
      const double rkVsFE  = std::abs( uRKPM - uFE ) / uFE;
      std::cout << "thick-walled cylinder, u_r(" << r << "): Lame " << exact << ", CX8RUL " << uFE << " (" << errorFE
                << "), particles " << uRKPM << " (" << errorRK << ", " << rkVsFE << " from CX8RUL)" << std::endl;
      throwExceptionOnFailure( errorFE < 1e-3,
                               Marmot::Utilities::MakeString()
                                 << "CX8RUL against Lame at r = " << r << ": " << errorFE );
      throwExceptionOnFailure( errorRK < 1e-2 && rkVsFE < 1e-2,
                               Marmot::Utilities::MakeString()
                                 << "the axisymmetric particles against Lame and CX8RUL at r = " << r << ": " << errorRK
                                 << ", " << rkVsFE );
    }
  }

  // A particle with a vertex at r < 0 is refused, one that touches the axis with a face is accepted; VCI and
  // CWFCORRECTION are refused; an increment that pulls the ring through the axis throws.
  void testRefusals( const ParticleType& type )
  {
    ReproducingKernelApproximation approximation{ 2, 1 };
    const auto                     create = [&]( const std::vector< double >& vertices ) {
      return std::unique_ptr< Particle >(
        Marmot::Registry::ParticleFactory::createParticle( type.name,
                                                           1,
                                                           vertices.data(),
                                                           8,
                                                           0.0,
                                                           type.material,
                                                           type.properties.data(),
                                                           static_cast< int >( type.properties.size() ),
                                                           approximation ) );
    };
    const auto throws = []( const std::function< void() >& f ) {
      try {
        f();
      }
      catch ( const std::exception& ) {
        return true;
      }
      return false;
    };
    throwExceptionOnFailure( throws( [&]() {
                               create( { -0.1, 0.0, 0.2, 0.0, 0.2, 0.3, -0.1, 0.3 } );
                             } ),
                             type.name + ": a vertex at r < 0 is refused" );
    throwExceptionOnFailure( !throws( [&]() {
      create( { 0.0, 0.0, 0.2, 0.0, 0.2, 0.3, 0.0, 0.3 } );
    } ),
                             type.name + ": a particle touching the axis is accepted" );

    Patch patch( type, 0.25, 0.0, 0.25, 2, 2 );
    auto& particle = *patch.particles[0];
    throwExceptionOnFailure( throws( [&]() { particle.vci_getNumberOfConstraints(); } ),
                             type.name + ": VCI is refused" );
    const double vciOrder = 1.0;
    throwExceptionOnFailure( throws( [&]() { particle.setProperty( "VCI order", &vciOrder ); } ),
                             type.name + ": a VCI order is refused" );
    throwExceptionOnFailure( particle.getSupportedDistributedLoadTypes().count( "CWFCORRECTION" ) == 0,
                             type.name + ": CWFCORRECTION is not offered" );
    throwExceptionOnFailure( throws( [&]() {
                               const int       n = static_cast< int >( patch.nodesOfParticle[0].size() ) * type.nBlock;
                               Eigen::VectorXd f = Eigen::VectorXd::Zero( n );
                               Eigen::MatrixXd k = Eigen::MatrixXd::Zero( n, n );
                               particle.computeDistributedLoad( 1, 1, nullptr, f.data(), k.data(), 1.0, 1.0 );
                             } ),
                             type.name + ": CWFCORRECTION is refused" );

    // u_r = -1.5 r pulls every ring through the axis
    const Eigen::VectorXd through = patch.field( []( double r, double ) -> std::array< double, 2 > {
      return { -1.5 * r, 0.0 };
    } );
    throwExceptionOnFailure( throws( [&]() { patch.forceAndTangent( through ); } ),
                             type.name + ": a ring pulled through the axis is refused" );
  }

  // The particle says its section (Particle::getSection()): axisymmetric, and Cartesian for its plane-strain twin.
  void testTheParticleSaysItsSection( const ParticleType& type )
  {
    Patch patch( type, 0.25, 0.0, 0.25, 1, 1 );
    throwExceptionOnFailure( patch.particles[0]->getSection() == Section::Axisymmetric,
                             type.name + ": the section is axisymmetric" );

    ParticleType planeStrainTwin = type;
    planeStrainTwin.name.replace( planeStrainTwin.name.find( "Axisymmetric" ), 12, "PlaneStrain" );
    Patch planeStrainPatch( planeStrainTwin, 0.25, 0.0, 0.25, 1, 1 );
    throwExceptionOnFailure( planeStrainPatch.particles[0]->getSection() == Section::Cartesian,
                             planeStrainTwin.name + ": the section is Cartesian" );
  }

} // namespace

int main()
{
  std::vector< std::function< void() > > testFunctions;
  for ( double r0 : { 2.0, 0.0 } ) {
    testFunctions.push_back( [r0]() { testUniformStatesArePatchTests( displacement, 0.02, -0.01, r0 ); } );
    testFunctions.push_back( [r0]() { testUniformStatesArePatchTests( gradientEnhanced, 2e-4, -1e-4, r0 ); } );
    testFunctions.push_back( [r0]() { testTractionPatchTest( displacement, 0.02, -0.01, r0 ); } );
    testFunctions.push_back( [r0]() { testTractionPatchTest( gradientEnhanced, 2e-4, -1e-4, r0 ); } );
    for ( const auto* type : { &displacement, &gradientEnhanced } )
      testFunctions.push_back( [type, r0]() { testLinearFieldsAreReproduced( *type, 0.01, r0 ); } );
  }
  for ( const auto* type : { &displacement, &gradientEnhanced } ) {
    testFunctions.push_back( [type]() { testRingVolumeAndLumpedMass( *type ); } );
    testFunctions.push_back( [type]() { testTangent( *type ); } );
    testFunctions.push_back( [type]() { testRefusals( *type ); } );
    testFunctions.push_back( [type]() { testTheParticleSaysItsSection( *type ); } );
  }
  testFunctions.push_back( []() { testStiffnessBound( displacement, false ); } );
  testFunctions.push_back( []() { testStiffnessBound( gradientEnhanced, true ); } );
  testFunctions.push_back( testThickWalledCylinder );

  executeTestsAndCollectExceptions( testFunctions );
  return 0;
}
