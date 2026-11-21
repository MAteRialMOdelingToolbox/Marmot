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
#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotElementProperty.h"
#include <Eigen/Dense>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

/**
 * @brief Helpers for the tests of finite elements.
 */
namespace Marmot::Testing {

  /**
   * @brief The node coordinates of a unit square or cube element, in the Abaqus node order: the corner nodes, then
   * the midside nodes, each at \f$(\xi + 1) / 2\f$ of its natural coordinates.
   * @param nDim The spatial dimension, 2 or 3.
   * @param nNodes The number of nodes: 4 or 8 in 2D, 8 or 20 in 3D.
   * @param lengthOfLastAxis The extent along the last axis, 1 for the square or cube, less for a slab.
   * @return The coordinates, node after node.
   */
  inline std::vector< double > unitElementCoordinates( int nDim, int nNodes, double lengthOfLastAxis )
  {
    std::vector< std::vector< double > > naturalCoordinates;
    if ( nDim == 2 ) {
      naturalCoordinates = { { -1, -1 }, { 1, -1 }, { 1, 1 }, { -1, 1 } };
      if ( nNodes == 8 )
        naturalCoordinates.insert( naturalCoordinates.end(), { { 0, -1 }, { 1, 0 }, { 0, 1 }, { -1, 0 } } );
    }
    else {
      naturalCoordinates = { { -1, -1, -1 },
                             { 1, -1, -1 },
                             { 1, 1, -1 },
                             { -1, 1, -1 },
                             { -1, -1, 1 },
                             { 1, -1, 1 },
                             { 1, 1, 1 },
                             { -1, 1, 1 } };
      if ( nNodes == 20 )
        naturalCoordinates.insert( naturalCoordinates.end(),
                                   { { 0, -1, -1 },
                                     { 1, 0, -1 },
                                     { 0, 1, -1 },
                                     { -1, 0, -1 },
                                     { 0, -1, 1 },
                                     { 1, 0, 1 },
                                     { 0, 1, 1 },
                                     { -1, 0, 1 },
                                     { -1, -1, 0 },
                                     { 1, -1, 0 },
                                     { 1, 1, 0 },
                                     { -1, 1, 0 } } );
    }
    if ( static_cast< int >( naturalCoordinates.size() ) != nNodes )
      throw std::invalid_argument( "unitElementCoordinates: no unit element with this number of nodes." );

    std::vector< double > coordinates;
    for ( const auto& xi : naturalCoordinates )
      for ( int d = 0; d < nDim; d++ )
        coordinates.push_back( 0.5 * ( xi[d] + 1.0 ) * ( d == nDim - 1 ? lengthOfLastAxis : 1.0 ) );
    return coordinates;
  }

  /**
   * @brief An element of the registry, with its coordinates, properties and state, initialized and ready to compute.
   */
  struct ElementUnderTest {

    /// The element.
    std::unique_ptr< Marmot::Discretizations::Elements::Element > element;
    /// Its node coordinates.
    std::vector< double > coordinates;
    /// Its element properties.
    std::vector< double > properties;
    /// Its material properties.
    std::vector< double > materialProperties;
    /// Its state variables.
    std::vector< double > stateVars;

    /**
     * @brief Create the element, assign everything and initialize the element and its material.
     * @param type The name of the element in the registry.
     * @param coordinates_ The node coordinates.
     * @param properties_ The element properties, e.g. the thickness of a plane element; may be empty.
     * @param material The name of the material.
     * @param materialProperties_ The material properties.
     */
    ElementUnderTest( const std::string&    type,
                      std::vector< double > coordinates_,
                      std::vector< double > properties_,
                      const std::string&    material,
                      std::vector< double > materialProperties_ )
      : element( Marmot::Registry::ElementFactory::createElement( type, 1 ) ),
        coordinates( std::move( coordinates_ ) ),
        properties( std::move( properties_ ) ),
        materialProperties( std::move( materialProperties_ ) )
    {
      element->assignNodeCoordinates( coordinates.data() );
      if ( !properties.empty() )
        element->assignProperty(
          Marmot::Discretizations::ElementProperties( properties.data(), static_cast< int >( properties.size() ) ) );
      element->assignProperty(
        Marmot::Discretizations::MaterialSection( material,
                                                  materialProperties.data(),
                                                  static_cast< int >( materialProperties.size() ) ) );
      stateVars.assign( element->getNumberOfRequiredStateVars(), 0.0 );
      element->assignStateVars( stateVars.data(), static_cast< int >( stateVars.size() ) );
      element->initializeYourself();
      element->setInitialConditions( Marmot::Discretizations::Elements::Element::MarmotMaterialInitialization,
                                     nullptr );
    }
  };

  /**
   * @brief The stable time increment \f$2 / \omega_{\max}\f$ of an element on its own, at rest: the largest
   * eigenvalue of \f$\mathbf{M}^{-1}\mathbf{K}\f$ of its displacement block, with the tangent computeKernels() returns
   * and the lumped mass computeLumpedInertia() assembles.
   * @param element The element; its state is evaluated at rest (its state variables change).
   * @param nDisplacementDof The size of the displacement block, the leading degrees of freedom of the element.
   * @return The stable time increment.
   */
  inline double stableTimeIncrementOfTheElement( Marmot::Discretizations::Elements::Element& element,
                                                 int                                         nDisplacementDof )
  {
    const int             nDof = element.getNDofPerElement();
    const Eigen::VectorXd zero = Eigen::VectorXd::Zero( nDof );
    Eigen::VectorXd       P    = Eigen::VectorXd::Zero( nDof );
    Eigen::MatrixXd       K    = Eigen::MatrixXd::Zero( nDof, nDof );
    element.computeKernels( zero.data(), zero.data(), P.data(), K.data(), 0.0, 1.0 );

    Eigen::VectorXd M = Eigen::VectorXd::Zero( nDof );
    element.computeLumpedInertia( M.data() );

    const Eigen::VectorXd scaling = M.head( nDisplacementDof ).cwiseInverse().cwiseSqrt();
    const Eigen::MatrixXd KUU     = K.topLeftCorner( nDisplacementDof, nDisplacementDof );
    const Eigen::MatrixXd scaled  = scaling.asDiagonal() * ( 0.5 * ( KUU + KUU.transpose() ) ) * scaling.asDiagonal();
    const Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > solver( scaled, Eigen::EigenvaluesOnly );
    return 2.0 / std::sqrt( solver.eigenvalues().maxCoeff() );
  }

} // namespace Marmot::Testing
