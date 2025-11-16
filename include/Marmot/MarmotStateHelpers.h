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
#include "Fastor/Fastor.h"
#include "Marmot/MarmotUtils.h"
#include <Eigen/Dense>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

/**
 * @struct StateMapper
 * @brief Template struct to map raw double pointers and sizes to specific views.
 * Specializations of this struct should provide a static `map` method that
 * takes a `double*` and a `std::size_t` size, returning the desired view type.
 */
template < class View, int... Args >
struct StateMapper;

/**
 * @struct StateMapper<double>
 * @brief Specialization of StateMapper for single double values.
 * This allows mapping a raw double pointer and size to a single double reference.
 */
template <>
struct StateMapper< double& > {
  static double& map( double* ptr, std::size_t n )
  {
    if ( n != 1 )
      throw std::runtime_error( "Size mismatch for double." );
    return *ptr;
  }
};

/**
 * @struct StateMapper<std::span<double>>
 * @brief Specialization of StateMapper for std::span<double>.
 * This allows mapping a raw double pointer and size to a mutable span.
 */
template <>
struct StateMapper< std::span< double > > {
  static std::span< double > map( double* ptr, std::size_t n ) { return { ptr, n }; }
};

/**
 * @struct StateMapper<std::span<const double>>
 * @brief Specialization of StateMapper for std::span<const double>.
 * This allows mapping a raw double pointer and size to a read-only span.
 */
template <>
struct StateMapper< std::span< const double > > {
  static std::span< const double > map( double* ptr, std::size_t n ) { return { ptr, n }; }
};

/**
 * @struct StateMapper<Fastor::Tensor<double,3,3>>
 * @brief Specialization of StateMapper for Fastor::Tensor<double,3,3>.
 * This allows mapping a raw double pointer and size to a Fastor 3x3 tensor.
 */
template <>
struct StateMapper< Fastor::Tensor< double, 3, 3 >& > {
  static Fastor::Tensor< double, 3, 3 >& map( double* ptr, std::size_t n )
  {
    if ( n != 9 )
      throw std::runtime_error( "Size mismatch for Fastor::Tensor<double,3,3>." );
    return *reinterpret_cast< Fastor::Tensor< double, 3, 3 >* >( ptr );
  }
};

/**
 * @struct StateMapper<Eigen::Map<Eigen::MatrixXd>>
 * @brief Specialization of StateMapper for Eigen::Map<Eigen::MatrixXd>.
 * This allows mapping a raw double pointer and size to an Eigen dynamic matrix map.
 */
template <>
struct StateMapper< Eigen::Map< Eigen::MatrixXd > > {
  template < class... Args >
  static Eigen::Map< Eigen::MatrixXd > map( double* ptr, std::size_t n, int Rows, int Cols )
  {
    // read size form Args
    if ( int( n ) != Rows * Cols )
      throw std::runtime_error( "Size mismatch for Eigen::Map." );
    return Eigen::Map< Eigen::MatrixXd >( ptr, Rows, n / Rows );
  }
};

class MarmotStateLayoutDynamic {
public:
  struct VarInfo {
    std::string name;
    std::size_t size;
    std::size_t offset;
  };

  void add( std::string name, std::size_t size )
  {
    if ( finalized )
      throw std::runtime_error( "Cannot add variables after finalize()." );

    if ( name_to_index.contains( name ) )
      throw std::runtime_error( "Duplicate state variable: " + name );

    name_to_index[name] = variables.size();
    variables.push_back( { std::move( name ), size, 0 } );
  }

  void finalize()
  {
    if ( finalized )
      throw std::runtime_error( "State layout already finalized." );
    std::size_t offset = 0;
    for ( auto& v : variables ) {
      v.offset = offset;
      offset += v.size;
    }
    total_sz  = offset;
    finalized = true;
  }

  const VarInfo& getInfo( const std::string& name ) const
  {
    auto it = name_to_index.find( name );
    if ( it == name_to_index.end() )
      throw std::runtime_error( "Unknown state variable: " + name );
    return variables[it->second];
  }

  std::pair< std::size_t, std::size_t > get( const std::string& name ) const
  {
    const auto& v = getInfo( name );
    return { v.offset, v.size };
  }

  double* getPtr( double* base, const std::string& name ) const
  {
    const auto& v = getInfo( name );
    return base + v.offset;
  }

  std::span< double > getSpan( double* base, const std::string& name ) const
  {
    const auto& v = getInfo( name );
    return { base + v.offset, v.size };
  }

  StateView getStateView( double* base, const std::string& name ) const
  {
    const auto& v = getInfo( name );
    return StateView{ base + v.offset, static_cast< int >( v.size ) };
  }

  template < class View, class... Args >
  View getAs( double* base, const std::string& name, Args&&... args ) const
  {
    const auto& v = getInfo( name );
    return StateMapper< View >::map( base + v.offset, v.size, std::forward< Args >( args )... );
  }

  int totalSize() const { return total_sz; }

  bool isFinalized() const { return finalized; }

private:
  std::vector< VarInfo >                         variables;
  std::unordered_map< std::string, std::size_t > name_to_index;

  bool finalized = false;
  int  total_sz  = 0;
};
