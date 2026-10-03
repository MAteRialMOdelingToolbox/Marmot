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

#include <Fastor/Fastor.h>
#include <algorithm>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

namespace nb = nanobind;

namespace nanobind::detail {
  template < typename T, size_t... Dims >
  struct type_caster< Fastor::Tensor< T, Dims... > > {
    using TensorType = Fastor::Tensor< T, Dims... >;
    NB_TYPE_CASTER( TensorType, const_name( "numpy.ndarray[" ) + make_caster< T >::Name + const_name( "]" ) );

    bool from_python( handle src, uint8_t flags, cleanup_list* cleanup ) noexcept
    {
      using NdArray = nb::ndarray< T, nb::shape< Dims... >, nb::c_contig, nb::device::cpu >;

      type_caster< NdArray > nd_caster;
      if ( !nd_caster.from_python( src, flags, cleanup ) ) {
        return false;
      }

      NdArray nd = nd_caster.value;
      std::copy( nd.data(), nd.data() + nd.size(), value.data() );
      return true;
    }

    // Not noexcept: the allocation below and the creation of the Python objects may throw,
    // and an exception escaping a noexcept function would call std::terminate.
    static handle from_cpp( const TensorType& src, rv_policy, cleanup_list* cleanup )
    {
      using NdArray = nb::ndarray< nb::numpy, T >;

      size_t shape[] = { Dims... };
      size_t ndim    = sizeof...( Dims );
      size_t size    = 1;
      for ( size_t i = 0; i < ndim; ++i )
        size *= shape[i];

      // The tensor is usually a temporary, so its data is copied into a buffer whose
      // lifetime is tied to the returned array by the capsule below.
      T* data = new T[size];
      std::copy( src.data(), src.data() + size, data );

      nb::capsule owner( data, []( void* p ) noexcept { delete[] (T*)p; } );

      // The copy above already decoupled the array from the C++ object, so the returned
      // array never aliases it and the caller's policy does not apply. rv_policy::reference
      // hands out the capsule-owned buffer without copying it a second time; forwarding the
      // caller's policy instead would make reference_internal fail on the capsule owner.
      return make_caster< NdArray >::from_cpp( NdArray( data, ndim, shape, owner ), rv_policy::reference, cleanup );
    }
  };
} // namespace nanobind::detail
