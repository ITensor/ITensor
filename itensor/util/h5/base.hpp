/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011-2014 by O. Parcollet
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once
#include "./base_public.hpp"
#include <triqs/utility/mini_vector.hpp>
#include <triqs/utility/is_complex.hpp>
#include <type_traits>
#include <hdf5.h>
#include <hdf5_hl.h>
/// This header is only for inclusion in the xxx.cpp files of the h5 lib,
/// not in general code.

namespace triqs {
  namespace h5 {

    using utility::mini_vector;

    // conversion of C type to HDF5 native
    inline hid_t native_type_from_C(char) { return H5T_NATIVE_CHAR; }
    inline hid_t native_type_from_C(signed char) { return H5T_NATIVE_SCHAR; }
    inline hid_t native_type_from_C(unsigned char) { return H5T_NATIVE_UCHAR; }
    inline hid_t native_type_from_C(short) { return H5T_NATIVE_SHORT; }
    inline hid_t native_type_from_C(unsigned short) { return H5T_NATIVE_USHORT; }
    inline hid_t native_type_from_C(int) { return H5T_NATIVE_INT; }
    inline hid_t native_type_from_C(unsigned) { return H5T_NATIVE_UINT; }
    inline hid_t native_type_from_C(long) { return H5T_NATIVE_LONG; }
    inline hid_t native_type_from_C(unsigned long) { return H5T_NATIVE_ULONG; }
    inline hid_t native_type_from_C(long long) { return H5T_NATIVE_LLONG; }
    inline hid_t native_type_from_C(unsigned long long) { return H5T_NATIVE_ULLONG; }
    inline hid_t native_type_from_C(float) { return H5T_NATIVE_FLOAT; }
    inline hid_t native_type_from_C(double) { return H5T_NATIVE_DOUBLE; }
    inline hid_t native_type_from_C(long double) { return H5T_NATIVE_LDOUBLE; }
    inline hid_t native_type_from_C(bool) { return H5T_NATIVE_SCHAR; }
    inline hid_t native_type_from_C(std::string) { return H5T_C_S1; }
    inline hid_t native_type_from_C(std::complex<double>) { return native_type_from_C(double()); }

    // conversion of C type to HDF5 native
    // We need to discuss which type we use in the file
    // NATIVE is only one possibility, like IEEE, etc...
    inline hid_t h5_type_from_C(char) { return H5T_NATIVE_CHAR; }
    inline hid_t h5_type_from_C(signed char) { return H5T_NATIVE_SCHAR; }
    inline hid_t h5_type_from_C(unsigned char) { return H5T_NATIVE_UCHAR; }
    inline hid_t h5_type_from_C(short) { return H5T_NATIVE_SHORT; }
    inline hid_t h5_type_from_C(unsigned short) { return H5T_NATIVE_USHORT; }
    inline hid_t h5_type_from_C(int) { return H5T_NATIVE_INT; }
    inline hid_t h5_type_from_C(unsigned) { return H5T_NATIVE_UINT; }
    inline hid_t h5_type_from_C(long) { return H5T_NATIVE_LONG; }
    inline hid_t h5_type_from_C(unsigned long) { return H5T_NATIVE_ULONG; }
    inline hid_t h5_type_from_C(long long) { return H5T_NATIVE_LLONG; }
    inline hid_t h5_type_from_C(unsigned long long) { return H5T_NATIVE_ULLONG; }
    inline hid_t h5_type_from_C(float) { return H5T_NATIVE_FLOAT; }
    inline hid_t h5_type_from_C(double) { return H5T_NATIVE_DOUBLE; }
    inline hid_t h5_type_from_C(long double) { return H5T_NATIVE_LDOUBLE; }
    inline hid_t h5_type_from_C(bool) { return H5T_NATIVE_SCHAR; }
    inline hid_t h5_type_from_C(std::string) { return H5T_C_S1; }
    inline hid_t h5_type_from_C(std::complex<double>) { return h5_type_from_C(double()); }

    //------------- compute the data type (removing complex) ----------------------------------

    // in memory
    template <typename T> hid_t data_type_memory() { return native_type_from_C(T{}); }

    // the type of data to put in the file_or_group
    template <typename T> hid_t data_type_file() { return h5_type_from_C(T{}); }

    //------------- compute void * pointer to the data ----------------------------------
    // 2 cases : complex or not. Complex are reinterpreted according to doc, as N+1 dim double array
    template <typename S>
    std14::enable_if_t<triqs::is_complex<S>::value, std14::conditional_t<std::is_const<S>::value, const void *, void *>> get_data_ptr(S *p) {
      using T = std14::conditional_t<std::is_const<S>::value, const typename S::value_type, typename S::value_type>;
      return reinterpret_cast<T *>(p);
    }

    template <typename S>
    std14::enable_if_t<!triqs::is_complex<S>::value, std14::conditional_t<std::is_const<S>::value, const void *, void *>> get_data_ptr(S *p) {
      return p;
    }

    // dataspace from lengths and strides. Correct for the complex. strides must be >0
    // used in array and vectors
    // implemented in base.cpp
    dataspace dataspace_from_LS(int R, bool is_complex, hsize_t const *Ltot, hsize_t const *L, hsize_t const *S, hsize_t const *offset = NULL);

  } // namespace h5
} // namespace triqs
