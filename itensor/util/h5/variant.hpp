
/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2017 by O. Parcollet
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
#include "./base.hpp"
#include "./group.hpp"
#include "./string.hpp"
#include "./generic.hpp"
#include <triqs/utility/variant.hpp>

namespace triqs::h5 {

  template <typename... T> struct hdf5_scheme_impl<std::variant<T...>> { static std::string invoke() = delete; };

  /**
   */
  template <typename... T> void h5_write(h5::group gr, std::string const &name, std::variant<T...> const &v) {
    visit([&](auto const &x) { h5_write(gr, name, x); }, v);
  }

  template <typename T> bool h5_is(datatype dt) {
    if constexpr (std::is_same<T, std::string>::value) {
      return H5Tget_class(dt) == H5T_STRING;
    } // H5T_INTEGER, H5T_FLOAT
    else
      return H5Tequal(dt, h5_type_from_C(T{}));
  }

  template <typename VT, typename U, typename... T> void h5_read_variant_helper(VT &v, datatype dt, h5::group gr, std::string const &name) {
    if (h5_is<U>(dt)) {
      v = VT{triqs::h5::h5_read<U>(gr, name)};
      return;
    }
    if constexpr (sizeof...(T) > 0)
      h5_read_variant_helper<VT, T...>(v, dt, gr, name);
    else
      TRIQS_RUNTIME_ERROR << " Error in h5_read: std::variant<...> not compatible with TRIQS_HDF5_data_scheme \n";
  }

  /**
   * Read vairant from the h5
   */
  template <typename... T> void h5_read(h5::group gr, std::string const &name, std::variant<T...> &v) {
    // name is a group --> triqs object
    // assume for the moment, name is a dataset.
    dataset ds  = gr.open_dataset(name);
    datatype dt = H5Dget_type(ds);
    h5_read_variant_helper<std::variant<T...>, T...>(v, dt, gr, name);
  }

} // namespace triqs::h5
