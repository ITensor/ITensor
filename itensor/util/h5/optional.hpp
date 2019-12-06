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
#include "./group.hpp"
#include "./string.hpp"
#include <optional>

namespace triqs::h5 {

  template <typename T> struct hdf5_scheme_impl<std::optional<T>> {
    static std::string invoke() { return hdf5_scheme_impl<T>::invoke(); }
  };

  /**
   * Optional : write if the value is set.
   */
  template <typename T> void h5_write(h5::group gr, std::string const &name, std::optional<T> const &v) {
    if (bool(v)) h5_write(gr, name, *v);
  }

  /**
   * Read optional from the h5
   */
  template <typename T> void h5_read(h5::group gr, std::string name, std::optional<T> &v) {
    v.reset();
    if (gr.has_key(name)) v.emplace(h5_read<T>(gr, name));
  }
} // namespace triqs::h5
