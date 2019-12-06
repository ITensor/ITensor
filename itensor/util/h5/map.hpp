/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011 by O. Parcollet
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
#include "group.hpp"
#include "string.hpp"
#include <map>

namespace triqs::h5 {

  template <typename T> struct hdf5_scheme_impl<std::map<std::string, T>> {
    static std::string invoke() { return "PythonDictWrap"; } //"std::map<string," + hdf5_scheme_impl<T>::invoke() + ">"; }
  };

  /**
   * Map of string and T as a subgroup with key_names
   */
  template <typename T> void h5_write(group f, std::string const &name, std::map<std::string, T> const &M) {
    auto gr = f.create_group(name);
    gr.write_hdf5_scheme(M);
    for (auto &pvp : M) h5_write(gr, pvp.first, pvp.second);
  }

  /**
   * Map of string and T
   */
  template <typename T> void h5_read(group f, std::string const &name, std::map<std::string, T> &M) {
    auto gr = f.open_group(name);
    M.clear();
    for (auto const &x : gr.get_all_subgroup_dataset_names()) {
      T value;
      h5_read(gr, x, value);
      M.emplace(x, std::move(value));
    }
  }

} // namespace triqs::h5
