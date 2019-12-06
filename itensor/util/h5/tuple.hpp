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
#include <tuple>

namespace triqs::h5 {

  template <typename... T> struct hdf5_scheme_impl<std::tuple<T...>> {
    static std::string invoke() { return "PythonTupleWrap"; }
  };

  template <typename... T, std::size_t... Is>
  void h5_write_impl(group f, std::string const &name, std::tuple<T...> const &tpl, std::index_sequence<Is...>) {
    auto gr = f.create_group(name);
    gr.write_hdf5_scheme(tpl);
    (h5_write(gr, std::to_string(Is), std::get<Is>(tpl)), ...);
  }

  /**
   * Tuple of T... as a subgroup with numbers
   */
  template <typename... T> void h5_write(group f, std::string const &name, std::tuple<T...> const &tpl) {
    h5_write_impl(f, name, tpl, std::index_sequence_for<T...>{});
  }

  template <typename... T, std::size_t... Is> void h5_read_impl(group f, std::string const &name, std::tuple<T...> &tpl, std::index_sequence<Is...>) {
    auto gr = f.open_group(name);
    if (gr.get_all_subgroup_dataset_names().size() != sizeof...(Is))
      TRIQS_RUNTIME_ERROR << "ERROR in std::tuple h5_read: Tuple size incompatible to number of group elements";
    (h5_read(gr, std::to_string(Is), std::get<Is>(tpl)), ...);
  }

  /**
   * Tuple of T...
   */
  template <typename... T> void h5_read(group f, std::string const &name, std::tuple<T...> &tpl) {
    h5_read_impl(f, name, tpl, std::index_sequence_for<T...>{});
  }
} // namespace triqs::h5
