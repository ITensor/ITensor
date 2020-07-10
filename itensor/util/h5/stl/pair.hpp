/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2018 by N. Wentzell, O. Parcollet
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
#include <utility>
#include "../group.hpp"
#include "./string.hpp"

namespace h5 {

  template <typename T1, typename T2>
  struct hdf5_format_impl<std::pair<T1, T2>> {
    static std::string invoke() { return "PythonTupleWrap"; }
  };

  /**
   * Write Pair of T1 and T2 as a subgroup with numbers
   */
  template <typename T1, typename T2>
  void h5_write(group f, std::string const &name, std::pair<T1, T2> const &p) {
    auto gr = f.create_group(name);
    gr.write_hdf5_format(p);
    h5_write(gr, "0", p.first);
    h5_write(gr, "1", p.second);
  }

  /**
   * Read Pair of T1 and T2 from group
   */
  template <typename T1, typename T2>
  void h5_read(group f, std::string const &name, std::pair<T1, T2> &p) {
    auto gr = f.open_group(name);
    if (gr.get_all_subgroup_dataset_names().size() != 2)
      throw std::runtime_error("ERROR in std::pair h5_read: Incompatible number of group elements");
    h5_read(gr, "0", p.first);
    h5_read(gr, "1", p.second);
  }
} // namespace h5
