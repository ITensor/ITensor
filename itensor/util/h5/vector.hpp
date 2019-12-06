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
#include "./group.hpp"
#include "./string.hpp"
#include <vector>
#include <complex>

namespace triqs::h5 {

  TRIQS_SPECIALIZE_HDF5_SCHEME2(std::vector<std::string>, vector<string>);

  template <typename T> struct hdf5_scheme_impl<std::vector<T>> {
    static std::string invoke() { return "PythonListWrap"; } //std::vector<" + hdf5_scheme_impl<T>::invoke() + ">"; }
    //static std::string invoke() { return "std::vector<" + hdf5_scheme_impl<T>::invoke() + ">"; }
  };

  void h5_write(group f, std::string const &name, std::vector<std::string> const &V);
  void h5_read(group f, std::string const &name, std::vector<std::string> &V);

  void h5_write_attribute(hid_t ob, std::string const &name, std::vector<std::vector<std::string>> const &V);
  void h5_read_attribute(hid_t ob, std::string const &name, std::vector<std::vector<std::string>> &V);

  //void h5_write_attribute (hid_t ob, std::string const & name, std::vector<std::string> const & V);
  //void h5_read_attribute (hid_t ob, std::string const & name, std::vector<std::string> & V);

  void h5_write(group f, std::string const &name, std::vector<int> const &V);
  void h5_write(group f, std::string const &name, std::vector<double> const &V);
  void h5_write(group f, std::string const &name, std::vector<std::complex<double>> const &V);
  void h5_write(group f, std::string const &name, std::vector<unsigned long> const &V);
  void h5_write(group f, std::string const &name, std::vector<unsigned long long> const &V);

  void h5_read(group f, std::string const &name, std::vector<int> &V);
  void h5_read(group f, std::string const &name, std::vector<double> &V);
  void h5_read(group f, std::string const &name, std::vector<std::complex<double>> &V);
  void h5_read(group f, std::string const &name, std::vector<unsigned long> &V);
  void h5_read(group f, std::string const &name, std::vector<unsigned long long> &V);

  // vector of non basic types

  /**
   * Vector of T as a subgroup with numbers
   */
  template <typename T> void h5_write(group f, std::string name, std::vector<T> const &v) {
    auto gr = f.create_group(name);
    gr.write_hdf5_scheme(v);
    for (int i = 0; i < v.size(); ++i) h5_write(gr, std::to_string(i), v[i]);
  }

  /**
   * Vector of T
   */
  template <typename T> void h5_read(group f, std::string name, std::vector<T> &v) {
    auto gr = f.open_group(name);
    v.resize(gr.get_all_dataset_names().size() + gr.get_all_subgroup_names().size());
    for (int i = 0; i < v.size(); ++i) { h5_read(gr, std::to_string(i), v[i]); }
  }
} // namespace triqs::h5
