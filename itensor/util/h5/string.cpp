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
#include "./string.hpp"
#include "./base.hpp"
namespace triqs {

  namespace h5 {

    void h5_write(group g, std::string const &name, std::string const &value) {

      datatype strdatatype = H5Tcopy(H5T_C_S1);
      // auto status = H5Tset_size (strdatatype, H5T_VARIABLE);
      //auto status = H5Tset_size(strdatatype, value.size() + 1);
      H5Tset_size(strdatatype, value.size() + 1);

      dataspace space = H5Screate(H5S_SCALAR);
      dataset ds      = g.create_dataset(name, strdatatype, space);

      auto err = H5Dwrite(ds, strdatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void *)(value.c_str()));
      if (err < 0) TRIQS_RUNTIME_ERROR << "Error writing the string named" << name << " in the group" << g.name();
    }

    // ------------

    void h5_read(group g, std::string const &name, std::string &value) {
      dataset ds            = g.open_dataset(name);
      h5::dataspace d_space = H5Dget_space(ds);
      int rank              = H5Sget_simple_extent_ndims(d_space);
      if (rank != 0) TRIQS_RUNTIME_ERROR << "Reading a string and got rank !=0";
      size_t size = H5Dget_storage_size(ds);

      datatype strdatatype = H5Tcopy(H5T_C_S1);
      H5Tset_size(strdatatype, size);
      //auto status = H5Tset_size(strdatatype, size);
      // auto status = H5Tset_size (strdatatype, H5T_VARIABLE);

      std::vector<char> buf(size + 1, 0x00);
      auto err = H5Dread(ds, strdatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &buf[0]);
      if (err < 0) TRIQS_RUNTIME_ERROR << "Error reading the string named" << name << " in the group" << g.name();

      value = "";
      value.append(&(buf.front()));
    }
  } // namespace h5
} // namespace triqs
