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
#include "./base.hpp"

namespace triqs {
  namespace h5 {

    //------------------------------------------------

    // dataspace from lengths and strides. Correct for the complex. strides must be >0
    dataspace dataspace_from_LS(int R, bool is_complex, hsize_t const *Ltot, hsize_t const *L, hsize_t const *S, hsize_t const *offset) {
      int rank = R + (is_complex ? 1 : 0);
      hsize_t totdimsf[rank], dimsf[rank], stridesf[rank], offsetf[rank]; // dataset dimensions
      for (size_t u = 0; u < R; ++u) {
        offsetf[u]  = (offset ? offset[u] : 0);
        dimsf[u]    = L[u];
        totdimsf[u] = Ltot[u];
        stridesf[u] = S[u];
      }
      if (is_complex) {
        offsetf[rank - 1]  = 0;
        dimsf[rank - 1]    = 2;
        totdimsf[rank - 1] = 2;
        stridesf[rank - 1] = 1;
      }

      dataspace ds = H5Screate_simple(rank, totdimsf, NULL);
      if (!ds.is_valid()) TRIQS_RUNTIME_ERROR << "Cannot create the dataset";

      herr_t err = H5Sselect_hyperslab(ds, H5S_SELECT_SET, offsetf, stridesf, dimsf, NULL);
      if (err < 0) TRIQS_RUNTIME_ERROR << "Cannot set hyperslab";

      return ds;
    }

    /****************** Write string attribute *********************************************/

    void h5_write_attribute(hid_t id, std::string const &name, std::string const &value) {

      datatype strdatatype = H5Tcopy(H5T_C_S1);
      auto status          = H5Tset_size(strdatatype, value.size() + 1);
      // auto status = H5Tset_size(strdatatype, H5T_VARIABLE);
      if (status < 0) TRIQS_RUNTIME_ERROR << "Internal error in H5Tset_size";

      dataspace space = H5Screate(H5S_SCALAR);

      attribute attr = H5Acreate2(id, name.c_str(), strdatatype, space, H5P_DEFAULT, H5P_DEFAULT);
      if (!attr.is_valid()) TRIQS_RUNTIME_ERROR << "Cannot create the attribute " << name;

      status = H5Awrite(attr, strdatatype, (void *)(value.c_str()));
      if (status < 0) TRIQS_RUNTIME_ERROR << "Cannot write the attribute " << name;
    }

    /****************** Read string attribute *********************************************/

    /// Return the attribute name of obj, and "" if the attribute does not exist.
    void h5_read_attribute(hid_t id, std::string const &name, std::string &s) {
      s = "";

      // if the attribute is not present, return 0
      if (H5LTfind_attribute(id, name.c_str()) == 0) return; // not present

      attribute attr = H5Aopen(id, name.c_str(), H5P_DEFAULT);
      if (!attr.is_valid()) TRIQS_RUNTIME_ERROR << "Cannot open the attribute " << name;

      dataspace space = H5Aget_space(attr);

      int rank = H5Sget_simple_extent_ndims(space);
      if (rank != 0) TRIQS_RUNTIME_ERROR << "Reading a string attribute and got rank !=0";

      datatype strdatatype = H5Aget_type(attr);

      std::vector<char> buf(H5Aget_storage_size(attr) + 1, 0x00);
      auto err = H5Aread(attr, strdatatype, (void *)(&buf[0]));
      if (err < 0) TRIQS_RUNTIME_ERROR << "Cannot read the attribute " << name;

      s.append(&(buf.front()));
    }
  } // namespace h5
} // namespace triqs
