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
#include "./scalar.hpp"
#include "./base.hpp"

namespace triqs {
  namespace h5 {

    namespace {

      // S must be scalar
      template <typename S> void h5_write_impl(group g, std::string const &key, S a) {

        g.unlink_key_if_exists(key);

        dataspace space = H5Screate(H5S_SCALAR);
        auto ds         = g.create_dataset(key, data_type_file<S>(), space);

        auto err = H5Dwrite(ds, data_type_memory<S>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, (void *)(&a));
        if (err < 0) TRIQS_RUNTIME_ERROR << "Error writing the scalar dataset " << key << " in the group" << g.name();
      }

      //-------------------------------------------------------------

      // attribute
      template <typename S> void h5_write_attribute_mpl(hid_t id, std::string const &name, S a) {

        if (H5LTfind_attribute(id, name.c_str()) != 0)
          TRIQS_RUNTIME_ERROR << "The attribute " << name << " is already present. Can not overwrite"; // not present
        dataspace space = H5Screate(H5S_SCALAR);

        attribute attr = H5Acreate2(id, name.c_str(), data_type_file<S>(), space, H5P_DEFAULT, H5P_DEFAULT);
        if (!attr.is_valid()) TRIQS_RUNTIME_ERROR << "Cannot create the attribute " << name;

        auto status = H5Awrite(attr, data_type_memory<S>(), (void *)(&a));
        if (status < 0) TRIQS_RUNTIME_ERROR << "Cannot write the attribute " << name;
      }

      //-------------------------------------------------------------

      template <typename S> void h5_read_impl(group g, std::string const &name, S &A) {

        dataset ds        = g.open_dataset(name);
        dataspace d_space = H5Dget_space(ds);

        // check that rank is 0, it is a scalar.
        int rank = H5Sget_simple_extent_ndims(d_space);
        if (rank != 0)
          TRIQS_RUNTIME_ERROR << "triqs::array::h5::read. Rank mismatch : expecting a scalar (rank =0)"
                              << " while the array stored in the hdf5 file has rank = " << rank;

        herr_t err = H5Dread(ds, data_type_memory<S>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, (void *)(&A));
        if (err < 0) TRIQS_RUNTIME_ERROR << "Error reading the scalar dataset " << name << " in the group" << g.name();
      }

      //-------------------------------------------------------------

      template <typename S> void h5_read_attribute_mpl(hid_t id, std::string const &name, S &a) {

        attribute attr = H5Aopen(id, name.c_str(), H5P_DEFAULT);
        if (!attr.is_valid()) TRIQS_RUNTIME_ERROR << "Cannot open the attribute " << name;

        dataspace space = H5Aget_space(attr);
        int rank        = H5Sget_simple_extent_ndims(space);
        if (rank != 0) TRIQS_RUNTIME_ERROR << "Reading a scalar attribute and got rank !=0";

        auto eq = H5Tequal(H5Aget_type(attr), data_type_memory<S>());
        if (eq < 0) TRIQS_RUNTIME_ERROR << "Type comparison failure in reading attribute";
        if (eq == 0) TRIQS_RUNTIME_ERROR << "Type mismatch in reading attribute";

        auto err = H5Aread(attr, data_type_memory<S>(), (void *)(&a));
        if (err < 0) TRIQS_RUNTIME_ERROR << "Cannot read the attribute " << name;
      }
    } // namespace

    // template <> void h5_write_impl(group g, std::string const &key, std::complex<double> a);
    // template <> void h5_read_impl(group g, std::string const &key, std::complex<double> &a);

    void h5_write(group g, std::string const &name, int const &x) { h5_write_impl(g, name, x); }
    void h5_write(group g, std::string const &name, long const &x) { h5_write_impl(g, name, x); }
    void h5_write(group g, std::string const &name, long long const &x) { h5_write_impl(g, name, x); }
    void h5_write(group g, std::string const &name, unsigned long long const &x) { h5_write_impl(g, name, x); }
    void h5_write(group g, std::string const &name, size_t const &x) { h5_write_impl(g, name, x); }
    void h5_write(group g, std::string const &name, char const &x) { h5_write_impl(g, name, x); }
    //void h5_write(group g, std::string const &name, bool const &x) { h5_write_impl(g, name, x); }
    void h5_write(group g, std::string const &name, double const &x) { h5_write_impl(g, name, x); }
    //void h5_write(group g, std::string const &name, std::complex<double> const &x) { h5_write_impl(g, name, x); }
    void h5_write(group g, std::string const &name, std::complex<double> const &x) {
      triqs::h5::group gr = g.create_group(name);
      h5_write(gr, "r", x.real());
      h5_write(gr, "i", x.imag());
    }

    void h5_read(group g, std::string const &name, int &x) { h5_read_impl(g, name, x); }
    void h5_read(group g, std::string const &name, long &x) { h5_read_impl(g, name, x); }
    void h5_read(group g, std::string const &name, long long &x) { h5_read_impl(g, name, x); }
    void h5_read(group g, std::string const &name, unsigned long long &x) { h5_read_impl(g, name, x); }
    void h5_read(group g, std::string const &name, size_t &x) { h5_read_impl(g, name, x); }
    void h5_read(group g, std::string const &name, char &x) { h5_read_impl(g, name, x); }
    //void h5_read(group g, std::string const &name, bool &x) { h5_read_impl(g, name, x); }
    void h5_read(group g, std::string const &name, double &x) { h5_read_impl(g, name, x); }
    //void h5_read(group g, std::string const &name, std::complex<double> &x) { h5_read_impl(g, name, x); }

    // special case : complex
    void h5_read(group g, std::string const &name, std::complex<double> &x) {
      triqs::h5::group gr = g.open_group(name);
      double r, i;
      h5_read(gr, "r", r);
      h5_read(gr, "i", i);
      x = std::complex<double>{r, i};
    }

    datatype bool_datatype() {
      datatype bool_enum_t = H5Tenum_create(H5T_NATIVE_CHAR);
      char val;
      H5Tenum_insert(bool_enum_t, "FALSE", (val = 0, &val));
      H5Tenum_insert(bool_enum_t, "TRUE", (val = 1, &val));
      return bool_enum_t;
    }

    // special case : bool. Same code as the general case
    void h5_write(group g, std::string const &key, bool const &x) {
      g.unlink_key_if_exists(key);
      datatype dt     = bool_datatype();
      dataspace space = H5Screate(H5S_SCALAR);
      auto ds         = g.create_dataset(key, dt, space);
      auto err        = H5Dwrite(ds, dt, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void *)(&x));
      if (err < 0) TRIQS_RUNTIME_ERROR << "Error writing the scalar dataset " << key << " in the group" << g.name();
    }

    // Seems to work with the default function (read as integer), but it is cleaner to implement a proper read function
    // with the enum.
    void h5_read(group g, std::string const &name, bool &x) {
      dataset ds        = g.open_dataset(name);
      dataspace d_space = H5Dget_space(ds);
      // check that rank is 0, it is a scalar.
      int rank = H5Sget_simple_extent_ndims(d_space);
      if (rank != 0)
        TRIQS_RUNTIME_ERROR << "triqs::array::h5::read. Rank mismatch : expecting a scalar (rank =0)"
                            << " while the array stored in the hdf5 file has rank = " << rank;
      datatype dt = bool_datatype();
      herr_t err  = H5Dread(ds, dt, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void *)(&x));
      if (err < 0) TRIQS_RUNTIME_ERROR << "Error reading the scalar dataset " << name << " in the group" << g.name();
    }

    void h5_write_attribute(hid_t id, std::string const &name, int x) { h5_write_attribute_mpl(id, name, x); }
    void h5_write_attribute(hid_t id, std::string const &name, long x) { h5_write_attribute_mpl(id, name, x); }
    void h5_write_attribute(hid_t id, std::string const &name, unsigned long x) { h5_write_attribute_mpl(id, name, x); }
    void h5_write_attribute(hid_t id, std::string const &name, long long x) { h5_write_attribute_mpl(id, name, x); }
    void h5_write_attribute(hid_t id, std::string const &name, unsigned long long x) { h5_write_attribute_mpl(id, name, x); }
    void h5_write_attribute(hid_t id, std::string const &name, double x) { h5_write_attribute_mpl(id, name, x); }

    void h5_read_attribute(hid_t id, std::string const &name, int &x) { h5_read_attribute_mpl(id, name, x); }
    void h5_read_attribute(hid_t id, std::string const &name, long &x) { h5_read_attribute_mpl(id, name, x); }
    void h5_read_attribute(hid_t id, std::string const &name, unsigned long &x) { h5_read_attribute_mpl(id, name, x); }
    void h5_read_attribute(hid_t id, std::string const &name, long long &x) { h5_read_attribute_mpl(id, name, x); }
    void h5_read_attribute(hid_t id, std::string const &name, unsigned long long &x) { h5_read_attribute_mpl(id, name, x); }
    void h5_read_attribute(hid_t id, std::string const &name, double &x) { h5_read_attribute_mpl(id, name, x); }
  } // namespace h5
} // namespace triqs
