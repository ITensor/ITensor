/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011-2013 by O. Parcollet
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
#include <complex>
namespace triqs {
  namespace h5 {

    // Issue several types are *implicitly* convertible to bool
    // it could be confusing. Better to use int in hdf5 files.
    void h5_write(group g, std::string const &name, int const &x);
    void h5_write(group g, std::string const &name, long const &x);
    void h5_write(group g, std::string const &name, long long const &x);
    void h5_write(group g, std::string const &name, unsigned long long const &x);
    void h5_write(group g, std::string const &name, size_t const &x);
    void h5_write(group g, std::string const &name, bool const &x);
    void h5_write(group g, std::string const &name, char const &x);
    void h5_write(group g, std::string const &name, double const &x);
    void h5_write(group g, std::string const &name, std::complex<double> const &x);

    void h5_read(group g, std::string const &name, int &x);
    void h5_read(group g, std::string const &name, long &x);
    void h5_read(group g, std::string const &name, long long &x);
    void h5_read(group g, std::string const &name, unsigned long long &x);
    void h5_read(group g, std::string const &name, size_t &x);
    void h5_read(group g, std::string const &name, bool &x);
    void h5_read(group g, std::string const &name, char &x);
    void h5_read(group g, std::string const &name, double &x);
    void h5_read(group g, std::string const &name, std::complex<double> &x);

    // attribute
    void h5_write_attribute(hid_t id, std::string const &name, int x);
    void h5_write_attribute(hid_t id, std::string const &name, long x);
    void h5_write_attribute(hid_t id, std::string const &name, unsigned long x);
    void h5_write_attribute(hid_t id, std::string const &name, long long const x);
    void h5_write_attribute(hid_t id, std::string const &name, unsigned long long x);
    void h5_write_attribute(hid_t id, std::string const &name, double x);

    void h5_read_attribute(hid_t id, std::string const &name, int &x);
    void h5_read_attribute(hid_t id, std::string const &name, long &x);
    void h5_read_attribute(hid_t id, std::string const &name, unsigned long &x);
    void h5_read_attribute(hid_t id, std::string const &name, long long &x);
    void h5_read_attribute(hid_t id, std::string const &name, unsigned long long &x);
    void h5_read_attribute(hid_t id, std::string const &name, double &x);
  } // namespace h5
} // namespace triqs
