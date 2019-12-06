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
#include <cstring>

namespace triqs::h5 {

  TRIQS_SPECIALIZE_HDF5_SCHEME2(std::string, string);

  /**
  * \brief Write a string  into an hdf5 file
  * \param f The h5 file or group
  * \param name The name of the hdf5 array in the file/group where the stack will be stored
  * \param value The string
  */
  void h5_write(group g, std::string const &name, std::string const &value);

  inline void h5_write(group g, std::string const &name, const char *s) { h5_write(g, name, std::string{s}); }

  /**
  * \brief Read a string from an hdf5 file
  * \param f The h5 file or group
  * \param name The name of the hdf5 array in the file/group where the stack will be stored
  * \param value The string to fill
  */
  void h5_read(group g, std::string const &name, std::string &value);

  inline void h5_read(group g, std::string const &name, char *s) = delete;

} // namespace triqs::h5
