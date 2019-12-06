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
#pragma once
#include "./base_public.hpp"

namespace triqs {
  namespace h5 {

    /**
  *  \brief A little handler for the file
  */
    class file : public h5_object {

      public:
      /**
   * Open the file name.
   * Flag char can be :
   *   - 'a' H5F_ACC_RDWR
   *   - 'r' H5F_ACC_RDONLY
   *   - 'w' H5F_ACC_TRUNC
   */
      file(const char *name, char flags);

      /**
   * Open the file name.
   * Flags can be :
   *   - H5F_ACC_RDWR
   *   - H5F_ACC_RDONLY
   *   - H5F_ACC_TRUNC
   *   - H5F_ACC_EXCL
   */
      file(const char *name, unsigned flags);

      /// Cf previous constructor
      file(std::string const &name, unsigned flags) : file(name.c_str(), flags) {}

      ///
      file(std::string const &name, char flags) : file(name.c_str(), flags) {}

      /// Internal : from an hdf5 id.
      file(hid_t id);
      file(h5_object obj);

      /// Name of the file
      std::string name() const;
    };
  } // namespace h5
} // namespace triqs
