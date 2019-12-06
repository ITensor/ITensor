/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011-2017 by O. Parcollet
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
#include "./base.hpp"
#include "../arrays.hpp"

#define CHECK_OR_THROW(Cond, Mess)                                                                                                                   \
  if (!(Cond)) TRIQS_RUNTIME_ERROR << "Error in h5 (de)serialization " << Mess;

namespace triqs {
  namespace h5 {

#if H5_VERSION_GE(1, 8, 9)

    using h5_serialization_char_t = unsigned char;

    template <typename T> arrays::array<h5_serialization_char_t, 1> serialize(T const &x) {

      proplist fapl = H5Pcreate(H5P_FILE_ACCESS);
      CHECK_OR_THROW((fapl >= 0), "creating fapl");

      auto err = H5Pset_fapl_core(fapl, (size_t)(64 * 1024), false);
      CHECK_OR_THROW((err >= 0), "setting core file driver in fapl.");

      h5::file f = H5Fcreate("serialization", 0, H5P_DEFAULT, fapl);
      CHECK_OR_THROW((f.is_valid()), "created core file");

      auto gr = triqs::h5::group(f);
      h5_write(gr, "object", x);

      err = H5Fflush(f, H5F_SCOPE_GLOBAL);
      CHECK_OR_THROW((err >= 0), "flushed core file.");

      ssize_t image_len = H5Fget_file_image(f, NULL, (size_t)0);
      CHECK_OR_THROW((image_len > 0), "got image file size");

      arrays::array<h5_serialization_char_t, 1> buf(image_len);
      buf() = 0;

      ssize_t bytes_read = H5Fget_file_image(f, (void *)buf.data_start(), (size_t)image_len);
      CHECK_OR_THROW(bytes_read == image_len, "wrote file into image buffer");

      return buf;
    }

    // -----------------------------

    template <typename T> T deserialize(arrays::array_const_view<h5_serialization_char_t, 1> buf) {

      proplist fapl = H5Pcreate(H5P_FILE_ACCESS);
      CHECK_OR_THROW((fapl >= 0), "creating fapl");

      auto err = H5Pset_fapl_core(fapl, (size_t)(64 * 1024), false);
      CHECK_OR_THROW((err >= 0), "setting core file driver in fapl.");

      err = H5Pset_file_image(fapl, (void *)buf.data_start(), first_dim(buf));
      CHECK_OR_THROW((err >= 0), "set file image in fapl.");

      h5::file f = H5Fopen("serialization", H5F_ACC_RDONLY, fapl);
      CHECK_OR_THROW((f.is_valid()), "opened received file image file");

      auto gr = triqs::h5::group(f);
      regular_type_if_view_else_type_t<T> res;
      h5_read(gr, "object", res);

      return res;
    }
  }

#endif

#undef CHECK_OR_THROW
} // namespace triqs
