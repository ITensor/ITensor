#include "./file.hpp"
#include "./base.hpp"

namespace triqs {
  namespace h5 {

    namespace {
      unsigned h5_char_to_int(char fl) {
        switch (fl) {
          case 'r': return H5F_ACC_RDONLY;
          case 'w': return H5F_ACC_TRUNC;
          case 'a': return H5F_ACC_RDWR;
        }
        TRIQS_RUNTIME_ERROR << " Internal error";
      }
    } // namespace

    file::file(const char *name, char flags) : file(name, h5_char_to_int(flags)) {}

    file::file(const char *name, unsigned flags) {

      if (flags == H5F_ACC_RDONLY) {
        id = H5Fopen(name, flags, H5P_DEFAULT);
        if (id < 0) TRIQS_RUNTIME_ERROR << "HDF5 : cannot open file " << name;
        return;
      }

      if (flags == H5F_ACC_RDWR) {
        id = H5Fopen(name, flags, H5P_DEFAULT);
        if (id < 0) {
          id = H5Fcreate(name, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
          if (id < 0) TRIQS_RUNTIME_ERROR << "HDF5 : cannot open file " << name;
        }
        return;
      }

      if (flags == H5F_ACC_TRUNC) {
        id = H5Fcreate(name, flags, H5P_DEFAULT, H5P_DEFAULT);
        if (id < 0) TRIQS_RUNTIME_ERROR << "HDF5 : cannot create file " << name;
        return;
      }

      if (flags == H5F_ACC_EXCL) {
        id = H5Fcreate(name, flags, H5P_DEFAULT, H5P_DEFAULT);
        if (id < 0) TRIQS_RUNTIME_ERROR << "HDF5 : cannot create file " << name << ". Does it exists ?";
        return;
      }

      TRIQS_RUNTIME_ERROR << "HDF5 file opening : flag not recognized";
    }

    //---------------------------------------------

    file::file(hid_t id_) : h5_object(h5_object(id_)) {}

    file::file(h5_object obj) : h5_object(std::move(obj)) {
      if (!H5Iis_valid(this->id)) TRIQS_RUNTIME_ERROR << "Invalid input in h5::file constructor from id";
      if (H5Iget_type(this->id) != H5I_FILE) TRIQS_RUNTIME_ERROR << "h5::file constructor must take the id of a file ";
    }

    //---------------------------------------------

    std::string file::name() const { // same function as for group
      char _n[1];
      ssize_t size = H5Fget_name(id, _n, 1); // first call, get the size only
      std::vector<char> buf(size + 1, 0x00);
      H5Fget_name(id, buf.data(), 1); // now get the name
      std::string res = "";
      res.append(&(buf.front()));
      return res;
    }
  } // namespace h5
} // namespace triqs
