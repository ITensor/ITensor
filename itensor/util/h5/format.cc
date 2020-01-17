#include <string>
#include "./group.hpp"
#include "./format.hpp"

namespace h5 {

// Code left
std::string group::read_hdf5_format() const {
  std::string s;
  h5_read_attribute(id, "TRIQS_HDF5_data_format", s);
  return s;
}

void group::assert_hdf5_format_as_string(const char *tag_expected, bool ignore_if_absent) const {
  auto tag_file = read_hdf5_format();
  if (ignore_if_absent and tag_file.empty()) return;
  if (tag_file != tag_expected)
    throw std::runtime_error("h5_read : mismatch of the tag TRIQS_HDF5_data_format tag in the h5 group : found " + tag_file + " while I expected "
                             + tag_expected);
}

} //namespace h5
