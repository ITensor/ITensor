#pragma once
#include "./base.hpp"

namespace h5 {

  template <typename T>
  std::vector<unsigned char> serialize(T const &x) {
    memory_file f;
    h5_write(f, "X", x);
    return f.as_buffer();
  }

  // -----------------------------

  template <typename T>
  T deserialize(std::vector<unsigned char> const &buf) {
    return h5_read<T>(memory_file f{buf}, "object");
  }
} // namespace h5
