#pragma once
#include <optional>
#include "../group.hpp"
#include "./string.hpp"

namespace h5 {

  template <typename T>
  struct hdf5_format_impl<std::optional<T>> {
    static std::string invoke() { return hdf5_format_impl<T>::invoke(); }
  };

  /**
   * Optional : write if the value is set.
   */
  template <typename T>
  void h5_write(h5::group gr, std::string const &name, std::optional<T> const &v) {
    if (bool(v)) h5_write(gr, name, *v);
  }

  /**
   * Read optional from the h5
   */
  template <typename T>
  void h5_read(h5::group gr, std::string name, std::optional<T> &v) {
    v.reset();
    if (gr.has_key(name)) v.emplace(h5_read<T>(gr, name));
  }
} // namespace h5
