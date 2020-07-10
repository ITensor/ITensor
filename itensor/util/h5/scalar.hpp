#pragma once
#include "./group.hpp"
#include "./array_interface.hpp"
namespace h5 {

  namespace array_interface {

    template <typename S>
    h5_array_view h5_array_view_from_scalar(S &&s) {
      return {hdf5_type<std::decay_t<S>>(), (void *)(&s), 0, h5::is_complex_v<std::decay_t<S>>};
    }
  } // namespace array_interface

  template <typename T>
  void h5_write(group g, std::string const &name, T const &x) H5_REQUIRES(std::is_arithmetic_v<T> or is_complex_v<T>) {
    array_interface::write(g, name, array_interface::h5_array_view_from_scalar(x), false);
  }

  template <typename T>
  void h5_read(group g, std::string const &name, T &x) H5_REQUIRES(std::is_arithmetic_v<T> or is_complex_v<T>) {
    array_interface::read(g, name, array_interface::h5_array_view_from_scalar(x), array_interface::get_h5_lengths_type(g, name));
  }

  template <typename T>
  void h5_write_attribute(hid_t id, std::string const &name, T const &x) H5_REQUIRES(std::is_arithmetic_v<T> or is_complex_v<T>) {
    array_interface::write_attribute(id, name, array_interface::h5_array_view_from_scalar(x));
  }

  template <typename T>
  void h5_read_attribute(hid_t id, std::string const &name, T &x) H5_REQUIRES(std::is_arithmetic_v<T> or is_complex_v<T>) {
    array_interface::read_attribute(id, name, array_interface::h5_array_view_from_scalar(x));
  }

} // namespace h5
