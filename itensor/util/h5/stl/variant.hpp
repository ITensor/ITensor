#pragma once

// FIXME : PB HERE
// LEAVE in TRIQS for the moment.

#include <triqs/utility/variant.hpp>
#include "../group.hpp"
#include "./string.hpp"
#include "../generic.hpp"

namespace h5 {

  // forward
  hid_t H5Dget_type(hid_t);

  template <typename... T>
  struct hdf5_format_impl<std::variant<T...>> {
    static std::string invoke() = delete;
  };

  /**
   */
  template <typename... T>
  void h5_write(h5::group gr, std::string const &name, std::variant<T...> const &v) {
    visit([&](auto const &x) { h5_write(gr, name, x); }, v);
  }

  template <typename VT, typename U, typename... T>
  void h5_read_variant_helper(VT &v, datatype dt, h5::group gr, std::string const &name) {
    if (types::is_same_type(typeid(U), dt)) {
      v = VT{h5::h5_read<U>(gr, name)};
      return;
    }
    if constexpr (sizeof...(T) > 0)
      h5_read_variant_helper<VT, T...>(v, dt, gr, name);
    else
      H5_ERROR << " Error in h5_read: std::variant<...> not compatible with TRIQS_HDF5_data_format \n";
  }

  /**
   * Read variant from the h5
   */
  template <typename... T>
  void h5_read(h5::group gr, std::string const &name, std::variant<T...> &v) {
    // name is a group --> triqs object
    // assume for the moment, name is a dataset.
    dataset ds  = gr.open_dataset(name);
    datatype dt = H5Dget_type(ds);
    h5_read_variant_helper<std::variant<T...>, T...>(v, dt, gr, name);
  }

} // namespace h5
