#pragma once
#include <map>
#include "../group.hpp"
#include "./string.hpp"

namespace h5 {

  template <typename T>
  struct hdf5_format_impl<std::map<std::string, T>> {
    static std::string invoke() { return "PythonDictWrap"; } //"std::map<string," + hdf5_format_impl<T>::invoke() + ">"; }
  };

  /**
   * Map of string and T as a subgroup with key_names
   */
  template <typename T>
  void h5_write(group f, std::string const &name, std::map<std::string, T> const &M) {
    auto gr = f.create_group(name);
    gr.write_hdf5_format(M);
    for (auto &pvp : M) h5_write(gr, pvp.first, pvp.second);
  }

  /**
   * Map of string and T
   */
  template <typename T>
  void h5_read(group f, std::string const &name, std::map<std::string, T> &M) {
    auto gr = f.open_group(name);
    M.clear();
    for (auto const &x : gr.get_all_subgroup_dataset_names()) {
      T value;
      h5_read(gr, x, value);
      M.emplace(x, std::move(value));
    }
  }

} // namespace h5
