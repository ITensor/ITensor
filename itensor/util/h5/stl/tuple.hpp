#pragma once
#include <tuple>
#include "../group.hpp"
#include "./string.hpp"

namespace h5 {

  template <typename... T>
  struct hdf5_format_impl<std::tuple<T...>> {
    static std::string invoke() { return "PythonTupleWrap"; }
  };

  namespace details {
    template <typename... T, std::size_t... Is>
    void h5_write_tuple_impl(group gr, std::string const &, std::tuple<T...> const &tpl, std::index_sequence<Is...>) {
      (h5_write(gr, std::to_string(Is), std::get<Is>(tpl)), ...);
    }

    template <typename... T, std::size_t... Is>
    void h5_read_tuple_impl(group gr, std::string const &, std::tuple<T...> &tpl, std::index_sequence<Is...>) {
      if (gr.get_all_subgroup_dataset_names().size() != sizeof...(Is))
        throw std::runtime_error("ERROR in std::tuple h5_read: Tuple size incompatible to number of group elements");
      (h5_read(gr, std::to_string(Is), std::get<Is>(tpl)), ...);
    }

  } // namespace details

  /**
   * Tuple of T... as a subgroup with numbers
   */
  template <typename... T>
  void h5_write(group f, std::string const &name, std::tuple<T...> const &tpl) {
    auto gr = f.create_group(name);
    gr.write_hdf5_format(tpl);
    h5_write_tuple_impl(gr, name, tpl, std::index_sequence_for<T...>{});
  }

  /**
   * Tuple of T...
   */
  template <typename... T>
  void h5_read(group f, std::string const &name, std::tuple<T...> &tpl) {
    auto gr = f.open_group(name);
    h5_read_tuple_impl(gr, name, tpl, std::index_sequence_for<T...>{});
  }
} // namespace h5
