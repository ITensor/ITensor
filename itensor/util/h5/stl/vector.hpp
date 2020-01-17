#pragma once
#include <vector>
#include <complex>
#include "../group.hpp"
#include "./string.hpp"
#include "../scalar.hpp"

namespace h5 {

  namespace array_interface {

    template <typename T>
    h5_array_view h5_array_view_from_vector(std::vector<T> const &v) {
      h5_array_view res{hdf5_type<T>(), (void *)v.data(), 1, is_complex_v<std::decay_t<T>>};
      res.slab.count[0] = v.size();
      res.L_tot[0]      = v.size();
      return res;
    }

  } // namespace array_interface

  // ----------------------------------------------------------------------------
  // FIXME : CLEAN THIS
  // Special case of vector < string >

  H5_SPECIALIZE_FORMAT2(std::vector<std::string>, vector<string>);

  template <typename T>
  struct hdf5_format_impl<std::vector<T>> {
    static std::string invoke() { return "PythonListWrap"; }
    //static std::string invoke() { return "std::vector<" + get_hdf5_format<T>() + ">"; }
  };

  // ----------------------------------------------------------------------------
  // details for string case
  char_buf to_char_buf(std::vector<std::string> const &v);
  void from_char_buf(char_buf const &cb, std::vector<std::string> &v);

  // ----------------------------------------------------------------------------

  /**
   * Writes std::vector into HDF5
   *
   * Format 
   *    * If T is a simple type (int, double, complex), it is a 1d array.
   *    * If T is std::string, it is a 2d array of char of dimensions (length of vector, max length of strings)
   *    * Otherwise, it opens a subgroup and writes each element as 0,1,2,3 ... in the subgroup
   *
   * @tparam T
   * @param g HDF5 group
   * @param name Name of the object in the HDF5 file
   * @param v Vector to save in the file
   */
  template <typename T>
  void h5_write(group g, std::string const &name, std::vector<T> const &v) {

    if constexpr (std::is_arithmetic_v<T> or is_complex_v<T>) {

      array_interface::write(g, name, array_interface::h5_array_view_from_vector(v), true);

    } else if constexpr (std::is_same_v<T, std::string>) {

      h5_write(g, name, to_char_buf(v));

    } else { // generic type

      auto gr = g.create_group(name);
      h5_write_attribute(gr, "format", get_hdf5_format(v));
      for (int i = 0; i < v.size(); ++i) h5_write(gr, std::to_string(i), v[i]);
    }
  }

  // ----------------------------------------------------------------------------

  /**
   * Reads std::vector from HDF5
   *
   * Format 
   *    * If T is a simple type (int, double, complex), it is a 1d array.
   *    * If T is std::string, it is a 2d array of char of dimensions (length of vector, max length of strings)
   *    * Otherwise, it opens a subgroup and writes each element as 0,1,2,3 ... in the subgroup
   *
   * @tparam T
   * @param g HDF5 group
   * @param name Name of the object in the HDF5 file
   * @param v Vector to read from the file
   */
  template <typename T>
  void h5_read(group g, std::string name, std::vector<T> &v) {

    if constexpr (std::is_arithmetic_v<T> or is_complex_v<T>) {

      auto lt = array_interface::get_h5_lengths_type(g, name);
      if (lt.rank() != 1 + is_complex_v<T>) throw make_runtime_error("h5 : reading a vector and I got an array of rank", lt.rank());
      v.resize(lt.lengths[0]);
      array_interface::read(g, name, array_interface::h5_array_view_from_vector(v), lt);

    } else if constexpr (std::is_same_v<T, std::string>) {

      char_buf cb;
      h5_read(g, name, cb);
      from_char_buf(cb, v);

    } else { // generic type

      auto g2 = g.open_group(name);
      v.resize(g2.get_all_dataset_names().size() + g2.get_all_subgroup_names().size());
      for (int i = 0; i < v.size(); ++i) { h5_read(g2, std::to_string(i), v[i]); }
    }
  }

  //void h5_write(group f, std::string const &name, std::vector<std::vector<std::string>> const &V);
  //void h5_read(group f, std::string const &name, std::vector<std::vector<std::string>> &V);

  // void h5_write_attribute(hid_t ob, std::string const &name, std::vector<std::vector<std::string>> const &V);
  // void h5_read_attribute(hid_t ob, std::string const &name, std::vector<std::vector<std::string>> &V);

  // void h5_write_attribute(hid_t ob, std::string const &name, std::vector<std::string> const &V);
  // void h5_read_attribute(hid_t ob, std::string const &name, std::vector<std::string> &V);

} // namespace h5
