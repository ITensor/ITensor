#pragma once
#include <type_traits>

namespace h5 {

  /**
   * A generic read
   * 
   * @tparam T 
   * @param g  HDF5 group
   * @param key Name of the element
   * @return The value read from the file
   */
  template <typename T>
  T h5_read(group g, std::string const &key) {
    if constexpr (std::is_default_constructible_v<T>) {
      T x;
      h5_read(g, key, x);
      return x;
    } else {
      return T::h5_read_construct(g, key);
    }
  }

  /**
   * A generic attribute read.
   * 
   * @tparam T 
   * @param obj The object to which the attribute is attached.
   * @param key Name of the element
   * @return The attribute name of obj, and "" if the attribute does not exist.
   */
  template <typename T>
  T h5_read_attribute(h5_object obj, std::string const &key) {
    using h5::h5_read_attribute;
    T x;
    h5_read_attribute(obj, key, x);
    return x;
  }

  /**
   * Try a read
   * 
   * @tparam T 
   * @param g  HDF5 group
   * @param key Name of the element
   * @param x Parameter to read.
   */
  template <typename T>
  inline bool h5_try_read(group g, std::string key, T &x) {
    if (g.has_key(key)) {
      h5_read(g, key, x);
      return true;
    }
    return false;
  }

} // namespace h5
