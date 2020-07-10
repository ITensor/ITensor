#pragma once
#include "./macros.hpp"

namespace h5 {

  // a class T has either :
  //  1- a static member hdf5_format -> std::string (or a constexpr char * ?)
  //  2- specializes hdf5_format_impl
  // user function is get_hdf5_format <T>() in all cases.
  // A claass which is not default constructible :
  //  -- 1 : implement static T h5_read_construct(gr, name) : rebuilt  a new T
  //  -- 2 : NOT IMPLEMENTED : if we want to make it non intrusive,
  //  specialize with a struct similarly to hdf5_format_impl
  // to be implemented if needed.

  template <typename T>
  struct hdf5_format_impl {
    static std::string invoke() { return T::hdf5_format(); }
  };

#define H5_SPECIALIZE_FORMAT2(X, Y)                                                                                                                  \
  template <>                                                                                                                                        \
  struct hdf5_format_impl<X> {                                                                                                                       \
    static std::string invoke() { return H5_AS_STRING(Y); }                                                                                          \
  }

#define H5_SPECIALIZE_FORMAT(X) H5_SPECIALIZE_FORMAT2(X, X)

  H5_SPECIALIZE_FORMAT(bool);
  H5_SPECIALIZE_FORMAT(int);
  H5_SPECIALIZE_FORMAT(long);
  H5_SPECIALIZE_FORMAT(long long);
  H5_SPECIALIZE_FORMAT(unsigned int);
  H5_SPECIALIZE_FORMAT(unsigned long);
  H5_SPECIALIZE_FORMAT(unsigned long long);
  H5_SPECIALIZE_FORMAT(float);
  H5_SPECIALIZE_FORMAT(double);
  H5_SPECIALIZE_FORMAT(long double);
  H5_SPECIALIZE_FORMAT2(std::complex<double>, complex);

  template <typename T>
  std::string get_hdf5_format() {
    return hdf5_format_impl<T>::invoke();
  }

  template <typename T>
  std::string get_hdf5_format(T const &) {
    return hdf5_format_impl<T>::invoke();
  }

  // A few helper functions

  /* /// Write the triqs tag*/
  //void write_hdf5_format_as_string(group g, const char *a) { h5_write_attribute(g, "TRIQS_HDF5_data_format", a); }

  ///// Write the triqs tag of the group if it is an object.
  //template <typename T> void write_hdf5_format(group g, T const &) { write_hdf5_format_as_string(g, ::h5::get_hdf5_format<T>().c_str()); }

  ///// Read the triqs tag of the group if it is an object. Returns the empty string "" if attribute is not present
  ////std::string read_hdf5_format() const;

  ///// Asserts that the tag of the group is the same as for T. Throws H5_ERROR if
  //void assert_hdf5_format_as_string(const char *tag_expected, bool ignore_if_absent = false) const;

  ///// Asserts that the tag of the group is the same as for T. Throws H5_ERROR if
  //template <typename T>
  //void assert_hdf5_format(T const &, bool ignore_if_absent = false) const {
  //assert_hdf5_format_as_string(get_hdf5_format<T>().c_str(), ignore_if_absent);
  /*}*/

} // namespace h5
