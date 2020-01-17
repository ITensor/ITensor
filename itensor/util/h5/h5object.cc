#include <type_traits>

#include <H5Ipublic.h>
#include <H5Fpublic.h>
#include <H5Gpublic.h>
#include <H5Ppublic.h>

// FIXME
//static_assert(std::is_same_v<hid_t, int64_t>, "Configuration error in HDF5. Check version.");

#include "./h5object.hpp"
#include <vector>
#include <algorithm>

namespace h5 {

  //static_assert(std::is_same<::hid_t, hid_t>::value, "Internal error");
  static_assert(std::is_same<::hsize_t, hsize_t>::value, "Internal error");

  // specializations for all basic types
  // clang-format off
  template <> hid_t hdf5_type<char>          (){return  H5T_NATIVE_CHAR;}
  template <> hid_t hdf5_type<signed char>   (){return  H5T_NATIVE_SCHAR;}
  template <> hid_t hdf5_type<unsigned char> (){return  H5T_NATIVE_UCHAR;}

  template <> hid_t hdf5_type<short>     (){return  H5T_NATIVE_SHORT;}
  template <> hid_t hdf5_type<int>       (){return  H5T_NATIVE_INT;}
  template <> hid_t hdf5_type<long>      (){return  H5T_NATIVE_LONG;}
  template <> hid_t hdf5_type<long long> (){return  H5T_NATIVE_LLONG;}

  template <> hid_t hdf5_type<unsigned short>     (){return  H5T_NATIVE_USHORT;}
  template <> hid_t hdf5_type<unsigned int>       (){return  H5T_NATIVE_UINT;}
  template <> hid_t hdf5_type<unsigned long>      (){return  H5T_NATIVE_ULONG;}
  template <> hid_t hdf5_type<unsigned long long> (){return  H5T_NATIVE_ULLONG;}

  template <> hid_t hdf5_type<float>       (){return  H5T_NATIVE_FLOAT;}
  template <> hid_t hdf5_type<double>      (){return  H5T_NATIVE_DOUBLE;}
  template <> hid_t hdf5_type<long double> (){return  H5T_NATIVE_LDOUBLE;}

  template <> hid_t hdf5_type<std::complex<double>>      (){return  H5T_NATIVE_DOUBLE;}
  template <> hid_t hdf5_type<std::complex<long double>> (){return  H5T_NATIVE_LDOUBLE;}
  template <> hid_t hdf5_type<std::complex<float>>       (){return  H5T_NATIVE_FLOAT;}
  // clang-format on

  // bool. Use a lambda to initialize it.
  template <>
  hid_t hdf5_type<bool>() {
    datatype bool_enum_h5type = H5Tenum_create(H5T_NATIVE_CHAR);
    char val;
    H5Tenum_insert(bool_enum_h5type, "FALSE", (val = 0, &val));
    H5Tenum_insert(bool_enum_h5type, "TRUE", (val = 1, &val));
    return bool_enum_h5type;
  }

  // -----------------------  name  ---------------------------

  struct h5_name_t {
    datatype hdf5_type; // type in hdf5
    std::string name;   // name of the type
  };

  //---------

  std::vector<h5_name_t> h5_name_table;

  //----------

  static void init_h5_name_table() {
    h5_name_table = std::vector<h5_name_t>{
       {hdf5_type<char>(), H5_AS_STRING(char)},
       {hdf5_type<signed char>(), H5_AS_STRING(signed char)},
       {hdf5_type<unsigned char>(), H5_AS_STRING(unsigned char)},
       {hdf5_type<short>(), H5_AS_STRING(short)},
       {hdf5_type<unsigned short>(), H5_AS_STRING(unsigned short)},
       {hdf5_type<int>(), H5_AS_STRING(int)},
       {hdf5_type<unsigned int>(), H5_AS_STRING(unsigned int)},
       {hdf5_type<long>(), H5_AS_STRING(long)},
       {hdf5_type<unsigned long>(), H5_AS_STRING(unsigned long)},
       {hdf5_type<long long>(), H5_AS_STRING(long long)},
       {hdf5_type<unsigned long long>(), H5_AS_STRING(unsigned long long)},
       {hdf5_type<float>(), H5_AS_STRING(float)},
       {hdf5_type<double>(), H5_AS_STRING(double)},
       {hdf5_type<long double>(), H5_AS_STRING(long double)},
       {hdf5_type<std::complex<float>>(), H5_AS_STRING(std::complex<float>)},
       {hdf5_type<std::complex<double>>(), H5_AS_STRING(std::complex<double>)},
       {hdf5_type<std::complex<long double>>(), H5_AS_STRING(std::complex<long double>)} //
    };
  }

  //--------

  std::string get_name_of_h5_type(datatype t) {

    if (h5_name_table.empty()) init_h5_name_table();
    auto _end = h5_name_table.end();
    auto pos  = std::find_if(h5_name_table.begin(), _end, [t](auto const &x) { return H5Tequal(x.hdf5_type, t) > 0; });
    if (pos == _end) throw std::logic_error("HDF5/Python : impossible error");
    return pos->name;
  }

  // -----------------------   Reference counting ---------------------------

  // xdecref, xincref manipulate the the ref, but ignore invalid (incl. 0) id.
  //  like XINC_REF and XDEC_REF in python
  inline void xdecref(hid_t id) {
    if (H5Iis_valid(id)) H5Idec_ref(id);
  }

  inline void xincref(hid_t id) {
    if (H5Iis_valid(id)) H5Iinc_ref(id);
  }

  // -----------------------  h5_object  ---------------------------

  h5_object::h5_object(h5_object const &x) : id(x.id) { xincref(id); } // a new copy, a new ref.

  // make an h5_object when the id is now owned (simply inc. the ref).
  h5_object h5_object::from_borrowed(hid_t id) {
    xincref(id);
    return h5_object(id);
  }

  h5_object &h5_object::operator=(h5_object &&x) noexcept { //steals the ref, after properly decref its own.
    xdecref(id);
    id   = x.id;
    x.id = 0;
    return *this;
  }

  void h5_object::close() {
    xdecref(id);
    id = 0;
  } // e.g. to close a file explicitely.

  bool h5_object::is_valid() const { return H5Iis_valid(id) == 1; }

} // namespace h5
