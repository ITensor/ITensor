#pragma once
#include <complex>
#include <vector>
#include <sstream>
#include "./macros.hpp"

namespace h5 {

  // We copy this from hdf5.h, and static_assert its validity in the cpp
  // in order to completely isolate our header from the hdf5 headers
  // Hence complex installation paths to hdf5 are only needed in the cpp file,
  // not by the users of the library.
  using hid_t   = int64_t;
  using hsize_t = unsigned long long;
  using v_t     = std::vector<hsize_t>;

  // Correspondance T -> hdf5 type
  template <typename T>
  hid_t hdf5_type();

  // impl trait to detect complex numbers
  template <typename T>
  struct _is_complex : std::false_type {};

  template <typename T>
  struct _is_complex<std::complex<T>> : std::true_type {};

  template <typename T>
  constexpr bool is_complex_v = _is_complex<T>::value;

  // impl
  template <typename... T>
  std::runtime_error make_runtime_error(T const &... x) {
    std::stringstream fs;
    (fs << ... << x);
    return std::runtime_error{fs.str()};
  }

  /*
   * A handle to the a general HDF5 object
   *
   * HDF5 uses a reference counting system, similar to python.
   * h5_object handles the proper reference couting (similar to pybind11::object e.g.)
   * using a RAII pattern (hence exception safety).
   */
  class h5_object {

    protected:
    hid_t id = 0; //NOLINT Ok, I want a protected variable ...

    public:
    /// make an h5_object from a simple borrowed ref (simply inc. the ref).
    static h5_object from_borrowed(hid_t id);

    /// Constructor from an owned id (or 0). It steals (take ownership) of the reference.
    h5_object(hid_t id = 0) : id(id) {}

    /// A new ref. No deep copy.
    h5_object(h5_object const &x);

    /// Steals the reference
    h5_object(h5_object &&x) noexcept : id(x.id) { x.id = 0; }

    /// Copy the reference and incref
    h5_object &operator=(h5_object const &x) {
      operator=(h5_object(x));
      return *this;
    }

    /// Steals the ref.
    h5_object &operator=(h5_object &&x) noexcept;

    ///
    ~h5_object() { close(); }

    /// Release the HDF5 handle and reset the object to default state (id =0).
    void close();

    /// cast operator to use it in the C function as its id
    operator hid_t() const { return id; }

    /// Ensure the id is valid (by H5Iis_valid).
    [[nodiscard]] bool is_valid() const;
  };

  //-----------------------------

  // simple but useful aliases. It does NOT check the h5 type of the object.
  // FIXME : derive and make a check ??
  using dataset   = h5_object;
  using datatype  = h5_object;
  using dataspace = h5_object;
  using proplist  = h5_object;
  using attribute = h5_object;

  // ------------------------------

  // A function to get the name of a datatype in clear (for error messages)
  std::string get_name_of_h5_type(datatype ty);

} // namespace h5
