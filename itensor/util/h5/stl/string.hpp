#pragma once
#include "../group.hpp"
#include <string>

namespace h5 {

  H5_SPECIALIZE_FORMAT2(std::string, string);

  /**
   * Write a string  into an hdf5 file
   *
   * Format : Fixed size string
   *
   * @tparam T
   * @param g HDF5 group
   * @param name Name of the object in the HDF5 file
   * @param s String to be saved.
   */
  void h5_write(group g, std::string const &name, std::string const &s, bool force_utf8 = false);

  /**
   * Write a string  into an hdf5 file
   *
   * Format : Fixed size string
   *
   * @tparam T
   * @param g HDF5 group
   * @param name Name of the object in the HDF5 file
   * @param s String to be saved.
   */
  inline void h5_write(group g, std::string const &name, const char *s, bool force_utf8 = false) { h5_write(g, name, std::string{s},force_utf8); }

  /**
   * Read a string from an hdf5 file
   *
   * @param f The h5 file or group
   * @param name The name of the hdf5 array in the file/group where the stack will be stored
   * @param value The string to read into
   */
  void h5_read(group g, std::string const &name, std::string &value);

  // Explicitly forbidden.
  inline void h5_read(group g, std::string const &name, char *s) = delete;

  /**
   * Write a string attribute
   *
   * @param f The h5 file or group
   * @param name The name of the hdf5 array in the file/group where the stack will be stored
   * @param s The string.
  */
  void h5_write_attribute(hid_t id, std::string const &name, std::string const &s, bool force_utf8 = false);

  /**
   * Write a string attribute
   *
   * @param f The h5 file or group
   * @param name The name of the hdf5 array in the file/group where the stack will be stored
   * @param s The string.
  */
  inline void h5_write_attribute(hid_t id, std::string const &name, const char *s) { h5_write_attribute(id, name, std::string{s}); }

  /**
  * Read a string attribute from id.
  * 
  * @param id  The object to which the attribute is attached
  * @param name The name of the attribute
  * @param value The string to fill
  */
  void h5_read_attribute(hid_t id, std::string const &name, std::string &s);

  // forbidden
  inline void h5_read_attribute(hid_t id, std::string const &name, char *s) = delete;

  // ---------------------   char_buf -----------------------

  // char_buf contains an n dimensional array of strings as fixed size strings, flatten in a 1d array of char.
  // the last dimension is the max length of the strings + 1, because of the ending 0 in C !
  struct char_buf {
    std::vector<char> buffer;
    v_t lengths;

    // the string datatype
    [[nodiscard]] datatype dtype() const;

    // the dataspace (without last dim, which is the string).
    [[nodiscard]] dataspace dspace() const;
  };

  // read/write for char_buf
  void h5_write(group g, std::string const &name, char_buf const &cb);
  void h5_write_attribute(hid_t id, std::string const &name, char_buf const &cb);
  void h5_read(group g, std::string const &name, char_buf &_cb);
  void h5_read_attribute(hid_t id, std::string const &name, char_buf &_cb);

} // namespace h5
