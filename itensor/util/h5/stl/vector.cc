#include "./vector.hpp"
#include <hdf5.h>
#include <hdf5_hl.h>
#include <array>
#include <cstring>

#include "./string.hpp"

namespace h5 {

  // FIXME : vector of vector string unused.
  // REMOVE CODE but KEEP IT in the history at final cleaning

  //------------------  to_char_buf ------------------------------

  // copy to the buffer, with each string having the same length
  char_buf to_char_buf(std::vector<std::string> const &v) {

    size_t s = 0;
    for (auto &x : v) s = std::max(s, x.size() + 1);
    auto len = v_t{v.size(), s};

    // copy to the buffer
    std::vector<char> buf;
    buf.resize(v.size() * s, 0x00);
    size_t i = 0;
    for (auto &x : v) {
      strcpy(&buf[i * s], x.c_str());
      ++i;
    }

    return {buf, len};
  }

  // copy to the buffer, with each string having the same length
  char_buf to_char_buf(std::vector<std::vector<std::string>> const &v) {

    size_t s = 0, lv = 0;
    for (auto &v1 : v) {
      lv = std::max(lv, v1.size());
      for (auto &x : v1) s = std::max(s, x.size() + 1);
    }
    auto len = v_t{v.size(), lv, s};

    // copy to the buffer
    std::vector<char> buf;
    buf.resize(v.size() * lv * s, 0x00);
    for (unsigned i = 0, k = 0; i < v.size(); i++)
      for (unsigned j = 0; j < lv; j++, k++) {
        if (j < v[i].size()) strcpy(&buf[k * s], v[i][j].c_str());
      }

    return {buf, len};
  }

  //------------------- from_char_buf-----------------------------

  void from_char_buf(char_buf const &cb, std::vector<std::string> &v) {
    v.clear();
    v.resize(cb.lengths[0]);
    auto len_string = cb.lengths[1];

    size_t i = 0;
    for (auto &x : v) {
      x = "";
      x.append(&cb.buffer[i * len_string]);
      ++i;
    }
  }

  //--------

  void from_char_buf(char_buf const &cb, std::vector<std::vector<std::string>> &v) {
    v.clear();
    v.resize(cb.lengths[0]);
    auto inner_vec_size = cb.lengths[1];
    auto len_string     = cb.lengths[2];

    long k = 0;
    for (auto &v_inner : v) {
      for (unsigned j = 0; j < inner_vec_size; ++j, ++k) {
        std::string s = "";
        s.append(&cb.buffer[k * len_string]);
        if (!s.empty()) v_inner.push_back(s);
      }
    }
  }

  // -----------   WRITE  ------------

  //  void h5_write(group g, std::string const &name, std::vector<std::string> const &v) { return h5_write(g, name, to_char_buf(v)); }
  //  void h5_write(group g, std::string const &name, std::vector<std::vector<std::string>> const &v) { return h5_write(g, name, to_char_buf(v)); }

  // -----------   WRITE  ATTRIBUTE ------------

  void h5_write_attribute(hid_t id, std::string const &name, std::vector<std::string> const &v) { h5_write_attribute(id, name, to_char_buf(v)); }

  //void h5_write_attribute(hid_t id, std::string const &name, std::vector<std::vector<std::string>> const &v) {
  //h5_write_attribute(id, name, to_char_buf(v));
  //}

  // -----------  READ  ------------

  /*
  void h5_read(group g, std::string const &name, std::vector<std::string> &v) {
    char_buf cb;
    h5_read(g, name, cb);
    from_char_buf(cb, v);
  }

  void h5_read(group g, std::string const &name, std::vector<std::vector<std::string>> &v) {
    char_buf cb;
    h5_read(g, name, cb);
    from_char_buf(cb, v);
  }
*/
  // -----------   READ  ATTRIBUTE ------------

  void h5_read_attribute(hid_t id, std::string const &name, std::vector<std::string> &v) {
    char_buf cb;
    h5_read_attribute(id, name, cb);
    from_char_buf(cb, v);
  }

  //void h5_read_attribute(hid_t id, std::string const &name, std::vector<std::vector<std::string>> &v) {
  //char_buf cb;
  //h5_read_attribute(id, name, cb);
  //from_char_buf(cb, v);
  //}

} // namespace h5
