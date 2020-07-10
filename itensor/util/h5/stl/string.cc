#include "./string.hpp"
#include <hdf5.h>
#include <hdf5_hl.h>
#include <vector>
#include <numeric>
#include <algorithm>

namespace h5 {

  static datatype str_datatype(long size, bool is_utf8 = false) {
    datatype dt = H5Tcopy(H5T_C_S1);
    auto err    = H5Tset_size(dt, size);
    if (is_utf8) H5Tset_cset(dt, H5T_CSET_UTF8);
    if (err < 0) throw std::runtime_error("Internal error in H5Tset_size");
    return dt;
  }

  static datatype str_datatype(std::string const &s) {
    auto is_nonascii = [](char c) { return static_cast<unsigned char>(c) > 127; };
    bool is_utf8 = std::any_of(begin(s), end(s), is_nonascii);
    return str_datatype(s.size() + 1, is_utf8);
  }

  // ------------------------------------------------------------------

  void h5_write(group g, std::string const &name, std::string const &value, bool force_utf8) {

    datatype dt     = force_utf8 ? str_datatype(value.size()+1,true) : str_datatype(value);
    dataspace space = H5Screate(H5S_SCALAR);

    // FIXME : remove create_dataset
    dataset ds = g.create_dataset(name, dt, space);

    auto err = H5Dwrite(ds, dt, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void *)(value.c_str()));
    if (err < 0) throw std::runtime_error("Error writing the string named" + name + " in the group" + g.name());
  }

  // ------------------------------------------------------------------

  void h5_write_attribute(hid_t id, std::string const &name, std::string const &value, bool force_utf8) {

    datatype dt     = force_utf8 ? str_datatype(value.size()+1,true) : str_datatype(value);
    dataspace space = H5Screate(H5S_SCALAR);

    attribute attr = H5Acreate2(id, name.c_str(), dt, space, H5P_DEFAULT, H5P_DEFAULT);
    if (!attr.is_valid()) throw std::runtime_error("Cannot create the attribute " + name);

    herr_t err = H5Awrite(attr, dt, (void *)(value.c_str()));
    if (err < 0) throw std::runtime_error("Cannot write the attribute " + name);
  }

  // -------------------- Read ----------------------------------------------

  void h5_read(group g, std::string const &name, std::string &value) {
    dataset ds            = g.open_dataset(name);
    h5::dataspace d_space = H5Dget_space(ds);
    int rank              = H5Sget_simple_extent_ndims(d_space);
    if (rank != 0) throw std::runtime_error("Reading a string and got rank !=0");
    size_t size = H5Dget_storage_size(ds);

    datatype dt = H5Dget_type(ds);
    H5_ASSERT(H5Tget_class(dt) == H5T_STRING);

    std::vector<char> buf(size + 1, 0x00);
    auto err = H5Dread(ds, dt, H5S_ALL, H5S_ALL, H5P_DEFAULT, &buf[0]);
    if (err < 0) throw std::runtime_error("Error reading the string named" + name + " in the group" + g.name());

    value = "";
    value.append(&(buf.front()));
  }

  // -------------------- Read ----------------------------------------------

  /// Return the attribute name of obj, and "" if the attribute does not exist.
  void h5_read_attribute(hid_t id, std::string const &name, std::string &s) {
    s = "";

    // if the attribute is not present, return 0
    if (H5LTfind_attribute(id, name.c_str()) == 0) return; // not present

    attribute attr = H5Aopen(id, name.c_str(), H5P_DEFAULT);
    if (!attr.is_valid()) throw std::runtime_error("Cannot open the attribute " + name);

    dataspace space = H5Aget_space(attr);

    int rank = H5Sget_simple_extent_ndims(space);
    if (rank != 0) throw std::runtime_error("Reading a string attribute and got rank !=0");

    datatype strdatatype = H5Aget_type(attr);
    H5_ASSERT(H5Tget_class(strdatatype) == H5T_STRING);

    std::vector<char> buf(H5Aget_storage_size(attr) + 1, 0x00);
    auto err = H5Aread(attr, strdatatype, (void *)(&buf[0]));
    if (err < 0) throw std::runtime_error("Cannot read the attribute " + name);

    s.append(&(buf.front()));
  }

  // -------------------------------------------------------------------
  // the string datatype
  datatype char_buf::dtype() const {
    return str_datatype(lengths.back());
  }

  // the dataspace (without last dim, which is the string).
  dataspace char_buf::dspace() const {
    dataspace ds = H5Screate_simple(lengths.size() - 1, lengths.data(), nullptr); // rank is size of length - 1
    if (!ds.is_valid()) throw make_runtime_error("Cannot create the dataset");
    return ds;
  }

  // -----------   WRITE  ------------

  void h5_write(group g, std::string const &name, char_buf const &cb) {
    auto dt     = cb.dtype();
    auto dspace = cb.dspace();

    h5::dataset ds = g.create_dataset(name, dt, dspace);

    auto err = H5Dwrite(ds, dt, dspace, H5S_ALL, H5P_DEFAULT, (void *)cb.buffer.data());
    if (err < 0) throw make_runtime_error("Error writing the vector<string> ", name, " in the group", g.name());
  }

  // -----------   WRITE  ATTRIBUTE ------------

  void h5_write_attribute(hid_t id, std::string const &name, char_buf const &cb) {
    auto dt     = cb.dtype();
    auto dspace = cb.dspace();

    attribute attr = H5Acreate2(id, name.c_str(), dt, dspace, H5P_DEFAULT, H5P_DEFAULT);
    if (!attr.is_valid()) throw make_runtime_error("Cannot create the attribute ", name);

    herr_t status = H5Awrite(attr, dt, (void *)cb.buffer.data());
    if (status < 0) throw make_runtime_error("Cannot write the attribute ", name);
  }

  // -----------  READ  ------------

  void h5_read(group g, std::string const &name, char_buf &_cb) {
    dataset ds            = g.open_dataset(name);
    h5::dataspace d_space = H5Dget_space(ds);
    datatype ty           = H5Dget_type(ds);

    char_buf cb_out;

    int dim = H5Sget_simple_extent_ndims(d_space);
    cb_out.lengths.resize(dim);
    H5Sget_simple_extent_dims(d_space, cb_out.lengths.data(), nullptr);

    size_t size = H5Tget_size(ty);

    cb_out.lengths.push_back(size); //  2 ?? last one is size of the string +1
    long ltot = std::accumulate(cb_out.lengths.begin(), cb_out.lengths.end(), 1, std::multiplies<>());
    cb_out.buffer.resize(ltot, 0x00);

    H5_PRINT(ltot);
    H5_PRINT(cb_out.lengths.size());
    H5_PRINT(cb_out.lengths[0]);
    H5_PRINT(cb_out.lengths[1]);
    H5_PRINT(cb_out.buffer.size());
    H5_PRINT(size);

    auto err = H5Dread(ds, cb_out.dtype(), cb_out.dspace(), H5S_ALL, H5P_DEFAULT, (void *)cb_out.buffer.data());
    if (err < 0) throw make_runtime_error("Error reading the vector<string> ", name, " in the group", g.name());

    _cb = std::move(cb_out);
  }

  // ----- read attribute -----

  void h5_read_attribute(hid_t id, std::string const &name, char_buf &_cb) {

    attribute attr = H5Aopen(id, name.c_str(), H5P_DEFAULT);
    if (!attr.is_valid()) throw make_runtime_error("Cannot open the attribute ", name);

    dataspace d_space = H5Aget_space(attr);

    char_buf cb_out;

    int dim = H5Sget_simple_extent_ndims(d_space);
    cb_out.lengths.resize(dim);
    H5Sget_simple_extent_dims(d_space, cb_out.lengths.data(), nullptr);

    size_t size = H5Aget_storage_size(attr);
    cb_out.lengths.push_back(size + 1); // last one is size of the string +1

    //long ltot = std::accumulate(cb_out.lengths.begin(), cb_out.lengths.end(), 1, std::multiplies<>());
    cb_out.buffer.resize(0x00);

    auto err = H5Aread(attr, cb_out.dtype(), (void *)cb_out.buffer.data());
    if (err < 0) throw make_runtime_error("Cannot read the attribute ", name);

    _cb = std::move(cb_out);
  }

} // namespace h5
