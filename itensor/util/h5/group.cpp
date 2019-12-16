#include <vector>
#include "./group.hpp"
#include "./stl/string.hpp"

#include <hdf5.h>
#include <hdf5_hl.h>

static_assert(std::is_same<hid_t, int64_t>::value or std::is_same<hid_t, int>::value, "Internal error");

using namespace std::string_literals;

namespace h5 {

  //static_assert(std::is_same<::hid_t, hid_t>::value, "Internal error");

  group::group(h5::file f) : h5_object(), parent_file(f) {
    id = H5Gopen2(f, "/", H5P_DEFAULT);
    if (id < 0) throw std::runtime_error("Cannot open the root group / in the file " + f.name());
  }

  //group::group(h5_object obj) : h5_object(std::move(obj)) {
  //if (!H5Iis_valid(this->id)) throw std::runtime_error("Invalid input in group constructor from id");
  //if (H5Iget_type(this->id) != H5I_GROUP) throw std::runtime_error("Group constructor must take the id of a group or a file ");
  //}

  //group::group(hid_t id_) : group(h5_object(id_)) {}

  std::string group::name() const {
    char _n[1];
    ssize_t size = H5Iget_name(id, _n, 1); // first call, get the size only
    std::vector<char> buf(size + 1, 0x00);
    H5Iget_name(id, buf.data(), size + 1); // now get the name
    std::string res = "";
    res.append(&(buf.front()));
    return res;
  }

  bool group::has_key(std::string const &key) const { return H5Lexists(id, key.c_str(), H5P_DEFAULT); }

  bool group::has_subgroup(std::string const &key) const {
    if (!has_key(key)) return false;
    hid_t id_node = H5Oopen(id, key.c_str(), H5P_DEFAULT);
    if (id_node <= 0) return false;
    bool r = (H5Iget_type(id_node) == H5I_GROUP);
    H5Oclose(id_node);
    return r;
  }

  bool group::has_dataset(std::string const &key) const {
    if (!has_key(key)) return false;
    hid_t id_node = H5Oopen(id, key.c_str(), H5P_DEFAULT);
    if (id_node <= 0) return false;
    bool r = (H5Iget_type(id_node) == H5I_DATASET);
    H5Oclose(id_node);
    return r;
  }

  void group::unlink(std::string const &key, bool error_if_absent) const {
    if (!has_key(key)) {
      if (error_if_absent) throw std::runtime_error("The key " + key + " is not present in the group " + name());
      return;
    }
    //auto err = H5Gunlink(id, key.c_str()); // deprecated function
    auto err = H5Ldelete(id, key.c_str(), H5P_DEFAULT);
    if (err < 0) throw std::runtime_error("Cannot unlink object " + key + " in group " + name());
  }

  group group::open_group(std::string const &key) const {
    if (key.empty()) return *this;
    if (!has_key(key)) throw std::runtime_error("no subgroup " + key + " in the group");
    hid_t sg = H5Gopen2(id, key.c_str(), H5P_DEFAULT);
    if (sg < 0) throw std::runtime_error("Error in opening the subgroup " + key);
    return group(h5_object{sg}, parent_file);
  }

  group group::create_group(std::string const &key, bool delete_if_exists) const {
    if (key.empty()) return *this;
    if (delete_if_exists) unlink(key);
    hid_t id_g = H5Gcreate2(id, key.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (id_g < 0) throw std::runtime_error("Cannot create the subgroup " + key + " of the group" + name());
    return group(h5_object{id_g}, parent_file);
  }

  /// Open an existing DataSet. Throw if it does not exist.
  dataset group::open_dataset(std::string const &key) const {
    if (!has_key(key)) throw std::runtime_error("no dataset " + key + " in the group");
    dataset ds = H5Dopen2(id, key.c_str(), H5P_DEFAULT);
    if (!ds.is_valid()) throw std::runtime_error("Cannot open dataset " + key + " in the group" + name());
    return ds;
  }

  /**
  * \brief Create a dataset.
  * \param key The name of the subgroup
  *
  * NB : It unlinks the dataset if it exists.
  */
  dataset group::create_dataset(std::string const &key, datatype ty, dataspace sp, hid_t pl) const {
    unlink(key);
    dataset ds = H5Dcreate2(id, key.c_str(), ty, sp, H5P_DEFAULT, pl, H5P_DEFAULT);
    if (!ds.is_valid()) throw std::runtime_error("Cannot create the dataset " + key + " in the group" + name());
    return ds;
  }

  dataset group::create_dataset(std::string const &key, datatype ty, dataspace sp) const { return create_dataset(key, ty, sp, H5P_DEFAULT); }

  //----------------------------------------------------------
  // Keep as an example of H5LTset_attribute_string
  /*
 void group::write_string_attribute (std::string const & obj_name, std::string const & attr_name, std::string const & value){
  herr_t err = H5LTset_attribute_string(id, obj_name.c_str(), attr_name.c_str(), value.c_str());
  if (err<0) throw std::runtime_error( "Error in setting attribute of "+ obj_name+" named "+ attr_name + " to " + value);
 }
*/
  //-----------------------------------------------------------------------

  // C callbacks for the next functions using H5Literate
  extern "C" {
  herr_t get_group_elements_name_ds(::hid_t loc_id, const char *name, const H5L_info_t *, void *opdata) {
    H5O_info_t object_info;
    herr_t err = H5Oget_info_by_name(loc_id, name, &object_info, H5P_DEFAULT);
    if (err < 0) throw std::runtime_error("get_group_elements_name_ds internal");
    if (object_info.type == H5O_TYPE_DATASET) static_cast<std::vector<std::string> *>(opdata)->push_back(name);
    return 0;
  }
  herr_t get_group_elements_name_grp(::hid_t loc_id, const char *name, const H5L_info_t *, void *opdata) {
    H5O_info_t object_info;
    herr_t err = H5Oget_info_by_name(loc_id, name, &object_info, H5P_DEFAULT);
    if (err < 0) throw std::runtime_error("get_group_elements_name_grp internal");
    if (object_info.type == H5O_TYPE_GROUP) static_cast<std::vector<std::string> *>(opdata)->push_back(name);
    return 0;
  }
  herr_t get_group_elements_name_ds_grp(::hid_t loc_id, const char *name, const H5L_info_t *, void *opdata) {
    H5O_info_t object_info;
    herr_t err = H5Oget_info_by_name(loc_id, name, &object_info, H5P_DEFAULT);
    if (err < 0) throw std::runtime_error("get_group_elements_name_grp internal");
    if ((object_info.type == H5O_TYPE_GROUP) or (object_info.type == H5O_TYPE_DATASET))
      static_cast<std::vector<std::string> *>(opdata)->push_back(name);
    return 0;
  }
  }
  //-----------------------------------------------------------------------

  std::vector<std::string> group::get_all_subgroup_names() const {
    std::vector<std::string> grp_name;
    int r = H5Literate(::hid_t(id), H5_INDEX_NAME, H5_ITER_NATIVE, nullptr, get_group_elements_name_grp, static_cast<void *>(&grp_name));
    if (r != 0) throw std::runtime_error("Iteration over subgroups of group " + name() + "failed");
    return grp_name;
  }

  std::vector<std::string> group::get_all_dataset_names() const {
    std::vector<std::string> ds_name;
    int r = H5Literate(::hid_t(id), H5_INDEX_NAME, H5_ITER_NATIVE, nullptr, get_group_elements_name_ds, static_cast<void *>(&ds_name));
    if (r != 0) throw std::runtime_error("Iteration over datasets of group " + name() + "failed");
    return ds_name;
  }

  std::vector<std::string> group::get_all_subgroup_dataset_names() const {
    std::vector<std::string> ds_name;
    int r = H5Literate(::hid_t(id), H5_INDEX_NAME, H5_ITER_NATIVE, nullptr, get_group_elements_name_ds_grp, static_cast<void *>(&ds_name));
    if (r != 0) throw std::runtime_error("Iteration over datasets of group " + name() + "failed");
    return ds_name;
  }

} // namespace h5
