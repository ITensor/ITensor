#include "./group.hpp"
#include "./base.hpp"

namespace triqs {
  namespace h5 {

    group::group(h5::file f) : h5_object() {
      id = H5Gopen2(f, "/", H5P_DEFAULT);
      if (id < 0) TRIQS_RUNTIME_ERROR << "Cannot open the root group / in the file " << f.name();
    }

    group::group(h5_object obj) : h5_object(std::move(obj)) {
      if (!H5Iis_valid(this->id)) TRIQS_RUNTIME_ERROR << "Invalid input in group constructor from id";
      if (H5Iget_type(this->id) != H5I_GROUP) TRIQS_RUNTIME_ERROR << "Group constructor must take the id of a group or a file ";
    }

    group::group(hid_t id_) : group(h5_object(id_)) {}

    void group::write_hdf5_scheme_as_string(const char *a) { h5_write_attribute(id, "TRIQS_HDF5_data_scheme", a); }

    std::string group::read_hdf5_scheme() const {
      std::string s;
      h5_read_attribute(id, "TRIQS_HDF5_data_scheme", s);
      return s;
    }

    void group::assert_hdf5_scheme_as_string(const char *tag_expected, bool ignore_if_absent) const {
      auto tag_file = read_hdf5_scheme();
      if (ignore_if_absent and tag_file.empty()) return;
      if (tag_file != tag_expected)
        TRIQS_RUNTIME_ERROR << "h5_read : mismatch of the tag TRIQS_HDF5_data_scheme tag in the h5 group : found " << tag_file << " while I expected "
                            << tag_expected;
    }

    std::string group::name() const {
      char _n[1];
      ssize_t size = H5Iget_name(id, _n, 1); // first call, get the size only
      std::vector<char> buf(size + 1, 0x00);
      H5Iget_name(id, buf.data(), 1); // now get the name
      std::string res = "";
      res.append(&(buf.front()));
      return res;
    }

    bool group::has_key(std::string const &key) const { return H5Lexists(id, key.c_str(), H5P_DEFAULT); }

    void group::unlink_key_if_exists(std::string const &key) const {
      if (!has_key(key)) return;
      //auto err = H5Gunlink(id, key.c_str()); // deprecated function
      auto err = H5Ldelete(id, key.c_str(), H5P_DEFAULT);
      if (err < 0) TRIQS_RUNTIME_ERROR << "Cannot unlink object " << key << " in group " << name();
    }

    group group::open_group(std::string const &key) const {
      if (key.empty()) return *this;
      if (!has_key(key)) TRIQS_RUNTIME_ERROR << "no subgroup " << key << " in the group";
      hid_t sg = H5Gopen2(id, key.c_str(), H5P_DEFAULT);
      if (sg < 0) TRIQS_RUNTIME_ERROR << "Error in opening the subgroup " << key;
      return group(sg);
    }

    /// Open an existing DataSet. Throw if it does not exist.
    dataset group::open_dataset(std::string const &key) const {
      if (!has_key(key)) TRIQS_RUNTIME_ERROR << "no dataset " << key << " in the group";
      dataset ds = H5Dopen2(id, key.c_str(), H5P_DEFAULT);
      if (!ds.is_valid()) TRIQS_RUNTIME_ERROR << "Cannot open dataset " << key << " in the group" << name();
      return ds;
    }

    group group::create_group(std::string const &key, bool delete_if_exists) const {
      if (key.empty()) return *this;
      unlink_key_if_exists(key);
      hid_t id_g = H5Gcreate2(id, key.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      if (id_g < 0) TRIQS_RUNTIME_ERROR << "Cannot create the subgroup " << key << " of the group" << name();
      return group(id_g);
    }

    /**
  * \brief Create a dataset.
  * \param key The name of the subgroup
  *
  * NB : It unlinks the dataset if it exists.
  */
    dataset group::create_dataset(std::string const &key, datatype ty, dataspace sp, hid_t pl) const {
      unlink_key_if_exists(key);
      dataset ds = H5Dcreate2(id, key.c_str(), ty, sp, H5P_DEFAULT, pl, H5P_DEFAULT);
      if (!ds.is_valid()) TRIQS_RUNTIME_ERROR << "Cannot create the dataset " << key << " in the group" << name();
      return ds;
    }

    //----------------------------------------------------------
    // Keep as an example of H5LTset_attribute_string
    /*
 void group::write_string_attribute (std::string const & obj_name, std::string const & attr_name, std::string const & value){
  herr_t err = H5LTset_attribute_string(id, obj_name.c_str(), attr_name.c_str(), value.c_str());
  if (err<0) TRIQS_RUNTIME_ERROR << "Error in setting attribute of "<< obj_name<<" named "<< attr_name << " to " << value;
 }
*/
    //-----------------------------------------------------------------------

    // C callbacks for the next functions using H5Literate
    extern "C" {
    herr_t get_group_elements_name_ds(hid_t loc_id, const char *name, const H5L_info_t *info, void *opdata) {
      H5O_info_t object_info;
      herr_t err = H5Oget_info_by_name(loc_id, name, &object_info, H5P_DEFAULT);
      if (err < 0) TRIQS_RUNTIME_ERROR << "get_group_elements_name_ds internal";
      if (object_info.type == H5O_TYPE_DATASET) static_cast<std::vector<std::string> *>(opdata)->push_back(name);
      return 0;
    }
    herr_t get_group_elements_name_grp(hid_t loc_id, const char *name, const H5L_info_t *info, void *opdata) {
      H5O_info_t object_info;
      herr_t err = H5Oget_info_by_name(loc_id, name, &object_info, H5P_DEFAULT);
      if (err < 0) TRIQS_RUNTIME_ERROR << "get_group_elements_name_grp internal";
      if (object_info.type == H5O_TYPE_GROUP) static_cast<std::vector<std::string> *>(opdata)->push_back(name);
      return 0;
    }
    herr_t get_group_elements_name_ds_grp(hid_t loc_id, const char *name, const H5L_info_t *info, void *opdata) {
      H5O_info_t object_info;
      herr_t err = H5Oget_info_by_name(loc_id, name, &object_info, H5P_DEFAULT);
      if (err < 0) TRIQS_RUNTIME_ERROR << "get_group_elements_name_grp internal";
      if ((object_info.type == H5O_TYPE_GROUP) or (object_info.type == H5O_TYPE_DATASET))
        static_cast<std::vector<std::string> *>(opdata)->push_back(name);
      return 0;
    }
    }
    //-----------------------------------------------------------------------

    std::vector<std::string> group::get_all_subgroup_names() const {
      std::vector<std::string> grp_name;
      int r = H5Literate(id, H5_INDEX_NAME, H5_ITER_NATIVE, NULL, get_group_elements_name_grp, static_cast<void *>(&grp_name));
      if (r != 0) TRIQS_RUNTIME_ERROR << "Iteration over subgroups of group " << name() << "failed";
      return grp_name;
    }

    std::vector<std::string> group::get_all_dataset_names() const {
      std::vector<std::string> ds_name;
      int r = H5Literate(id, H5_INDEX_NAME, H5_ITER_NATIVE, NULL, get_group_elements_name_ds, static_cast<void *>(&ds_name));
      if (r != 0) TRIQS_RUNTIME_ERROR << "Iteration over datasets of group " << name() << "failed";
      return ds_name;
    }

    std::vector<std::string> group::get_all_subgroup_dataset_names() const {
      std::vector<std::string> ds_name;
      int r = H5Literate(id, H5_INDEX_NAME, H5_ITER_NATIVE, NULL, get_group_elements_name_ds_grp, static_cast<void *>(&ds_name));
      if (r != 0) TRIQS_RUNTIME_ERROR << "Iteration over datasets of group " << name() << "failed";
      return ds_name;
    }

  } // namespace h5
} // namespace triqs
