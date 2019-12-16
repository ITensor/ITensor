#pragma once
#include <utility>

#include "./file.hpp"
#include "./format.hpp"

namespace h5 {

  /**
   *  HDF5 group
   */
  class group : public h5_object {

    h5::file parent_file;

    public:
    /// Takes the "/" group at the top of the file
    group(h5::file f);

    ///
    group(group const &) = default;

    private:
    // construct from the bare object and the parent
    // internal use only for open/create subgroup
    group(h5_object obj, h5::file _parent_file) : h5_object{obj}, parent_file(std::move(_parent_file)) {}

    public:
    /// Name of the group
    [[nodiscard]] std::string name() const;

    /// Access to the parent file
    [[nodiscard]] h5::file get_file() const { return parent_file; }

    /**
     * True iff key is an object in the group
     *
     * @param key
     */
    [[nodiscard]] bool has_key(std::string const &key) const;

    /**
     * True iff key is a subgroup of this.
     *
     * @param key
     */
    [[nodiscard]] bool has_subgroup(std::string const &key) const;

    /**
     * True iff key is a dataset of this.
     *
     * @param key
     */
    [[nodiscard]] bool has_dataset(std::string const &key) const;

    /**
     * Unlinks the subgroup key if it exists
     * No error is thrown if key does not exists
     * NB : unlink is almost a remove, but it does not really remove from the file (memory is not freed).
     * After unlinking a large object, a h5repack may be needed. Cf HDF5 documentation.
     *
     * @param key The name of the subgroup to be removed.
     * @param error_if_absent If True, throws an error if the key is missing.
     */
    void unlink(std::string const &key, bool error_if_absent = false) const;

    /**
     * Open a subgroup.
     * Throws std::runtime_error if it does not exist.
     *
     * @param key  The name of the subgroup. If empty, return this group
     */
    [[nodiscard]] group open_group(std::string const &key) const;

    /**
     * Create a subgroup in this group
     * 
     * @param key  The name of the subgroup. If empty, return this group.
     * @param delete_if_exists  Unlink the group if it exists
     */
    group create_group(std::string const &key, bool delete_if_exists = true) const; // NOLINT

    /**
     * Open a existing DataSet in the group.
     * Throws std::runtime_error if it does not exist.
     *
     * @param key  The name of the subgroup. If empty, return this group
     */
    [[nodiscard]] dataset open_dataset(std::string const &key) const;

    /**
     * Create a dataset in this group
     * 
     * @param key  The name of the dataset. 
     * @param ty Datatype
     * @param sp Dataspace
     */
    [[nodiscard]] dataset create_dataset(std::string const &key, datatype ty, dataspace sp) const;

    /**
     * Create a dataset in this group
     * 
     * @param key  The name of the dataset. 
     * @param ty Datatype
     * @param sp Dataspace
     * @param pl Property list
     */
    [[nodiscard]] dataset create_dataset(std::string const &key, datatype ty, dataspace sp, hid_t pl) const;

    /// Returns all names of subgroup of  G
    [[nodiscard]] std::vector<std::string> get_all_subgroup_names() const;

    /// Returns all names of dataset of G
    [[nodiscard]] std::vector<std::string> get_all_dataset_names() const;

    /// Returns all names of dataset of G
    [[nodiscard]] std::vector<std::string> get_all_subgroup_dataset_names() const;
  };

} // namespace h5
