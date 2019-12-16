#include "./array_interface.hpp"
#include "./stl/string.hpp"

#include <hdf5.h>
#include <hdf5_hl.h>

#include <numeric>
#include <iostream> // DEBUG

namespace h5::array_interface {

  //------------------------------------------------
  // the dataspace corresponding to the array. Contiguous data only...
  dataspace make_mem_dpace(h5_array_view const &v) {

    if (v.rank() == 0) return H5Screate(H5S_SCALAR);

    dataspace ds = H5Screate_simple(v.slab.rank(), v.L_tot.data(), nullptr);
    if (!ds.is_valid()) throw std::runtime_error("Cannot create the dataset");

    herr_t err = H5Sselect_hyperslab(ds, H5S_SELECT_SET, v.slab.offset.data(), v.slab.stride.data(), v.slab.count.data(),
                                     (v.slab.block.empty() ? nullptr : v.slab.block.data()));
    if (err < 0) throw std::runtime_error("Cannot set hyperslab");

    return ds;
  }

  //------------------------------------------------

  std::pair<v_t, v_t> get_L_tot_and_strides_h5(long const *stri, int rank, long total_size) {
    v_t Ltot(rank), strides_h5(rank);
    for (int u = 0; u < rank; ++u) strides_h5[u] = stri[u];
    Ltot[0] = total_size;

    for (int u = rank - 2; u >= 0; --u) {
      // L[u+1] as  gcd of size and stride[u] ... stride[0]
      long L = strides_h5[u];
      // L becomes the  gcd
      for (int v = u - 1; v >= 0; --v) { L = std::gcd(L, strides_h5[v]); }
      // divides
      for (int v = u; v >= 0; --v) { strides_h5[v] /= L; }
      Ltot[u + 1] = L;
    }

    return {Ltot, strides_h5};
  }

  //-------------------------------------------------------
  //                    write
  //--------------------------------------------------------

  void write(group g, std::string const &name, h5_array_view const &v, bool compress) {

    g.unlink(name);

    // Some properties for the dataset : add compression
    proplist cparms = H5P_DEFAULT;
    if (compress and (v.rank() != 0)) {
      int n_dims = v.rank();
      std::vector<hsize_t> chunk_dims(n_dims);
      for (int i = 0; i < v.rank(); ++i) chunk_dims[i] = std::max(v.slab.count[i], hsize_t{1});
      cparms = H5Pcreate(H5P_DATASET_CREATE);
      H5Pset_chunk(cparms, n_dims, chunk_dims.data());
      // FIXME : OLD COMMIT
      H5Pset_deflate(cparms, 8);
    }

    // dataspace for the dataset in the file
    dataspace file_dspace = H5Screate_simple(v.slab.rank(), v.slab.count.data(), nullptr);

    // create the dataset in the file
    dataset ds = H5Dcreate2(g, name.c_str(), v.ty, file_dspace, H5P_DEFAULT, cparms, H5P_DEFAULT);
    if (!ds.is_valid()) throw std::runtime_error("Cannot create the dataset " + name + " in the group" + g.name());

    // memory data space
    dataspace mem_d_space = make_mem_dpace(v);
    if (H5Sget_simple_extent_npoints(mem_d_space) > 0) { // avoid writing empty arrays
      herr_t err = H5Dwrite(ds, v.ty, mem_d_space, H5S_ALL, H5P_DEFAULT, v.start);
      if (err < 0) throw std::runtime_error("Error writing the scalar dataset " + name + " in the group" + g.name());
    }

    // If we are dealing with complex, we had the attribute
    if (v.is_complex) h5_write_attribute(ds, "__complex__", "1");
  }

  //-------------------------------------------------------------

  void write_attribute(hid_t id, std::string const &name, h5_array_view v) {

    if (H5LTfind_attribute(id, name.c_str()) != 0) throw std::runtime_error("The attribute " + name + " is already present. Can not overwrite");

    dataspace mem_d_space = make_mem_dpace(v);

    attribute attr = H5Acreate2(id, name.c_str(), v.ty, mem_d_space, H5P_DEFAULT, H5P_DEFAULT);
    if (!attr.is_valid()) throw std::runtime_error("Cannot create the attribute " + name);

    herr_t err = H5Awrite(attr, v.ty, v.start);
    if (err < 0) throw std::runtime_error("Cannot write the attribute " + name);
  }

  //-------------------------------------------------------
  //                    READ
  //--------------------------------------------------------

  h5_lengths_type get_h5_lengths_type(group g, std::string const &name) {

    dataset ds = g.open_dataset(name);

    bool has_complex_attribute = H5LTfind_attribute(ds, "__complex__"); // the array in file should be interpreted as a complex
    dataspace d_space          = H5Dget_space(ds);
    int rank                   = H5Sget_simple_extent_ndims(d_space);

    // need to use hsize_t here and the vector is truncated at rank
    v_t dims_out(rank);
    H5Sget_simple_extent_dims(d_space, dims_out.data(), nullptr);

    //  get the type from the file
    datatype ty = H5Dget_type(ds);
    return {std::move(dims_out), ty, has_complex_attribute};
  }

  //--------------------------------------------------------

  void read(group g, std::string const &name, h5_array_view v, h5_lengths_type lt) {

    dataset ds             = g.open_dataset(name);
    dataspace file_d_space = H5Dget_space(ds);

    // Checks
    if (H5Tequal(v.ty, lt.ty) <= 0)
      throw std::runtime_error("h5 read. Type mismatch : expecting a " + get_name_of_h5_type(v.ty)
                               + " while the array stored in the hdf5 file has type " + get_name_of_h5_type(lt.ty));

    if (lt.rank() != v.rank())
      throw std::runtime_error("h5 read. Rank mismatch : expecting in file a rank " + std::to_string(v.rank())
                               + " while the array stored in the hdf5 file has rank " + std::to_string(lt.rank()));

    //if (lt.lengths != v.slab.count)
    //throw std::runtime_error("h5 read. Lengths mismatch : expecting a rank " + std::to_string(v.rank())
    //+ " while the array stored in the hdf5 file has rank = " + std::to_string(lt.rank()));

    dataspace mem_d_space = make_mem_dpace(v);
    if (H5Sget_simple_extent_npoints(file_d_space) > 0) {
      herr_t err = H5Dread(ds, v.ty, mem_d_space, file_d_space, H5P_DEFAULT, v.start);
      if (err < 0) throw std::runtime_error("Error reading the scalar dataset " + name + " in the group" + g.name());
    }
  }

  //-------------------------------------------------------------

  void read_attribute(hid_t id, std::string const &name, h5_array_view v) {

    //if (v.rank() != 0) throw std::runtime_error("Non scalar attribute not implemented");

    attribute attr = H5Aopen(id, name.c_str(), H5P_DEFAULT);
    if (!attr.is_valid()) throw std::runtime_error("Cannot open the attribute " + name);

    dataspace space = H5Aget_space(attr);
    int rank        = H5Sget_simple_extent_ndims(space);
    if (rank != 0) throw std::runtime_error("Reading a scalar attribute and got rank !=0");

    auto eq = H5Tequal(H5Aget_type(attr), v.ty);
    if (eq < 0) throw std::runtime_error("Type comparison failure in reading attribute");
    if (eq == 0) throw std::runtime_error("Type mismatch in reading attribute");

    auto err = H5Aread(attr, v.ty, v.start);
    if (err < 0) throw std::runtime_error("Cannot read the attribute " + name);
  }

} // namespace h5::array_interface
