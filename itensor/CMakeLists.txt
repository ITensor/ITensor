#
#  This file is part of the ITensor library
#
# Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#ITensor headers
set(ITENSOR_HEADER_FILES
    all.h
    all_basic.h
    all_mps.h
    arrow.h
    core.h
    decomp.h
    decomp_impl.h
    detail/algs.h
    detail/call_rewrite.h
    detail/gcounter.h
    detail/skip_iterator.h
    global.h
    index.h
    index_impl.h
    indexset.h
    indexset_impl.h
    itdata/applyfunc.h
    itdata/combiner.h
    itdata/dense.h
    itdata/diag.h
    itdata/dotask.h
    itdata/itdata.h
    itdata/itlazy.h
    itdata/qcombiner.h
    itdata/qdense.h
    itensor.h
    itensor_impl.h
    iterativesolvers.h
    qn.h
    real.h
    smallstring.h
    spectrum.h
    tagset.h
    types.h
    )

set(ITENSOR_SOURCE_FILES
        util/args.cc
        util/input.cc
        util/cputime.cc
        tensor/cmake_lapack_wrap.cc
        tensor/vec.cc
        tensor/mat.cc
        tensor/gemm.cc
        tensor/algs.cc
        tensor/contract.cc
        itdata/dense.cc
        itdata/combiner.cc
        itdata/diag.cc
        itdata/qdense.cc
        itdata/qcombiner.cc
        itdata/qdiag.cc
        itdata/scalar.cc
        qn.cc
        tagset.cc
        index.cc
        indexset.cc
        itensor.cc
        spectrum.cc
        decomp.cc
        hermitian.cc
        svd.cc
        global.cc
        mps/mps.cc
        mps/mpsalgs.cc
        mps/mpo.cc
        mps/mpoalgs.cc
        mps/autompo.cc)

if(ENABLE_HDF5)
  set(H5_SOURCE_FILES
          util/h5/array_interface.cc
          util/h5/file.cc
          util/h5/group.cc
          util/h5/h5object.cc
          util/h5/stl/string.cc
          util/h5/stl/vector.cc)

  set(ITENSOR_SOURCE_FILES
          ${ITENSOR_SOURCE_FILES}
          ${H5_SOURCE_FILES})
endif(ENABLE_HDF5)

# Create the ITensor library
add_library(itensor ${ITENSOR_SOURCE_FILES} ${ITENSOR_HEADER_FILES})

set(targetname itensor)

target_include_directories(${targetname} INTERFACE
                           $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
                           $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
                           $<INSTALL_INTERFACE:${ITENSOR_INSTALL_INCLUDEDIR}>
                           )

set(ITENSOR_LINK_LIBRARIES blaspp lapackpp ${blaspp_headers} CACHE STRING "List of libraries which ITensor is dependent on")
target_link_libraries(${targetname} LINK_PRIVATE ${ITENSOR_LINK_LIBRARIES} blaspp lapackpp)
if(ENABLE_HDF5)
  target_link_libraries(${targetname} PUBLIC ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} hdf5 hdf5_hl)
endif(ENABLE_HDF5)

# append current CMAKE_CXX_FLAGS
string(REPLACE " " ";" CMAKE_CXX_FLAG_LIST "${CMAKE_CXX_FLAGS}")
target_compile_options(${targetname} PUBLIC ${CMAKE_CXX_FLAG_LIST})
target_compile_features(${targetname} PUBLIC "cxx_std_${CMAKE_CXX_STANDARD}")

if (LAPACK_INCLUDE_DIRS)
  target_include_directories(${targetname} PUBLIC ${LAPACK_INCLUDE_DIRS})
endif (LAPACK_INCLUDE_DIRS)
if (LAPACK_COMPILE_OPTIONS)
  target_compile_options(${targetname} PUBLIC ${LAPACK_COMPILE_OPTIONS})
endif (LAPACK_COMPILE_OPTIONS)
if (LAPACK_COMPILE_DEFINITIONS)
  target_compile_definitions(${targetname} PUBLIC ${LAPACK_COMPILE_DEFINITIONS})
endif (LAPACK_COMPILE_DEFINITIONS)

# Add library to the list of installed components
install(TARGETS itensor EXPORT itensor COMPONENT itensor
        LIBRARY DESTINATION "${ITENSOR_INSTALL_LIBDIR}"
        ARCHIVE DESTINATION "${ITENSOR_INSTALL_LIBDIR}")
# Install header files
install(
    DIRECTORY
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_BINARY_DIR}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.h"
    PATTERN "CMakeFiles" EXCLUDE
)
