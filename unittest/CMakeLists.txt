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

set(executable itensor_check)
#ifdef test
#SOURCES+= $(test)_test.cc
#else
set(itensor_test_src_files
        util_test.cc
        algorithm_test.cc
        real_test.cc
        args_test.cc
        matrix_test.cc
        tensor_test.cc
        contract_test.cc
        sparse_contract_test.cc
        index_test.cc
        indexset_test.cc
        itensor_test.cc
        qn_test.cc
        decomp_test.cc
        mps_test.cc
        mpo_test.cc
        autompo_test.cc
        iterativesolvers_test.cc
        regression_test.cc
        localop_test.cc
        siteset_test.cc
        test.cc
        )

if(ENABLE_HDF5)
    set(h5_test_source_files
            hdf5_test.cc)
    set(itensor_test_src_files
            ${itensor_test_src_files}
            ${h5_test_source_files}
            )
endif(ENABLE_HDF5)

add_executable(${executable} EXCLUDE_FROM_ALL ${itensor_test_src_files})
# Add include directories and compiler flags for ta_test
target_include_directories(${executable} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_BINARY_DIR}
        )
target_link_libraries(${executable} itensor)
if(ENABLE_HDF5)
    target_link_libraries(${executable} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} HDF5::HDF5)
endif(ENABLE_HDF5)

# Add test suite
add_test(itensor/unit/build "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --target ${executable})
set_tests_properties(itensor/unit/build PROPERTIES FIXTURES_SETUP ITENSOR_TEST_EXEC)
add_test(NAME itensor/unit/run
        COMMAND ${executable})
set_tests_properties(itensor/unit/run
        PROPERTIES
        FIXTURES_REQUIRED ITENSOR_TEST_EXEC
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/unittest"
        )
target_compile_options(${executable} PUBLIC ${CMAKE_CXX_FLAG_LIST})
target_compile_features(${executable} PUBLIC "cxx_std_${CMAKE_CXX_STANDARD}")
