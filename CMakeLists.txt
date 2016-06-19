cmake_minimum_required (VERSION 2.8)

project (ITensor CXX)

# Disable build in source (not to overlap with existing Makefile)
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

if(NOT DEFINED ITENSOR_BUILD_TYPE)
    set(ITENSOR_BUILD_TYPE STATIC)
endif()

# Select default build type : Release = -03 -DNDEBUG
if(NOT CMAKE_BUILD_TYPE:String AND NOT CMAKE_CONFIGURATION_TYPES)
   SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
else()
   SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# rpath fix on a mac
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
 set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
else()
 set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

# enable C++11 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs")
    if (NOT (GCC_VERSION VERSION_GREATER 4.6 OR GCC_VERSION VERSION_EQUAL 4.6))
        message(FATAL_ERROR "${PROJECT_NAME} requires g++ 4.6 or greater.")
    endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
    if ("${CMAKE_SYSTEM_NAME}" MATCHES "Linux")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++")
    elseif ("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
        message( STATUS "Adding -stdlib=libc++ flag")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    else()
        message(FATAL_ERROR "Platform undefined")
    endif()
else ()
    message(FATAL_ERROR "Your C++ compiler does not support C++11.")
endif ()


# LAPACK
# Check if LAPACK vendor is not defined and try to set it.
if (NOT DEFINED BLA_VENDOR)
    if (${APPLE})
        message(STATUS "Setting lapack platform : Apple")
        set(BLA_VENDOR Apple)
    else()
        set(BLA_VENDOR All)
    endif()
endif()
# define platform flag inside of include/lapack_wrap.h
message(STATUS "Searching for LAPACK vendor : ${BLA_VENDOR}")
if (${BLA_VENDOR} MATCHES "Intel" OR "mkl")
    add_definitions(-DPLATFORM_mkl)
elseif(${BLA_VENDOR} MATCHES ACML)
    add_definitions(-PLATFORM_mkl)
elseif(${BLA_VENDOR} MATCHES Apple)
    add_definitions(-DPLATFORM_macos)
else()
    add_definitions(-DPLATFORM_lapack)
endif()
find_package(LAPACK REQUIRED)
list(APPEND CMAKE_LINKER_FLAGS ${LAPACK_LINKER_FLAGS})
link_libraries(${LAPACK_LIBRARIES})

# build itensor 
include_directories(".")
file(GLOB_RECURSE itensor_sources RELATIVE ${CMAKE_SOURCE_DIR} itensor/*.cc) 
list(REMOVE_ITEM itensor_sources "itensor/itdata/itlazy.cc")
add_library(itensor ${ITENSOR_BUILD_TYPE} ${itensor_sources})

# Enable testing
option(Testing "Enable testing" ON)
if (Testing)
    enable_testing()
    list(APPEND itensor_test_sources 
        unittest/util_test.cc
        unittest/algorithm_test.cc
        unittest/real_test.cc
        unittest/args_test.cc
        unittest/matrix_test.cc
        unittest/tensor_test.cc
        unittest/sparse_contract_test.cc
        unittest/index_test.cc
        unittest/indexset_test.cc
        unittest/itensor_test.cc
        unittest/qn_test.cc
        unittest/iqindex_test.cc
        unittest/iqtensor_test.cc
        unittest/decomp_test.cc
        unittest/autompo_test.cc
        unittest/regression_test.cc
        unittest/localop_test.cc
        )

    message(STATUS ${itensor_test_sources})
    foreach(test_src ${itensor_test_sources})
        string(REPLACE ".cc" "" test_name ${test_src})
        string(REPLACE "unittest/" "" test_name ${test_name})
        add_executable(${test_name} ${test_src} ./unittest/test.cc)
        target_link_libraries(${test_name} itensor)
        add_test(${test_name} ${test_name})
    endforeach()
    
    
#    add_subdirectory(unittest)
endif (Testing)

# Build example
option(Tutorial "Build tutorial" ON)
if (Tutorial)
    message(STATUS "Building tutorial")
 #   add_subdirectory(tutorial)
endif (Tutorial)

# Build sample
option(Sample "Build samples" ON)
if (Sample)
    message(STATUS "Building samples")
  #  add_subdirectory(sample)
endif (Sample)

# Install pkg-config file
#configure_file(ITensor.pc.in ITensor.pc @ONLY)
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/ITensor.pc DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/pkgconfig)
