#; -*-CMake-*-

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

cmake_minimum_required(VERSION 3.17.0) # decouples C++ and CUDA standards, see https://gitlab.kitware.com/cmake/cmake/issues/19123
set(ITENSORS_REQUIRED_CMAKE_VERSION 3.17.0 CACHE INTERNAL "")

###############################################################################
# CMake defaults to address key pain points
###############################################################################

# safety net for dev workflow: accidental install will not affect FindOrFetch*
if (NOT DEFINED CACHE{CMAKE_FIND_NO_INSTALL_PREFIX})
  set(CMAKE_FIND_NO_INSTALL_PREFIX ON CACHE BOOL "Whether find_* commands will search CMAKE_INSTALL_PREFIX and CMAKE_STAGING_PREFIX; see https://cmake.org/cmake/help/latest/variable/CMAKE_FIND_NO_INSTALL_PREFIX.html#variable:CMAKE_FIND_NO_INSTALL_PREFIX")
endif()

###############################################################################
# Bring ValeevGroup cmake toolkit to autofind linalg
###############################################################################
include(FetchContent)
if (DEFINED PROJECT_BINARY_DIR)
  set(VG_CMAKE_KIT_PREFIX_DIR PROJECT_BINARY_DIR)
else ()
  set(VG_CMAKE_KIT_PREFIX_DIR CMAKE_CURRENT_BINARY_DIR)
endif()
FetchContent_Declare(
    vg_cmake_kit
    QUIET
    GIT_REPOSITORY      https://github.com/ValeevGroup/kit-cmake.git
    SOURCE_DIR ${${VG_CMAKE_KIT_PREFIX_DIR}}/cmake/vg
    BINARY_DIR ${${VG_CMAKE_KIT_PREFIX_DIR}}/cmake/vg-build
    SUBBUILD_DIR ${${VG_CMAKE_KIT_PREFIX_DIR}}/cmake/vg-subbuild
)
FetchContent_MakeAvailable(vg_cmake_kit)
list(APPEND CMAKE_MODULE_PATH "${vg_cmake_kit_SOURCE_DIR}/modules")

###############################################################################
# Announce ourselves
###############################################################################
# see https://semver.org/
set (ITENSORS_MAJOR_VERSION 4)
set (ITENSORS_MINOR_VERSION 0)
set (ITENSORS_PATCH_VERSION 0)
set (ITENSORS_PRERELEASE_ID beta.1)
set (ITENSORS_BUILD_ID )

set(ITENSORS_VERSION "${ITENSORS_MAJOR_VERSION}.${ITENSORS_MINOR_VERSION}.${ITENSORS_PATCH_VERSION}")
if (ITENSORS_PRERELEASE_ID)
  set(ITENSORS_EXT_VERSION "${ITENSORS_VERSION}-${ITENSORS_PRERELEASE_ID}")
else(ITENSORS_PRERELEASE_ID)
  set(ITENSORS_EXT_VERSION "${ITENSORS_VERSION}")
endif(ITENSORS_PRERELEASE_ID)
if (ITENSORS_BUILD_ID)
  set(ITENSORS_EXT_VERSION "${ITENSORS_EXT_VERSION}+${ITENSORS_BUILD_ID}")
endif(ITENSORS_BUILD_ID)

# extra cmake files are shipped with ITENSORS
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/modules")

include(AppendFlags)

# Extract the git revision tag information
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
  find_package(Git REQUIRED)
  execute_process(
      COMMAND ${GIT_EXECUTABLE} rev-parse -q HEAD
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      OUTPUT_VARIABLE ITENSORS_REVISION )
  string(REGEX MATCH "[0-9a-f]*"
         ITENSORS_REVISION "${ITENSORS_REVISION}")
else()
  set(ITENSORS_REVISION "unknown")
endif()

project(ITENSORS
        VERSION ${ITENSORS_VERSION}
        DESCRIPTION "ITensor: High-Performance Tensor Software Inspired By Tensor Diagrams"
        LANGUAGES CXX
        HOMEPAGE_URL "https://itensor.org/")
enable_language(C)  # C needed even for basic platform introspection

###############################################################################
# Preconditions
###############################################################################
set(TARGET_ARCH "${CMAKE_SYSTEM}-${CMAKE_SYSTEM_PROCESSOR}")

##########################
# Standard build variables
##########################
set(ITENSORS_INSTALL_BINDIR "bin"
    CACHE PATH "ITENSORS BIN install directory")
set(ITENSORS_INSTALL_INCLUDEDIR "include"
    CACHE PATH "ITENSORS INCLUDE install directory")
set(ITENSORS_INSTALL_LIBDIR "lib"
    CACHE PATH "ITENSORS LIB install directory")
set(ITENSORS_INSTALL_SHAREDIR "share/ITENSORS/${ITENSORS_MAJOR_VERSION}.${ITENSORS_MINOR_VERSION}.${ITENSORS_MICRO_VERSION}"
    CACHE PATH "ITENSORS SHARE install directory")
set(ITENSORS_INSTALL_DATADIR "${ITENSORS_INSTALL_SHAREDIR}/data"
    CACHE PATH "ITENSORS DATA install directory")
set(ITENSORS_INSTALL_DOCDIR "${ITENSORS_INSTALL_SHAREDIR}/doc"
    CACHE PATH "ITENSORS DOC install directory")
set(ITENSORS_INSTALL_CMAKEDIR "lib/cmake/ITENSORS"
    CACHE PATH "ITENSORS CMAKE install directory")

# Get standard build variables from the environment if they have not already been set
if(NOT CMAKE_C_FLAGS OR NOT DEFINED CMAKE_C_FLAGS)
  set(CMAKE_C_FLAGS "$ENV{CPPFLAGS}")
  append_flags(CMAKE_C_FLAGS "$ENV{CFLAGS}")
endif()
if(NOT CMAKE_CXX_FLAGS OR NOT DEFINED CMAKE_CXX_FLAGS)
  set(CMAKE_CXX_FLAGS "$ENV{CPPFLAGS}")
  append_flags(CMAKE_CXX_FLAGS "$ENV{CXXFLAGS}")
endif()
if(NOT CMAKE_EXE_LINKER_FLAGS OR NOT DEFINED CMAKE_EXE_LINKER_FLAGS)
  set(CMAKE_EXE_LINKER_FLAGS "$ENV{LDFLAGS}")
endif()
if (NOT CMAKE_CXX_COMPILER)
  message(FATAL_ERROR "C++ compiler not found")
endif()

set(CMAKE_SKIP_RPATH FALSE)

##########################
# We use C++17 features
##########################
# but insist on strict standard
set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ ISO Standard version")
if (NOT(CMAKE_CXX_STANDARD EQUAL 17 OR CMAKE_CXX_STANDARD EQUAL 20))
  message(FATAL_ERROR "C++ 2017 ISO Standard or higher is required to compile ITENSORS")
endif()
# C++20 is only configurable via compile features with cmake 3.12 and older
if (CMAKE_CXX_STANDARD EQUAL 20 AND CMAKE_VERSION VERSION_LESS 3.12.0)
  cmake_minimum_required (VERSION 3.12.0)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF CACHE BOOL  "Whether to use extensions of C++ ISO Standard version")

##########################
# Check type support
include(CheckTypeSize)
##########################

##########################
# Load extra CMake features
##########################

include(GNUInstallDirs)
include(CMakeDependentOption)
include(CMakePackageConfigHelpers)
include(FeatureSummary)

set(MPI_CXX_SKIP_MPICXX TRUE CACHE BOOL "MPI_CXX_SKIP_MPICXX")

##########################
# Configure options
##########################
option(ENABLE_MPI "Enable MPI" ON)
add_feature_info(MPI ENABLE_MPI "Message-Passing Interface supports distributed-memory parallel programs")

option(ENABLE_CUDA "Enable use of CUDA with ITensor" OFF)
add_feature_info(CUDA ENABLE_CUDA "NVIDIA CUDA support for GPU")

option(ENABLE_BLAS_LAPACK "Enable BLAS and LAPACK Libraries" ON)
add_feature_info(BLAS_LAPACK ENABLE_BLAS_LAPACK "Fast standardized linear algebra definitions")

option(ENABLE_HDF5 "Compile ITensor with HDF5 option" OFF)
add_feature_info(HDF5 ENABLE_HDF5 "Storage of data via the HDF5 file protocol")

option(ENABLE_OMP "Compile ITensor with OpenMP option" OFF)
add_feature_info(OMP ENABLE_OMP "Message passing interface which supports unified-memory threadwise-parallel programs")

option(ENABLE_WFN91_LINALG_DISCOVERY_KIT "Use linear algebra discovery kit from github.com/wavefunction91 [recommended]" ON)
add_feature_info(WFN91LinearAlgebraDiscoveryKit ENABLE_WFN91_LINALG_DISCOVERY_KIT "Linear algebra discovery kit from github.com/wavefunction91 supports many more corner cases than the default CMake modules and/or ICL's BLAS++/LAPACK++ modules")

#if(ENABLE_CUDA)
#    option(ENABLE_CUDA_ERROR_CHECK "ITensor will always check errors in CUDA calls" ON)
#    add_feature_info(CUDA_ERROR_CHECK ENABLE_CUDA_ERROR_CHECK "Checks CUDA Error")
#endif()

##########################
# miscellaneous cmake platform-neutral and platform-specific configuration
##########################
set(CMAKE_FIND_NO_INSTALL_PREFIX TRUE)  # do not search in CMAKE_INSTALL_PREFIX
set(CMAKE_SKIP_RPATH FALSE)
set(CMAKE_SKIP_BUILD_RPATH  FALSE)
set(CMAKE_SKIP_INSTALL_RPATH FALSE)
set(CMAKE_NO_SYSTEM_FROM_IMPORTED TRUE)  # do not use -isystem by default to avoid include dir ordering issues as well as other issues like https://gcc.gnu.org/onlinedocs/cpp/System-Headers.html

if(CMAKE_SYSTEM_NAME MATCHES "Darwin")
  # look for frameworks and appbundles last
  set(CMAKE_FIND_FRAMEWORK LAST)
  set(CMAKE_FIND_APPBUNDLE LAST)
endif()

##########################
# Include source dirctories
##########################
include_directories(${PROJECT_SOURCE_DIR}/itensor ${PROJECT_BINARY_DIR}/itensor)

##########################
# external dependencies
##########################
add_library(ITENSORS INTERFACE)
# required deps:
# 1. CUDA first since others may depend on it
if(ENABLE_CUDA)
    include(external/cuda.cmake)
endif(ENABLE_CUDA)

##########################
# discover linear algebra
##########################

include(external/linalgpp.cmake)

##########################
# sources
##########################

add_subdirectory(itensor)

##########################
# pkg-config variables
##########################
foreach(_inc ${ITENSORS_CONFIG_INCLUDE_DIRS})
  append_flags(ITENSORS_PC_CFLAGS "-I${_inc}")
endforeach()
foreach(_lib ${ITENSORS_CONFIG_LIBRARIES})
  append_flags(ITENSORS_PC_LIBS "${_lib}")
endforeach()

##########################
# wrap up
##########################

# Force cache refresh for compile flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}" CACHE STRING "C compile flags" FORCE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" CACHE STRING "C++ compile flags" FORCE)

#CONFIGURE_FILE(
#    ${PROJECT_SOURCE_DIR}/cmake/config.h.in
#    ${PROJECT_BINARY_DIR}/cmake/config.h
#)

CONFIGURE_FILE(
    ${PROJECT_SOURCE_DIR}/itensors.pc.in
    ${PROJECT_BINARY_DIR}/itensors.pc
)

# install config files
install(FILES ${PROJECT_BINARY_DIR}/itensor.pc
        DESTINATION lib/pkgconfig)

## include extra cmake files
#install(FILES
#        "${PROJECT_SOURCE_DIR}/cmake/modules/SanitizeCUDAImplicitDirectories.cmake"
#        DESTINATION "${TILEDARRAY_INSTALL_CMAKEDIR}")

# Create the version file
write_basic_package_version_file(itensor-config-version.cmake
                                 VERSION ${ITENSORS_VERSION} COMPATIBILITY AnyNewerVersion)

# Create the targets file
export(EXPORT itensor
       FILE "${PROJECT_BINARY_DIR}/itensor-targets.cmake")

# Create the configure file
configure_package_config_file(cmake/itensor-config.cmake.in
                              "${PROJECT_BINARY_DIR}/itensor-config.cmake"
                              INSTALL_DESTINATION "${ITENSORS_INSTALL_CMAKEDIR}"
                              PATH_VARS CMAKE_INSTALL_PREFIX ITENSORS_INSTALL_BINDIR
                              ITENSORS_INSTALL_INCLUDEDIR ITENSORS_INSTALL_LIBDIR
                              ITENSORS_INSTALL_DOCDIR ITENSORS_INSTALL_CMAKEDIR)

# Install config, version, and target files
install(EXPORT itensor
        FILE "itensor-targets.cmake"
        DESTINATION "${ITENSORS_INSTALL_CMAKEDIR}"
        COMPONENT ITENSORS)

install(FILES
        "${PROJECT_BINARY_DIR}/itensor-config.cmake"
        "${PROJECT_BINARY_DIR}/itensor-config-version.cmake"
        DESTINATION "${ITENSORS_INSTALL_CMAKEDIR}"
        COMPONENT ITENSORS)


# Add target to allow on-the-fly switching of build type

ADD_CUSTOM_TARGET(debug
                  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_CURRENT_SOURCE_DIR}
                  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_CURRENT_BINARY_DIR} --target all
                  COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
                  )

ADD_CUSTOM_TARGET(release
                  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_CURRENT_SOURCE_DIR}
                  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_CURRENT_BINARY_DIR} --target all
                  COMMENT "Switch CMAKE_BUILD_TYPE to Release"
                  )

feature_summary(WHAT ALL
                DESCRIPTION "=== ITensor Package/Feature Info ===")