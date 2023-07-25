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
else(ITENSOR_PRERELEASE_ID)
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
        DESCRIPTION "ITENSORS: Tensor network library"
        LANGUAGES CXX
        HOMEPAGE_URL "https://itensors.org")
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

# Check type support
include(CheckTypeSize)