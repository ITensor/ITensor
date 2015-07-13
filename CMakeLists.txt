cmake_minimum_required (VERSION 2.8)
project (ITensor)

# Disable build in source (not to overlap with existing Makefile)
if (${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
    message(FATAL_ERROR "In source builds are disabled. Please use a separate build directory.")
endif()
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}")

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

option (CXX11 OFF)
if (CXX11) 
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
endif (CXX11)


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
    message(STATUS "Setting lapack platform : Intel")
    add_definitions(-DPLATFORM_mkl)
    find_package(MKL REQUIRED)
    list(APPEND CMAKE_LINKER_FLAGS ${MKL_LINKER_FLAGS})
    include_directories(${MKL_INCLUDE_DIRS})
    link_libraries(${MKL_LIBRARIES})
else()
  if(${BLA_VENDOR} MATCHES ACML)
      message(STATUS "Setting lapack platform : ACML")
      add_definitions(-PLATFORM_mkl)
  elseif(${BLA_VENDOR} MATCHES Apple)
      message(STATUS "Setting lapack platform : Apple")
      add_definitions(-DPLATFORM_macos)
  else()
      message(STATUS "Setting lapack platform : Generic")
      add_definitions(-DPLATFORM_lapack)
  endif()
  find_package(LAPACK REQUIRED)
  list(APPEND CMAKE_LINKER_FLAGS ${LAPACK_LINKER_FLAGS})
  link_libraries(${LAPACK_LIBRARIES})
endif()

# build itensor parts
add_subdirectory(utilities)
add_subdirectory(matrix)
add_subdirectory(itensor)

# Enable testing
option(Testing "Enable testing" ON)
if (Testing)
    enable_testing()
    add_subdirectory(unittest)
endif (Testing)

# Build example
option(Tutorial "Build tutorial" ON)
if (Tutorial)
    message(STATUS "Building tutorial")
    add_subdirectory(tutorial)
endif (Tutorial)

# Build sample
option(Sample "Build samples" ON)
if (Sample)
    message(STATUS "Building samples")
    add_subdirectory(sample)
endif (Sample)

# Install pkg-config file
#configure_file(ITensor.pc.in ITensor.pc @ONLY)
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/ITensor.pc DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/pkgconfig)
