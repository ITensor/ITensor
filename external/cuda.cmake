# cmake 3.17 decouples C++ and CUDA standards, see https://gitlab.kitware.com/cmake/cmake/issues/19123
# cmake 3.18 knows that CUDA 11 provides cuda_std_17
cmake_minimum_required(VERSION 3.18.0)
set(CMAKE_CUDA_STANDARD 20)
set(CMAKE_CUDA_EXTENSIONS OFF)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
# N.B. need relaxed constexpr for std::complex
# see https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#constexpr-functions%5B/url%5D:
if (DEFINED CMAKE_CUDA_FLAGS)
  set(CMAKE_CUDA_FLAGS "--expt-relaxed-constexpr ${CMAKE_CUDA_FLAGS}")
else()
  set(CMAKE_CUDA_FLAGS "--expt-relaxed-constexpr")
endif()
enable_language(CUDA)

set(CUDA_FOUND TRUE)
set(ITENSOR_HAS_CUDA 1 CACHE BOOL "Whether ITensor has CUDA support")

if(ENABLE_CUDA_ERROR_CHECK)
  set (ITENSOR_CHECK_CUDA_ERROR 1)
endif(ENABLE_CUDA_ERROR_CHECK)

# find CUDA toolkit
# NB CUDAToolkit does NOT have COMPONENTS
find_package(CUDAToolkit REQUIRED)

foreach (library cublas;nvToolsExt)
  if (NOT TARGET CUDA::${library})
    message(FATAL_ERROR "CUDA::${library} not found")
  endif()
endforeach()

if (NOT DEFINED CUDAToolkit_ROOT)
  get_filename_component(CUDAToolkit_ROOT "${CUDAToolkit_INCLUDE_DIR}/../" ABSOLUTE CACHE)
endif(NOT DEFINED CUDAToolkit_ROOT)

# sanitize implicit dirs if CUDA host compiler != C++ compiler
message(STATUS "CMAKE Implicit Include Directories: ${CMAKE_CUDA_IMPLICIT_INCLUDE_DIRECTORIES}")
message(STATUS "CMAKE Implicit Link Directories: ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES}")
include(SanitizeCUDA)
sanitize_cuda_implicit_directories()
message(STATUS "CMAKE Implicit Include Directories: ${CMAKE_CUDA_IMPLICIT_INCLUDE_DIRECTORIES}")
message(STATUS "CMAKE Implicit Link Directories: ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES}")

