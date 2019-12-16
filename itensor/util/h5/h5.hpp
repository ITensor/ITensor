#pragma once

#include "std_addons/complex.hpp"

#include "./file.hpp"
#include "./group.hpp"
#include "./format.hpp"
#include "./scalar.hpp"
#include "./generic.hpp"
#include "./stl/string.hpp"
#include "./stl/vector.hpp"
//#include "./stl/map.hpp"
//#include "./stl/pair.hpp"
//#include "./stl/tuple.hpp"
//#include "./stl/optional.hpp"
//#include "./stl/variant.hpp"

// FIXME : Still needed ?
// for python code generator, we need to know what has to been included.
//#define TRIQS_INCLUDED_H5

// in some old version of hdf5 (Ubuntu 12.04 e.g.), the macro is not yet defined.
#ifndef H5_VERSION_GE

#define H5_VERSION_GE(Maj, Min, Rel)                                                                                                                 \
  (((H5_VERS_MAJOR == Maj) && (H5_VERS_MINOR == Min) && (H5_VERS_RELEASE >= Rel)) || ((H5_VERS_MAJOR == Maj) && (H5_VERS_MINOR > Min))               \
   || (H5_VERS_MAJOR > Maj))

#endif
