//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
#ifndef __ITENSOR_WRAP_H5_H
#define __ITENSOR_WRAP_H5_H

// The purpose of this file is to 
// do any ITensor specific customizations
// of the h5 library, such as putting
// using declarations of certain functions
// to bring them into the itensor namespace

#ifdef ITENSOR_USE_HDF5

#include "itensor/util/h5/h5.hpp"

namespace itensor {

using h5::h5_write;
using h5::h5_write_attribute;
using h5::h5_read;
using h5::h5_read_attribute;


h5::file inline
h5_open(const char* name, char mode)
    {
    return h5::file(name,mode);
    }

h5::file inline
h5_open(std::string const& name, char mode)
    {
    return h5::file(name.c_str(),mode);
    }

}

#endif //ITENSOR_USE_HDF5

#endif
