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
#ifndef __ITENSOR_ALL_MPS_H
#define __ITENSOR_ALL_MPS_H

//
// all_mps.h - convenience header file that
//            includes all headers related
//            to MPS and MPO such as 
//            DMRG, methods for time-evolving
//            MPS, 2D lattice helpers, etc.
//
//          (the headers explicitly included
//          here are not meant to be an 
//          exhaustive list, but are 
//          a minimal set which pull in
//          the key features)
//

#include "all_basic.h"

#include "util/stats.h"

#include "mps/dmrg.h"
#include "mps/tevol.h"
#include "mps/autompo.h"

#include "mps/lattice/square.h"
#include "mps/lattice/triangular.h"

#include "mps/sites/spinhalf.h"
#include "mps/sites/spinone.h"
#include "mps/sites/spintwo.h"
#include "mps/sites/customspin.h"
#include "mps/sites/boson.h"
#include "mps/sites/electron.h"
#include "mps/sites/fermion.h"
#include "mps/sites/tj.h"
#include "mps/sites/Z3.h"

#endif
