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

#include "itensor/all_basic.h"

#include "itensor/util/stats.h"

#include "itensor/mps/dmrg.h"
#include "itensor/mps/tevol.h"
#include "itensor/mps/autompo.h"

#include "itensor/mps/lattice/square.h"
#include "itensor/mps/lattice/triangular.h"

#include "itensor/mps/sites/spinhalf.h"
#include "itensor/mps/sites/spinone.h"
#include "itensor/mps/sites/spintwo.h"
#include "itensor/mps/sites/customspin.h"
#include "itensor/mps/sites/boson.h"
#include "itensor/mps/sites/electron.h"
#include "itensor/mps/sites/fermion.h"
#include "itensor/mps/sites/tj.h"
#include "itensor/mps/sites/Z3.h"

#endif
