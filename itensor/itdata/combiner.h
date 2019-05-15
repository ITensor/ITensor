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
#ifndef __ITENSOR_ITCOMBINER_H
#define __ITENSOR_ITCOMBINER_H

#include "itensor/itdata/dense.h"

namespace itensor {

class QN;

class Combiner
    {
    public:

    Combiner() { }

    };

const char*
typeNameOf(Combiner const& d);

void inline
read(std::istream& s, Combiner& dat) { }

void inline
write(std::ostream& s, Combiner const& dat) { }

Cplx
doTask(GetElt const& g, Combiner const& c);

Real
doTask(NormNoScale, Combiner const& d);

void
doTask(Conj, Combiner const& d);

template<typename V>
void
doTask(Contract & C,
       Dense<V>   const& d,
       Combiner const& cmb,
       ManageStore     & m);

template<typename V>
void
doTask(Contract & C,
       Combiner const& cmb,
       Dense<V>   const& d,
       ManageStore     & m);

bool
doTask(CheckComplex, Combiner const& d);

void
doTask(PrintIT& P, Combiner const& d);

auto inline
doTask(StorageType const& S, Combiner const& d) 
    ->StorageType::Type 
    { return StorageType::Combiner; }

QN 
doTask(CalcDiv const& C, Combiner const& d);

} //namespace itensor

#endif

