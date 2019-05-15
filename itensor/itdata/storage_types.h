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
#include "itensor/types.h"
#include "itensor/util/typelist.h"

namespace itensor 
{

//
// To register a new storage type:
//
// (1) Forward declare the storage type
//
// (2) Add a line to the definition of StorageTypes
//     below, following the same format
//
// (3) Include the header file defining
//     the new type
//

//(1) Forward declare storage types

template<typename T>
class Dense;

template<typename T>
class Diag;

class Combiner;

template<typename T>
class QDense;

class QCombiner;

template<typename T>
class QDiag;

template<typename T>
class Scalar;



using 
StorageTypes = TypeList< 
//-----------
//(2) Register storage type names
Dense<Real>,
Dense<Cplx>,
Combiner,
Diag<Real>,
Diag<Cplx>,
QDense<Real>,
QDense<Cplx>,
QCombiner,
QDiag<Real>,
QDiag<Cplx>,
Scalar<Real>,
Scalar<Cplx>
//-----------
>;

}

//(3) Register header file names
#ifdef REGISTER_ITDATA_HEADER_FILES
#include "itensor/itdata/dense.h"
#include "itensor/itdata/combiner.h"
#include "itensor/itdata/diag.h"
#include "itensor/itdata/qdense.h"
#include "itensor/itdata/qcombiner.h"
#include "itensor/itdata/qdiag.h"
#include "itensor/itdata/scalar.h"
#endif
