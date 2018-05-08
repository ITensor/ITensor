//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
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
class QMixed;

template<typename T>
class Scalar;

//class ITLazy;


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
QMixed<Real>,
QMixed<Cplx>,
Scalar<Real>,
Scalar<Cplx>
//ITLazy
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
#include "itensor/itdata/qmixed.h"
#include "itensor/itdata/scalar.h"
////#include "itensor/itdata/itlazy.h"
#endif
