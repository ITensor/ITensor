//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
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

class ITReal;

class ITCplx;

template<typename T>
class ITDiag;

class ITCombiner;

class IQTReal;

class IQTCombiner;

class IQTDiag;

//class ITLazy;


using 
StorageTypes = TypeList< 
//-----------
//(2) Register storage type names
ITReal,
ITCplx,
ITDiag<Real>,
ITDiag<Cplx>,
ITCombiner,
IQTReal,
IQTCombiner,
IQTDiag
//ITLazy
//-----------
>;

}

//(3) Register header file names
#ifdef REGISTER_ITDATA_HEADER_FILES
#include "itensor/itdata/itreal.h"
#include "itensor/itdata/itcplx.h"
#include "itensor/itdata/itdiag.h"
#include "itensor/itdata/itcombiner.h"
#include "itensor/itdata/iqtreal.h"
#include "itensor/itdata/iqtcombiner.h"
#include "itensor/itdata/iqtdiag.h"
//#include "itensor/itdata/itlazy.h"
#endif
