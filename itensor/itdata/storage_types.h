//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//

//
// To register a new storage type:
//
// (1) Forward declare the storage type
//
// (2) Add a line to the REGISTER_ITDATA_TYPE
//     macro below, following the same format
//
// (3) Include the header file defining
//     the new type
//


//(1) Forward declare storage types
namespace itensor 
{

class ITReal;

class ITCplx;

template<typename T>
class ITDiag;

class ITCombiner;

class IQTData;

}


//(2) Register storage type names
#define REGISTER_ITDATA_TYPES(NewType,EndType) \
        NewType     ITReal             EndType \
        NewType     ITCplx             EndType \
        NewType     ITCombiner         EndType \
        NewType     ITDiag<Real>       EndType \
        NewType     ITDiag<Cplx>       EndType \
        NewType     IQTData            EndType \


//(3) Register header file names
#ifdef REGISTER_ITDATA_HEADER_FILES
#include "itensor/itdata/itreal.h"
#include "itensor/itdata/itcplx.h"
#include "itensor/itdata/itdiag.h"
#include "itensor/itdata/itcombiner.h"
#include "itensor/itdata/iqtdata.h"
#endif
