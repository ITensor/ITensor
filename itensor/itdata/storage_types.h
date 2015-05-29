//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_STORAGE_TYPES_H_
#define __ITENSOR_STORAGE_TYPES_H_

//
// To register a new storage type:
//
// (1) Forward declare the storage type
//
// (2) Add a line to the REGISTER macro
//     below, following the same format
//
// (3) Include the header file in
//     which the new type is defined
//

//(1)
namespace itensor {

class ITReal;
class ITCplx;
template<typename T>
class ITDiag;
class ITCombiner;
class IQTData;

} //namespace itensor


//(2)
#define REGISTER_TYPES(NewType,EndType)    \
        NewType   ITReal           EndType \
        NewType   ITCplx           EndType \
        NewType   ITCombiner       EndType \
        NewType   ITDiag<Real>     EndType \
        NewType   ITDiag<Cplx>     EndType \
        NewType   IQTData          EndType  


//(3)
//#ifdef STORAGE_INCLUDES
//
//#include "itensor/itdata/itreal.h"
//#include "itensor/itdata/itcplx.h"
//#include "itensor/itdata/itdiag.h"
//#include "itensor/itdata/itcombiner.h"
//#include "itensor/itdata/iqtdata.h"
//
//#endif//STORAGE_INCLUDES

#endif
