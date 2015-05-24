//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_STORAGE_TYPES_H_
#define __ITENSOR_STORAGE_TYPES_H_

//
// To register a new storage type:
//
// (1) Include the header defining the storage type below
//
// (2) Add a new line to the REGISTER macro below, following the same format
//     and no trailing \ on the last line.
//

//(1)
#include "itensor/itdata/itreal.h"
#include "itensor/itdata/itcplx.h"
#include "itensor/itdata/itdiag.h"
#include "itensor/itdata/itcombiner.h"
#include "itensor/itdata/iqtdata.h"

//(2)
#define REGISTER_TYPES(Type_,EndType_)   \
        Type_   ITReal           EndType_ \
        Type_   ITCplx           EndType_ \
        Type_   ITCombiner       EndType_ \
        Type_   ITDiag<Real>     EndType_ \
        Type_   ITDiag<Complex>  EndType_ \
        Type_   IQTData<Real>    EndType_

#endif
