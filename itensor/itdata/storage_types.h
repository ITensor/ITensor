//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_STORAGE_TYPES_H_
#define __ITENSOR_STORAGE_TYPES_H_

//
// To register a new storage type:
//
// (1) Include the header file defining the new storage type
//
// (2) Add a new line to the REGISTER macro below, 
//     following the same format
//

//(1)
#include "itensor/itdata/itreal.h"
#include "itensor/itdata/itcplx.h"
#include "itensor/itdata/itdiag.h"
#include "itensor/itdata/itcombiner.h"
#include "itensor/itdata/iqtdata.h"

//(2)
#define REGISTER_TYPES(NewType,EndType)    \
        NewType   ITReal           EndType \
        NewType   ITCplx           EndType \
        NewType   ITCombiner       EndType \
        NewType   ITDiag<Real>     EndType \
        NewType   ITDiag<Complex>  EndType \
        NewType   IQTData          EndType \

#endif
