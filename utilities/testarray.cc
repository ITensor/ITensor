//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//

// test.cc -- Test the matrix package
#define THIS_IS_MAIN

#include "indent.h"
#include "tarray1.h"
#include "intarray2.h"

void main()
    {
    IntArray1 ia(10);
    ia[4] = 7;
    cout << "ia->Numref() is " << ia->Numref() << iendl;
    IntArray1 ib(ia),id,ie,ig;
    cout << "ia->Numref() is " << ia->Numref() << iendl;
    cout << "ib->Numref() is " << ib->Numref() << iendl;
    cout << "ib = " << iendl << ib << iendl;
    ib->EnLarge(5);
    cout << "ib = " << iendl << ib << iendl;
    cout << "ia = " << iendl << ia << iendl;
    cout << "ib->Size() is " << ib->Size() << iendl;
    cout << "ib->StoreSize() is " << ib->StoreSize() << iendl;
    ia[5] = ib(4);
    cout << "ib = " << iendl << ib << iendl;
    cout << "ia = " << iendl << ia << iendl;
    ib = ia;
    cout << "ib = " << iendl << ib << iendl;
    ie = ig = id = ia;
    cout << "ie = " << iendl << ie << iendl;
    cout << "ie->Numref() is " << ie->Numref() << iendl;
    ie->ReduceDimension(7);
    ie = 4;
    cout << "after ReduceDimension, ie = " << iendl << ie << iendl;
    cout << "ie->Numref() is " << ie->Numref() << iendl;
    cout << "id->Numref() is " << id->Numref() << iendl;
    id[2] = 7;
    cout << "id = " << iendl << id << iendl;
    cout << "id->Numref() is " << id->Numref() << iendl;
    cout << "ia->Numref() is " << ia->Numref() << iendl;
    ie->ReduceDimension(12);
    ie = 2;
    cout << "after ReduceDimension, ie = " << iendl << ie << iendl;
    cout << "ie->Numref() is " << ie->Numref() << iendl;
    id = ie;
    ie->ReSize(10);
    ie = 3;
    cout << "after ReSize, ie = " << iendl << ie << iendl;
    cout << "ie->Numref() is " << ie->Numref() << iendl;
    ie->SetSize(7);
    cout << "after SetSize, ie = " << iendl << ie << iendl;
    cout << "ie->Numref() is " << ie->Numref() << iendl;
    cout << "ia->Size() is " << ia->Size() << iendl;
    cout << "ie->Size() is " << ie->Size() << iendl;
    cout << "ia->Size() + ie->Size() = " << ia->Size() + ie->Size() << iendl;
    ie->ReSize(ia->Size());
    cout << "ia->Size() is " << ia->Size() << iendl;
    cout << "ie->Size() is " << ie->Size() << iendl;
    IntArray2 IA(5,10);
    IA = 0;
    IntArray2 IB(IA);
    IA = 5;
    cout << "IB(2,2) is " << IB(2,2) << iendl;
    cout << "IA(2,2) is " << IA(2,2) << iendl;
    IA(2,2) = 7;
    cout << "IB(2,2) is " << IB(2,2) << iendl;
    cout << "IA(2,2) is " << IA(2,2) << iendl;
    }
