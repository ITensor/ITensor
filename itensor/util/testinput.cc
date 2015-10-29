//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//

// test.cc -- Test the matrix package
#define THIS_IS_MAIN

#include "input.h"
#include <Accelerate/Accelerate.h>

int main(int argc, char* argv[])
    {
    if(argc != 2) 
	cout << "Need to input argument on command line\n";
    string infilename(argv[1]);
    InputFile infile(infilename);
    cout << infile;
    InputGroup lattice(infile,"lattice", "lattice parameters");
    int nx,ny;
    lattice.GetIntM("nx",nx);
    lattice.GetIntM("ny",ny);
    InputGroup pbc(lattice,"pbc","group for periodic boundaries");
    int dopbc; 
    if(!pbc.GetYesNo("dopbc",dopbc,"do periodic boundaries"))
	cout << "didnt get dopbc " << endl;
    cout << "sizeof(long int) is " << sizeof(long int) << endl;
    cout << "sizeof(__CLPK_integer) is " << sizeof(__CLPK_integer) << endl;

    }
