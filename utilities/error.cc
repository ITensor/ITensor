//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//

// error.cc -- Print out an error message and abort for debugging

#include "error.h"
#include <cstdlib>
using namespace std;

void 
error(const string& s)
    {
    cerr << endl << s << endl;
    cout << endl << s << endl;
    cout.flush();
    abort();
    }

void 
error(const string& s, int line, const char* file = 0)
    {
    cerr << "From line " << line;
    if(file != 0) cerr << ", file " << file;
    cerr << endl;

    cerr << endl << s << endl;
    cout << endl << s << endl;
    cout.flush();
    cerr.flush();
    abort();
    }

