// error.cc -- Print out an error message and abort for debugging

#include "error.h"
using std::cout;
using std::cerr;
using std::endl;

void error(const char* s)
{
    cerr << endl << s << endl;
    cout << endl << s << endl;
    cout.flush();
    abort();
}

void error(const char* s, int line, const char* file = 0)
{
    cerr << "From line " << line;
    if(file != 0) cerr << ", file " << file;
    cerr << endl;

    cerr << endl << s << endl;
    cout << endl << s << endl;
    cout.flush();
    throw ITError(s);
}

