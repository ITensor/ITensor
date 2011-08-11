// input.h -- classes for reading from input files

#ifndef _input_h
#define _input_h
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
//#include "String.h"
#include <cctype>
#define SP << " " <<
typedef double Real;
void error(const char* s);

using std::cout;
using std::cerr;
using std::endl;
using std::ostream;
using std::ostringstream;
using std::setprecision;
using std::string;
using std::istream;
using std::ifstream;

class InputFile
    {
public:
    string filename;
    ifstream file;
    int opened;
    InputFile(string fname) : filename(fname), opened(0) {}
    void open();
    void close();
    };
ostream & operator << (ostream &s, InputFile &a);

typedef long int lint;

class InputGroup
    {
public:
    InputFile & infile;
    InputGroup * parent;
    string name;
    bool quiet;
    InputGroup(InputFile& inf, string nam,const char* c = 0)
		: infile(inf), parent(0), name(nam), quiet(false)
	{
	cout << "Making input group " << name;
	if(c) cout << ": " << c;
	cout << endl;
	}
    InputGroup(InputGroup& par, string nam,const char* c = 0)
		: infile(par.infile), parent(&par), name(nam), quiet(false)
	{
	cout << "Making input group " << parent->name << "." << name;
	if(c) cout << ": " << c;
	cout << endl;
	}

    int GotoGroup();		// Goes to group, then eats "{" + whitespace
    int GotoToken(string s);	// Goes to Token, then eats "=" + whitespace

// The following go to s, and read into i,r,t, or yes, printing c.

    int GetInt(string s, int& i,const char* c = 0);
    int GetLong(string s,lint& i,const char* c = 0);
    int GetReal(string s, Real& r,const char* c = 0);	
    int GetString(string s, string& t,const char* c = 0);
    int GetYesNo(string s, int& yes,const char* c = 0);	 // understands yes/no

// The following are mandatory versions; if they doesn't get it, we quit
    void GetIntM(string s, int& i,const char* c = 0);	
    void GetLongM(string s, lint& i,const char* c = 0);	
    void GetRealM(string s, Real& r,const char* c = 0);
    void GetStringM(string s, string& t,const char* c = 0);
    void GetYesNoM(string s, int& yes,const char* c = 0);

    void SkipLine();
    };

    /*
InputGroup(InputFile& inf, String nam,const char* c)
	    : infile(inf), name(nam), parent(0) 
    {
    cout << "Making input group " << name;
    if(c) cout << ", " << c;
    cout << endl;
    }
    */

/* 
To read in a table:

in input file:

tablename
    {
    a	x	y	z
    1	4.0	5.0	7.0
    2	4.0	2.0	7.0
    4	3.0	5.0	7.0
    5	4.0	5.0	1.0
    }

Then in program:
    InputGroup table(parent,"tablename");
    if(table.GotoGroup())
	{
	table.SkipLine();
	for(int i = 1; i <= n; i++)
	    table.infile.file >> a[i] >> x[i] >> y[i] >> z[i];
	}
*/

int gettoken(istream& is, string& s);

#endif
