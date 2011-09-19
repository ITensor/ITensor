#include <iostream>
#include <sstream>

#ifndef _indent_h
#define _indent_h 1

// indent.h -- Steve's Indent and PrintLevel classes

class Indent	// new-line and indent
    {
private:
    static int indentlev;
public:
    static int indentlevel() {return indentlev;}
    static void indent() {indentlev++;}
    static void unindent() 
	    {indentlev--; if(indentlev < 0)indentlev = 0;}
    };

// std::ostream io manipulator to do newline and indent

inline std::ostream&	iendl(std::ostream &s)
    {
    s << std::endl;
    int i;
    for(i=0; i < Indent::indentlevel()-1; i += 2)
	    s << "	";
    if(i == Indent::indentlevel()-1)
	s << "    ";
    return s;
    }

// std::istream io manipulator to ignore rest of line

inline std::istream& nextline(std::istream &is)
    {
    char c;
    while(is.get(c) && c != '\n') ;
    return is;
    }

class PrintLevel	// what to print
    {
private:
    static int printlev;
public:
    static void setprint(int i)
	{ printlev = i; }
    static int print(int i)
	{ return i <= printlev; }
    static int printlevel()
	{ return printlev; }
    };

inline int print(int i) {return PrintLevel::print(i);}
inline int printlevel() {return PrintLevel::printlevel();}

#ifdef THIS_IS_MAIN
int Indent::indentlev = 0;
int PrintLevel::printlev = 0;
#endif

#endif
