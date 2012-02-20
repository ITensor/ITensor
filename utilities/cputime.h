//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//

#ifndef _CPUTIME_h
#define _CPUTIME_h
#include <iostream>
#include <iomanip>

double mytime();

class init_time
    {
public:
    double dummy;
    init_time() { dummy = mytime(); }
    friend class cpu_time;
    };

class cpu_time
    {
public:
    double time;		// in seconds
    friend std::ostream & operator << (std::ostream & s, const cpu_time & t);
    cpu_time()
	{ time = mytime(); }
    void mark() 
	{ time = mytime(); }
    cpu_time sincemark();

    static init_time& init()
        {
        static init_time init_;
        return init_;
        }
    };

#endif
