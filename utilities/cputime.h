//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//

#ifndef _CPUTIME_h
#define _CPUTIME_h
#include <iostream>
#include <iomanip>

double cpu_mytime();
double cpu_mywall();

class cpu_time
    {
public:
    double time;		// in seconds
    double wall;
    cpu_time()
	{ time = cpu_mytime(); wall = cpu_mywall(); }
    void mark() 
	{ time = cpu_mytime(); wall = cpu_mywall(); }
    cpu_time sincemark() const;
    };

std::ostream & operator << (std::ostream & s, const cpu_time & t);
#endif
