//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef _CPUTIME_h
#define _CPUTIME_h

#include <iostream>

namespace itensor {

double cpu_mytime();
double cpu_mywall();

struct cpu_time
    {
    double time = 0; // in seconds
    double wall = 0;

    cpu_time() { mark(); }

    void 
    mark() { time = cpu_mytime(); wall = cpu_mywall(); }

    cpu_time 
    sincemark() const;
    };

std::ostream& 
operator<<(std::ostream & s, const cpu_time& t);

const double&
firstwall();

std::string 
showtime(double time);

std::ostream& 
operator<<(std::ostream & s, const cpu_time& t);

} //namespace itensor

#endif
