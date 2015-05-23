//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef _CPUTIME_h
#define _CPUTIME_h

#include <iostream>
#include <iomanip>
#include <sstream>
#include <chrono>
//#include <sys/time.h>
//#include <sys/resource.h>

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

std::ostream& 
operator<<(std::ostream & s, const cpu_time& t);

//
// Implementations
//

double inline
cpu_mytime() // cpu time in seconds used by program
    {
    struct rusage result;
    getrusage(RUSAGE_SELF,&result);
    return result.ru_utime.tv_sec + 1e-6 * result.ru_utime.tv_usec;
    }

double inline
cpu_mywall() // wall time since program started
    {
    static auto initial = std::chrono::steady_clock::now();
    auto current = std::chrono::steady_clock::now();
    auto dur = std::chrono::duration_cast<std::chrono::microseconds>(current-initial);
    return dur.count() * 1.0e-6;
    }

inline const double&
firstwall()
    {
    static auto firstwall_ = cpu_mywall();
    return firstwall_;
    }

std::string inline
showtime(double time)
    {
    if(time < -1.0e-5) 
        {
        std::cout << "CPU time is negative!" << time << std::endl;
        std::cerr << "CPU time is negative!" << time << std::endl;
        }
    if(time < 0.0) time = 0.0;
    int hours = (int) time/3600;
    time -= hours * 3600;
    int minutes = (int) time/60;
    time -= minutes * 60;
    double seconds = time;
    std::ostringstream oh;
    //if(hours != 0) oh << hours << (hours == 1 ? " hour, " : " hours, ");
    if(hours != 0) oh << std::fixed << hours << "h, ";
    if(minutes != 0) oh << std::fixed << minutes << 
    		(minutes == 1 ? "m, " : "m, ");
    int prec = 6;
    if(seconds > 0.001) prec = 5;
    if(seconds > 0.01) prec = 4;
    if(seconds > 0.1) prec = 3;
    if(seconds >= 10) prec = 2;
    if(minutes > 0) prec = 1;
    if(minutes >= 10) prec = 0;
    if(hours > 0) prec = 0;
    oh << std::setprecision(prec) << std::fixed << seconds << "s";
    return oh.str();
    }

inline std::ostream& 
operator<<(std::ostream & s, const cpu_time& t)
    {
    s << "[CPU: " << showtime(t.time) << 
	", Wall: " << showtime(t.wall) << "]";
    return s;
    }

cpu_time inline
cpu_time::sincemark() const
    {
    cpu_time res;
    res.time -= time;
    res.wall -= wall;
    return res;
    }

#endif
