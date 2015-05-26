//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//

#include "cputime.h"
#include <string>
#include <sstream>
#include <iomanip>
#include <chrono>

using namespace std;
using namespace std::chrono;

#if defined(_MSC_VER)
#include <windows.h>

double cpu_mytime()
    {
    FILETIME ftCreation, ftExit, ftKernel, ftUser;
    SYSTEMTIME stUser;

    HANDLE hProcess = GetCurrentProcess();
    if (GetProcessTimes(hProcess, &ftCreation, &ftExit, &ftKernel, &ftUser) != -1) 
        {
        if (FileTimeToSystemTime(&ftUser, &stUser) != -1)
            {
            return stUser.wHour * 3600.0 + stUser.wMinute * 60.0 + stUser.wSecond + stUser.wMilliseconds * 1E-3;
            }
        }
    return 0.0;
    }
#else
#include <sys/time.h>
#include <sys/resource.h>

double cpu_mytime()		// cpu time in seconds used by program
    {
    struct rusage result;
    getrusage(RUSAGE_SELF,&result);
    return result.ru_utime.tv_sec + 1e-6 * result.ru_utime.tv_usec;
    }
#endif

double cpu_mywall()		// wall time since program started
    {
    static auto initial = steady_clock::now();
    auto current = steady_clock::now();
    auto dur = duration_cast<microseconds>(current-initial);
    return dur.count() * 1.0e-6;
    }

double firstwall = cpu_mywall();	// So cpu_mywall gets called at begining of program

static string showtime(double time)
    {
    if(time < -1.0e-5) 
	{
	cout << "CPU time is negative!" << time << endl;
	cerr << "CPU time is negative!" << time << endl;
	}
    if(time < 0.0) time = 0.0;
    int hours = (int) time/3600;
    time -= hours * 3600;
    int minutes = (int) time/60;
    time -= minutes * 60;
    double seconds = time;
    ostringstream oh;
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
    oh << setprecision(prec) << std::fixed << seconds << "s";
    return oh.str();
    }

ostream & operator << (ostream & s, const cpu_time & t)
    {
    s << "[CPU: " << showtime(t.time) << 
	", Wall: " << showtime(t.wall) << "]";
    return s;
    }

cpu_time cpu_time::sincemark() const
    {
    cpu_time res;
    res.time -= time;
    res.wall -= wall;
    return res;
    }
