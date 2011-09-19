#include "cputime.h"
using std::cout;
using std::cerr;
using std::endl;
using std::ostream;
using std::ostringstream;
using std::setprecision;

#define __alpha
#ifdef __alpha
#include <sys/time.h>
#include <sys/resource.h>
#include <string>
#include <sstream>

double mytime()
    {
    struct rusage result;
    getrusage(RUSAGE_SELF,&result);
    return result.ru_utime.tv_sec + 1e-6 * result.ru_utime.tv_usec;
    }

#else

#include <time.h>

double mytime()
    { return clock() * (1.0 / CLOCKS_PER_SEC); }

#endif

ostream & operator << (ostream & s, const cpu_time & t)
    {
    double time = t.time;
    if(time < -1.0e-5) 
	{
	cout << "CPU time is negative!" << time << endl;
	s << "CPU time is negative!" << time << endl;
	cerr << "CPU time is negative!" << time << endl;
	}
    if(time < 0.0) time = 0.0;
    int hours = (int) time/3600;
    time -= hours * 3600;
    int minutes = (int) time/60;
    time -= minutes * 60;
    double seconds = time;
    if(hours != 0) s << hours << (hours == 1 ? " hour, " : " hours, ");
    if(minutes != 0) s << minutes << 
    		(minutes == 1 ? " minute, " : " minutes, ");
    ostringstream oh;
    oh << setprecision(3) << seconds << " seconds";
    s << oh.str();
    // s << "(" << t.time << ")";
    return s;
    }

cpu_time cpu_time::sincemark()
    {
    cpu_time res;
    res.time -= time;
    return res;
    }
    
    


