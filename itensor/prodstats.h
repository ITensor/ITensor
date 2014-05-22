//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_PRODSTATS_H
#define __ITENSOR_PRODSTATS_H


//
// Usage of Timers feature:
// 1. Uncomment #define COLLECT_PRODSTATS line below and recompile library.
// 2. In either library or client code, start and stop timers around sections
//    of code using START_TIMER(#); ... STOP_TIMER(#); where # is any integer
//    from 1 to 100.
// 3. At end of main function of client code, add the line PRINT_TIMERS to see
//    a printout of the times of each section.
//


//#define COLLECT_PRODSTATS

#ifdef COLLECT_PRODSTATS
#define DO_IF_PS(x) { x }
#else
#define DO_IF_PS(x) { }
#endif
#ifdef COLLECT_PRODSTATS
#define START_TIMER(x) { Timers::global().start_section(x); }
#else
#define START_TIMER(x) { }
#endif
#ifdef COLLECT_PRODSTATS
#define STOP_TIMER(x) { Timers::global().finish_section(x); }
#else
#define STOP_TIMER(x) { }
#endif

#ifdef COLLECT_PRODSTATS
#define PRINT_TIMERS { Timers::global().print(); }
#else
#define PRINT_TIMERS { }
#endif

#ifdef COLLECT_PRODSTATS

#include <chrono>
#include <vector>
#include <map>
#include <iostream>
#include "error.h"

#define NTIMERS 100

namespace itensor {

class Timers
        {
        typedef std::chrono::high_resolution_clock
        Clock;

        typedef Clock::time_point
        TimePoint;

        typedef std::chrono::microseconds
        TimeUnits;

        ////
        std::vector<double> time;
        std::vector<int>  tcount;
        std::vector<TimePoint> marker;
        std::vector<bool> timer_running;
        ////

        public:

        Timers()
            :
            time(NTIMERS,0.0),
            tcount(NTIMERS,0),
            marker(NTIMERS),
            timer_running(NTIMERS,false)
            { }

        void 
        start_section(int j) 
            { 
            if(timer_running.at(j)) Error("Timer already running.");
            timer_running.at(j) = true;
            marker.at(j) = Clock::now();
            }

        void 
        finish_section(int j)
            {
            if(!timer_running.at(j)) Error("No timer running.");
            timer_running.at(j) = false; 
            auto end_time = Clock::now();
            time.at(j) += 1E-6*std::chrono::duration_cast<TimeUnits>(end_time-marker.at(j)).count();
            tcount.at(j) += 1; 
            }

        void 
        print() const
            {
            std::cout << "Timers ----------------- " << std::endl;
            for(int j = 0; j < (int) time.size(); ++j)
                {
                const double count = tcount.at(j);
                if(count > 0) 
                    {
                    printfln("Section %d, Average Time = %.6f\n",j,(time.at(j)/count));
                    }
                }

            for(int j = 0; j < (int) time.size(); ++j)
                {
                if(tcount.at(j) > 0) 
                    {
                    printfln("Section %d, Total Time = %.6f\n",j,time.at(j));
                    }
                }
            std::cout << std::endl;
            }

    static Timers& global()
        {
        static Timers global_;
        return global_;
        }

    };

/*
class Prodstats
        {
        public:
        typedef std::map<std::pair<int,int>,int>
        CountMap;
        CountMap global;
        CountMap ps32;
        std::vector<int> perms_of_3;
        std::vector<int> perms_of_4;
        std::vector<int> perms_of_5;
        std::vector<int> perms_of_6;
        int total, did_matrix;
        int c1,c2,c3,c4;

        Prodstats()
            {
            //for(int i = 0; i <= 20; ++i)
            //for(int j = i; j <= 20; ++j)
                //global[std::make_pair(j,i)] = 0;
            //total = 0;
            //did_matrix = 0;
            //c1 = c2 = c3 = c4 = 0;
            //perms_of_3 = std::vector<int>(81,0);
            //perms_of_4 = std::vector<int>(256,0);
            //perms_of_5 = std::vector<int>(3125,0);
            //perms_of_6 = std::vector<int>(46656,0);

            }


        void 
        print() const
            {
            //std::cout << "\n-------- Product Statistics ----------\n";
            //std::cout << "Global Count: " << std::endl;
            //for(CountMap::const_iterator pp = global.begin(); pp != global.end(); ++pp)
            //    printfln("(%d,%d) = %d",pp->first.first,pp->first.second,pp->second);
            //std::cout << "Total = " << total << std::endl;
            //printfln("# Matrices = %d (%.2f%%)",(did_matrix),(total == 0 ? 0 : (100.0*(1.*did_matrix/(2*total)))));

            //std::cout << "# Case 1 = " << c1 << std::endl;
            //std::cout << "# Case 2 = " << c2 << std::endl;
            //std::cout << "# Case 3 = " << c3 << std::endl;
            //std::cout << "# Case 4 = " << c4 << std::endl;

            //std::cout << "Permutations of 3 Count: " << std::endl;
            //for(int j = 0; j < (int) perms_of_3.size(); ++j)
            //    {
            //    if(perms_of_3[j] == 0) continue;
            //    int c = j;
            //    int i3 = (c%3 == 0 ? 3 : c%3);
            //    c = (c-i3)/3+1;
            //    int i2 = (c%3 == 0 ? 3 : c%3);
            //    c = (c-i2)/3+1;
            //    int i1 = (c%3 == 0 ? 3 : c%3);
            //    int idx = ((i1-1)*3+i2-1)*3+i3;
            //    if(idx != j) std::cout << "Incorrect idx val (perms of 3)." << std::endl;
            //    printfln("(%02d) %d, %d, %d = %d",j,i1,i2,i3,perms_of_3[j]);
            //    }

            //std::cout << "Permutations of 4 Count: " << std::endl;
            //for(int j = 0; j < (int) perms_of_4.size(); ++j)
            //    {
            //    if(perms_of_4[j] == 0) continue;
            //    int c = j;
            //    int i4 = (c%4 == 0 ? 4 : c%4);
            //    c = (c-i4)/4+1;
            //    int i3 = (c%4 == 0 ? 4 : c%4);
            //    c = (c-i3)/4+1;
            //    int i2 = (c%4 == 0 ? 4 : c%4);
            //    c = (c-i2)/4+1;
            //    int i1 = (c%4 == 0 ? 4 : c%4);
            //    int idx = (((i1-1)*4+i2-1)*4+i3-1)*4+i4;
            //    if(idx != j) std::cout << "Incorrect idx val (perms of 4)." << std::endl;
            //    printfln("(%02d) %d, %d, %d, %d = %d",j,i1,i2,i3,i4,perms_of_4[j]);
            //    }

            //std::cout << "Permutations of 5 Count: " << std::endl;
            //for(int j = 0; j < (int) perms_of_5.size(); ++j)
            //    {
            //    if(perms_of_5[j] == 0) continue;
            //    int c = j;
            //    int i5 = (c%5 == 0 ? 5 : c%5);
            //    c = (c-i5)/5+1;
            //    int i4 = (c%5 == 0 ? 5 : c%5);
            //    c = (c-i4)/5+1;
            //    int i3 = (c%5 == 0 ? 5 : c%5);
            //    c = (c-i3)/5+1;
            //    int i2 = (c%5 == 0 ? 5 : c%5);
            //    c = (c-i2)/5+1;
            //    int i1 = (c%5 == 0 ? 5 : c%5);
            //    int idx = ((((i1-1)*5+i2-1)*5+i3-1)*5+i4-1)*5+i5;
            //    if(idx != j) std::cout << "Incorrect idx val (perms of 5)." << std::endl;
            //    printfln("(%02d) %d, %d, %d, %d, %d = %d",j,i1,i2,i3,i4,i5,perms_of_5[j]);
            //    }

            //std::cout << "Permutations of 6 Count: " << std::endl;
            //for(int j = 0; j < (int) perms_of_6.size(); ++j)
            //    {
            //    //printfln("po6[%d] = %d",j,perms_of_6[j]);
            //    if(perms_of_6[j] == 0) continue;
            //    int c = j;
            //    int i6 = (c%6 == 0 ? 6 : c%6);
            //    c = (c-i6)/6+1;
            //    int i5 = (c%6 == 0 ? 6 : c%6);
            //    c = (c-i5)/6+1;
            //    int i4 = (c%6 == 0 ? 6 : c%6);
            //    c = (c-i4)/6+1;
            //    int i3 = (c%6 == 0 ? 6 : c%6);
            //    c = (c-i3)/6+1;
            //    int i2 = (c%6 == 0 ? 6 : c%6);
            //    c = (c-i2)/6+1;
            //    int i1 = (c%6 == 0 ? 6 : c%6);
            //    int idx = (((((i1-1)*6+i2-1)*6+i3-1)*6+i4-1)*6+i5-1)*6+i6;
            //    if(idx != j) std::cout << "Incorrect idx val (perms of 6)." << std::endl;
            //    printfln("(%02d) %d, %d, %d, %d, %d, %d = %d",j,i1,i2,i3,i4,i5,i6,perms_of_6[j]);
            //    }

            }

    static Prodstats& stats()
        {
        static Prodstats stats_;
        return stats_;
        }

    };
*/

}; //namespace itensor

#endif //COLLECT_PRODSTATS

#endif
