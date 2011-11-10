#ifndef __ITENSOR_PRODSTATS_H
#define __ITENSOR_PRODSTATS_H

//#define COLLECT_PRODSTATS

#ifdef COLLECT_PRODSTATS
#define DO_IF_PS(x) { x }
#else
#define DO_IF_PS(x) { }
#endif
#ifdef COLLECT_PRODSTATS
#define START_TIMER(x) { prodstats.start_section(x); }
#else
#define START_TIMER(x) { }
#endif
#ifdef COLLECT_PRODSTATS
#define STOP_TIMER(x) { prodstats.finish_section(x); }
#else
#define STOP_TIMER(x) { }
#endif

#ifdef COLLECT_PRODSTATS
#include <cputime.h>

#include <map>
#define NTIMERS 70
class Prodstats
{
    std::vector<Real> time;
    std::vector<int>  tcount;
    std::vector<cpu_time> cpu;
    std::vector<bool> timer_running;
public:
    typedef std::pair<std::pair<int,int>,int> gitertype;
    std::map<std::pair<int,int>,int> global;
    std::map<std::pair<int,int>,int> ps32;
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
        total = 0;
        did_matrix = 0;
        c1 = c2 = c3 = c4 = 0;
        perms_of_3 = std::vector<int>(81,0);
        perms_of_4 = std::vector<int>(256,0);
        perms_of_5 = std::vector<int>(3125,0);
        perms_of_6 = std::vector<int>(46656,0);

        time = std::vector<Real>(NTIMERS,0);
        tcount = std::vector<int>(NTIMERS,0);
        timer_running = std::vector<bool>(NTIMERS,false);
        cpu = std::vector<cpu_time>(NTIMERS);
    }

    void start_section(int j) 
    { 
        if(timer_running.at(j)) Error("Timer already running.");
        timer_running.at(j) = true;
        cpu.at(j) = cpu_time();
    }

    void finish_section(int j)
    {
        if(!timer_running.at(j)) Error("No timer running.");
        timer_running.at(j) = false; 
        cpu_time since(cpu.at(j).sincemark());
        time.at(j) += since.time; 
        tcount.at(j) += 1; 
    }

    void print() const
    {
        std::cerr << "\n-------- Product Statistics ----------\n";
        std::cerr << "Global Count: " << std::endl;
        Foreach(gitertype pp, global) std::cerr << boost::format("(%d,%d) = %d\n")%pp.first.first%pp.first.second%pp.second;
        std::cerr << "Total = " << total << std::endl;
        std::cerr << boost::format("# Matrices = %d (%.2f%%)\n") % (did_matrix) % (total == 0 ? 0 : (100.0*(1.*did_matrix/(2*total))));

        std::cerr << "# Case 1 = " << c1 << std::endl;
        std::cerr << "# Case 2 = " << c2 << std::endl;
        std::cerr << "# Case 3 = " << c3 << std::endl;
        std::cerr << "# Case 4 = " << c4 << std::endl;

        std::cerr << "Permutations of 3 Count: " << std::endl;
        for(int j = 0; j < (int) perms_of_3.size(); ++j)
        {
            if(perms_of_3[j] == 0) continue;
            int c = j;
            int i3 = (c%3 == 0 ? 3 : c%3);
            c = (c-i3)/3+1;
            int i2 = (c%3 == 0 ? 3 : c%3);
            c = (c-i2)/3+1;
            int i1 = (c%3 == 0 ? 3 : c%3);
            int idx = ((i1-1)*3+i2-1)*3+i3;
            if(idx != j) std::cerr << "Incorrect idx val (perms of 3)." << std::endl;
            std::cerr << boost::format("(%02d) %d, %d, %d = %d\n") % j % i1 % i2 % i3 % perms_of_3[j];
        }

        std::cerr << "Permutations of 4 Count: " << std::endl;
        for(int j = 0; j < (int) perms_of_4.size(); ++j)
        {
            if(perms_of_4[j] == 0) continue;
            int c = j;
            int i4 = (c%4 == 0 ? 4 : c%4);
            c = (c-i4)/4+1;
            int i3 = (c%4 == 0 ? 4 : c%4);
            c = (c-i3)/4+1;
            int i2 = (c%4 == 0 ? 4 : c%4);
            c = (c-i2)/4+1;
            int i1 = (c%4 == 0 ? 4 : c%4);
            int idx = (((i1-1)*4+i2-1)*4+i3-1)*4+i4;
            if(idx != j) std::cerr << "Incorrect idx val (perms of 4)." << std::endl;
            std::cerr << boost::format("(%02d) %d, %d, %d, %d = %d\n") % j % i1 % i2 % i3 % i4 % perms_of_4[j];
        }

        std::cerr << "Permutations of 5 Count: " << std::endl;
        for(int j = 0; j < (int) perms_of_5.size(); ++j)
        {
            if(perms_of_5[j] == 0) continue;
            int c = j;
            int i5 = (c%5 == 0 ? 5 : c%5);
            c = (c-i5)/5+1;
            int i4 = (c%5 == 0 ? 5 : c%5);
            c = (c-i4)/5+1;
            int i3 = (c%5 == 0 ? 5 : c%5);
            c = (c-i3)/5+1;
            int i2 = (c%5 == 0 ? 5 : c%5);
            c = (c-i2)/5+1;
            int i1 = (c%5 == 0 ? 5 : c%5);
            int idx = ((((i1-1)*5+i2-1)*5+i3-1)*5+i4-1)*5+i5;
            if(idx != j) std::cerr << "Incorrect idx val (perms of 5)." << std::endl;
            std::cerr << boost::format("(%02d) %d, %d, %d, %d, %d = %d\n") % j % i1 % i2 % i3 % i4 % i5 % perms_of_5[j];
        }

        std::cerr << "Permutations of 6 Count: " << std::endl;
        for(int j = 0; j < (int) perms_of_6.size(); ++j)
        {
            //std::cerr << boost::format("po6[%d] = %d\n") % j % perms_of_6[j];
            if(perms_of_6[j] == 0) continue;
            int c = j;
            int i6 = (c%6 == 0 ? 6 : c%6);
            c = (c-i6)/6+1;
            int i5 = (c%6 == 0 ? 6 : c%6);
            c = (c-i5)/6+1;
            int i4 = (c%6 == 0 ? 6 : c%6);
            c = (c-i4)/6+1;
            int i3 = (c%6 == 0 ? 6 : c%6);
            c = (c-i3)/6+1;
            int i2 = (c%6 == 0 ? 6 : c%6);
            c = (c-i2)/6+1;
            int i1 = (c%6 == 0 ? 6 : c%6);
            int idx = (((((i1-1)*6+i2-1)*6+i3-1)*6+i4-1)*6+i5-1)*6+i6;
            if(idx != j) std::cerr << "Incorrect idx val (perms of 6)." << std::endl;
            std::cerr << boost::format("(%02d) %d, %d, %d, %d, %d, %d = %d\n") % j % i1 % i2 % i3 % i4 % i5 % i6 % perms_of_6[j];
        }

        for(int j = 0; j < (int) time.size(); ++j)
        {
            Real count = tcount.at(j);
            if(time.at(j) > 0) std::cerr << boost::format("Section %d, Average CPU Time = %.2E\n") % j % (time.at(j)/count);
        }

        for(int j = 0; j < (int) time.size(); ++j)
        {
            if(time.at(j) > 0) std::cerr << boost::format("Section %d, Total CPU Time = %f\n") % j % time.at(j);
        }
    }
};


extern Prodstats prodstats;
#ifdef THIS_IS_MAIN
Prodstats prodstats;
#endif

#endif //COLLECT_PRODSTATS

#endif
