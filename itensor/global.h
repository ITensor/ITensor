//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_GLOBAL_H
#define __ITENSOR_GLOBAL_H

#include <cstdlib>
#include <string>
#include <cstring>
#include <random>
#include <unistd.h>
#include "itensor/util/count.h"
#include "itensor/util/error.h"
#include "itensor/util/args.h"
#include "itensor/real.h"
#include "itensor/util/timers.h"

namespace itensor {

enum Direction { Fromright, Fromleft, Both, None };

static const Real MIN_CUT = 1E-15;
static const int MAX_M = 5000;

// The PAUSE macro is useful for debugging. 
// Prints the current line number and pauses
// execution until the enter key is pressed.
#ifndef PAUSE
#define PAUSE { std::cout << "(Paused, Line " << __LINE__ << ")"; std::cin.get(); }
#endif

#ifndef EXIT
#define EXIT exit(0);
#endif

#ifndef DEBUG
#ifndef NDEBUG
#define NDEBUG //turn off asserts
#endif
#endif

#ifdef DEBUG
#define DO_IF_DEBUG(X) X
#else
#define DO_IF_DEBUG(X)
#endif

#ifdef DEBUG
#define GET(container,j) (container.at(j))
#else
#define GET(container,j) (container[j])
#endif	

enum Printdat { ShowData, HideData };


class Global
    {
    public:
    static bool& 
    checkArrows()
        {
        static bool checkArrows_ = true;
        return checkArrows_;
        }
    static bool& 
    debug1()
        {
        static bool debug1_ = false;
        return debug1_;
        }
    static bool& 
    debug2()
        {
        static bool debug2_ = false;
        return debug2_;
        }
    static bool& 
    debug3()
        {
        static bool debug3_ = false;
        return debug3_;
        }
    static bool& 
    debug4()
        {
        static bool debug4_ = false;
        return debug4_;
        }
    //Global named arguments
    static Args&
    args()
        {
        return Args::global();
        }
    void static
    args(const Args::Name& name, bool bval)
        {
        Args::global().add(name,bval);
        }
    void static
    args(const Args::Name& name, int ival)
        {
        Args::global().add(name,ival);
        }
    void static
    args(const Args::Name& name, Real rval)
        {
        Args::global().add(name,rval);
        }
    void static
    args(const Args::Name& name, const std::string& sval)
        {
        Args::global().add(name,sval);
        }
    static bool& 
    printdat()
        {
        static bool printdat_ = false;
        return printdat_;
        }
    static Real& 
    printScale()
        {
        static Real printScale_ = 1E-10;
        return printScale_;
        }
    static Real
    random()//int seed = 0)
        {
        int seed = 0;
        using Generator = std::mt19937;
        using Distribution = std::uniform_real_distribution<Real>;

        static Generator rng(std::time(NULL)+getpid());
        static Distribution dist(0,1);

        if(seed != 0)  //reseed rng
            {
            rng = Generator(seed);
            }

        return dist(rng);
        }
    void static
    warnDeprecated(const std::string& message)
        {
        static int depcount = 1;
        if(depcount <= 10)
            {
            println("\n\n",message,"\n");
            ++depcount;
            }
        }
    };

template<typename T>
void
PrintEither(bool pdat,
            const char* tok,
            T const& X)
    {
    bool savep = Global::printdat();
    Global::printdat() = pdat;
    auto pre = format("%s = ",tok);
    auto str = format("%s",X);
    std::cout << pre;
    bool has_newline = false;
    for(auto c : str)
        if(c == '\n')
            {
            has_newline = true;
            break;
            }
    if(has_newline || (pre.size() + str.size() > 60)) std::cout << "\n";
    std::cout << str << std::endl;
    Global::printdat() = savep;
    }
#define Print(X)     PrintEither(false,#X,X)
#define PrintDat(X)  PrintEither(true,#X,X)
#define PrintData(X) PrintEither(true,#X,X)

} //namespace itensor

#endif
