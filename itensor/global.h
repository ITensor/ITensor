//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_GLOBAL_H
#define __ITENSOR_GLOBAL_H

#include <cstdlib>
#include <string>
#include <cstring>
#include <memory>
#include <random>
#include <unistd.h>
#include <memory>
#include "itensor/util/range.h"
#include "itensor/util/error.h"
#include "itensor/util/args.h"
#include "itensor/real.h"
#include "itensor/util/timers.h"
#include "itensor/detail/algs.h"

namespace itensor {

enum Direction { Fromright, Fromleft, BothDir, NoDir };

Direction inline
toDirection(int i)
    {
    if(static_cast<int>(Fromright) == i) return Fromright;
    if(static_cast<int>(Fromleft) == i) return Fromleft;
    if(static_cast<int>(BothDir) == i) return BothDir;
    return NoDir;
    }

Direction inline
getDirection(Args const& args,
             Args::Name const& name)
    {
    return toDirection(args.getInt(name));
    }

Direction inline
getDirection(Args const& args,
             Args::Name const& name,
             Direction default_val)
    {
    return toDirection(args.getInt(name,default_val));
    }

const Real MIN_CUT = 1E-15;
const int MAX_M = 5000;

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

#ifdef USESCALE
#define IF_USESCALE(X) X
#else
#define IF_USESCALE(X) 
#endif

enum Printdat { ShowData, HideData };


class Global
    {
    public:
    static bool& checkArrows();
    static bool& debug1();
    static bool& debug2();
    static bool& debug3();
    static bool& debug4();

    //Global named arguments
    static Args& args();

    void static args(const Args::Name& name, bool bval);
    void static args(const Args::Name& name, int ival);
    void static args(const Args::Name& name, Real rval);
    void static args(const Args::Name& name, const std::string& sval);
    static bool& printdat();
    static Real& printScale();
    static bool& showIDs();
    static Real random(int seed = 0);
    void static warnDeprecated(const std::string& message);
    static bool& read32BitIDs();
    };

#define PrintData(X) PrintEither(true,#X,X)
#define PrintDat(X)  PrintEither(true,#X,X)

template<typename T>
void
PrintEither(bool pdat,
            const char* tok,
            T const& X)
    {
    auto savep = Global::printdat();
    Global::printdat() = pdat;
    PrintNice(tok,X);
    Global::printdat() = savep;
    }

void inline
seedRNG(int seed)
    {
    Global::random(seed);
    detail::seed_quickran(seed);
    }

} //namespace itensor

#endif
