//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
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
#include "itensor/util/iterate.h"
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
const int MAX_DIM = 5000;
const int MAX_TAGS = 4;

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
