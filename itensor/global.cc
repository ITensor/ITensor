#include "itensor/global.h"

namespace itensor {

bool&
Global::checkArrows()
    {
    static bool checkArrows_ = true;
    return checkArrows_;
    }
bool&
Global::debug1()
    {
    static bool debug1_ = false;
    return debug1_;
    }
bool&
Global::debug2()
    {
    static bool debug2_ = false;
    return debug2_;
    }
bool&
Global::debug3()
    {
    static bool debug3_ = false;
    return debug3_;
    }
bool&
Global::debug4()
    {
    static bool debug4_ = false;
    return debug4_;
    }
//Global named arguments
Args&
Global::args()
    {
    return Args::global();
    }
void
Global::args(const Args::Name& name, bool bval)
    {
    Args::global().add(name,bval);
    }
void
Global::args(const Args::Name& name, int ival)
    {
    Args::global().add(name,ival);
    }
void
Global::args(const Args::Name& name, Real rval)
    {
    Args::global().add(name,rval);
    }
void
Global::args(const Args::Name& name, const std::string& sval)
    {
    Args::global().add(name,sval);
    }
bool&
Global::printdat()
    {
    static bool printdat_ = false;
    return printdat_;
    }
Real&
Global::printScale()
    {
    static Real printScale_ = 1E-10;
    return printScale_;
    }
bool&
Global::showIDs()
    {
    static bool showIDs_ = true;
    return showIDs_;
    }
Real
Global::random(int seed)
    {
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
void
Global::warnDeprecated(const std::string& message)
    {
    static int depcount = 1;
    if(depcount <= 10)
        {
        println("\n\n",message,"\n");
        ++depcount;
        }
    }
bool&
Global::read32BitIDs()
    {
    static bool read32_ = false;
    return read32_;
    }
} //namespace itensor
