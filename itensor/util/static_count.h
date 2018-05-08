#ifndef __ITENSOR_STATIC_COUNT_H__
#define __ITENSOR_STATIC_COUNT_H__

#include "itensor/util/print.h"

namespace itensor {

//
// Usage:
//
//    static auto sc = StaticCount("Incremented %d times");
//    ++sc;
//    //will automatically print "Incremented ## times"
//    //when program ends
//

struct StaticCount
    {
    long count = 0;
    const char* fstring = "";

    StaticCount(const char* fstring_) : fstring(fstring_) { }

    ~StaticCount()
        {
        printfln(fstring,count);
        }

    void
    operator++() { ++count; }
    };

} //namespace itensor

#endif
