//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_PRINT_H
#define __ITENSOR_PRINT_H

#include "tinyformat.h"

namespace itensor {

using tinyformat::format;

template <typename... Args>
void
printf(const char* fmt_string, const Args&... args)
    {
    tinyformat::printf(fmt_string,args...);
    std::cout.flush();
    }

template <typename... Args>
void
printfln(const char* fmt_string, const Args&... args)
    {
    tinyformat::printf(fmt_string,args...);
    std::cout << std::endl;
    }

void inline
println()
    {
    std::cout << std::endl;
    }

template <typename T>
void
println(const T& arg)
    {
    std::cout << arg << std::endl;
    }

template <typename T, typename... Args>
void
println(const T& arg1, const Args&... args)
    {
    std::cout << arg1;
    println(args...);
    }

template <typename T>
void
print(const T& arg)
    {
    std::cout << arg;
    std::cout.flush();
    }

template <typename T, typename... Args>
void
print(const T& arg1, const Args&... args)
    {
    std::cout << arg1;
    print(args...);
    }



} //namespace itensor

#endif
