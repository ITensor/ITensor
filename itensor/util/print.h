//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_PRINT_H
#define __ITENSOR_PRINT_H

#include "tinyformat.h"

namespace itensor {

using tinyformat::format;

template <typename... VArgs>
void
printf(const char* fmt_string, VArgs&&... vargs)
    {
    tinyformat::printf(fmt_string,std::forward<VArgs>(vargs)...);
    std::cout.flush();
    }

template <typename... VArgs>
void
printfln(const char* fmt_string, VArgs&&... vargs)
    {
    tinyformat::printf(fmt_string,std::forward<VArgs>(vargs)...);
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

template <typename T, typename... VArgs>
void
println(const T& arg1, VArgs&&... vargs)
    {
    std::cout << arg1;
    println(std::forward<VArgs>(vargs)...);
    }

template <typename T>
void
print(const T& arg)
    {
    std::cout << arg;
    std::cout.flush();
    }

template <typename T, typename... VArgs>
void
print(const T& arg1, VArgs&&... vargs)
    {
    std::cout << arg1;
    print(std::forward<VArgs>(vargs)...);
    }



} //namespace itensor

#endif
