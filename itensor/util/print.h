//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_PRINT_H
#define __ITENSOR_PRINT_H

#include "itensor/util/tinyformat.h"

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

#define Print(X)  itensor::PrintNice(#X,X)

template<typename T>
void
PrintNice(const char* tok,
          T const& X)
    {
    auto pre = format("%s = ",tok);
    auto str = format("%s",X);

    //Put a newline after '=' if
    //output is large or output contains
    //newline character
    bool put_newline = false;
    if(pre.size() + str.size() > 60)
        {
        put_newline = true;
        }
    else
        {
        for(auto c : str)
            if(c == '\n')
                {
                put_newline = true;
                break;
                }
        }
    if(put_newline) println(pre);
    else            print(pre);
    println(str);
    }


} //namespace itensor

#endif
