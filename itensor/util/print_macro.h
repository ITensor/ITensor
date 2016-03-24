//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_PRINT_MACRO_H
#define __ITENSOR_PRINT_MACRO_H

#include "itensor/util/print.h"

#ifndef Print
#define Print(X)  itensor::PrintNice(#X,X)
#endif

namespace itensor {

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
