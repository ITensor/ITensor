//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_PRINT_H
#define __ITENSOR_PRINT_H

#include "tinyformat.h"

namespace itensor {

using tinyformat::printf;
using tinyformat::format;

#ifdef USE_CPP11

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

#else //if not USE_CPP11

template <typename T1>
void
printfln(const char* fmt_string, 
         const T1& arg1)
    {
    tinyformat::printf(fmt_string,arg1);
    std::cout << std::endl;
    }

template <typename T1,typename T2>
void
printfln(const char* fmt_string, 
         const T1& arg1,const T2& arg2)
    {
    tinyformat::printf(fmt_string,arg1,arg2);
    std::cout << std::endl;
    }

template <typename T1,typename T2,
          typename T3>
void
printfln(const char* fmt_string, 
         const T1& arg1,const T2& arg2,
         const T3& arg3)
    {
    tinyformat::printf(fmt_string,arg1,arg2,arg3);
    std::cout << std::endl;
    }

template <typename T1,typename T2,
          typename T3,typename T4>
void
printfln(const char* fmt_string, 
         const T1& arg1,const T2& arg2,
         const T3& arg3,const T4& arg4)
    {
    tinyformat::printf(fmt_string,arg1,arg2,arg3,arg4);
    std::cout << std::endl;
    }

template <typename T1,typename T2,
          typename T3,typename T4,
          typename T5>
void
printfln(const char* fmt_string, 
         const T1& arg1,const T2& arg2,
         const T3& arg3,const T4& arg4,
         const T5& arg5)
    {
    tinyformat::printf(fmt_string,arg1,arg2,arg3,arg4,
                                  arg5);
    std::cout << std::endl;
    }

template <typename T1,typename T2,
          typename T3,typename T4,
          typename T5,typename T6>
void
printfln(const char* fmt_string, 
         const T1& arg1,const T2& arg2,
         const T3& arg3,const T4& arg4,
         const T5& arg5,const T6& arg6)
    {
    tinyformat::printf(fmt_string,arg1,arg2,arg3,arg4,
                                  arg5,arg6);
    std::cout << std::endl;
    }

template <typename T1,typename T2,
          typename T3,typename T4,
          typename T5,typename T6,
          typename T7>
void
printfln(const char* fmt_string, 
         const T1& arg1,const T2& arg2,
         const T3& arg3,const T4& arg4,
         const T5& arg5,const T6& arg6,
         const T7& arg7)
    {
    tinyformat::printf(fmt_string,arg1,arg2,arg3,arg4,
                                  arg5,arg6,arg7);
    std::cout << std::endl;
    }

template <typename T1,typename T2,
          typename T3,typename T4,
          typename T5,typename T6,
          typename T7,typename T8>
void
printfln(const char* fmt_string, 
         const T1& arg1,const T2& arg2,
         const T3& arg3,const T4& arg4,
         const T5& arg5,const T6& arg6,
         const T7& arg7,const T8& arg8)
    {
    tinyformat::printf(fmt_string,arg1,arg2,arg3,arg4,
                                  arg5,arg6,arg7,arg8);
    std::cout << std::endl;
    }

void inline
println()
    {
    std::cout << std::endl;
    }

template <typename T1>
void
println(const T1& arg1)
    {
    std::cout << arg1 << std::endl;
    }

template <typename T1,typename T2>
void
println(const T1& arg1,const T2& arg2)
    {
    std::cout << arg1 << arg2 << std::endl;
    }

template <typename T1,typename T2,
          typename T3>
void
println(const T1& arg1,const T2& arg2,
        const T3& arg3)
    {
    std::cout << arg1;
    println(arg2,arg3);
    }

template <typename T1,typename T2,
          typename T3,typename T4>
void
println(const T1& arg1,const T2& arg2,
        const T3& arg3,const T4& arg4)
    {
    std::cout << arg1;
    println(arg2,arg3,arg4);
    }

template <typename T1,typename T2,
          typename T3,typename T4,
          typename T5>
void
println(const T1& arg1,const T2& arg2,
        const T3& arg3,const T4& arg4,
        const T5& arg5)
    {
    std::cout << arg1;
    println(arg2,arg3,arg4,arg5);
    }

template <typename T1,typename T2,
          typename T3,typename T4,
          typename T5,typename T6>
void
println(const T1& arg1,const T2& arg2,
        const T3& arg3,const T4& arg4,
        const T5& arg5,const T6& arg6)
    {
    std::cout << arg1;
    println(arg2,arg3,arg4,arg5,arg6);
    }

template <typename T1,typename T2,
          typename T3,typename T4,
          typename T5,typename T6,
          typename T7>
void
println(const T1& arg1,const T2& arg2,
        const T3& arg3,const T4& arg4,
        const T5& arg5,const T6& arg6,
        const T7& arg7)
    {
    std::cout << arg1;
    println(arg2,arg3,arg4,arg5,arg6,arg7);
    }

template <typename T1,typename T2,
          typename T3,typename T4,
          typename T5,typename T6,
          typename T7,typename T8>
void
println(const T1& arg1,const T2& arg2,
        const T3& arg3,const T4& arg4,
        const T5& arg5,const T6& arg6,
        const T7& arg7,const T8& arg8)
    {
    std::cout << arg1;
    println(arg2,arg3,arg4,arg5,arg6,arg7,arg8);
    }

template <typename T1,typename T2,
          typename T3,typename T4,
          typename T5,typename T6,
          typename T7,typename T8,
          typename T9>
void
println(const T1& arg1,const T2& arg2,
        const T3& arg3,const T4& arg4,
        const T5& arg5,const T6& arg6,
        const T7& arg7,const T8& arg8,
        const T9& arg9)
    {
    std::cout << arg1;
    println(arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9);
    }

template <typename T1,typename T2,
          typename T3,typename T4,
          typename T5,typename T6,
          typename T7,typename T8,
          typename T9,typename T10>
void
println(const T1& arg1,const T2& arg2,
        const T3& arg3,const T4& arg4,
        const T5& arg5,const T6& arg6,
        const T7& arg7,const T8& arg8,
        const T9& arg9,const T10& arg10)
    {
    std::cout << arg1;
    println(arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10);
    }

template <typename T1,typename T2,
          typename T3,typename T4,
          typename T5,typename T6,
          typename T7,typename T8,
          typename T9,typename T10,
          typename T11>
void
println(const T1& arg1,const T2& arg2,
        const T3& arg3,const T4& arg4,
        const T5& arg5,const T6& arg6,
        const T7& arg7,const T8& arg8,
        const T9& arg9,const T10& arg10,
        const T11& arg11)
    {
    std::cout << arg1;
    println(arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11);
    }

template <typename T1,typename T2,
          typename T3,typename T4,
          typename T5,typename T6,
          typename T7,typename T8,
          typename T9,typename T10,
          typename T11,typename T12>
void
println(const T1& arg1,const T2& arg2,
        const T3& arg3,const T4& arg4,
        const T5& arg5,const T6& arg6,
        const T7& arg7,const T8& arg8,
        const T9& arg9,const T10& arg10,
        const T11& arg11,const T12& arg12)
    {
    std::cout << arg1;
    println(arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12);
    }

template <typename T1>
void
print(const T1& arg1)
    {
    std::cout << arg1;
    std::cout.flush();
    }

template <typename T1,typename T2>
void
print(const T1& arg1,const T2& arg2)
    {
    std::cout << arg1 << arg2;
    std::cout.flush();
    }

template <typename T1,typename T2,
          typename T3>
void
print(const T1& arg1,const T2& arg2,
      const T3& arg3)
    {
    std::cout << arg1;
    print(arg2,arg3);
    }

template <typename T1,typename T2,
          typename T3,typename T4>
void
print(const T1& arg1,const T2& arg2,
      const T3& arg3,const T4& arg4)
    {
    std::cout << arg1;
    print(arg2,arg3,arg4);
    }

template <typename T1,typename T2,
          typename T3,typename T4,
          typename T5>
void
print(const T1& arg1,const T2& arg2,
      const T3& arg3,const T4& arg4,
      const T5& arg5)
    {
    std::cout << arg1;
    print(arg2,arg3,arg4,arg5);
    }

template <typename T1,typename T2,
          typename T3,typename T4,
          typename T5,typename T6>
void
print(const T1& arg1,const T2& arg2,
      const T3& arg3,const T4& arg4,
      const T5& arg5,const T6& arg6)
    {
    std::cout << arg1;
    print(arg2,arg3,arg4,arg5,arg6);
    }

template <typename T1,typename T2,
          typename T3,typename T4,
          typename T5,typename T6,
          typename T7>
void
print(const T1& arg1,const T2& arg2,
      const T3& arg3,const T4& arg4,
      const T5& arg5,const T6& arg6,
      const T7& arg7)
    {
    std::cout << arg1;
    print(arg2,arg3,arg4,arg5,arg6,arg7);
    }

template <typename T1,typename T2,
          typename T3,typename T4,
          typename T5,typename T6,
          typename T7,typename T8>
void
print(const T1& arg1,const T2& arg2,
      const T3& arg3,const T4& arg4,
      const T5& arg5,const T6& arg6,
      const T7& arg7,const T8& arg8)
    {
    std::cout << arg1;
    print(arg2,arg3,arg4,arg5,arg6,arg7,arg8);
    }

template <typename T1,typename T2,
          typename T3,typename T4,
          typename T5,typename T6,
          typename T7,typename T8,
          typename T9>
void
print(const T1& arg1,const T2& arg2,
      const T3& arg3,const T4& arg4,
      const T5& arg5,const T6& arg6,
      const T7& arg7,const T8& arg8,
      const T9& arg9)
    {
    std::cout << arg1;
    print(arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9);
    }

template <typename T1,typename T2,
          typename T3,typename T4,
          typename T5,typename T6,
          typename T7,typename T8,
          typename T9,typename T10>
void
print(const T1& arg1,const T2& arg2,
      const T3& arg3,const T4& arg4,
      const T5& arg5,const T6& arg6,
      const T7& arg7,const T8& arg8,
      const T9& arg9,const T10& arg10)
    {
    std::cout << arg1;
    print(arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10);
    }

template <typename T1,typename T2,
          typename T3,typename T4,
          typename T5,typename T6,
          typename T7,typename T8,
          typename T9,typename T10,
          typename T11>
void
print(const T1& arg1,const T2& arg2,
      const T3& arg3,const T4& arg4,
      const T5& arg5,const T6& arg6,
      const T7& arg7,const T8& arg8,
      const T9& arg9,const T10& arg10,
      const T11& arg11)
    {
    std::cout << arg1;
    print(arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11);
    }

template <typename T1,typename T2,
          typename T3,typename T4,
          typename T5,typename T6,
          typename T7,typename T8,
          typename T9,typename T10,
          typename T11,typename T12>
void
print(const T1& arg1,const T2& arg2,
      const T3& arg3,const T4& arg4,
      const T5& arg5,const T6& arg6,
      const T7& arg7,const T8& arg8,
      const T9& arg9,const T10& arg10,
      const T11& arg11,const T12& arg12)
    {
    std::cout << arg1;
    print(arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12);
    }

#endif //USE_CPP11

}; //namespace itensor

#endif
