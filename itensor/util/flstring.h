//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//

#ifndef _ITENSOR_STRING_H_
#define _ITENSOR_STRING_H_
#include <string>
#include "itensor/tensor/types.h"
#include "itensor/util/error.h"

namespace itensor {


//
// Fixed length string.
// Allocates memory on the stack
// for efficiency purposes.
//

template <size_t Maxlen>
class FLString
    {
    public:

    FLString();

    FLString(const char* str);
     
    FLString(const std::string& str);

    size_t 
    length() const;

    size_t 
    size() const { return length(); }

    bool
    empty() const { return (data_[0] == '\0'); }

    operator bool() const { return !empty(); }

    std::string
    toString() const { return std::string(data_); }

    char
    operator[](size_t n) const;

    const char*
    data() const { return &(data_[0]); }

    private:
    char data_[Maxlen+1];
    };


template <size_t Maxlen>
std::ostream& 
operator<<(std::ostream &s, const FLString<Maxlen>& a);

template <size_t Maxlen>
bool
operator==(const FLString<Maxlen>& s1, const FLString<Maxlen>& s2);

template <size_t Maxlen>
bool
operator!=(const FLString<Maxlen>& s1, const FLString<Maxlen>& s2) { return !operator==(s1,s2); }

template <size_t Maxlen>
bool
operator==(const FLString<Maxlen>& s1, const std::string& s2);

template <size_t Maxlen>
bool
operator==(const std::string& s1, const FLString<Maxlen>& s2) { return operator==(s2,s1); }

template <size_t Maxlen>
bool
operator!=(const FLString<Maxlen>& s1, const std::string& s2) { return !operator==(s1,s2); }

template <size_t Maxlen>
bool
operator!=(const std::string& s1, const FLString<Maxlen>& s2) { return !operator==(s2,s1); }

template <size_t Maxlen>
bool
operator==(const FLString<Maxlen>& s1, const char* s2);

template <size_t Maxlen>
bool
operator==(const char* s1, const FLString<Maxlen>& s2) { return operator==(s2,s1); }

template <size_t Maxlen>
bool
operator!=(const FLString<Maxlen>& s1, const char* s2) { return !operator==(s1,s2); }

template <size_t Maxlen>
bool
operator!=(const char* s1, const FLString<Maxlen>& s2) { return !operator==(s2,s1); }

//
// Implementations
//

template <size_t Maxlen>
FLString<Maxlen>::
FLString()
    { 
    for(size_t n = 0; n <= Maxlen; ++n)
        {
        data_[n] = '\0';
        }
    }

template <size_t Maxlen>
FLString<Maxlen>::
FLString(const char* str)
    {
    size_t n = 0;
    for(;str[n] != '\0' && n < Maxlen; ++n)
        {
        data_[n] = str[n];
        }
#ifdef DEBUG
    if(n == Maxlen && str[n] != '\0')
        {
        std::cout << "Warning: input to FLString exceeds " << Maxlen << " characters, truncating." << std::endl;
        }
#endif
    for(;n <= Maxlen; ++n) data_[n] = '\0';
    }

template <size_t Maxlen>
FLString<Maxlen>::
FLString(const std::string& str)
    {
    size_t ll = std::min(str.size(),Maxlen);
    size_t n = 0;
    for(; n < ll; ++n)
        {
        data_[n] = str[n];
        }
#ifdef DEBUG
    if(str.size() > Maxlen)
        {
        std::cout << "Warning: input to FLString exceeds " << Maxlen << " characters, truncating." << std::endl;
        }
#endif
    for(;n <= Maxlen; ++n) data_[n] = '\0';
    }

template <size_t Maxlen>
char FLString<Maxlen>::
operator[](size_t n) const
    {
#ifdef DEBUG
    if(n >= Maxlen) Error("FLString index out of range");
#endif
    return data_[n];
    }

template <size_t Maxlen>
size_t FLString<Maxlen>::
length() const 
    {
    size_t l = 0;
    while(data_[l] != '\0' && l < Maxlen)
        {
        ++l;
        }
    return l;
    }

template <size_t Maxlen>
std::ostream&
operator<<(std::ostream& s, const FLString<Maxlen>& str)
    {
    size_t n = 0;
    const char* d = str.data();
    for(; d[n] != '\0' && n < Maxlen; ++n)
        {
        s << d[n];
        }
    return s;
    }

template <size_t Maxlen>
bool
operator==(const FLString<Maxlen>& s1, const FLString<Maxlen>& s2)
    {
    for(size_t n = 0; n < Maxlen; ++n)
        {
        if(s1[n] != s2[n]) return false;
        }
    return true;
    }

template <size_t Maxlen>
bool
operator==(const FLString<Maxlen>& s1, 
           const std::string& s2)
    {
    const size_t ll = std::min(s2.size(),Maxlen);
    for(size_t n = 0; n < ll; ++n)
        {
        if(s1[n] != s2[n]) return false;
        }
    return true;
    }

template <size_t Maxlen>
bool
operator==(const FLString<Maxlen>& s1, const char* s2)
    {
    for(size_t n = 0; n < Maxlen; ++n)
        {
        if(s1[n] != s2[n]) return false;
        if(s2[n] == '\0') break;
        }
    return true;
    }

//template <size_t Maxlen>
//bool
//operator==(const FLString<Maxlen>& s1, const char s2[7])
//    {
//    size_t N = 7;
//    const size_t ll = std::min(N,Maxlen);
//    for(size_t n = 0; n < ll; ++n)
//        {
//        if(s1[n] != s2[n]) return false;
//        }
//    return true;
//    }
//
//template <size_t Maxlen, size_t N>
//bool
//operator==(const char s1[7], const FLString<Maxlen>& s2)
//    {
//    size_t N = 7;
//    const size_t ll = std::min(N,Maxlen);
//    for(size_t n = 0; n < ll; ++n)
//        {
//        if(s1[n] != s2[n]) return false;
//        }
//    return true;
//    }

} //namespace itensor

#endif
