//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_STR_H
#define __ITENSOR_STR_H

namespace itensor {

using std::string;

// Take an int and turn it into a string
string 
str(int n)
    { 
    return std::to_string(n);
    }

string  
str(string s, int n)
    {  
    return s+str(n);
    }

string  
str(int n, string s)
    {  
    return str(n)+s;
    }

template<typename... VarArgs>
string
str(string s1, int n1, string s2, VarArgs&&... vargs)
    {  
    return str(s1,n1)+str(s2,std::forward<VarArgs>(vargs)...);
    }

template<typename... VarArgs>
string
str(int n1, string s1, int n2, VarArgs&&... vargs)
    {  
    return str(n1,s1)+str(n2,std::forward<VarArgs>(vargs)...);
    }

} //namespace itensor
#endif
