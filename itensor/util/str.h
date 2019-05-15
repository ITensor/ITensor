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
#ifndef __ITENSOR_STR_H
#define __ITENSOR_STR_H

namespace itensor {

using std::string;

// Take an int and turn it into a string
string  inline
str(int n)
    { 
    return std::to_string(n);
    }

string inline
str(string s, int n)
    {  
    return s+str(n);
    }

string inline
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
