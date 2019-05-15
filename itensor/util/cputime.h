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
#ifndef _CPUTIME_h
#define _CPUTIME_h

#include <iostream>

namespace itensor {

double cpu_mytime();
double cpu_mywall();

struct cpu_time
    {
    double time = 0; // in seconds
    double wall = 0;

    cpu_time() { mark(); }

    void 
    mark() { time = cpu_mytime(); wall = cpu_mywall(); }

    cpu_time 
    sincemark() const;
    };

std::ostream& 
operator<<(std::ostream & s, const cpu_time& t);

const double&
firstwall();

std::string 
showtime(double time);

std::ostream& 
operator<<(std::ostream & s, const cpu_time& t);

} //namespace itensor

#endif
