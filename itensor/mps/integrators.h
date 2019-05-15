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
#ifndef __ITENSOR_INTEGRATORS_H
#define __ITENSOR_INTEGRATORS_H

#include "itensor/global.h"

namespace itensor {

//
// 4th Order Runge-Kutta.
// Assumes a time-independent Force.
// The vector v should be 1 indexed.
//
template <class Tensor, typename Deriv>
void
rungeKutta4(const Deriv& D, Real tstep, std::vector<Tensor>& v, 
            const Args& args = Args::global())
    {
    int N = int(v.size())-1;
    while(!v.at(N) && N > 1) --N;

    if(N <= 0) Error("Empty vector v (v should be 1-indexed)");

    std::vector<Tensor> k1,k2,k3,k4;

    k1 = D(v);

    //d = v + (tstep/2)*k1
    std::vector<Tensor> d(v);
    for(int j = 1; j <= N; ++j)
        {
        d.at(j) += (tstep/2.)*k1.at(j);
        }
    k2 = D(d);

    //d = v + (tstep/2)*k2
    d = v;
    for(int j = 1; j <= N; ++j)
        {
        d.at(j) += (tstep/2.)*k2.at(j);
        }
    k3 = D(d);

    //d = v + (tstep)*k3
    d = v;
    for(int j = 1; j <= N; ++j)
        {
        d.at(j) += (tstep)*k3.at(j);
        }
    k4 = D(d);


    for(int j = 1; j <= N; ++j)
        {
        v.at(j) += (tstep/6.)*(k1.at(j) + 2*k2.at(j) + 2*k3.at(j) + k4.at(j));
        }
    }

template <class Tensor, typename Deriv>
void
midpointMethod(const Deriv& D, Real tstep, std::vector<Tensor>& v, 
               const Args& args = Args::global())
    {
    int N = int(v.size())-1;
    while(!v.at(N) && N > 1) --N;

    if(N <= 0) Error("Empty vector v (v should be 1-indexed)");

    std::vector<Tensor> k1,k2,k3,k4;

    k1 = D(v);

    std::vector<Tensor> d(v);
    for(int j = 1; j <= N; ++j)
        {
        d.at(j) += (tstep/2.)*k1.at(j);
        }
    k2 = D(d);

    for(int j = 1; j <= N; ++j)
        {
        v.at(j) += tstep*k2.at(j);
        }
    }

} //namespace itensor

#endif
