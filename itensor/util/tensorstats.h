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
#ifndef __ITENSOR_TENSORSTATS_H
#define __ITENSOR_TENSORSTATS_H

#include <cmath>
#include "itensor/util/stdx.h"
#include "itensor/util/print.h"
#include "itensor/itensor.h"

//#define COLLECT_TSTATS

#ifdef COLLECT_TSTATS
namespace itensor {

struct TStats
    {
    using iarray = std::array<int,8>;

    int Ar = 0;
    int Br = 0;
    int Cr = 0;
    iarray Adims = iarray{};
    iarray Alabs = iarray{};
    iarray Bdims = iarray{};
    iarray Blabs = iarray{};
    iarray Cdims = iarray{};
    iarray Clabs = iarray{};

    TStats() { }

    template<typename RangeT, typename VA, typename VB>
    TStats(TenRefc<RangeT,VA> A, Labels const& ai, 
           TenRefc<RangeT,VB> B, Labels const& bi, 
           TenRef<RangeT,common_type<VA,VB>>  C, 
           Labels const& ci)
        {
        Ar = A.order();
        for(auto n : range(Ar))
            {
            Adims[n] = A.extent(n);
            Alabs[n] = ai[n];
            }
        Br = B.order();
        for(auto n : range(Br))
            {
            Bdims[n] = B.extent(n);
            Blabs[n] = bi[n];
            }
        Cr = C.order();
        for(auto n : range(Cr))
            {
            Cdims[n] = C.extent(n);
            Clabs[n] = ci[n];
            }
        }

    void
    write(std::ostream & s) const;
    void
    read(std::istream & s);

    };

void inline TStats::
write(std::ostream & s) const
    {
    itensor::write(s,Ar);
    itensor::write(s,Br);
    itensor::write(s,Cr);
    itensor::write(s,Adims);
    itensor::write(s,Bdims);
    itensor::write(s,Cdims);
    itensor::write(s,Alabs);
    itensor::write(s,Blabs);
    itensor::write(s,Clabs);
    }

void inline TStats::
read(std::istream & s)
    {
    itensor::read(s,Ar);
    itensor::read(s,Br);
    itensor::read(s,Cr);
    itensor::read(s,Adims);
    itensor::read(s,Bdims);
    itensor::read(s,Cdims);
    itensor::read(s,Alabs);
    itensor::read(s,Blabs);
    itensor::read(s,Clabs);
    }

inline std::vector<TStats>&
global_tstats()
    {
    static std::vector<TStats> gts;
    return gts;
    }

template<typename... VArgs>
void
tstats(VArgs&&... vargs)
    {
    global_tstats().emplace_back(std::forward<VArgs&&>(vargs)...);
    }

inline std::ostream&
operator<<(std::ostream& s, TStats const& T)
    {
    s << format("A (%d) [ ",T.Ar);
    for(auto n : range(T.Ar)) s << T.Adims[n] << " ";
    s << "] { ";
    for(auto n : range(T.Ar)) s << T.Alabs[n] << " ";
    s << "}\n";

    s << format("B (%d) [ ",T.Br);
    for(auto n : range(T.Br)) s << T.Bdims[n] << " ";
    s << "] { ";
    for(auto n : range(T.Br)) s << T.Blabs[n] << " ";
    s << "}\n";

    s << format("C (%d) [ ",T.Cr);
    for(auto n : range(T.Cr)) s << T.Cdims[n] << " ";
    s << "] { ";
    for(auto n : range(T.Cr)) s << T.Clabs[n] << " ";
    s << "}\n";

    return s;
    }

void inline
doHistogram1(Real bucket_size = 0.5)
    {
    Real max_value = 10.;
	int nbucket = (int)std::ceil(max_value/bucket_size);
	auto histogram = std::vector<int>(nbucket);
    //printfln("nbucket = %d",nbucket);

    for(auto& t : global_tstats())
        {
        Real logAdim = 0.;
        for(auto n : range(t.Ar)) logAdim += std::log10(t.Adims[n]);
		int bucket = (int)floor(logAdim/bucket_size);
		//histogram[bucket] += 1;
		histogram.at(bucket) += 1;
        }

    std::ofstream f("hist.dat");
    for(auto b : range(histogram))
        {
        printfln(f,"%.12f %.12f",b*bucket_size,histogram[b]);
        }
    f.close();
    }

void inline
makeBenchmark(std::string name)
    {
    auto charOf = [](int n) -> char
        {
        char i = 'a'-1;
        return i+abs(n);
        };
    auto printLabels = [charOf](std::ofstream & f,
                                TStats::iarray labels,
                                int r)

        {
        for(auto l : range(r)) 
            {
            if(l != 0) print(f,",");
            print(f,charOf(labels[l]));
            }
        };
    auto printSizes = [charOf](std::ofstream & f,
                               int j,
                               TStats::iarray labels,
                               TStats::iarray dims,
                               int r)

        {
        for(auto i : range(r))
            {
            if(abs(labels[i]) == j)
                {
                print(f," ",charOf(j),":",dims[i]);
                return true;
                }
            }
        return false;
        };
    std::ofstream f(name+".dat");
    for(auto& t : global_tstats())
        {
        print(f,"C[");
        printLabels(f,t.Clabs,t.Cr);
        print(f,"] = A[");
        printLabels(f,t.Alabs,t.Ar);
        print(f,"] * B[");
        printLabels(f,t.Blabs,t.Br);
        print(f,"] &");

        auto Nind = (t.Ar+t.Br+t.Cr)/2;
        for(int j : range1(Nind))
            {
            bool found = printSizes(f,j,t.Clabs,t.Cdims,t.Cr);
            if(!found) found = printSizes(f,j,t.Alabs,t.Adims,t.Ar);
            if(!found) printSizes(f,j,t.Blabs,t.Bdims,t.Br);
            }
        print(f,";\n");
        }
    f.close();
    }

} //namespace itensor
#endif //COLLECT_TSTATS

#endif
