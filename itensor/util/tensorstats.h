//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TENSORSTATS_H
#define __ITENSOR_TENSORSTATS_H

#include <cmath>
#include "itensor/util/stdx.h"
#include "itensor/util/print.h"
#include "itensor/itensor_interface.h"

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
        Ar = A.r();
        for(auto n : range(Ar))
            {
            Adims[n] = A.extent(n);
            Alabs[n] = ai[n];
            }
        Br = B.r();
        for(auto n : range(Br))
            {
            Bdims[n] = B.extent(n);
            Blabs[n] = bi[n];
            }
        Cr = C.r();
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

} //namespace itensor
#endif //COLLECT_TSTATS

#endif
