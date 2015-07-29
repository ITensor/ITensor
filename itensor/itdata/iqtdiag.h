//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_IQTDIAG_H
#define __ITENSOR_IQTDIAG_H

#include "itensor/itdata/iqtreal.h"

namespace itensor {

class ManageStore;
class ITCombiner;
class IQTReal;


class IQTDiag
    {
    public:

    //////////////
    std::vector<IQTReal::BlOf> offsets;
        //^ Block index / data offset pairs.
        //  Assumed that block indices are
        //  in increasing order.

    std::vector<Real> data;
        //^ *diagonal* tensor data stored contiguously
    //////////////

    IQTDiag() { }

    IQTDiag(IQIndexSet const& is, 
            QN const& div_);

    explicit operator bool() const { return !data.empty(); }

    long
    updateOffsets(IQIndexSet const& is,
                  QN const& div);

    };

void inline
write(std::ostream & s, IQTDiag const& dat)
    {
    itensor::write(s,dat.offsets);
    itensor::write(s,dat.data);
    }

void inline
read(std::istream & s, IQTDiag & dat)
    {
    itensor::read(s,dat.offsets);
    itensor::read(s,dat.data);
    }

QN
doTask(CalcDiv const& C, IQTDiag const& d);

template<typename F>
void
doTask(ApplyIT<F> & A, IQTDiag & d)
    {
    for(auto& elt : d.data)
        elt = A.f(elt);
    }

template<typename F>
void
doTask(VisitIT<F> & V, IQTDiag const& d)
    {
    for(auto& elt : d.data)
        V.f(elt*V.scale_fac);
    }

template<typename F>
void
doTask(GenerateIT<F,Real> & G, IQTDiag & d)
    {
    std::generate(d.data.begin(),d.data.end(),G.f);
    }

template<typename F>
void
doTask(GenerateIT<F,Cplx> & G, IQTDiag const& cd, ManageStore & mp)
    {
    Error("Complex version of IQTensor generate not yet supported");
    }

void
doTask(MultReal & M, IQTDiag & d);

void
doTask(Contract<IQIndex>& Con,
       IQTReal const& A,
       IQTDiag const& B,
       ManageStore & mp);

void
doTask(Contract<IQIndex>& Con,
       IQTDiag const& A,
       IQTReal const& B,
       ManageStore & mp);

void
doTask(Conj, IQTDiag const& d);

bool inline
doTask(CheckComplex, IQTDiag const& d) { return false; }

Real
doTask(NormNoScale, IQTDiag const& d);

void
doTask(PrintIT<IQIndex> & P, IQTDiag const& d);

void
doTask(Write & W, IQTReal const& d);

} //namespace itensor

#endif

