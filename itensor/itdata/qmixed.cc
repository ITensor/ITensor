//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/itdata/qmixed.h"
#include "itensor/detail/gcounter.h"

using std::vector;
using std::move;

namespace itensor {

const char*
typeNameOf(QMixed<Real> const& d) { return "QMixed<Real>"; }
const char*
typeNameOf(QMixed<Cplx> const& d) { return "QMixed<Cplx>"; }

template<typename V>
Cplx 
doTask(GetElt<IQIndex> const& g, QMixed<V> const& d)
    {
    return d[offset(g.is,g.inds)];
    }
template Cplx 
doTask(GetElt<IQIndex> const& g, QMixed<Real> const& d);
template Cplx 
doTask(GetElt<IQIndex> const& g, QMixed<Cplx> const& d);

namespace detail {
    template<typename E, typename T>
    struct SetEltHelper
        {
        void static
        set(SetElt<E,IQIndex> const& S, QMixed<T> const& D, ManageStore& m)
            {
            auto& Dnc = *m.modifyData(D);
            Dnc[offset(S.is,S.inds)] = S.elt;
            }
        };
    template<>
    struct SetEltHelper<Cplx,Real>
        {
        void static
        set(SetElt<Cplx,IQIndex> const& S, QMixed<Real> const& D, ManageStore & m)
            {
            auto& nd = *m.makeNewData<QMixed<Cplx>>(D.begin(),D.end());
            nd[offset(S.is,S.inds)] = S.elt;
            }
        };
} //namespace detail

template<typename E, typename T>
void
doTask(SetElt<E,IQIndex> const& S, QMixed<T> const& d, ManageStore & m)
    {
    detail::SetEltHelper<E,T>::set(S,d,m);
    }
template void
doTask(SetElt<Real,IQIndex> const& S, QMixed<Real> const& d, ManageStore & m);
template void
doTask(SetElt<Real,IQIndex> const& S, QMixed<Cplx> const& d, ManageStore & m);
template void
doTask(SetElt<Cplx,IQIndex> const& S, QMixed<Real> const& d, ManageStore & m);
template void
doTask(SetElt<Cplx,IQIndex> const& S, QMixed<Cplx> const& d, ManageStore & m);

template<typename T>
void
doTask(PrintIT<IQIndex>& P, QMixed<T> const& D)
    {
    auto name = std::is_same<T,Real>::value ? "QMixed Real"
                                            : "QMixed Cplx";
    P.printInfo(D,name);
     
    auto r = rank(P.is);
    if(r == 0) 
        {
        P.s << "  ";
        P.s << formatVal(P.scalefac*D.store.front()) << "\n";
        return;
        }

    if(!P.print_data) return;

    auto gc = detail::GCounter(r);
    for(auto i : range(r))
        gc.setRange(i,0,P.is.extent(i)-1);

    for(; gc.notDone(); ++gc)
        {
        auto val = P.scalefac*D[offset(P.is,gc.i)];
        if(std::norm(val) >= Global::printScale())
            {
            P.s << "(";
            for(auto ii : range1(gc.i.mini(),gc.i.maxi()))
                {
                P.s << (1+gc[ii]);
                if(ii < gc.i.maxi()) P.s << ",";
                }
            P.s << ") ";

            P.s << formatVal(val) << "\n";
            }
        }
    }
template void
doTask(PrintIT<IQIndex>& P, QMixed<Real> const& D);
template void
doTask(PrintIT<IQIndex>& P, QMixed<Cplx> const& D);

} //namespace itensor

