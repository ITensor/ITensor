//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_QDIAG_H
#define __ITENSOR_QDIAG_H

#include "itensor/itdata/qdense.h"

namespace itensor {

template<typename T>
class QDiag;

using QDiagReal = QDiag<Real>;
using QDiagCplx = QDiag<Cplx>;

template<typename T>
class QDiag
    {
    static_assert(not std::is_const<T>::value,
                  "Template argument of QDense must be non-const");
    public:
    using value_type = T;
    using storage_type = std::vector<value_type>;
    using iterator = typename storage_type::iterator;
    using const_iterator = typename storage_type::const_iterator;

    //////////////
    std::vector<BlOf> offsets;
        //^ Block index / data offset pairs.
        //  Assumed that block indices are
        //  in increasing order.

    storage_type store;
        //^ *diagonal* tensor data stored contiguously
    //////////////

    QDiag() { }

    QDiag(IQIndexSet const& is, 
          QN const& div_);

    template<typename... SArgs>
    QDiag(std::vector<BlOf> const& off,
          SArgs&&... sargs)
         : offsets(off),
           store(std::forward<SArgs>(sargs)...)
           { }

    explicit operator bool() const { return !store.empty(); }

    value_type*
    data() { return store.data(); }

    value_type const*
    data() const { return store.data(); }

    size_t
    size() const { return store.size(); }

    long
    updateOffsets(IQIndexSet const& is,
                  QN const& div);

    iterator
    begin() { return store.begin(); }

    iterator
    end() { return store.end(); }

    const_iterator
    begin() const { return store.begin(); }

    const_iterator
    end() const { return store.end(); }

    };

template<typename T>
bool constexpr
isReal(QDiag<T> const& t) { return std::is_same<T,Real>::value; }

template<typename T>
bool constexpr
isCplx(QDiag<T> const& t) { return std::is_same<T,Cplx>::value; }

Data inline
realData(QDiagReal & d) { return Data(d.data(),d.size()); }

Datac inline
realData(QDiagReal const& d) { return Datac(d.data(),d.size()); }

Data inline
realData(QDiagCplx & d) { return Data(reinterpret_cast<Real*>(d.data()),2*d.size()); }

Datac inline
realData(QDiagCplx const& d) { return Datac(reinterpret_cast<const Real*>(d.data()),2*d.size()); }

template<typename T, typename Indexable>
T const*
getElt(QDiag<T> const& D,
       IQIndexSet const& is,
       Indexable const& ind)
    {
    auto r = is.r();
#ifdef DEBUG
    if(is.r() != decltype(r)(ind.size())) 
        {
        printfln("is.r() = %d, ind.size() = %d",is.r(),ind.size());
        Error("Wrong number of indices passed to getElt");
        }
#endif
    if(r == 0) return D.data();
    long bind = 0, //block index (total)
         bstr = 1; //block stride so far
    auto last_elt_subind = ind[0];
    for(decltype(r) i = 0; i < r; ++i)
        {
        auto& I = is[i];
        long block_subind = 0,
             elt_subind = ind[i];
        while(elt_subind >= I[block_subind].m()) //elt_subind 0-indexed
            {
            elt_subind -= I[block_subind].m();
            ++block_subind;
            }
        if(i != 0 && elt_subind != last_elt_subind) return nullptr;

        last_elt_subind = elt_subind;
        bind += block_subind*bstr;
        bstr *= I.nindex();
        }
    //Do a binary search (equal_range) to see
    //if there is a block with block index "bind"
    auto boff = offsetOf(D.offsets,bind);
    if(boff != -1)
        {
        auto eoff = last_elt_subind;
#ifdef DEBUG
        if(size_t(boff+eoff) >= D.store.size()) Error("get_elt out of range");
#endif
        return D.data()+boff+eoff;
        }
    return nullptr;
    }

//template<typename Indexable>
//Real*
//getElt(IQTDiag & D,
//       IndexSetT<IQIndex> const& is,
//       Indexable const& ind)
//    {
//    return const_cast<Real*>(getElt(D,is,ind));
//    }

template<typename T>
void
write(std::ostream & s, QDiag<T> const& dat)
    {
    itensor::write(s,dat.offsets);
    itensor::write(s,dat.store);
    }

template<typename T>
void
read(std::istream & s, QDiag<T> & dat)
    {
    itensor::read(s,dat.offsets);
    itensor::read(s,dat.store);
    }
 
template<typename T>
Cplx
doTask(GetElt<IQIndex>& G, QDiag<T> const& D);

template<typename T>
QN
doTask(CalcDiv const& C, QDiag<T> const& d);

template<typename F, typename T>
void
doTask(ApplyIT<F> & A, QDiag<T> & d)
    {
    for(auto& el : d.store) A(el);
    }

template<typename F, typename T>
void
doTask(VisitIT<F> & V, QDiag<T> const& d)
    {
    for(auto& elt : d.store) detail::call<void>(V.f,elt*V.scale_fac);
    }

template<typename F>
void
doTask(GenerateIT<F,Real>& G, QDiagReal & D)
    {
    stdx::generate(D,G.f);
    }

template<typename F>
void
doTask(GenerateIT<F,Real>& G, QDiagCplx const& D, ManageStore & m)
    {
    auto *nD = m.makeNewData<QDiagReal>(D.offsets,D.size());
    stdx::generate(*nD,G.f);
    }

template<typename F>
void
doTask(GenerateIT<F,Cplx>& G, QDiagReal const& D, ManageStore & m)
    {
    auto *nD = m.makeNewData<QDiagCplx>(D.offsets,D.size());
    stdx::generate(*nD,G.f);
    }

template<typename F>
void
doTask(GenerateIT<F,Cplx>& G, QDiagCplx & D)
    {
    stdx::generate(D,G.f);
    }

template<typename T>
Cplx
doTask(SumEls<IQIndex>, QDiag<T> const& d);

template<typename T>
void
doTask(Mult<Real> const& M, QDiag<T>& D);

void
doTask(Mult<Cplx> const& M, QDiag<Cplx> & d);

void
doTask(Mult<Cplx> const& M, QDiag<Real> const& d, ManageStore & m);


void inline
doTask(Conj, QDiagReal const& d) { }

void
doTask(Conj, QDiagCplx & d);

template<typename T>
bool constexpr
doTask(CheckComplex, QDiag<T> const& d) { return isCplx(d); }

template<typename T>
Real
doTask(NormNoScale, QDiag<T> const& d);

template<typename T>
void
doTask(PrintIT<IQIndex>& P, QDiag<T> const& d);

auto inline constexpr
doTask(StorageType const& S, QDiagReal const& d) ->StorageType::Type { return StorageType::QDiagReal; }

auto inline constexpr
doTask(StorageType const& S, QDiagCplx const& d) ->StorageType::Type { return StorageType::QDiagCplx; }

template<typename VA, typename VB>
void
doTask(Contract<IQIndex>& Con,
       QDiag<VA> const& A,
       QDense<VB> const& B,
       ManageStore& m);

template<typename VA, typename VB>
void
doTask(Contract<IQIndex>& Con,
       QDense<VA> const& A,
       QDiag<VB> const& B,
       ManageStore& m);

} //namespace itensor

#endif

