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
                  "Template argument of QDiag must be non-const");
    public:
    using value_type = T;
    using storage_type = vector_no_init<value_type>;
    using iterator = typename storage_type::iterator;
    using const_iterator = typename storage_type::const_iterator;

    //////////////
    storage_type store;
        //^ *diagonal* tensor elements stored contiguously

    T val = 0;
    size_t length = 0ul;
    //////////////

    QDiag() { }

    QDiag(IndexSet const& is);

    //Special "allSame" mode where non-zero
    //elements assumed to have the same value "val"
    QDiag(IndexSet const& is, 
          T val_);

    template<typename V>
    explicit
    QDiag(QDiag<V> const& D)
      : store(D.begin(),D.end()),
        val(D.val),
        length(D.length)
        { }

    QDiag(size_t size)
      : store(size,0.),
        length(size)
        { }

    explicit operator bool() const { return !store.empty(); }

    bool
    allSame() const { return store.empty(); }

    value_type*
    data() { return store.data(); }

    value_type const*
    data() const { return store.data(); }

    size_t
    size() const { return length; }

    iterator
    begin() { return store.begin(); }

    iterator
    end() { return store.end(); }

    const_iterator
    begin() const { return store.begin(); }

    const_iterator
    end() const { return store.end(); }

    };

const char*
typeNameOf(QDiagReal const& d);
const char*
typeNameOf(QDiagCplx const& d);

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

template<typename T>
void
write(std::ostream & s, QDiag<T> const& dat)
    {
    itensor::write(s,dat.val);
    itensor::write(s,dat.length);
    itensor::write(s,dat.store);
    }

template<typename T>
void
read(std::istream & s, QDiag<T> & dat)
    {
    itensor::read(s,dat.val);
    itensor::read(s,dat.length);
    itensor::read(s,dat.store);
    }
 
template<typename T>
Cplx
doTask(GetElt& G, QDiag<T> const& D);

template<typename T>
QN
doTask(CalcDiv const& C, QDiag<T> const& d);

template<typename F, typename T>
void
doTask(ApplyIT<F> & A, QDiag<T> const& d, ManageStore & m)
    {
    using new_type = ApplyIT_result_of<T,F>;
    if(switchesType<T>(A))
        {
        auto *nd = m.makeNewData<QDiag<new_type>>(d.size());
        assert(nd->store.size() == d.store.size());
        A(d.val,nd->val);
        for(auto n : range(d.store.size()))
            {
            A(d.store[n],nd->store[n]);
            }
        }
    else
        {
        auto *md = m.modifyData(d);
        A(md->val);
        for(auto& el : md->store) A(el);
        }
    }

template<typename F, typename T>
void
doTask(VisitIT<F> & V, QDiag<T> const& d)
    {
    if(d.allSame()) 
        {
        for(decltype(d.length) j = 0; j < d.length; ++j) 
            {
            detail::call<void>(V.f,V.scale_fac * d.val);
            }
        }
    else
        {
        for(auto& elt : d.store) 
            {
            detail::call<void>(V.f,elt*V.scale_fac);
            }
        }
    }

template<typename F>
void
doTask(GenerateIT<F,Real>& G, QDiagReal & D)
    {
    if(D.allSame())
        {
        D.val = 0;
        D.store.resize(D.length);
        }
    stdx::generate(D,G.f);
    }

template<typename F>
void
doTask(GenerateIT<F,Real>& G, QDiagCplx const& D, ManageStore & m)
    {
    auto *nD = m.makeNewData<QDiagReal>();
    nD->store.resize(D.length);
    stdx::generate(*nD,G.f);
    }

template<typename F>
void
doTask(GenerateIT<F,Cplx>& G, QDiagReal const& D, ManageStore & m)
    {
    auto *nD = m.makeNewData<QDiagCplx>();
    nD->store.resize(D.length);
    stdx::generate(*nD,G.f);
    }

template<typename F>
void
doTask(GenerateIT<F,Cplx>& G, QDiagCplx & D)
    {
    if(D.allSame())
        {
        D.val = 0;
        D.store.resize(D.length);
        }
    stdx::generate(D,G.f);
    }

template<typename T>
Cplx
doTask(SumEls, QDiag<T> const& d);

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
doTask(PrintIT& P, QDiag<T> const& d);

auto inline constexpr
doTask(StorageType const& S, QDiagReal const& d) ->StorageType::Type { return StorageType::QDiagReal; }

auto inline constexpr
doTask(StorageType const& S, QDiagCplx const& d) ->StorageType::Type { return StorageType::QDiagCplx; }

template<typename VA, typename VB>
void
doTask(Contract& Con,
       QDiag<VA> const& A,
       QDense<VB> const& B,
       ManageStore& m);

template<typename VA, typename VB>
void
doTask(Contract& Con,
       QDense<VA> const& A,
       QDiag<VB> const& B,
       ManageStore& m);

template<typename T>
void
doTask(Order const& P,
       QDiag<T> & dA) { }

template<typename Indexable>
std::tuple<size_t,size_t,IntArray>
diagBlockBounds(IndexSet const& is,
                Indexable const& block_ind)
    {
    long nb = -1;
    long ne = std::numeric_limits<long>::max();
    auto starts = IntArray(order(is),0);
    for(auto n : range(is))
        {
        for(auto j : range(block_ind[n])) starts[n] += is[n].blocksize0(j);
        nb = std::max(nb,starts[n]);
        ne = std::min(ne,starts[n]+is[n].blocksize0(block_ind[n]));
        }
    for(auto n : range(is))
        {
        starts[n] = nb-starts[n];
        }
    return std::make_tuple(nb,ne,starts);
    }

template<typename V, typename Indexable>
DataRange<const V>
getBlock(QDiag<V> const& D,
         IndexSet const& is,
         Indexable const& block_ind)
    {
    long nb = -1, ne = -1;
    auto starts = IntArray{};

    if(block_ind.size()==0 && order(is)==0)
        {
        nb = 0;
        ne = 1;
        }
    else
        {
        //print("block_ind:"); for(auto& el : block_ind) print(" ",el); println();
        std::tie(nb,ne,starts) = diagBlockBounds(is,block_ind);
        if(nb >= ne) return DataRange<const V>{};
        }

    if(D.allSame())
        {
        return DataRange<const V>(&D.val,1ul);
        }
    //printfln("nb=%d ne=%d",nb,ne);
    return sliceData(makeDataRange(D.data(),D.size()),nb,ne);
    }

template<typename V, typename Indexable>
DataRange<V>
getBlock(QDiag<V> & D,
         IndexSet const& is,
         Indexable const& block_ind)
    {
    auto const& cD = D;
    auto cdr = getBlock(cD,is,block_ind);
    //const_cast safe here because we know
    //original QDiag d is non-const
    auto ncd = const_cast<V*>(cdr.data());
    return DataRange<V>{ncd,cdr.size()};
    }

//template<typename V>
//ITensor
//doTask(ToITensor & T, QDiag<V> const& d);

template<typename V>
bool
doTask(IsEmpty, QDiag<V> const& d) { return (d.length == 0ul); }

template<typename T>
bool
doTask(IsDense,
       QDiag<T> const& d);

template<typename V>
void
doTask(RemoveQNs &, QDiag<V> const&, ManageStore &);

template<typename T>
void
doTask(ToDense &, QDiag<T> const&, ManageStore &);

} //namespace itensor

#endif

