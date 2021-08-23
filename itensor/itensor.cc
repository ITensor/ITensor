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
#include "itensor/util/print_macro.h"
//#include "itensor/util/iterate.h"
#include "itensor/util/safe_ptr.h"
#include "itensor/itensor.h"
#include "itensor/tensor/lapack_wrap.h"
#include "itensor/tensor/contract.h"

using std::array;
using std::ostream;
using std::vector;
using std::move;
using std::string;

namespace itensor {

//
// ITensor Constructors
//

ITensor::
ITensor(IndexSet const& is)
  : is_(is)
    { 
    IF_USESCALE(scale_ = LogNum(1.);)
    }

ITensor::
ITensor(IndexSet iset,
        storage_ptr&& pdat,
        LogNum const& scale)
    :
    is_(std::move(iset)),
    store_(std::move(pdat))
    { 
    IF_USESCALE(scale_ = scale;)
    }

ITensor::
ITensor(std::initializer_list<Index> inds)
  : is_(inds)
    { 
    IF_USESCALE(scale_ = LogNum(1.);)
    }

    
ITensor::
ITensor(Cplx val) 
    { 
    IF_USESCALE(scale_ = LogNum(1.);)
    if(val.imag() == 0)
        {
        store_ = newITData<ScalarReal>(val.real());
        }
    else
        {
        store_ = newITData<ScalarCplx>(val);
        }
    }

ITensor::
ITensor(QN q, IndexSet const& is)
  :
  is_(std::move(is))
  {
  store_ = newITData<QDenseReal>(is,q); 
  }

Cplx ITensor::
eltC() const
    {
    if(inds().order() != 0)
        {
        PrintData(inds());
        Error(format("Wrong number of IndexVals passed to elt/eltC (expected 0, got %d)",inds().order()));
        }
    constexpr size_t size = 0;
    auto inds = IntArray(size);
    auto z = itensor::doTask(GetElt{is_,inds},store_);
#ifndef USESCALE
    return z;
#else
    try {
        return z*scale_.real0(); 
        }
    catch(TooBigForReal const& e)
        {
        println("too big for real in eltC(...), scale = ",scale());
        throw e;
        }
    catch(TooSmallForReal const&)
        {
        println("warning: too small for real in eltC(...)");
        return Cplx(0.,0.);
        }
    return Cplx(NAN,NAN);
#endif
    }

Cplx ITensor::
eltC(std::vector<IndexVal> const& ivs) const
    {
    if(!store()) Error("tensor storage unallocated");

    auto size = ivs.size();
    if(size != size_t(inds().order()))
        {
        println("---------------------------------------------");
        println("Tensor indices = \n",inds(),"\n");
        println("---------------------------------------------");
        println("Indices provided = ");
        for(auto& iv : ivs) println(iv.index);
        println("---------------------------------------------");
        Error(format("Wrong number of IndexVals passed to elt/eltC (expected %d, got %d)",inds().order(),size));
        }

    auto ints = IntArray(size);
    detail::permute_map(inds(),ivs,ints,
                [](IndexVal const& iv) { return iv.val-1; });
    auto z = itensor::doTask(GetElt{inds(),ints},store_);
#ifndef USESCALE
    return z;
#else
    try {
        return z*scale().real0();
        }
    catch(TooBigForReal const& e)
        {
        println("too big for real in eltC(...), scale = ",scale());
        throw e;
        }
    catch(TooSmallForReal const&)
        {
        println("warning: too small for real in eltC(...)");
        return Cplx(0.,0.);
        }
    return Cplx(NAN,NAN);
#endif
    }



void ITensor::
set(std::vector<IndexVal> const& ivals,
    Cplx val)
    {
    auto size = ivals.size();
    if(size != size_t(inds().order())) 
        {
        println("---------------------------------------------");
        println("Tensor indices = \n",inds(),"\n");
        println("---------------------------------------------");
        println("Indices provided = ");
        for(auto& iv : ivals) println(iv.index);
        println("---------------------------------------------");
        Error(format("Wrong number of IndexVals passed to set (expected %d, got %d)",
                     inds().order(),size));
        }
    auto inds = IntArray(is_.order(),0);
    detail::permute_map(is_,ivals,inds,
                        [](IndexVal const& iv) { return iv.val-1; });
    if(!store_) detail::allocReal(*this,inds); 
    scaleTo(1.);
    if(val.imag()==0.0)
        {
        doTask(SetElt<Real>{val.real(),is_,inds},store_);
        }
    else
        {
        doTask(SetElt<Cplx>{val,is_,inds},store_);
        }
    }

void ITensor::
set(Cplx val)
    {
    if(0 != size_t(inds().order())) 
        {
        Error(format("Wrong number of IndexVals passed to set (expected %d, got 0)",
                     inds().order()));
        }
    auto inds = IntArray(0,1);
    if(!store_) detail::allocReal(*this,inds); 
    scaleTo(1.);
    if(val.imag()==0.)
        {
        doTask(SetElt<Real>{val.real(),is_,inds},store_);
        }
    else
        {
        doTask(SetElt<Cplx>{val,is_,inds},store_);
        }
    }

void ITensor::
set(std::vector<int> const& ints,
    Cplx val)
    {
    auto size = ints.size();
    if(size != size_t(inds().order())) 
        {
        println("---------------------------------------------");
        println("Tensor indices = \n",inds(),"\n");
        println("---------------------------------------------");
        println("Indices provided = ");
        for(auto& iv : ints) println(iv);
        println("---------------------------------------------");
        Error(format("Wrong number of IndexVals passed to set (expected %d, got %d)",
                     inds().order(),size));
        }
    auto inds = IntArray(is_.order(),0);
    for(auto i : range(size))
        inds[i] = ints[i]-1;
    //TODO: if !store_ and !is_real, call allocCplx instead
    //detail::permute_map(is_,ivals,inds,
    //                    [](IndexVal const& iv) { return iv.val-1; });
    if(!store_) detail::allocReal(*this,inds); 
    scaleTo(1.);
    if(val.imag()==0.0)
        {
        doTask(SetElt<Real>{val.real(),is_,inds},store_);
        }
    else
        {
        doTask(SetElt<Cplx>{val,is_,inds},store_);
        }
    }

ITensor& ITensor::
conj()
    {
    doTask(Conj{},store_);
    return *this;
    }


ITensor& ITensor::
dag() 
    { 
    if(hasQNs(is_)) is_.dag();
    return conj(); 
    }

ITensor& ITensor::
takeReal()
    {
    doTask(TakeReal{},store_);
    return *this;
    }

ITensor& ITensor::
takeImag()
    {
    doTask(TakeImag{},store_);
    return *this;
    }

ITensor& ITensor::
makeCplx()
    {
    doTask(MakeCplx{},store_);
    return *this;
    }

ITensor& ITensor::
randomize(Args const& args)
    {
    if(!this->store()) detail::allocReal(*this);
#ifdef DEBUG
    if(!(*this)) Error("default initialized tensor in randomize");
#endif
    auto cplx = args.getBool("Complex",false);
    if(cplx) this->generate(detail::quickranCplx);
    else     this->generate(detail::quickran);
    return *this;
    }

size_t
nnzblocks(ITensor const& A)
    {
    if(hasQNs(A)) return doTask(NNZBlocks{},A.store());
    return 1;
    }

long
nnz(ITensor const& A)
    {
    return doTask(NNZ{},A.store());
    }

ITensor& ITensor::
fixBlockDeficient()
    {
    if(itensor::hasQNs(*this) && itensor::isDense(*this))
        {
        // If the ITensor has QNs, we may need to
        // expand the storage since it may be
        // block deficient
        auto itflux = itensor::flux(*this);
        auto [bofs,size] = getBlockOffsets(inds(),itflux);
        if(bofs.size() != itensor::nnzblocks(*this))
            {
            // Make a copy of the original ITensor
            auto Torig = *this;
            // The QDense storage is block deficient,
            // need to allocate new memory.
            if(isReal(*this))
              store_ = newITData<QDense<Real>>(bofs,size);
            else
              store_ = newITData<QDense<Cplx>>(bofs,size);
            *this += Torig;
            }
        }
    return *this;
    }

ITensor& ITensor::
fill(Cplx z)
    {
    if(!store_) 
        {
        if(is_) detail::allocReal(*this);
        else Error("Can't fill default-constructed tensor");
        }
    IF_USESCALE(scale_ = scale_type(1.);)
    if(itensor::hasQNs(*this))
        {
        // If the ITensor has QNs, we may need to
        // expand the storage since it may be
        // block deficient
        auto itflux = itensor::flux(*this);
        auto [bofs,size] = getBlockOffsets(inds(),itflux);
        if(bofs.size() != itensor::nnzblocks(*this))
            {
            // The QDense storage is block deficient,
            // need to allocate new memory.
            // Make the new memory undefined since it
            // will be overwritten anyway.
            if(z.imag() == 0)
              store_ = newITData<QDense<Real>>(undef,bofs,size);
            else
              store_ = newITData<QDense<Cplx>>(undef,bofs,size);
            }
        }
    if(z.imag() == 0)
        doTask(Fill<Real>{z.real()},store_);
    else
        doTask(Fill<Cplx>{z},store_);
    return *this;
    }

#ifdef USESCALE
void ITensor::
scaleTo(scale_type const& newscale)
    {
    if(scale_ == newscale) return;
    if(newscale.sign() == 0) Error("Trying to scale an ITensor to a 0 scale");
    scale_ /= newscale;
    doTask(Mult<Real>{scale_.real0()},store_);
    scale_ = newscale;
    }

void ITensor::
scaleTo(Real newscale) { scaleTo(LogNum{newscale}); }
#endif

void ITensor::
swap(ITensor & other)
    {
    is_.swap(other.is_);
    store_.swap(other.store_);
    IF_USESCALE(scale_.swap(other.scale_);)
    }

ITensor
conj(ITensor T)
    {
    T.conj();
    return T;
    }

ITensor
dag(ITensor T)
    {
    T.dag();
    return T;
    }

bool
hasQNs(ITensor const& T) { return hasQNs(inds(T)); }

long
order(ITensor const& T) { return order(inds(T)); }

IndexSet const&
inds(ITensor const& A) { return A.inds(); }

long
minDim(ITensor const& T) { return minDim(inds(T)); }

long
maxDim(ITensor const& T) { return maxDim(inds(T)); }

std::vector<IndexSet>
inds(std::vector<ITensor> const& A)
  {
  auto is = std::vector<IndexSet>(A.size());
  for( auto i : range(A.size()) )
      is[i] = inds(A[i]);
  return is;
  }

Index const&
index(ITensor const& A, RangeT<Index>::size_type I) { return A.index(I); }

Index
findIndex(ITensor const& T,
          TagSet const& tsmatch)
    {
    return findIndex(inds(T),tsmatch);
    }

IndexSet
findInds(ITensor const& T,
         TagSet const& tsmatch)
    {
    return findInds(inds(T),tsmatch);
    }

IndexSet
commonInds(ITensor const& A,
           ITensor const& B)
    {
    return commonInds(inds(A),inds(B));
    }

IndexSet
commonInds(ITensor const& A,
           ITensor const& B,
           TagSet const& tsmatch)
    {
    return findInds(commonInds(inds(A),inds(B)),tsmatch);
    }

Index
commonIndex(ITensor const& A, 
            ITensor const& B)
    {
    return findIndex(commonInds(inds(A),inds(B)));
    }

Index
commonIndex(ITensor const& A, 
            ITensor const& B, 
            TagSet const& tsmatch)
    {
    return findIndex(commonInds(inds(A),inds(B)),tsmatch);
    }

IndexSet
uniqueInds(ITensor const& A,
           ITensor const& B)
    {
    return uniqueInds(inds(A),inds(B));
    }

IndexSet
uniqueInds(ITensor const& A,
           std::vector<ITensor> const& B)
    {
    return uniqueInds(inds(A),inds(B));
    }

IndexSet
uniqueInds(ITensor const& A,
           std::initializer_list<ITensor> B)
    {
    return uniqueInds(A,std::vector<ITensor>(B));
    }

Index
uniqueIndex(ITensor const& A, 
            ITensor const& B)
    {
    return findIndex(uniqueInds(A,B));
    }

Index
uniqueIndex(ITensor const& A, 
            ITensor const& B, 
            TagSet const& tsmatch)
    {
    return findIndex(uniqueInds(A,B),tsmatch);
    }

Index
uniqueIndex(ITensor const& A,
            std::vector<ITensor> const& B)
    {
    return findIndex(uniqueInds(A,B));
    }

Index
uniqueIndex(ITensor const& A,
            std::vector<ITensor> const& B,
            TagSet const& tsmatch)
    {
    return findIndex(uniqueInds(A,B),tsmatch);
    }

Index
uniqueIndex(ITensor const& A,
            std::initializer_list<ITensor> B)
    {
    return uniqueIndex(A,std::vector<ITensor>(B));
    }

Index
uniqueIndex(ITensor const& A,
            std::initializer_list<ITensor> B,
            TagSet const& tsmatch)
    {
    return uniqueIndex(A,std::vector<ITensor>(B),tsmatch);
    }

ITensor
setTags(ITensor A,
        TagSet const& ts,
        IndexSet const& is)
  {
  A.setTags(ts,is);
  return A;
  }

ITensor
noTags(ITensor A,
       IndexSet const& is)
  {
  A.noTags(is);
  return A;
  }

ITensor
addTags(ITensor A,
        TagSet const& ts,
        IndexSet const& is)
  {
  A.addTags(ts,is);
  return A;
  }

ITensor
removeTags(ITensor A,
           TagSet const& ts,
           IndexSet const& is)
  {
  A.removeTags(ts,is);
  return A;
  }

ITensor
replaceTags(ITensor A,
            TagSet const& ts1,
            TagSet const& ts2,
            IndexSet const& is)
  {
  A.replaceTags(ts1,ts2,is);
  return A;
  }

ITensor
swapTags(ITensor A,
         TagSet const& ts1,
         TagSet const& ts2,
         IndexSet const& is)
  {
  A.swapTags(ts1,ts2,is);
  return A;
  }

ITensor
prime(ITensor A,
      int plev,
      IndexSet const& is)
  {
  A.prime(plev,is);
  return A;
  }

ITensor
prime(ITensor A,
      IndexSet const& is)
  {
  A.prime(is);
  return A;
  }

ITensor
setPrime(ITensor A,
         int plev,
         IndexSet const& is)
  {
  A.setPrime(plev,is);
  return A;
  }

ITensor
noPrime(ITensor A,
        IndexSet const& is)
  {
  A.noPrime(is);
  return A;
  }

ITensor& ITensor::
permute(IndexSet const& iset)
    {
    auto& A = *this;
    auto Ais = A.inds();
    auto r = Ais.order();

    if(size_t(r) != size_t(iset.order()))
        {
        println("---------------------------------------------");
        println("Tensor indices = \n",Ais,"\n");
        println("---------------------------------------------");
        println("Indices provided = \n",iset,"\n");
        println("---------------------------------------------");
        Error(format("Wrong number of Indexes passed to permute (expected %d, got %d)",r,iset.order()));
        }

    // Get permutation
    auto P = Permutation(r);
    calcPerm(Ais,iset,P);
    if(isTrivial(P))
        {
        return A;
        }
    // If not trivial, use permutation to get new index set
    // This is necessary to preserve the proper arrow direction of QN Index
    auto bind = RangeBuilderT<IndexSet>(r);
    for(auto i : range(r))
        {
        bind.setIndex(P.dest(i),Ais[i]);
        }
    auto Bis = bind.build();

    auto O = Order{P,Ais,Bis};
    if(A.store())
        {
        doTask(O, A.store());
        }

    A.is_.swap(Bis);

    return A;
    }

ITensor
permute(ITensor A,
        IndexSet const& is)
    {
    A.permute(is);
    return A;
    }

ITensor& ITensor::
replaceInds(IndexSet const& is1,
            IndexSet const& is2)
    {
#ifdef DEBUG
    if( itensor::order(is1) != itensor::order(is2) ) Error("In replaceInds, must replace with equal number of Indices");
#endif
    auto& T = *this;
    // Add a random prime to account for possible
    // Index swaps
    auto plev_temp = 43218432;
    auto is2p = itensor::prime(is2,plev_temp);
    auto isT = itensor::inds(T);
    for(auto& J : isT)
        {
        for(auto i : range(itensor::order(is1)))
            {
            if( is1[i] && (J == is1[i]) )
                {
                if( dim(J) != dim(is2[i]) )
                    {
                    printfln("Old dim = %d",dim(J));
                    printfln("New dim would be = %d",dim(is2[i]));
                    throw ITError("Mismatch of index dimension in replaceInds");
                    }
                // Make the arrow directions correct
                J.dag();
                auto& Jnew = is2p[i];
                if( dir(J)==dir(Jnew) ) Jnew.dag();
                T *= delta(J,Jnew);
                break;
                }
            }
        }

    // Bring the prime levels back down to the original
    // desired ones
    T.prime(-plev_temp,is2p);

    return T;
    }

ITensor
replaceInds(ITensor T,
            IndexSet const& is1,
            IndexSet const& is2)
    {
    T.replaceInds(is1,is2);
    return T;
    }

ITensor& ITensor::
swapInds(IndexSet const& is1,
         IndexSet const& is2)
    {
#ifdef DEBUG
    if( itensor::order(is1) != itensor::order(is2) ) Error("In swapInds, must swap equal numbers of Indices");
#endif
    auto& T = *this;
    T.replaceInds({is1,is2},{is2,is1});
    return T;
    }

ITensor
swapInds(ITensor T,
         IndexSet const& is1,
         IndexSet const& is2)
    {
    T.swapInds(is1,is2);
    return T;
    }

Real
norm(ITensor const& T)
    {
#ifdef DEBUG
    if(!T) Error("Default initialized tensor in norm(ITensor)");
#endif

#ifndef USESCALE
    return doTask(NormNoScale{},T.store());
#else
    auto fac = std::fabs(T.scale().real0());
    return fac * doTask(NormNoScale{},T.store());
#endif
    }

void ITensor::
write(std::ostream& s) const
    {
    itensor::write(s,inds());
    itensor::write(s,scale());
    auto type = StorageType::Null;
    if(store()) 
        {
        type = doTask(StorageType{},store());
        }
    itensor::write(s,type);
    if(store()) 
        {
        doTask(Write{s},store());
        }
    }

void ITensor::
read(std::istream& s)
    {
    itensor::read(s,is_);
    LogNum scale;
    itensor::read(s,scale);
    IF_USESCALE(scale_ = scale;)
    auto type = StorageType::Null;
    itensor::read(s,type);
    if(type==StorageType::Null) { /*intentionally left blank*/  }
    else if(type==StorageType::DenseReal) { store_ = readType<DenseReal>(s); }
    else if(type==StorageType::DenseCplx) { store_ = readType<DenseCplx>(s); }
    else if(type==StorageType::Combiner) { store_ = readType<Combiner>(s); }
    else if(type==StorageType::DiagReal) { store_ = readType<Diag<Real>>(s); }
    else if(type==StorageType::DiagCplx) { store_ = readType<Diag<Cplx>>(s); }
    else if(type==StorageType::QDenseReal) { store_ = readType<QDense<Real>>(s); }
    else if(type==StorageType::QDenseCplx) { store_ = readType<QDense<Cplx>>(s); }
    else if(type==StorageType::QDiagReal) { store_ = readType<QDiag<Real>>(s); }
    else if(type==StorageType::QDiagCplx) { store_ = readType<QDiag<Cplx>>(s); }
    else if(type==StorageType::QCombiner) { store_ = readType<QCombiner>(s); }
    else if(type==StorageType::ScalarReal) { store_ = readType<ScalarReal>(s); }
    else if(type==StorageType::ScalarCplx) { store_ = readType<ScalarCplx>(s); }
    else
        {
        Error("Unrecognized type when reading tensor from istream");
        }
    }

#ifdef ITENSOR_USE_HDF5

void
h5_write(h5::group parent, std::string const& name, ITensor const& T)
    {
    auto g = parent.create_group(name);
    h5_write_attribute(g,"type","ITensor",true);
    h5_write_attribute(g,"version",long(1));
    h5_write(g,"inds",T.inds());
    doTask(H5Write(g,"storage"),T.store());
    }

void
h5_read(h5::group parent, std::string const& name, ITensor & I)
    {
    auto g = parent.open_group(name);
    auto type = h5_read_attribute<string>(g,"type");
    if(type != "ITensor") Error("Group does not contain ITensor data in HDF5 file");

    auto is = h5_read<IndexSet>(g,"inds");

    std::string store_name;
    if(g.has_subgroup("storage"))    store_name = "storage";
    else if(g.has_subgroup("store")) store_name = "store";
    else error("Expected ITensor HDF5 data to have group named \"storage\" or \"store\"");

    auto sg = g.open_group(store_name);
    auto s_type = h5_read_attribute<string>(sg,"type");
    ITensor::storage_ptr store;
    if(s_type == "Dense{Float64}") store = h5_readStore<DenseReal>(g,store_name); 
    else if(s_type == "Dense{ComplexF64}") store = h5_readStore<DenseCplx>(g,store_name); 
    else if(s_type == "BlockSparse{Float64}") store = h5_readStore<QDenseReal>(g,store_name); 
    else if(s_type == "BlockSparse{ComplexF64}") store = h5_readStore<QDenseCplx>(g,store_name); 
    else error(format("Reading of ITensor storage type %s not yet supported",s_type));

    I = ITensor(is,std::move(store));
    }

#endif //ITENSOR_USE_HDF5


namespace detail {

void
allocReal(ITensor& T)
    {
    if(hasQNs(T)) Error("Can't allocate quantum ITensor with undefined divergence");
    T.store() = newITData<DenseReal>(dim(inds(T)),0);
    }

void
allocReal(ITensor& T, IntArray const& ints)
    {
    if(not hasQNs(T))
        {
        T.store() = newITData<DenseReal>(dim(inds(T)),0);
        }
    else
        {
        QN div;
        for(auto i : range(inds(T)))
            {
            auto iv = (inds(T)[i])(1+ints[i]);
            div += qn(iv)*dir(iv);
            }
        T.store() = newITData<QDenseReal>(inds(T),div);
        }
    }

void
allocCplx(ITensor& T)
    {
    if(hasQNs(T)) Error("Can't allocate quantum ITensor with undefined divergence");
    T.store() = newITData<DenseCplx>(dim(inds(T)),0);
    }


void
checkArrows(IndexSet const& is1,
            IndexSet const& is2,
            bool shouldMatch = false)
    {
    if(hasQNs(is1) && hasQNs(is2))
        {
        for(auto I1 : is1)
        for(auto I2 : is2)
            {
            if(I1 == I2)
                {
                auto cond = shouldMatch ^ (I1.dir() == I2.dir());
                if(cond)
                    {
                    println("----------------------------------------");
                    println("IndexSet 1 = \n",is1);
                    println("----------------------------------------");
                    println("IndexSet 2 = \n",is2);
                    println("----------------------------------------");
                    printfln("Mismatched QN Index from set 1 %s",I1);
                    printfln("Mismatched QN Index from set 2 %s",I2);
                    Error("Mismatched QN Index arrows");
                    }
                }
            }
        }
    }

void
checkSameDiv(ITensor const& T1,
             ITensor const& T2)
    {
    if(hasQNs(T1) && hasQNs(T2))
        {
        if(div(T1) != div(T2)) 
            {
            Error(format("div(T1)=%s must equal div(T2)=%s when adding T1+T2",div(T1),div(T2)));
            }
        }
    }

} //namespace detail

//TODO: implement proper Dense*QDense to avoid conversion cost
ITensor& ITensor::
operator*=(ITensor const& R)
    {
    auto& L = *this;

    if(!L || !R) Error("Default constructed ITensor in product");

    if(L.order() == 0)
        {
        auto z = L.eltC();
        *this = R*z;
        return *this;
        }
    else if(R.order()==0)
        {
        auto z = R.eltC();
        *this *= z;
        return *this;
        }

    if(Global::checkArrows()) detail::checkArrows(L.inds(),R.inds());

    //TODO: create a proper doTask(Contract,Dense,QDense)
    auto hqL = hasQNs(L);
    auto hqR = hasQNs(R);
    auto Rdense = R;
    if(hqL && !hqR) L = removeQNs(L);
    else if(!hqL && hqR) Rdense = removeQNs(Rdense);

    auto C = doTask(Contract{L.inds(),Rdense.inds()},
                    L.store(),
                    Rdense.store());

#ifdef USESCALE
    L.scale_ *= Rdense.scale();
    if(!std::isnan(C.scalefac)) L.scale_ *= C.scalefac;
#endif

#ifdef DEBUG
    //Check for duplicate indices
    checkIndexSet(C.Nis);
#endif

    L.is_.swap(C.Nis);

    return L;
    }

#ifndef USESCALE

//for Diag and QDiag
//Diag  -> Dense
//QDiag -> QDense
ITensor
toDense(ITensor T)
    {
    if(T.store()) doTask(ToDense{T.inds()},T.store());
    return ITensor{move(T.inds()),move(T.store()),T.scale()};
    }

bool
isDense(ITensor const& T)
    {
    return doTask(IsDense{},T.store());
    }

//TODO: make this use a RemoveQNs task type that does:
//QDense -> Dense
//QDiag  -> Diag
ITensor
removeQNs(ITensor T)
    {
    if(not hasQNs(T)) return T;
    if(T.store()) doTask(RemoveQNs{inds(T)},T.store());
    auto nis = inds(T);
    nis.removeQNs();
    return ITensor{move(nis),move(T.store()),T.scale()};
    }

ITensor& ITensor::
operator*=(Real r)
    {
    doTask(Mult<Real>{r},store_);
    return *this;
    }

ITensor& ITensor::
operator/=(Real r)
    {
    auto fac = 1./r;
    doTask(Mult<Real>{fac},store_);
    return *this;
    }

#endif

ITensor& ITensor::
operator*=(Cplx z)
    {
    if(z.imag() == 0) return operator*=(z.real());
    doTask(Mult<Cplx>{z},store_);
    return *this;
    }

//Non-contracting product
ITensor& ITensor::
operator/=(ITensor const& R)
    {
    auto& L = *this;

    if(!L || !R) Error("Default constructed ITensor in product");

    if(Global::checkArrows()) detail::checkArrows(L.inds(),R.inds(),true);

    auto C = doTask(NCProd{L.inds(),R.inds()},
                    L.store(),
                    R.store());

#ifdef USESCALE
    L.scale_ *= R.scale();
    if(!std::isnan(C.scalefac)) L.scale_ *= C.scalefac;
#endif

#ifdef DEBUG
    //Check for duplicate indices
    checkIndexSet(C.Nis);
#endif

    L.is_.swap(C.Nis);

    return L;
    }

#ifndef USESCALE

ITensor ITensor::
operator-() const
    {
    auto res = *this;
    doTask(Mult<Real>(-1.),res.store());
    return res;
    }

#else

ITensor ITensor::
operator-() const
    {
    auto res = *this;
    res.scale_.negate();
    return res;
    }

#endif

ITensor
operator*(ITensor A, ITensor const& B) { A *= B; return A; }
ITensor
operator*(ITensor const& A, ITensor&& B) { B *= A; return std::move(B); }
ITensor
operator*(ITensor T, Real fac) { T *= fac; return T; }
ITensor
operator*(Real fac, ITensor T) { T *= fac; return T; }
ITensor
operator*(ITensor T, Complex fac) { T *= fac; return T; }
ITensor
operator*(Complex fac, ITensor T) { T *= fac; return T; }
ITensor
operator/(ITensor T, Real fac) { T /= fac; return T; }
ITensor
operator/(ITensor T, Complex fac) { T /= fac; return T; }
ITensor
operator+(ITensor A, ITensor const& B) { A += B; return A; }
ITensor
operator+(ITensor const& A, ITensor&& B) { B += A; return std::move(B); }
ITensor
operator-(ITensor A, ITensor const& B) { A -= B; return A; }
ITensor
operator-(ITensor const& A, ITensor&& B) { B -= A; B *= -1; return std::move(B); }
ITensor
operator/(ITensor A, ITensor const& B) { A /= B; return A; }
ITensor
operator/(ITensor const& A, ITensor && B) { B /= A; return std::move(B); }

// Create some sparse tensors to help with
// a partial direct sum
Index
directSumITensors(Index const& i,
                  Index const& j,
                  ITensor& D1,
                  ITensor& D2,
                  Args const& args = Args::global())
  {
  auto ij = directSum(i,j,args);
  if(not hasQNs(i) && not hasQNs(j))
    {
    D1 = delta(i,ij);
    auto S = Matrix(dim(j),dim(ij));
    for(auto jj : range(dim(j)))
      {
      S(jj,dim(i)+jj) = 1;
      }
    D2 = matrixITensor(std::move(S),j,ij);
    }
  else
    {
    D1 = ITensor(dag(i),ij);
    int n = 1;
    auto nblock_i = nblock(i);
    for(auto iq1 : range1(nblock_i))
        {
        auto blocksize_i = blocksize(i,iq1);
        auto blocksize_ij = blocksize(ij,n);
        auto D = Matrix(blocksize_i,blocksize_ij);
        auto minsize = std::min(blocksize_i,blocksize_ij);
        for(auto ii : range(minsize)) D(ii,ii) = 1.0;
        getBlock<Real>(D1,{iq1,n}) &= D;
        ++n;
        }
    auto nblock_j = nblock(j);
    D2 = ITensor(dag(j),ij);
    for(auto iq2 : range1(nblock_j))
        {
        auto blocksize_j = blocksize(j,iq2);
        auto blocksize_ij = blocksize(ij,n);
        auto D = Matrix(blocksize_j,blocksize_ij);
        auto minsize = std::min(blocksize_j,blocksize_ij);
        for(auto ii : range(minsize)) D(ii,ii) = 1.0;
        getBlock<Real>(D2,{iq2,n}) &= D;
        ++n;
        }
    }
  return ij;
  }

std::tuple<ITensor,IndexSet>
directSum(ITensor const& A, ITensor const& B,
          IndexSet const& I, IndexSet const& J,
          Args const& args)
  {
  if( order(I) != order(J) ) Error("In directSum(ITensor, ITensor, ...), must sum equal number of indices");
  auto AD = A;
  auto BD = B;
  auto newinds = IndexSetBuilder(I.size());
  for( auto n : range1(order(I)) )
    {
    auto In = I(n);
    auto Jn = J(n);
    if( dir(A,In) != dir(In) ) In.dag();
    if( dir(B,Jn) != dir(Jn) ) Jn.dag();
    ITensor D1, D2;
    auto IJn = directSumITensors(In,Jn,D1,D2,args);
    newinds.nextIndex(IJn);
    AD *= D1;
    BD *= D2;
    }
  auto IJ = newinds.build();
  auto C = AD+BD;
  return std::make_tuple(C,IJ);
  }

// Direct sum A and B over indices i on A and j on B
std::tuple<ITensor,Index>
directSum(ITensor const& A, ITensor const& B,
          Index const& i, Index const& j,
          Args const& args)
  {
  auto [C,IJ] = directSum(A,B,IndexSet(i),IndexSet(j),args);
  auto ij = IJ.index(1);
  return std::make_tuple(C,ij);
  }

void
daxpy(ITensor & L,
      ITensor const& R,
      Real alpha)
    {
    if(L.order() != R.order()) Error("ITensor::operator+=: different number of indices");
    if(nnzblocks(R) == 0) return;
    detail::checkSameDiv(L,R);

    using permutation = typename PlusEQ::permutation;

    auto P = permutation(inds(L).size());

    try {
        calcPerm(inds(R),inds(L),P);
        }
    catch(std::exception const& e)
        {
        println("L = ",L);
        println("R = ",R);
        Error("ITensoITensor::operator+=: different index structure");
        }

    if(Global::checkArrows()) 
        {
        auto shouldMatch = true;
        detail::checkArrows(inds(L),inds(R),shouldMatch);
        }

    if(!L.store()) Error("L not initialized in daxpy");

#ifdef USESCALE
    if(L.scale().magnitudeLessThan(R.scale())) 
        {
        L.scaleTo(R.scale()); 
        }
    else
        {
        alpha *= (R.scale()/L.scale()).real();
        }
#endif

    auto PEq = PlusEQ{P,inds(L),inds(R),alpha};
    doTask(PEq,L.store(),R.store());
    } 

ITensor& ITensor::
operator+=(ITensor const& R)
    {
    auto& L = *this;
    if(!L || !L.store()) { return (L=R); } //special case when this (L) is not initialized
    if(!R) Error("Right-hand-side of ITensor += is default constructed");
    if(&L == &R) return operator*=(2.);

    daxpy(L,R,1.);
    
    return L;
    } 

ITensor& ITensor::
operator-=(ITensor const& R)
    {
    auto& L = *this;
    if(!L || !L.store()) { return (L = -R); } //special case when this (L) is not initialized
    if(!R) Error("Right-hand-side of ITensor -= is default constructed");
    if(&L == &R) 
        { 
        L *= 0.;
        return L;
        }
    daxpy(L,R,-1.);
    return L;
    }

detail::IndexValIter
iterInds(ITensor const& T)
    {
    return iterInds(inds(T));
    }

ITensor
multSiteOps(ITensor A, ITensor const& B) 
    {
    A.prime("Site");
    A *= B;
    A.replaceTags("2","1","Site");
    return A;
    }

bool
hasIndex(ITensor const& T,
         Index const& imatch)
    {
    return hasIndex(inds(T),imatch);
    }

Arrow
dir(ITensor const& T,
    Index const& i)
    {
    return dir(inds(T),i);
    }

bool
hasInds(ITensor const& T,
        IndexSet const& ismatch)
    {
    return hasInds(inds(T),ismatch);
    }

bool
isComplex(ITensor const& T)
    {
    return doTask(CheckComplex{},T.store());
    }

bool
isReal(ITensor const& T)
    {
    return not isComplex(T);
    }


ITensor 
random(ITensor T, const Args& args)
    {
    T.randomize(args);
    return T;
    }

Real 
sumels(ITensor const& t)
    {
    auto z = sumelsC(t);
    if(z.imag() != 0) Error("ITensor has non-zero imaginary part, use sumelsC");
    return z.real();
    }

Cplx 
sumelsC(ITensor const& t)
    {
    auto z = doTask(SumEls{inds(t)},t.store());
#ifndef USESCALE
    return z;
#else
    return t.scale().real0()*z;
#endif
    }

ostream& 
operator<<(ostream & s, ITensor const& t)
    {
    s << "ITensor ord=" << order(t) << ": "; 
    if(hasQNs(t)) 
        {
        if(t.order() > 0) s << "\n";
        for(auto& I : inds(t)) s << I << "\n";
        }
    else
        {
        s << inds(t);
        }
    if(not hasQNs(t)) s << "\n";
    if(not t.store()) 
        {
        s << "{Zero / Not yet allocated}\n";
        }
    else
        {
        //Checking whether std::ios::floatfield is set enables 
        //printing the contents of an ITensor when using the printf
        //format string %f (or another float-related format string)
        bool ff_set = (std::ios::floatfield & s.flags()) != 0;
        bool print_data = (ff_set || Global::printdat());
        doTask(PrintIT{s,t.scale(),inds(t),print_data},t.store());
        }
    return s;
    }

ITensor
matrixITensor(Matrix&& M,
              IndexSet const& is)
    {
#ifdef DEBUG
    if( order(is) != 2 )
        Error("matrixITensor(Matrix,...) constructor only accepts 2 indices");
#endif
    auto res = ITensor(is,DenseReal{std::move(M.storage())});
    M.clear();
    return res;
    }

ITensor
matrixITensor(Matrix const& M,
              IndexSet const& is)
    {
    return matrixITensor(Matrix(M),is);
    }

ITensor
matrixITensor(CMatrix&& M,
              IndexSet const& is)
    {
#ifdef DEBUG
    if( order(is) != 2 )
        Error("matrixITensor(Matrix,...) constructor only accepts 2 indices");
#endif
    bool isReal = true;
    for(auto& el : M)
    if(std::fabs(el.imag()) > 1E-14)
        {
        isReal = false;
        break;
        }
    ITensor res;
    if(isReal)
        {
        auto store = vector<Real>(M.size());
        for(auto n : range(M.size())) store[n] = M.store()[n].real();
        res = ITensor(is,DenseReal{std::move(store)});
        }
    else
        {
        res = ITensor(is,DenseCplx{std::move(M.storage())});
        }
    M.clear();
    return res;
    }

ITensor
matrixITensor(CMatrix const& M,
              IndexSet const& is)
    {
    return matrixITensor(CMatrix(M),is);
    }

ITensor
matrixITensor(Matrix && M,
              Index const& i1, Index const& i2)
  {
  return matrixITensor(M,IndexSet(i1,i2));
  }
ITensor
matrixITensor(Matrix const& M,
              Index const& i1, Index const& i2)
  {
  return matrixITensor(M,IndexSet(i1,i2));
  }
ITensor
matrixITensor(CMatrix && M,
              Index const& i1, Index const& i2)
  {
  return matrixITensor(M,IndexSet(i1,i2));
  }
ITensor
matrixITensor(CMatrix const& M,
              Index const& i1, Index const& i2)
  {
  return matrixITensor(M,IndexSet(i1,i2));
  }

std::tuple<ITensor,Index>
combiner(IndexSet const& inds, Args const& args)
    {
    auto itagset = getTagSet(args,"Tags","Link,CMB");

    if(not hasQNs(inds))
        {
        if(inds.empty()) Error("No indices passed to combiner");
        long rm = 1;
        for(const auto& i : inds) rm *= dim(i);
        //create combined index
        auto cind = Index(rm,itagset);
        //create new IndexSet with combined index in front
        auto newind = IndexSetBuilder(1+inds.order());
        newind.nextIndex(std::move(cind));
        for(auto& I : inds)
            newind.nextIndex(std::move(I));
        return std::tuple<ITensor,Index>(ITensor(newind.build(),Combiner()),cind);
        }
    else if(hasQNs(inds))
        {
        if(inds.empty()) Error("No indices passed to combiner");

        auto cdir = Out;
        if(args.defined("IndexDir"))
            {
            cdir = toArrow(args.getInt("IndexDir"));
            }
        else
            {
            //If not specified by user, make combined Index
            //point Out unless all combined indices have In arrows.
            auto allin = true;
            for(auto& i : inds) 
                if(i.dir() != In)
                    {
                    allin = false;
                    break;
                    }
            if(allin) cdir = In;
            }


        auto C = QCombiner(inds);

        //Build the combined Index,
        //saving information about
        //how we're doing it in C
        struct QNm { QN q; long m = 1; };
        auto qms = vector<QNm>{};

        //Loop over all possible QN sectors that
        //can be formed by the combined indices
        for(auto I : C.range())
            {
            QNm qm;
            //For this sector, figure out the total QN (qm.q)
            //and combined sector size (qm.m)
            for(auto j : range(inds))
                {
                qm.q += inds[j].qn(1+I[j]) * inds[j].dir() * cdir;
                qm.m *= inds[j].blocksize0(I[j]);
                }

            size_t n = 0;
            for(; n < qms.size(); ++n) if(qms[n].q==qm.q) break;

            if(n < qms.size())
                {
                C.setBlockRange(I,n,qms[n].m,qm.m);
                qms[n].m += qm.m;
                }
            else 
                {
                C.setBlockRange(I,n,0,qm.m);
                qms.push_back(qm);
                }
            }

        auto cstore = stdx::reserve_vector<QNInt>(qms.size());
        for(auto n : range(qms)) 
            {
            cstore.emplace_back(qms[n].q,qms[n].m);
            }
        auto cind = Index(std::move(cstore),cdir,itagset);

        // TODO: it seems like the next steps remove the QNs,
        // why? Does std::move(cind) destroy the QN storage?
        // To output the correct Index, we need to make a copy
        auto c = cind;

        auto newind = IndexSetBuilder(1+inds.size());
        newind.nextIndex(std::move(cind));
        for(auto& I : inds) 
            {
            newind.nextIndex(dag(I));
            }
        return std::tuple<ITensor,Index>(ITensor(newind.build(),std::move(C)),c);
        }
    return std::tuple<ITensor,Index>(ITensor(),Index());
    }

struct IsCombiner
    {
    template<typename D>
    bool 
    operator()(D const& d) { return false; }
    bool
    operator()(Combiner const& d) { return true; }
    };

//TODO: this doesn't work if the combiner is permuted
Index
combinedIndex(ITensor const& C)
    {
#ifdef DEBUG
    auto iscombiner = applyFunc(IsCombiner{},C.store());
    if(not iscombiner)
        {
        throw ITError("Called combinedIndex on ITensor that is not a combiner");
        }
#endif
    return inds(C).front();
    }

ITensor
delta(IndexSet const& is)
    {
    if(hasQNs(is))
        {
        return ITensor(std::move(is),QDiagReal(is,1.));
        }
    auto len = minDim(is);
    return ITensor(std::move(is),DiagReal(len,1.));
    }

ITensor
randomITensor(IndexSet const& inds)
    {
    return random(ITensor(inds));
    }
ITensor
randomITensorC(IndexSet const& inds)
    {
    return random(ITensor(inds),{"Complex=",true});
    }

ITensor
randomITensor(QN q, IndexSet const& is)
    {
#ifdef DEBUG
    if(not hasQNs(is)) 
        Error("Cannot use randomITensor(QN,...) to create non-QN-conserving ITensor");
#endif
    ITensor T;
    auto dat = QDenseReal{is,q};
    T = ITensor(std::move(is),std::move(dat));
    if(nnz(T) == 0) Error("Requested QN for random ITensor resulted in zero allowed blocks (QN not satisfiable by any settings of the indices)");
    T.generate(detail::quickran);
    return T;
    }

ITensor
randomITensorC(QN q, IndexSet const& is)
    {
#ifdef DEBUG
    if(not hasQNs(is)) 
        Error("Cannot use randomITensor(QN,...) to create non-QN-conserving ITensor");
#endif
    ITensor T;
    auto dat = QDenseCplx{is,q};
    T = ITensor(std::move(is),std::move(dat));
    if(nnz(T) == 0) Error("Requested QN for random ITensor resulted in zero allowed blocks (QN not satisfiable by any settings of the indices)");
    T.generate(detail::quickranCplx);
    return T;
    }

QN
div(ITensor const& T) 
    { 
    if(not hasQNs(T)) Error("div(ITensor) not defined for non QN conserving ITensor");
    if(!T) Error("div(ITensor) not defined for unallocated IQTensor");
    return doTask(CalcDiv{inds(T)},T.store());
    }

QN
flux(ITensor const& T) 
    {
    return div(T);
    }


IndexSet
moveToFront(IndexSet const& isf, IndexSet const& is)
    {
    auto rf = isf.order();
    auto r = is.order();

    if(rf >= r)
        {
        println("---------------------------------------------");
        println("Tensor indices = \n",is,"\n");
        println("---------------------------------------------");
        println("Indices provided = \n",isf," '...'\n");
        println("---------------------------------------------");
        Error(format("Wrong number of indices passed to permute (expected < %d, got %d)",r,rf));
        }

    auto iso = IndexSet(r);

    auto i = 0;
    for(auto& I : isf) 
        {
        if(!hasIndex(is,I))
            {
            println("---------------------------------------------");
            println("Tensor indices = \n",is,"\n");
            println("---------------------------------------------");
            println("Indices provided = \n",isf," '...'\n");
            println("---------------------------------------------");
            Error(format("Bad index passed to permute"));
            }
        iso[i] = I;
        i++;
        }

    auto j = rf;
    for(auto& J : is)
        {
        if(!hasIndex(isf,J))
            {
            iso[j] = J;
            j++;
            }
        }

    return iso;
    }

IndexSet 
moveToBack(IndexSet const& isb, IndexSet const& is)
    {
    auto rb = isb.order();
    auto r = is.order();

    if(rb >= r)
        {
        println("---------------------------------------------");
        println("Tensor indices = \n",is,"\n");
        println("---------------------------------------------");
        println("Indices provided = \n'...' ",isb,"\n");
        println("---------------------------------------------");
        Error(format("Wrong number of indices passed to permute (expected < %d, got %d)",r,rb));
        }

    auto iso = IndexSet(r);

    auto i = r-rb;
    for(auto& I : isb) 
        {
        if(!hasIndex(is,I))
            {
            println("---------------------------------------------");
            println("Tensor indices = \n",is,"\n");
            println("---------------------------------------------");
            println("Indices provided = \n'...' ",isb,"\n");
            println("---------------------------------------------");
            Error(format("Bad index passed to permute"));
            }
        iso[i] = I;
        i++;
        }

    auto j = 0;
    for(auto& J : is)
        {
        if(!hasIndex(isb,J))
            {
            iso[j] = J;
            j++;
            }
        }

    return iso;
    }

//
//Deprecated
//

#ifdef USESCALE
void ITensor::
scaleTo(scale_type const& newscale)
    {
    if(scale_ == newscale) return;
    if(newscale.sign() == 0) Error("Trying to scale an ITensor to a 0 scale");
    scale_ /= newscale;
    doTask(Mult<Real>{scale_.real0()},store_);
    scale_ = newscale;
    }

void ITensor::
scaleTo(Real newscale) { scaleTo(LogNum{newscale}); }
#endif

long
rank(ITensor const& T)
  {
  Global::warnDeprecated("rank(ITensor) is deprecated in favor of order(ITensor)");
  return order(T);
  }

void
randomize(ITensor & T, Args const& args)
    {
    Global::warnDeprecated("randomize(ITensor,args) is deprecated in favor of .randomize(args)");
    T.randomize(args);
    }

ITensor
matrixTensor(Matrix const& M, Index const& i1, Index const& i2)
    {
    Global::warnDeprecated("matrixTensor(Matrix,Index,Index) is deprecated in favor of matrixITensor(Matrix,Index,Index)");
    return matrixITensor(M,i1,i2);
    }

ITensor
matrixTensor(CMatrix&& M, Index const& i1, Index const& i2)
    {
    Global::warnDeprecated("matrixTensor(CMatrix,Index,Index) is deprecated in favor of matrixITensor(CMatrix,Index,Index)");
    return matrixITensor(M,i1,i2);
    }

ITensor
matrixTensor(Matrix&& M, Index const& i1, Index const& i2)
    {
    Global::warnDeprecated("matrixTensor(Matrix,Index,Index) is deprecated in favor of matrixITensor(Matrix,Index,Index)");
    return matrixITensor(M,i1,i2);
    }

ITensor
matrixTensor(CMatrix const& M, Index const& i1, Index const& i2)
    {
    Global::warnDeprecated("matrixTensor(CMatrix,Index,Index) is deprecated in favor of matrixITensor(CMatrix,Index,Index)");
    return matrixITensor(M,i1,i2);
    }


} //namespace itensor
