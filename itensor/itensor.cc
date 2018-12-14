//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/util/print_macro.h"
#include "itensor/util/range.h"
#include "itensor/util/safe_ptr.h"
#include "itensor/itensor.h"
#include "itensor/tensor/lapack_wrap.h"
#include "itensor/tensor/contract.h"

using std::array;
using std::ostream;
using std::vector;

namespace itensor {

//
// ITensor Constructors
//

    
template<>
ITensor::
ITensorT(Cplx val) 
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
    //if(val.imag() == 0)
    //    store_ = newITData<Diag<Real>>(1,val.real());
    //else
    //    store_ = newITData<Diag<Cplx>>(1,val);
    }


//template<>
//void ITensor::
//scaleOutNorm()
//    {
//    auto nrm = doTask(NormNoScale<Index>{is_},store_);
//    //If norm already 1 return so
//    //we don't have to call MultReal
//    if(fabs(nrm-1.) < 1E-12) return;
//    if(nrm == 0)
//        {
//        scale_ = LogNumber(1.);
//        return;
//        }
//    doTask(MultReal{1./nrm},store_);
//    scale_ *= nrm;
//    }

//template<>
//void ITensor::
//equalizeScales(ITensor& other)
//    {
//    if(scale_.sign() != 0)
//        {
//        other.scaleTo(scale_);
//        }
//    else //*this is equivalent to zero
//        {
//        fill(0);
//        scale_ = other.scale_;
//        }
//    }

ostream& 
operator<<(ostream & s, const ITensor& t)
    {
    s << "ITensor r=" << t.r() << ": " << t.inds() << "\n";
    if(!t.store()) 
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
        doTask(PrintIT<Index>{s,t.scale(),t.inds(),print_data},t.store());
        }
    return s;
    }

ITensor
matrixTensor(Matrix&& M, Index const& i1, Index const& i2)
    {
    auto res = ITensor({i1,i2},DenseReal{std::move(M.storage())});
    M.clear();
    return res;
    }

ITensor
matrixTensor(Matrix const& M, const Index& i1, const Index& i2)
    {
    return matrixTensor(Matrix(M),i1,i2);
    }

ITensor
matrixTensor(CMatrix&& M, Index const& i1, Index const& i2)
    {
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
        res = ITensor({i1,i2},DenseReal{std::move(store)});
        }
    else
        {
        res = ITensor({i1,i2},DenseCplx{std::move(M.storage())});
        }
    M.clear();
    return res;
    }

ITensor
matrixTensor(CMatrix const& M, const Index& i1, const Index& i2)
    {
    return matrixTensor(CMatrix(M),i1,i2);
    }


ITensor
combiner(IndexSet const& inds, Args const& args)
    {
    if(inds.empty()) Error("No indices passed to combiner");
    long rm = 1;
    for(const auto& i : inds) rm *= i.m();
    //create combined index
    auto cname = args.getString("IndexName","cmb");
    auto itype = getIndexType(args,"IndexType",Link);
    auto cind = Index(cname,rm,itype);
    //create new IndexSet with combined index in front
    auto newind = IndexSetBuilder(1+inds.r());
    newind.nextIndex(std::move(cind));
    for(auto& I : inds)
        {
        newind.nextIndex(std::move(I));
        }
    return ITensor(newind.build(),Combiner{});
    }

ITensor
combiner(std::vector<Index> const& inds, Args const& args)
    {
    return combiner(IndexSet(inds),args);
    }

ITensor
combiner(std::initializer_list<Index> inds, Args const& args)
    {
    return combiner(IndexSet(inds),args);
    }

struct IsCombiner
    {
    template<typename D>
    bool 
    operator()(D const& d) { return false; }
    bool
    operator()(Combiner const& d) { return true; }
    };

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
    return C.inds().front();
    }

} //namespace itensor
