//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/matrix/lapack_wrap.h"
#include "itensor/itensor.h"
#include "itensor/tensor/contract.h"
#include "itensor/util/count.h"
#include "itensor/util/safe_ptr.h"

using std::array;
using std::ostream;
using std::vector;

namespace itensor {

//
// ITensor Constructors
//


template<>
ITensor::
ITensorT(const Index& i1) 
  : is_(i1),
    scale_(1.)
	{ }


template<>
ITensor::
ITensorT(const Index& i1,const Index& i2) 
  : is_(i1,i2),
    scale_(1.)
	{ }
    
template<>
ITensor::
ITensorT(Cplx val) 
  : scale_(1.)
    { 
    if(val.imag() == 0)
        store_ = newITData<ITReal>(1,val.real());
    else
        store_ = newITData<ITCplx>(1,val);
    //if(val.imag() == 0)
    //    store_ = newITData<ITDiag<Real>>(1,val.real());
    //else
    //    store_ = newITData<ITDiag<Cplx>>(1,val);
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
matrixTensor(Mat&& M, const Index& i1, const Index& i2)
    {
    auto res = ITensor({i1,i2},ITReal{std::move(M.store())});
    M.clear();
    return res;
    }


ITensor
combiner(std::vector<Index> inds, const Args& args)
    {
    if(inds.empty()) Error("No indices passed to combiner");
    long rm = 1;
    for(const auto& i : inds)
        {
        rm *= i.m();
        }
    //increase size by 1
    inds.push_back(Index());
    //shuffle contents to the end
    for(size_t j = inds.size()-1; j > 0; --j)
        {
        inds[j] = inds[j-1];
        }
    //create combined index
    auto cname = args.getString("IndexName","cmb");
    auto itype = getIndexType(args,"IndexType",Link);
    inds.front() = Index(cname,rm,itype);
    return ITensor(IndexSet(std::move(inds)),ITCombiner());
    }

ITensor
deltaTensor(const Index& i1, const Index& i2)
    {
#ifdef DEBUG
    if(i1.m() != i2.m()) Error("delta: indices must have same dimension");
#endif
    return ITensor({i1,i2},ITCombiner());
    }


} //namespace itensor
