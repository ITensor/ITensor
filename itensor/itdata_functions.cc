//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itdata_functions.h"
#include "detail/gcounter.h"
#include "lapack_wrap.h"

using std::vector;

namespace itensor {

void Contract::
computeNis(SortOption sort)
    {
    long ncont = 0;
    for(const auto& i : Lind_) if(i < 0) ++ncont;
    long nuniq = Lis_.r()+Ris_.r()-2*ncont;
    vector<Index> newind(nuniq);
    long nn = 0;
    for(int j = 0; j < Lis_.r(); ++j)
        {
        if(Lind_[j] > 0) newind[nn++] = Lis_[j];
        }
    for(int j = 0; j < Ris_.r(); ++j)
        {
        if(Rind_[j] > 0) newind[nn++] = Ris_[j];
        }
    if(sort == Sort)
        {
        auto comp = [](const Index& i1, const Index& i2) { return i1 > i2; };
        std::sort(newind.begin(),newind.end(),comp);
        }
    Nis_ = IndexSet(std::move(newind));
    }

ITResult Contract::
operator()(const ITDense<Real>& a1,
           const ITDense<Real>& a2)
    {
    computeNis(Sort);
    
    Label Nind(Nis_.r());
    for(size_t i = 0; i < Nis_.r(); ++i)
        {
        auto j = findindex(Lis_,Nis_[i]);
        if(j >= 0)
            {
            Nind[i] = Lind_[j];
            }
        else
            {
            j = findindex(Ris_,Nis_[i]);
            Nind[i] = Rind_[j];
            }
        }

    //PRI(Lind);
    //PRI(Rind);
    //PRI(Nind);

    auto res = make_newdata<ITDense<Real>>(area(Nis_),0.);
    auto t1 = make_tensorref(a1.data.data(),Lis_),
         t2 = make_tensorref(a2.data.data(),Ris_),
         tr = make_tensorref(res->data.data(),Nis_);
    contractloop(t1,Lind_,t2,Rind_,tr,Nind);
    scalefac_ = 0;
    for(auto elt : res->data)
        {
        scalefac_ += elt*elt;
        }
    scalefac_ = std::sqrt(scalefac_);
    for(auto& elt : res->data)
        {
        elt /= scalefac_;
        }
    return std::move(res);
    }

ITResult Contract::
diagDense(const ITDiag<Real>& d,
          const IndexSet& dis,
          const Label& dind,
          const ITDense<Real>& t,
          const IndexSet& tis,
          const Label& tind)
    {
    computeNis(Sort);

    long t_cstride = 0; //total t-stride of contracted inds of t
    size_t ntu = 0; //number uncontracted inds of t
    assert(int(tind.size()) == tis.size());
    for(size_t j = 0; j < tind.size(); ++j)
        {
        //if index j is contracted, add its stride to t_cstride:
        if(tind[j] < 0) t_cstride += tis.stride(j);
        else            ++ntu;
        }

    long d_ustride = 0; //total result-stride of uncontracted inds of d
    for(size_t i = 0; i < Nis_.r(); ++i)
        {
        auto j = findindex(dis,Nis_[i]);
        if(j >= 0) d_ustride += Nis_.stride(i);
        }

    if(ntu > 0)
        {
        vector<long> tstride(ntu,0),
                     rstride(ntu,0);
        detail::GCounter C(0,ntu,0);
        size_t n = 0;
        for(size_t j = 0; j < tind.size(); ++j)
            {
            if(tind[j] > 0)
                {
#ifdef DEBUG
                if(n >= ntu) Error("n out of range");
#endif
                C.setInd(n,0,tis.dim(j)-1);
                tstride.at(n) = tis.stride(j);
                auto k = findindex(Nis_,tis[j]);
#ifdef DEBUG
                if(k < 0) Error("Index not found");
#endif
                rstride.at(n) = Nis_.stride(k);
                ++n;
                }
            }
        auto res = make_newdata<ITDense<Real>>(area(Nis_),0.);
        auto *pr = res->data.data();
        const auto *pt = t.data.data();

        if(d.allSame())
            {
            auto size = minM(dis);
            for(;C.notDone();++C)
                {
                size_t roffset = 0,
                       toffset = 0;
                for(size_t i = 0; i < ntu; ++i)
                    {
                    auto ii = C.i.fast(i);
                    toffset += ii*tstride[i];
                    roffset += ii*rstride[i];
                    }
                for(long J = 0; J < size; ++J)
                    {
                    pr[J*d_ustride+roffset] += d.val*pt[J*t_cstride+toffset];
                    }
                }
            }
        else
            {
            auto* pd = d.data.data();
            auto Md = d.data.size();
            for(;C.notDone();++C)
                {
                size_t roffset = 0,
                       toffset = 0;
                for(size_t i = 0; i < ntu; ++i)
                    {
                    auto ii = C.i.fast(i);
                    toffset += ii*tstride[i];
                    roffset += ii*rstride[i];
                    }
                for(size_t J = 0; J < Md; ++J)
                    {
                    pr[J*d_ustride+roffset] += pd[J]*pt[J*t_cstride+toffset];
                    }
                }
            }
        return std::move(res);
        }
    else
        {
        //all of t's indices contracted with d
        //result will be diagonal:
        // o scalar if all of d's inds contracted also
        // o dot product of d's data and t's diagonal otherwise
        //auto res = make_newdata<ITDiag<Real>>(???,0.);
        Error("Case not handled");
        }
    Error("Case not handled");
    return ITResult();
    }

NewData Contract::
combine(const ITDense<Real>& d,
        const IndexSet& dis,
        const IndexSet& Cis)
    {
    //TODO: try to make use of Lind,Rind label vectors
    //      to simplify combine logic
    const auto& cind = Cis[0];
    int jc = findindex(dis,cind);
    if(jc >= 0) //has cind
        {
        //dis has cind, replace with other inds
        vector<Index> newind;
        newind.reserve(dis.r()+Cis.r()-2);
        for(int j = 0; j < dis.r(); ++j)
            if(j == jc)
                {
                for(int k = 1; k < Cis.size(); ++k)
                    newind.push_back(Cis[k]);
                }
            else
                {
                newind.push_back(dis[j]);
                }
        Nis_ = IndexSet(std::move(newind));
        return NewData();
        }
    else
        {
        //dis doesn't have cind, replace
        //Cis[1], Cis[2], ... with cind
        //may need to reshape
        int J1 = findindex(dis,Cis[1]);
        if(J1 < 0) 
            {
            println("IndexSet of dense tensor = \n",dis);
            println("IndexSet of combiner/delta = \n",Cis);
            Error("No contracted indices in combiner-tensor product");
            }
        //Check if Cis[1],Cis[2],... are grouped together (contiguous)
        bool contig = true;
        for(int j = J1+1, c = 2; c < Cis.r() && j < dis.r(); ++j,++c)
            if(dis[j] != Cis[c])
                {
                contig = false;
                break;
                }
        if(contig)
            {
            vector<Index> newind;
            newind.reserve(dis.r()-Cis.r()+1);
            for(int j = 0; j < J1; ++j) 
                newind.push_back(dis[j]);
            newind.push_back(cind);
            for(int j = J1+Cis.r()-1; j < dis.r(); ++j) 
                newind.push_back(dis[j]);
            Nis_ = IndexSet(std::move(newind));
            return NewData();
            }
        else
            {
            Permutation P(dis.r());
            //Set P destination values to -1 to mark
            //indices that need to be assigned destinations:
            for(int i = 0; i < P.size(); ++i) P.setFromTo(i,-1);

            //permute combined indices to the front, in same
            //order as in Cis:
            long ni = 0;
            for(int c = 1; c < Cis.r(); ++c)
                {
                int j = findindex(dis,Cis[c]);
                if(j < 0) 
                    {
                    println("IndexSet of dense tensor =\n  ",dis);
                    println("IndexSet of combiner/delta =\n  ",Cis);
                    println("Missing index: ",Cis[c]);
                    Error("Combiner: missing index");
                    }
                P.setFromTo(j,ni++);
                }
            //permute uncombined indices to back, keeping relative order:
            vector<Index> newind;
            vector<long> pdims(dis.r(),-1);
            newind.reserve(dis.r()-Cis.r()+1);
            newind.push_back(cind);
            for(int j = 0; j < dis.r(); ++j)
                {
                if(P.dest(j) == -1) 
                    {
                    P.setFromTo(j,ni++);
                    newind.push_back(dis[j]);
                    }
                pdims[j] = dis[P.dest(j)].m();
                }
            Range rr(pdims);
            Nis_ = IndexSet(std::move(newind));
            auto res = make_newdata<ITDense<Real>>(area(Nis_));
            auto td = make_tensorref(d.data.data(),dis);
            auto tr = make_tensorref(res->data.data(),rr);
            reshape(td,P,tr);
            return std::move(res);
            }
        }
    return NewData();
    }


ITResult FillReal::
operator()(ITDense<Real>& d) const
    {
    std::fill(d.data.begin(),d.data.end(),r_);
    return ITResult();
    }

ITResult FillReal::
operator()(const ITDense<Complex>& d) const
    {
    auto nd = make_newdata<ITDense<Real>>(d.data.size());
    operator()(*nd);
    return std::move(nd);
    }

ITResult FillReal::
operator()(ITDiag<Real>& d) const
    {
    return make_newdata<ITDiag<Real>>(r_);
    }

ITResult FillReal::
operator()(const ITDiag<Complex>& d) const
    {
    return make_newdata<ITDiag<Real>>(r_);
    }

ITResult FillCplx::
operator()(const ITDense<Real>& d) const
    {
    return make_newdata<ITDense<Complex>>(d.data.size(),z_);
    }

ITResult FillCplx::
operator()(ITDense<Complex>& d) const
    {
    std::fill(d.data.begin(),d.data.end(),z_);
    return ITResult();
    }

ITResult MultComplex::
operator()(const ITDense<Real>& d) const
    {
    auto nd = make_newdata<ITDense<Complex>>(d.data.begin(),d.data.end());
    operator()(*nd);
    return std::move(nd);
    }

ITResult MultComplex::
operator()(ITDense<Complex>& d) const
    {
    //TODO: use BLAS algorithm
    for(auto& elt : d.data)
        elt *= z_;
    return ITResult();
    }

ITResult MultReal::
operator()(ITDense<Real>& d) const
    {
    //TODO: use BLAS algorithm
    for(auto& elt : d.data)
        elt *= r_;
    return ITResult();
    }

ITResult MultReal::
operator()(ITDense<Complex>& d) const
    {
    //TODO: use BLAS algorithm
    for(auto& elt : d.data)
        elt *= r_;
    return ITResult();
    }

void
plusEqData(Real fac, Real *d1, const Real *d2, LAPACK_INT size)
    {
    LAPACK_INT inc = 1;
    daxpy_wrapper(&size,&fac,d2,&inc,d1,&inc);
    }

ITResult PlusEQ::
operator()(ITDense<Real>& a1,
           const ITDense<Real>& a2)
    {
#ifdef DEBUG
    if(a1.data.size() != a2.data.size()) Error("Mismatched sizes in plusEq");
#endif
    if(permute_)
        {
        auto ref1 = tensorref<Real,IndexSet>(a1.data.data(),*is1_),
             ref2 = tensorref<Real,IndexSet>(a2.data.data(),*is2_);
        auto f = fac_;
        auto add = [f](Real& r1, Real r2) { r1 += f*r2; };
        reshape(ref2,*P_,ref1,add);
        }
    else
        {
        plusEqData(fac_,a1.data.data(),a2.data.data(),a1.data.size());
        }
    return ITResult();
    }

ITResult PlusEQ::
operator()(ITDiag<Real>& a1,
           const ITDiag<Real>& a2)
    {
#ifdef DEBUG
    if(a1.data.size() != a2.data.size()) Error("Mismatched sizes in plusEq");
#endif
    if(a1.allSame() || a2.allSame()) Error("ITDiag plusEq allSame case not implemented");
    plusEqData(fac_,a1.data.data(),a2.data.data(),a1.data.size());
    return ITResult();
    }

void
printVal(std::ostream& s,
         Real val)
    {
    if(std::fabs(val) > 1E-10)
        s << val << "\n";
    else
        s << format("%.8E\n",val);
    }

void
printVal(std::ostream& s,
         const Complex& val)
    {
    if(std::norm(val) > 1E-10)
        {
        auto sgn = (val.imag() < 0 ? '-' : '+');
        s << val.real() << sgn << std::fabs(val.imag()) << "i\n";
        }
    else
        {
        s << format("%.8E\n",val);
        }
    }

template<typename T>
ITResult PrintIT::
operator()(const ITDense<T>& d) const
    {
    s_ << "}\n";
    Real scalefac = 1.0;
    if(!x_.isTooBigForReal()) scalefac = x_.real0();
    else s_ << "  (omitting too large scale factor)\n";

    auto rank = is_.r();
    if(rank == 0) return ITResult();

    auto gc = detail::GCounter(0,rank-1,0);
    for(int i = 0; i < rank; ++i)
        gc.setInd(i,0,is_.dim(i)-1);

    for(; gc.notDone(); ++gc)
        {
        auto val = scalefac*d.data[ind(is_,gc.i)];
        if(std::norm(val) > Global::printScale())
            {
            s_ << "  (";
            for(auto ii = gc.i.mini(); ii <= gc.i.maxi(); ++ii)
                {
                s_ << (1+gc.i(ii));
                if(ii < gc.i.maxi()) s_ << ",";
                }
            s_ << ") ";

            printVal(s_,val);
            }
        }
    return ITResult();
    }
template ITResult PrintIT::operator()(const ITDense<Real>& d) const;
template ITResult PrintIT::operator()(const ITDense<Complex>& d) const;

template<typename T>
ITResult PrintIT::
operator()(const ITDiag<T>& d) const
    {
    auto allsame = d.allSame();
    s_ << " Diag" << (allsame ? "(all same)" : "") << "}\n";
    Real scalefac = 1.0;
    if(!x_.isTooBigForReal()) scalefac = x_.real0();
    else s_ << "  (omitting too large scale factor)\n";

    auto size = minM(is_);
    for(size_t i = 0; i < size; ++i)
        {
        auto val = scalefac*(allsame ? d.val : d.data[i]);
        if(std::norm(val) > Global::printScale())
            {
            s_ << "  (";
            for(int j = 1; j < is_.size(); ++j)
                {
                s_ << (1+i) << ",";
                }
            s_ << (1+i) << ") ";
            printVal(s_,val);
            }
        }
    return ITResult();
    }
template ITResult PrintIT::operator()(const ITDiag<Real>& d) const;
template ITResult PrintIT::operator()(const ITDiag<Complex>& d) const;

}; //namespace itensor
