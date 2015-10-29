//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/util/count.h"
#include "itensor/itdata/itcombiner.h"
#include "itensor/itdata/itdata.h"
//#include "itensor/itdata/itcplx.h"
#include "itensor/tensor/contract.h"
#include "itensor/tensor/sliceten.h"
#include "itensor/iqindex.h"

using std::vector;

namespace itensor {

Cplx
doTask(GetElt<Index> const& g, ITCombiner const& c)
    {
    if(g.inds.size()!=0) Error("GetElt not defined for non-scalar ITCombiner storage");
    return Cplx(1.,0.);
    }

Real
doTask(NormNoScale, ITCombiner const& d) { return 0; }

void
doTask(Conj,ITCombiner const& d) { }

template<typename T>
void
permuteStore(Dense<T>   const& d,
             IndexSet    const& dis,
             Permutation const& P,
             ManageStore      & m)
    {
    auto tfrom = makeTenRef(d.data(),d.size(),&dis);
    auto to = Ten<Range,T>(permute(tfrom,P));
    m.makeNewData<Dense<T>>(move(to.storage()));
    }

//void
//permuteStore(ITCplx      const& d,
//             IndexSet    const& dis,
//             Permutation const& P,
//             ManageStore      & m)
//    {
//    auto *nd = m.makeNewData<ITCplx>(d.size());
//    auto csize = d.csize();
//    auto fromre = makeTenRef(d.rstart(),csize,&dis);
//    auto tore = makeTenRef(nd->rstart(),csize,&dis);
//    auto fromim = makeTenRef(d.istart(),csize,&dis);
//    auto toim = makeTenRef(nd->istart(),csize,&dis);
//    tore &= permute(fromre,P);
//    toim &= permute(fromim,P);
//    }

template<typename Storage>
void
combine(Storage  const& d,
        IndexSet const& dis,
        IndexSet const& Cis,
        IndexSet      & Nis,
        ManageStore   & m)
    {
    //TODO: try to make use of Lind,Rind label vectors
    //      to simplify combine logic
    auto const& cind = Cis[0];
    auto jc = findindex(dis,cind);
    if(jc >= 0) //has cind, uncombining
        {
        //dis has cind, replace with other inds
        auto newind = IndexSetBuilder(dis.r()+Cis.r()-2);
        long i = 0;
        for(auto j : count(dis.r()))
            if(j == jc)
                {
                for(size_t k = 1; k < Cis.size(); ++k)
                    newind.setIndex(i++,Cis[k]);
                }
            else
                {
                newind.setIndex(i++,dis[j]);
                }
        Nis = newind.build();
        }
    else //combining
        {
        //dis doesn't have cind, replace
        //Cis[1], Cis[2], ... with cind
        //may need to permute
        auto J1 = findindex(dis,Cis[1]);
        if(J1 < 0) 
            {
            println("IndexSet of dense tensor = \n",dis);
            println("IndexSet of combiner/delta = \n",Cis);
            Error("No contracted indices in combiner-tensor product");
            }
        //Check if Cis[1],Cis[2],... are grouped together (contiguous)
        //and in same order as on combiner
        bool contig_sameord = true;
        decltype(Cis.r()) c = 2;
        for(auto j = J1+1; c < Cis.r() && j < dis.r(); ++j,++c)
            if(dis[j] != Cis[c])
                {
                contig_sameord = false;
                break;
                }
        if(c != Cis.r()) contig_sameord = false;

        if(contig_sameord)
            {
            auto newind = IndexSetBuilder(dis.r()+2-Cis.r());
            long i = 0;
            for(auto j : count(J1))
                newind.setIndex(i++,dis[j]);
            newind.setIndex(i++,cind);
            for(auto j : count(J1+Cis.r()-1, dis.r()))
                newind.setIndex(i++,dis[j]);
            Nis = newind.build();
            }
        else
            {
            auto P = Permutation(dis.r());
            //Set P destination values to -1 to mark
            //indices that need to be assigned destinations:
            for(auto i : index(P)) P.setFromTo(i,-1);

            //permute combined indices to the front, in same
            //order as in Cis:
            long ni = 0;
            for(auto c : count(1,Cis.r()))
                {
                auto j = findindex(dis,Cis[c]);
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
            auto newind = IndexSetBuilder(dis.r()+2-Cis.r());
            long i = 0;
            newind.setIndex(i++,cind);
            for(auto j : count(dis.r()))
                if(P.dest(j) == -1) 
                    {
                    P.setFromTo(j,ni++);
                    newind.setIndex(i++,dis[j]);
                    }
            Nis = newind.build();
            permuteStore(d,dis,P,m);
            }
        }
    }

void
doTask(Contract<Index>& C,
       DenseReal const& d,
       ITCombiner const& cmb,
       ManageStore& m)
    {
    combine(d,C.Lis,C.Ris,C.Nis,m);
    }
void
doTask(Contract<Index>& C,
       ITCombiner const& cmb,
       DenseReal const& d,
       ManageStore& m)
    { 
    combine(d,C.Ris,C.Lis,C.Nis,m);
    if(!m.newData()) m.assignPointerRtoL();
    }

void
doTask(Contract<Index> & C,
       DenseCplx  const& d,
       ITCombiner const& cmb,
       ManageStore     & m)
    {
    combine(d,C.Lis,C.Ris,C.Nis,m);
    }
void
doTask(Contract<Index> & C,
       ITCombiner const& cmb,
       DenseCplx  const& d,
       ManageStore     & m)
    { 
    combine(d,C.Ris,C.Lis,C.Nis,m);
    if(!m.newData()) m.assignPointerRtoL();
    }

bool
doTask(CheckComplex, ITCombiner const& d) { return false; }

void
doTask(PrintIT<Index>& P, ITCombiner const& d)
    {
    P.printInfo(d,"Combiner");
    }

void
doTask(PrintIT<IQIndex>& P, ITCombiner const& d)
    {
    P.s << "ITCombiner";
    }

void
doTask(Write& W, const ITCombiner& d) 
    { 
    W.writeType(StorageType::ITCombiner,d); 
    }

QN 
doTask(CalcDiv const& C, 
       ITCombiner const& d)
    {
    return QN{};
    }

} //namespace itensor
