//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/util/count.h"
#include "itensor/tensor/sliceten.h"
#include "itensor/itdata/iqtcombiner.h"
#include "itensor/itdata/itdata.h"
#include "itensor/iqindex.h"

using std::vector;

namespace itensor {

Cplx
doTask(GetElt<IQIndex> const& g, IQTCombiner const& c)
    {
    if(g.inds.size()!=0) Error("GetElt not defined for non-scalar IQTCombiner storage");
    return Cplx(1.,0.);
    }

void inline
doTask(Write& W, IQTCombiner const& d)
    { W.writeType(StorageType::IQTCombiner,d);}

void
permuteIQ(const Permutation& P,
          const IQIndexSet& Ais,
          const IQTReal& dA,
          IQIndexSet& Bis,
          IQTReal& dB)
    {
#ifdef DEBUG
    if(isTrivial(P)) Error("Calling permuteIQ for trivial Permutation");
#endif
    auto r = Ais.r();
    auto bind = IQIndexSetBuilder(r);
    for(auto i : count(r))
        {
        bind.setExtent(P.dest(i),Ais[i]);
        }
    Bis = bind.build();
    dB = IQTReal(Bis,doTask(CalcDiv{Ais},dA));

    Label Ablock(r,-1),
          Bblock(r,-1);
    Range Arange,
          Brange;
    for(auto aio : dA.offsets)
        {
        //Compute bi, new block index of blk
        computeBlockInd(aio.block,Ais,Ablock);
        for(auto j : index(Ablock)) 
            Bblock.at(P.dest(j)) = Ablock[j];
        Arange.init(make_indexdim(Ais,Ablock));
        Brange.init(make_indexdim(Bis,Bblock));

        auto aref = makeTenRef(dA.data(),aio.offset,dA.size(),&Arange);

        auto bblock = getBlock(dB,Bis,Bblock);
        auto bref = makeTenRef(bblock,&Brange);

        do_permute(aref,P,bref);
        }
    }

IQIndexSet
replaceInd(IQIndexSet const& is,
           long              loc,
           IQIndex    const& replacement)
    {
    auto newind = IQIndexSetBuilder(is.r());
    long i = 0;
    for(long j = 0; j < loc; ++j)
        newind.setExtent(i++,is[j]);
    newind.setExtent(i++,replacement);
    for(decltype(is.r()) j = loc+1; j < is.r(); ++j)
        newind.setExtent(i++,is[j]);
    return newind.build();
    }

void
combine(IQTReal     const& d,
        IQTCombiner const& C,
        IQIndexSet  const& dis,
        IQIndexSet  const& Cis,
        IQIndexSet       & Nis,
        ManageStore      & m)
    {
#ifdef DEBUG
    for(auto i : count(1,Cis.r()))
        {
        auto jc = findindex(dis,Cis[i]);
        if(jc == -1)
            {
            printfln("Indices of tensor = \n%s\n------",dis);
            printfln("Extra index = \n%s",Cis[i]);
            Error("Combiner has extra index not found on other tensor");
            }
        }
#endif

    auto& cind = Cis[0];
    auto cr = Cis.r();
    auto dr = dis.r();
    auto ncomb = cr-1;
    auto nr = dr-ncomb+1;

    auto dtoC = Label(dr,0);
    auto to_comb = Label(ncomb,-1);
    for(auto i : count(dr)) 
        {
        auto jc = findindex(Cis,dis[i]);
        if(jc >= 0)  //combined index
            {
            dtoC[i] = jc;
            to_comb[jc-1] = i;
            }
        }

    //Create new IQIndexSet
    auto newind = IQIndexSetBuilder(nr);
    newind.nextExtent(cind);
    for(auto i : count(dr)) if(!dtoC[i]) newind.nextExtent(dis[i]);
    Nis = newind.build();

    //Allocate new data
    auto& nd = *m.makeNewData<IQTReal>(Nis,doTask(CalcDiv{dis},d));

    auto drange = Range(dr),
         nrange = Range(nr);
    auto dblock = Label(dr),
         nblock = Label(nr),
         cblock = Label(ncomb);
    size_t start = 0,
           end   = 0;
    for(auto io : d.offsets) //loop over non-zero blocks
        {
        computeBlockInd(io.block,dis,dblock);
        drange.init(make_indexdim(dis,dblock));
        auto dref = makeTenRef(d.data(),io.offset,d.size(),&drange);

        //QN of combined indices for this block
        //total dimension of combined indices for this block
        size_t nu = 1;
        for(auto i : count(dr)) 
            {
            if(!dtoC[i]) //uncombined
                nblock[nu++] = dblock[i];
            else            //combined
                cblock[dtoC[i]-1] = dblock[i];
            }

        std::tie(nblock[0],start,end) = C.getBlockRange(cblock);

        //Get full block of new storage
        nrange.init(make_indexdim(Nis,nblock));
        auto nref = makeTenRef(getBlock(nd,Nis,nblock),&nrange);
        //Do tensor slicing to get subblock where data will go
        auto nsub = subIndex(nref,0,start,end);

        //Call groupInds on current block permutes combined
        //inds to front, then groups them into a single ind
        nsub &= groupInds(dref,to_comb);
        }
    }

void
combReplaceIndex(IQIndexSet  const& dis,
                 IQIndexSet  const& Cis,
                 IQIndexSet       & Nis)
    {
    auto jc = findindex(dis,Cis[0]);

    if(jc >= 0) //uncombining
        {
        //Has Cis[0], replace with Cis[1]
        Nis = replaceInd(dis,jc,Cis[1]);
        }
    else //combining
        {
        //Has Cis[1], replace with cind
        auto ju = findindex(dis,Cis[1]);
        if(ju < 0)
            {
            println("IQIndexSet of regular IQTensor =\n",dis);
            println("IQIndexSet of combiner/delta =\n",Cis);
            println("Missing IQIndex: ",Cis[1]);
            Error("IQCombiner: missing IQIndex");
            }
        Nis = replaceInd(dis,ju,Cis[0]);
        }
    }

void
uncombine(IQTReal     const& d,
          IQIndexSet  const& dis,
          IQIndexSet  const& Cis,
          IQIndexSet       & Nis,
          ManageStore      & m,
          bool              own_data)
    {
    //cind is special "combined index"
    auto const& cind = Cis[0];
    auto jc = findindex(dis,cind);

    IQTReal* pd = nullptr;
    //Call this if necessary to modify the data
    auto copyDataSetpd = [&]()
        {
        if(pd) return;
        if(own_data) pd = m.modifyData(d);
        else         pd = m.makeNewData<IQTReal>(d);
        };

    Permutation P;

    if(jc != 0) //cind not at front
        {
        P = Permutation(dis.r());
        P.setFromTo(jc,0);
        long ni = 1;
        for(auto j : count(dis.r()))
            if(j != jc) P.setFromTo(j,ni++);
        jc = 0;
        }

    if(P) 
        {
        copyDataSetpd(); //sets pd
        permuteIQ(P,dis,d,Nis,*pd);
        }
    //Pis means 'permuted' index set
    auto& Pis = (P ? Nis : dis);

    auto newr = Pis.r()+Cis.r()-2;
    auto offset = Cis.r()-1;
    auto newind = IQIndexSetBuilder(newr);
    for(auto j : count(offset))
        newind.setExtent(j,Cis[1+j]);
    for(auto j : count(Pis.r()-1))
        newind.setExtent(offset+j,Pis[1+j]);
    Nis = newind.build();

    copyDataSetpd();
    auto div = doTask(CalcDiv{dis},d);
    pd->updateOffsets(Nis,div);
    }

void
doTask(Contract<IQIndex>      & C,
       IQTReal           const& d,
       IQTCombiner        const& cmb,
       ManageStore            & m)
    {
    if(C.Ris.r()==2)
        combReplaceIndex(C.Lis,C.Ris,C.Nis);
    else if(hasindex(C.Lis,C.Ris[0]))
        uncombine(d,C.Lis,C.Ris,C.Nis,m,true);
    else
        combine(d,cmb,C.Lis,C.Ris,C.Nis,m);
    }

void
doTask(Contract<IQIndex>      & C,
       IQTCombiner       const& cmb,
       IQTReal           const& d,
       ManageStore            & m)
    { 
    if(C.Lis.r()==2)
        {
        combReplaceIndex(C.Ris,C.Lis,C.Nis);
        m.assignPointerRtoL();
        }
    else if(hasindex(C.Ris,C.Lis[0]))
        {
        uncombine(d,C.Ris,C.Lis,C.Nis,m,false);
        }
    else
        {
        combine(d,cmb,C.Ris,C.Lis,C.Nis,m);
        }
    }


} //namespace itensor
