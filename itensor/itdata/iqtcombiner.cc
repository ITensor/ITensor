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

void 
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
        bind.setIndex(P.dest(i),Ais[i]);
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

        auto bblock = getBlock(dB,Bis,Bblock);
        auto bref = TensorRef(bblock,&Brange);

        auto aref = makeTenRef(dA.data(),aio.offset,dA.size(),&Arange);

        bref += permute(aref,P);
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
        newind.setIndex(i++,is[j]);
    newind.setIndex(i++,replacement);
    for(decltype(is.r()) j = loc+1; j < is.r(); ++j)
        newind.setIndex(i++,is[j]);
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

    using size_type = decltype(rank(dis));
    auto dr = rank(dis);
    auto ncomb = rank(Cis)-1;
    auto nr = dr-ncomb+1;

    auto dperm = Label(dr,-1);
    auto uncomb_dest = ncomb;
    for(auto i : count(dr)) 
        {
        auto jc = findindex(Cis,dis[i]);
        if(jc >= 0) dperm[i] = jc-1;
        else        dperm[i] = uncomb_dest++;
        }

    auto combined = [&dperm,ncomb](size_type i) { return dperm[i] < long(ncomb); };

    //Create new IQIndexSet
    auto newind = IQIndexSetBuilder(nr);
    newind.nextIndex(Cis[0]);
    for(auto i : count(dr)) if(!combined(i)) newind.nextIndex(dis[i]);
    Nis = newind.build();

    //Allocate new data
    auto& nd = *m.makeNewData<IQTReal>(Nis,doTask(CalcDiv{dis},d));

    auto drange = Range(dr), //block range of current storage
         nrange = Range(nr); //block range of new storage
    auto dblock = Label(dr), //block index of current storage
         nblock = Label(nr), //block index of new storage
         cblock = Label(ncomb); //corresponding subblock of combiner
    size_t start = 0, //offsets within sector of combined
           end   = 0; //IQIndex where block will go
    for(auto io : d.offsets) //loop over non-zero blocks
        {
        //Figure out this block's "block index"
        computeBlockInd(io.block,dis,dblock);

        //Make TensorRef for this block of d
        drange.init(make_indexdim(dis,dblock));
        auto dref = makeTenRef(d.data(),io.offset,d.size(),&drange);

        //Permute combined indices to front, then
        //group the first ncomb indices into one
        auto Pdref = Tensor{permute(dref,dperm)};
        auto gPdref = groupInds(Pdref,0,ncomb);

        //Figure out "block index" where this block will
        //go in new storage (nblock) and which sector of
        //combined indices maps to new combined index (cblock)
        size_t nu = 1;
        for(auto i : count(dr)) 
            {
            if(combined(i)) cblock[dperm[i]] = dblock[i];
            else            nblock[nu++] = dblock[i];
            }

        //Use cblock to recover info about structure of combined IQIndex,
        //which sector to map to, where this subsector starts, and ends
        std::tie(nblock[0],start,end) = C.getBlockRange(cblock);

        //Get full block of new storage
        nrange.init(make_indexdim(Nis,nblock));
        auto nref = TensorRef(getBlock(nd,Nis,nblock),&nrange);

        //Slice this new-storage block to get subblock where data will go
        auto nsub = subIndex(nref,0,start,end);
        nsub &= gPdref;
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
          IQTCombiner const& C,
          IQIndexSet  const& dis,
          IQIndexSet  const& Cis,
          IQIndexSet       & Nis,
          ManageStore      & m,
          bool              own_data)
    {
    using size_type = decltype(rank(dis));
    auto& cind = Cis[0];
    auto dr = rank(dis);
    auto cr = rank(Cis);
    auto ncomb = cr-1;
    auto nr = dr-1+ncomb;

    decltype(dr) jc = 0;
    auto newind = IQIndexSetBuilder(nr);
    for(auto n : count(dr)) 
        {
        if(dis[n] == cind)
            {
            jc = n;
            for(auto n : count(1,cr)) newind.nextIndex(Cis[n]);
            }
        else
            {
            newind.nextIndex(dis[n]);
            }
        }
    Nis = newind.build();

    //Allocate new data
    auto& nd = *m.makeNewData<IQTReal>(Nis,doTask(CalcDiv{dis},d));

    auto drange = Range(dr), //block range of current storage
         nrange = Range(nr); //block range of new storage
    auto dblock = Label(dr), //block index of current storage
         nblock = Label(nr); //block index of new storage
    for(auto io : d.offsets) //loop over non-zero blocks
        {
        //Figure out this block's "block index"
        computeBlockInd(io.block,dis,dblock);

        //Make TensorRef for this block of d
        drange.init(make_indexdim(dis,dblock));
        auto dref = makeTenRef(d.data(),io.offset,d.size(),&drange);

        auto n = dblock[jc];

        for(auto o : index(C.store_))
            {
            auto& br = C.store_[o];

            //Only loop over subblocks of combined
            //indices compatible with current sector (==n)
            //of combined index cind (==dis[jc])
            if(br.block != size_type(n)) continue;

            //"invert" offset o into 
            //indices of nblock corresponding to
            //newly restored uncombined indices
            //using algorithm
            //similar to computeBlockInd
            for(auto m : count(ncomb-1))
                {
                nblock[jc+m] = o % Cis[1+m].nindex();
                o = (o-nblock[jc+m])/Cis[1+m].nindex();
                }
            nblock[jc+ncomb-1] = o;

            //fill out rest of nblock
            for(auto m : count(jc)) nblock[m] = dblock[m];
            for(auto m : count(1+jc,dr)) nblock[ncomb+m-1] = dblock[m];

            //Get subblock of d data
            auto dsub = subIndex(dref,jc,br.start,br.start+br.extent);

            nrange.init(make_indexdim(Nis,nblock));
            auto nb = getBlock(nd,Nis,nblock);
            assert(nb.data() != nullptr);
            auto nref = TensorRef(nb,&nrange);
            auto nslice = groupInds(nref,jc,jc+ncomb);

            nslice &= dsub;

            } //for br in C storage
        } //for blocks of d
    }//uncombine

void
doTask(Contract<IQIndex>      & C,
       IQTReal           const& d,
       IQTCombiner        const& cmb,
       ManageStore            & m)
    {
    if(C.Ris.r()==2)
        {
        combReplaceIndex(C.Lis,C.Ris,C.Nis);
        }
    else if(hasindex(C.Lis,C.Ris[0]))
        {
        uncombine(d,cmb,C.Lis,C.Ris,C.Nis,m,true);
        }
    else
        {
        combine(d,cmb,C.Lis,C.Ris,C.Nis,m);
        }
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
        uncombine(d,cmb,C.Ris,C.Lis,C.Nis,m,false);
        }
    else
        {
        combine(d,cmb,C.Ris,C.Lis,C.Nis,m);
        }
    }


} //namespace itensor
