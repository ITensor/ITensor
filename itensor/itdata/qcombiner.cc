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
#include "itensor/util/iterate.h"
#include "itensor/tensor/sliceten.h"
#include "itensor/itdata/qcombiner.h"
#include "itensor/itdata/itdata.h"
#include "itensor/itdata/qutil.h"

using std::vector;
using std::tie;

namespace itensor {

const char*
typeNameOf(QCombiner const& d) { return "QCombiner"; }

void 
read(std::istream& s, QCombiner & dat)
    {
    Range range;
    itensor::read(s,range);

    using storage_type = typename QCombiner::storage_type;
    storage_type store;
    itensor::read(s,store);

    dat = QCombiner(std::move(range),std::move(store));
    }

void
write(std::ostream& s, QCombiner const& dat)
    {
    itensor::write(s,dat.range());
    itensor::write(s,dat.store());
    }

std::ostream&
operator<<(std::ostream & s, QCombiner const& dat)
    {
    for(auto const& br : dat.store_) br.print(s);
    s << dat.range();
    return s;
    }

Cplx
doTask(GetElt const& g, QCombiner const& c)
    {
    if(g.inds.size()!=0) Error("GetElt not defined for non-scalar QCombiner storage");
    return Cplx(1.,0.);
    }

void
permuteIQ(Permutation const& P,
          IndexSet const& Ais,
          QDenseReal const& dA,
          IndexSet & Bis,
          QDenseReal & dB)
    {
#ifdef DEBUG
    if(isTrivial(P)) Error("Calling permuteIQ for trivial Permutation");
#endif
    auto r = Ais.order();
    auto bind = IndexSetBuilder(r);
    for(auto i : range(r))
        {
        bind.setIndex(P.dest(i),Ais[i]);
        }
    Bis = bind.build();
    dB = QDenseReal(Bis,doTask(CalcDiv{Ais},dA));

    auto Bblock = Block(r,-1);
    Range Arange,
          Brange;
    for(auto aio : dA.offsets)
        {
        //Compute bi, new block index of blk
        for(auto j : range(aio.block)) 
            {
            Bblock.at(P.dest(j)) = aio.block[j];
            }
        Arange.init(make_indexdim(Ais,aio.block));
        Brange.init(make_indexdim(Bis,Bblock));

        auto bblock = getBlock(dB,Bis,Bblock);
        auto bref = TensorRef(bblock,&Brange);

        auto aref = makeTenRef(dA.data(),aio.offset,dA.size(),&Arange);

        bref += permute(aref,P);
        }
    }

IndexSet
replaceInd(IndexSet const& is,
           long            loc,
           Index    const& replacement)
    {
    auto newind = IndexSetBuilder(is.order());
    long i = 0;
    for(auto j : range(loc))
        {
        newind.setIndex(i++,is[j]);
        }
    newind.setIndex(i++,replacement);
    for(decltype(is.order()) j = loc+1; j < is.order(); ++j)
        {
        newind.setIndex(i++,is[j]);
        }
    return newind.build();
    }

template<typename T>
void
combine(QDense<T> const& d,
        QCombiner const& C,
        IndexSet  const& dis,
        IndexSet  const& Cis,
        IndexSet       & Nis,
        ManageStore    & m)
    {
#ifdef DEBUG
    for(auto i : range(1,Cis.order()))
        {
        auto jc = indexPosition(dis,Cis[i]);
        if(jc == -1)
            {
            printfln("Indices of tensor = \n%s\n------",dis);
            printfln("Extra index = \n%s",Cis[i]);
            Error("Combiner has extra index not found on other tensor");
            }
        }
#endif

    using size_type = decltype(order(dis));
    auto dr = order(dis);
    auto ncomb = order(Cis)-1;
    auto nr = dr-ncomb+1;

    auto dperm = Labels(dr,-1);
    auto uncomb_dest = ncomb;
    for(auto i : range(dr)) 
        {
        auto jc = indexPosition(Cis,dis[i]);
        if(jc >= 0) dperm[i] = jc-1;
        else        dperm[i] = uncomb_dest++;
        }

    auto combined = [&dperm,ncomb](size_type i) { return dperm[i] < long(ncomb); };

    //Create new IndexSet
    auto newind = IndexSetBuilder(nr);
    newind.nextIndex(Cis[0]);
    for(auto i : range(dr)) if(!combined(i)) newind.nextIndex(dis[i]);
    Nis = newind.build();

    //Allocate new data
    auto& nd = *m.makeNewData<QDense<T>>(Nis,doTask(CalcDiv{dis},d));

    auto drange = Range(dr), //block range of current storage
         nrange = Range(nr); //block range of new storage
    auto nblock = Block(nr), //block index of new storage
         cblock = Block(ncomb); //corresponding subblock of combiner
    size_t start = 0, //offsets within sector of combined
           end   = 0; //Index where block will go
    for(auto io : d.offsets) //loop over non-zero blocks
        {
        //Make TensorRef for this block of d
        drange.init(make_indexdim(dis,io.block));
        auto dref = makeTenRef(d.data(),io.offset,d.size(),&drange);

        //Permute combined indices to front, then
        //group the first ncomb indices into one
        auto Pdref = Ten<Range,T>{permute(dref,dperm)};
        auto gPdref = groupInds(Pdref,0,ncomb);

        //Figure out "block index" where this block will
        //go in new storage (nblock) and which sector of
        //combined indices maps to new combined index (cblock)
        size_t nu = 1;
        for(auto i : range(dr)) 
            {
            if(combined(i)) cblock[dperm[i]] = io.block[i];
            else            nblock[nu++] = io.block[i];
            }

        //Use cblock to recover info about structure of combined Index,
        //which sector to map to, where this subsector starts, and ends
        tie(nblock[0],start,end) = C.getBlockRange(cblock);

        //Get full block of new storage
        nrange.init(make_indexdim(Nis,nblock));
        auto nb = getBlock(nd,Nis,nblock);
        assert(nb.data() != nullptr);
        auto nref = makeRef(nb,&nrange);

        //Slice this new-storage block to get subblock where data will go
        auto nsub = subIndex(nref,0,start,end);
        nsub &= gPdref;
        }
    }


template<typename T>
void
uncombine(QDense<T> const& d,
          QCombiner const& C,
          IndexSet  const& dis,
          IndexSet  const& Cis,
          IndexSet       & Nis,
          ManageStore    & m,
          bool              own_data)
    {
    using size_type = decltype(order(dis));
    auto& cind = Cis[0];
    auto dr = order(dis);
    auto cr = order(Cis);
    auto ncomb = cr-1;
    auto nr = dr-1+ncomb;

    decltype(dr) jc = 0;
    auto newind = IndexSetBuilder(nr);
    for(auto n : range(dr)) 
        {
        if(dis[n] == cind)
            {
            jc = n;
            for(auto n : range(1,cr)) newind.nextIndex(Cis[n]);
            }
        else
            {
            newind.nextIndex(dis[n]);
            }
        }
    Nis = newind.build();

    //Allocate new data
    auto& nd = *m.makeNewData<QDense<T>>(Nis,doTask(CalcDiv{dis},d));

    auto drange = Range(dr), //block range of current storage
         nrange = Range(nr); //block range of new storage
    auto nblock = Block(nr); //block index of new storage
    for(auto io : d.offsets) //loop over non-zero blocks
        {
        //Make TensorRef for this block of d
        drange.init(make_indexdim(dis,io.block));
        auto dref = makeTenRef(d.data(),io.offset,d.size(),&drange);

        auto n = io.block[jc];

        for(auto o : range(C.store_))
            {
            auto& br = C.store_[o];

            //Only loop over subblocks of combined
            //indices compatible with current sector (==n)
            //of combined index cind (==dis[jc])
            if(size_type(br.block) != size_type(n)) continue;

            //"invert" offset o into 
            //indices of nblock corresponding to
            //newly restored uncombined indices
            //using algorithm
            //similar to computeBlockInd
            for(auto m : range(ncomb-1))
                {
                nblock[jc+m] = o % Cis[1+m].nblock();
                o = (o-nblock[jc+m])/Cis[1+m].nblock();
                }
            nblock[jc+ncomb-1] = o;

            //fill out rest of nblock
            for(auto m : range(jc)) nblock[m] = io.block[m];
            for(auto m : range(1+jc,dr)) nblock[ncomb+m-1] = io.block[m];

            //Get subblock of d data
            auto dsub = subIndex(dref,jc,br.start,br.start+br.extent);

            nrange.init(make_indexdim(Nis,nblock));
            auto nb = getBlock(nd,Nis,nblock);
            assert(nb.data() != nullptr);
            auto nref = makeRef(nb,&nrange);
            auto nslice = groupInds(nref,jc,jc+ncomb);

            nslice &= dsub;

            } //for br in C storage
        } //for blocks of d
    }//uncombine

template<typename T>
void
doTask(Contract & C,
       QDense<T>  const& d,
       QCombiner  const& cmb,
       ManageStore     & m)
    {
    if(hasIndex(C.Lis,C.Ris[0]))
        {
        uncombine(d,cmb,C.Lis,C.Ris,C.Nis,m,true);
        }
    else
        {
        combine(d,cmb,C.Lis,C.Ris,C.Nis,m);
        }
    }
template void doTask(Contract &,QDense<Real> const&,QCombiner const&,ManageStore &);
template void doTask(Contract &,QDense<Cplx> const&,QCombiner const&,ManageStore &);

template<typename T>
void
doTask(Contract & C,
       QCombiner  const& cmb,
       QDense<T>    const& d,
       ManageStore       & m)
    { 
    if(hasIndex(C.Ris,C.Lis[0]))
        {
        uncombine(d,cmb,C.Ris,C.Lis,C.Nis,m,false);
        }
    else
        {
        combine(d,cmb,C.Ris,C.Lis,C.Nis,m);
        }
    }
template void doTask(Contract &,QCombiner const&,QDense<Real> const&,ManageStore &);
template void doTask(Contract &,QCombiner const&,QDense<Cplx> const&,ManageStore &);


} //namespace itensor
