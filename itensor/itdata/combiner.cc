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
#include "itensor/itdata/combiner.h"
#include "itensor/itdata/itdata.h"
//#include "itensor/itdata/itcplx.h"
#include "itensor/tensor/contract.h"
#include "itensor/tensor/sliceten.h"

using std::vector;

namespace itensor {

const char*
typeNameOf(Combiner const& d) { return "Combiner"; }

Cplx
doTask(GetElt const& g, Combiner const& c)
    {
    if(g.inds.size()!=0) Error("GetElt not defined for non-scalar Combiner storage");
    return Cplx(1.,0.);
    }

Real
doTask(NormNoScale, Combiner const& d) { return 0; }

void
doTask(Conj,Combiner const& d) { }

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
    auto jc = indexPosition(dis,cind);
    if(jc >= 0) //has cind, uncombining
        {
        //dis has cind, replace with other inds
        auto newind = IndexSetBuilder(dis.order()+Cis.order()-2);
        long i = 0;
        for(auto j : range(dis.order()))
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
        auto J1 = indexPosition(dis,Cis[1]);
        if(J1 < 0) 
            {
            println("IndexSet of dense tensor = \n",dis);
            println("IndexSet of combiner/delta = \n",Cis);
            Error("No contracted indices in combiner-tensor product");
            }
        //Check if Cis[1],Cis[2],... are grouped together (contiguous)
        //and in same order as on combiner
        bool contig_sameord = true;
        decltype(Cis.order()) c = 2;
        for(auto j = J1+1; c < Cis.order() && j < dis.order(); ++j,++c)
            if(dis[j] != Cis[c])
                {
                contig_sameord = false;
                break;
                }
        if(c != Cis.order()) contig_sameord = false;

        if(contig_sameord)
            {
            auto newind = IndexSetBuilder(dis.order()+2-Cis.order());
            long i = 0;
            for(auto j : range(J1))
                newind.setIndex(i++,dis[j]);
            newind.setIndex(i++,cind);
            for(auto j : range(J1+Cis.order()-1, dis.order()))
                newind.setIndex(i++,dis[j]);
            Nis = newind.build();
            }
        else
            {
            auto P = Permutation(dis.order());
            //Set P destination values to -1 to mark
            //indices that need to be assigned destinations:
            for(auto i : range(P)) P.setFromTo(i,-1);

            //permute combined indices to the front, in same
            //order as in Cis:
            long ni = 0;
            for(auto c : range(1,Cis.order()))
                {
                auto j = indexPosition(dis,Cis[c]);
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
            auto newind = IndexSetBuilder(dis.order()+2-Cis.order());
            long i = 0;
            newind.setIndex(i++,cind);
            for(auto j : range(dis.order()))
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

template<typename V>
void
doTask(Contract & C,
       Dense<V>   const& d,
       Combiner const& cmb,
       ManageStore     & m)
    {
    combine(d,C.Lis,C.Ris,C.Nis,m);
    }
template void doTask(Contract &,DenseReal const&,Combiner const&,ManageStore&);
template void doTask(Contract &,DenseCplx const&,Combiner const&,ManageStore&);

template<typename V>
void
doTask(Contract & C,
       Combiner const& cmb,
       Dense<V>   const& d,
       ManageStore     & m)
    { 
    combine(d,C.Ris,C.Lis,C.Nis,m);
    if(!m.newData()) m.assignPointerRtoL();
    }
template void doTask(Contract &,Combiner const&,DenseReal const&,ManageStore&);
template void doTask(Contract &,Combiner const&,DenseCplx const&,ManageStore&);

bool
doTask(CheckComplex, Combiner const& d) { return false; }

void
doTask(PrintIT& P, Combiner const& d)
    {
    P.printInfo(d,"Combiner");
    }

QN 
doTask(CalcDiv const& C, 
       Combiner const& d)
    {
    return QN{};
    }

} //namespace itensor
