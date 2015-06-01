//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/itdata/itcombiner.h"
#include "itensor/itdata/itdata.h"
#include "itensor/tensor/contract.h"
#include "itensor/util/count.h"

using std::vector;

namespace itensor {

Cplx
doTask(const GetElt<Index>& g, const ITCombiner& c)
    {
    if(g.inds.size()!=0) Error("GetElt not defined for non-scalar ITCombiner storage");
    return Cplx(1.,0.);
    }

Real
doTask(const NormNoScale<Index>& N, const ITCombiner& d) { return 0; }

void
doTask(Conj,const ITCombiner& d) { }

void
combine(const ITReal& d,
        const IndexSet& dis,
        const IndexSet& Cis,
        IndexSet& Nis,
        ManagePtr& mp)
    {
    //TODO: try to make use of Lind,Rind label vectors
    //      to simplify combine logic
    const auto& cind = Cis[0];
    auto jc = findindex(dis,cind);
    if(jc >= 0) //has cind
        {
        //dis has cind, replace with other inds
        IndexSet::storage_type newind(dis.r()+Cis.r()-2);
        long i = 0;
        for(auto j : count(dis.r()))
            if(j == jc)
                {
                for(size_t k = 1; k < Cis.size(); ++k)
                    newind.at(i++).ext = Cis[k];
                }
            else
                {
                newind.at(i++).ext = dis[j];
                }
        Nis = IndexSet(move(newind));
        }
    else
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
        int c = 2;
        for(int j = J1+1; c < Cis.r() && j < dis.r(); ++j,++c)
            if(dis[j] != Cis[c])
                {
                contig_sameord = false;
                break;
                }
        if(c != Cis.r()) contig_sameord = false;

        //printfln("%s:",contig_sameord?"Contig":"Not Contig");
        //println("  dis = ",dis);
        //println("  Cis = ",Cis);

        if(contig_sameord)
            {
            IndexSet::storage_type newind(dis.r()+2-Cis.r());
            long i = 0;
            for(int j = 0; j < J1; ++j) 
                newind.at(i++).ext = dis[j];
            newind.at(i++).ext = cind;
            for(int j = J1+Cis.r()-1; j < dis.r(); ++j) 
                newind.at(i++).ext = dis[j];
            Nis = IndexSet(move(newind));
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
            Range::storage_type pdims(dis.r());
            IndexSet::storage_type newind(dis.r()+2-Cis.r());
            long i = 0;
            newind[i++].ext = cind;
            for(auto j : count(dis.r()))
                {
                if(P.dest(j) == -1) 
                    {
                    P.setFromTo(j,ni++);
                    newind.at(i++).ext = dis[j];
                    }
                pdims.at(P.dest(j)).ext = dis[j].m();
                }
            assert(i==dis.r()+2-Cis.r());
            Range rr(move(pdims));
            Nis = IndexSet(move(newind));
            auto nd = mp.makeNewData<ITReal>(area(Nis));
            auto tr = makeTensorRef(nd->data(),rr);
            auto td = makeTensorRef(d.data(),dis);
            permute(td,P,tr);
            }
        }
    }

void
doTask(Contract<Index>& C,
       const ITReal& d,
       const ITCombiner& cmb,
       ManagePtr& mp)
    {
    combine(d,C.Lis,C.Ris,C.Nis,mp);
    }
void
doTask(Contract<Index>& C,
       const ITCombiner& cmb,
       const ITReal& d,
       ManagePtr& mp)
    { 
    combine(d,C.Ris,C.Lis,C.Nis,mp);
    if(!mp.newData()) mp.assignPointerRtoL();
    }

bool
doTask(CheckComplex, const ITCombiner& d) { return false; }

void
doTask(PrintIT<Index>& P, const ITCombiner& d)
    {
    P.printInfo(d,"Combiner");
    }

void
doTask(Write& W, const ITCombiner& d) 
    { 
    W.writeType(StorageType::ITCombiner,d); 
    }

} //namespace itensor
