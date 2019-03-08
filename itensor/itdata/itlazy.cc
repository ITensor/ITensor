#include "itensor/itdata/itlazy.h"
#include "itensor/itdata/dotask.h"
#define REGISTER_ITDATA_HEADER_FILES
#include "itensor/itdata/storage_types.h"
#include "itensor/tensor/contract.h"

namespace itensor {

struct PrintStore { };

void
doTask(PrintStore,const ITLazy& L) { println("Storage type of P is ITLazy"); }
void
doTask(PrintStore,const ITReal& R) { println("Storage type of P is ITReal"); }

bool
hasResult(const ITLazy& Z) { return Z.hasResult(); }

PData
evaluate(ITLazy& Z)
    {
    if(Z.hasResult()) return Z.result();

    long maxr = 0;
    for(auto& t : Z.todo())
        {
        maxr = std::max(maxr,t.i.order());
        }

    //Convention for contA, contB
    //- negative means contracted
    //0 zero means uncontracted, not on C
    //+ positive means uncontracted, on C
    InfArray<int,12ul> contA(maxr,0),
                       contB(maxr,0);

    auto P = Z.store(0);
    IndexSet Ais = Z.iset(0);
    IndexSet Cis;
    IndexSetBuilder cis;

    auto N = long(Z.todo().size());
    for(decltype(N) n = 0; n < N-2; ++n)
        {
        contA.fill(0);
        contB.fill(0);

        auto& Bis = Z.iset(n+1);
        auto& Nis = Z.iset(n+2);

        long ncont = 0;
        for(decltype(Ais.size()) na = 0; na < Ais.size(); ++na)
        for(decltype(Bis.size()) nb = 0; nb < Bis.size(); ++nb)
            if(Ais[na] == Bis[nb])
                {
                ++ncont;
                contA[na] = -(nb+1);
                contB[nb] = -(na+1);
                }

        auto nuncont = Ais.order() + Bis.order() - 2*ncont;
        cis.resize(nuncont);

        //Indices we want in Cis:
        //o Uncontracted inds of A & B
        //o *If* in Nis, in same order
        //Possible algorithm:
        //1. Loop over inds of Nis in order,
        //   putting any inds shared with A
        //   or B into Cis
        //2. Put rest of uncontracted inds of 
        //   A & B into Cis (loop over Cis to
        //   see if already put in?)

        size_t cn = 0;
        for(const auto& In : Nis)
            {
            for(size_t na = 0; na < Ais.size(); ++na)
                if(In == Ais[na])
                    {
                    cis.setExtent(cn++,In);
                    contA[na] = cn;
                    break;
                    }
            for(decltype(Bis.size()) nb = 0; nb < Bis.size(); ++nb)
                if(In == Bis[nb])
                    {
                    cis.setExtent(cn++,In);
                    contB[nb] = cn;
                    break;
                    }
            }

        for(decltype(Ais.size()) na = 0; na < Ais.size(); ++na)
            if(contA[na] == 0)
                {
                cis.setExtent(cn++,Ais[na]);
                }

        for(decltype(Bis.size()) nb = 0; nb < Bis.size(); ++nb)
            if(contB[nb] == 0)
                {
                cis.setExtent(cn++,Bis[nb]);
                }

        Cis = cis.build();

        println("Ais = ",Ais);
        println("Bis = ",Bis);
        println("-> Cis = ",Cis);
        println("----------------");
        println("Nis = ",Nis);

        //Global::debug3() = true;
        doTask(Contract<Index>{Ais,Bis,Cis,true},P,Z.store(n+1));

        Ais = Cis;
        }

    //println("Setting Z result");
    //doTask(PrintStore{},P);

    Z.setResult(std::move(P));
    return Z.result();
    }

void
doTask(Contract<Index>& C,
       ITLazy& L,
       const ITLazy& R)
    {
    //Eventually handle these cases into dotask.h?
    if(L.hasResult() && R.hasResult())
        {
        throw ITError("Not yet implemented 1");
        }
    else if(L.hasResult())
        {
        throw ITError("Not yet implemented 2");
        }
    else if(R.hasResult())
        {
        throw ITError("Not yet implemented 3");
        }
    else
        {
        L.addStore(R);
        }
    }

void
doTask(Contract<Index>& C,
       ITLazy& L,
       const PData& R)
    {
    println("In doTask(Contract,ITLazy, const CPData)");
    if(L.hasResult())
        {
        //Eventually handle this in dotask.h?
        //--> call doTask(Contract,L.result(),R);
        }
    contractIS(C.Lis,C.Ris,C.Nis);
    L.addStore(C.Nis,R);
    }

void
doTask(Contract<Index>& C,
       const PData& R,
       const ITLazy& L,
       ManageStore& m)
    {
    if(L.hasResult())
        {
        //Eventually handle this in dotask.h?
        //--> call doTask(Contract,R,L.result());
        }
    contractIS(C.Lis,C.Ris,C.Nis);
    auto* nL = m.makeNewData<ITLazy>(L);
    nL->addStore(C.Nis,R);
    }

} // namespace itensor
