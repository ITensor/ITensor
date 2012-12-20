//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_HAMS_H
#define __ITENSOR_HAMS_H
#include "mpo.h"

namespace Internal {

template <class Tensor>
class HamBuilder
    {
    public:

    typedef Tensor TensorT;
    typedef typename Tensor::IndexT IndexT;

    private:

    const Model& mod;
    const int N_;
    std::vector<IndexT> currentlinks;

    public:

    int N() const { return N_; }
    IndexT si(int i) const { return mod.si(i); }

    HamBuilder(const Model& mod_) : mod(mod_), N_(mod.N()), currentlinks(N_) { }

    void 
    newlinks(std::vector<Index>& currentlinks)
        {
        currentlinks.resize(N_);
        static int ver = 0; ++ver;
        for(int i = 1; i < N_; i++)
            {
            std::stringstream ss;
            ss << "h" << ver << "-" << i;
            currentlinks[i] = IndexT(ss.str(),1);
            }
        }

    void 
    newlinks(std::vector<IQIndex>& currentlinks)
        {
        //
        // Needs to be updated to work for IQIndex
        //
        currentlinks.resize(N_);
        static int ver = 0; ++ver;
        for(int i = 1; i < N_; i++)
            {
            std::stringstream ss;
            ss << "h" << ver << "-" << i;
            currentlinks[i] = IndexT(ss.str());
            }
        }

    void 
    getidentity(Real factor, MPOt<ITensor>& res)
        {
        newlinks(currentlinks);
        res = MPOt<ITensor>(mod,res.maxm(),res.cutoff());
        res.Aref(1) = mod.id(1); 
        res.Aref(1).addindex1(GET(currentlinks,1));
        res.Aref(1) *= factor;
        res.Aref(N_) = mod.id(N_); 
        res.Aref(N_).addindex1(GET(currentlinks,N_-1));
        for(int i = 2; i < N_; ++i)
            {
            res.Aref(i) = mod.id(i);
            res.Aref(i).addindex1(GET(currentlinks,i-1));
            res.Aref(i).addindex1(GET(currentlinks,i));
            }
        }

    void 
    getMPO(Real factor, int i, Tensor op, MPOt<ITensor>& res)
        {
        getidentity(1,res);
        res.Aref(i) = op;
        if(i > 1) res.Aref(i).addindex1(GET(currentlinks,i-1));
        if(i < N_) res.Aref(i).addindex1(GET(currentlinks,i));
        res *= factor;
        }

    void 
    getMPO(Real factor, int i1, Tensor op1, int i2, Tensor op2, MPOt<ITensor>& res)
        {
        if(i1 == i2) Error("HamBuilder::getMPO: i1 cannot equal i2.");
        getMPO(1,i2,op2,res);
        res.Aref(i1) = op1;
        if(i1 > 1) res.Aref(i1).addindex1(GET(currentlinks,i1-1));
        if(i1 < N_) res.Aref(i1).addindex1(GET(currentlinks,i1));
        res *= factor;
        }

    template <typename Iterable1, typename Iterable2>
    void getMPO(Real factor, Iterable1 sites, Iterable2 ops, MPOt<ITensor>& res)
        {
        for(int i = 0; i < (int) sites.size(); ++i)
        for(int j = 0; j < (int) sites.size(); ++j)
            {
            if(i == j) continue;
            if(sites[i] == sites[j]) Error("HamBuilder::getMPO: all sites should be unique.");
            }

        if(sites.size() != ops.size()) Error("HamBuilder::getMPO: need same number of sites as ops.");

        getidentity(1,res);
        for(int i = 0; i < (int) sites.size(); ++i)
            {
            const int s = sites[i];
            res.Aref(s) = ops.at(i);
            //assert(GET(ops,i).hasindex(si(sites[i])));
            if(s > 1) res.Aref(s).addindex1(GET(currentlinks,s-1));
            if(s < N_) res.Aref(s).addindex1(GET(currentlinks,s));
            }
        res *= factor;
        }

    };

} //namespace Internal
//typedef Internal::HamBuilder<ITensor> HamBuilder;
//typedef Internal::HamBuilder<IQTensor> IQHamBuilder;


/*
namespace SpinOne {
class BBChain : public MPOBuilder
{
public:
    const int Ns; 

    BBChain(IndexSiteSet& iss) : MPOBuilder(iss), Ns(iss.NN()) { }

    void getIMPO(Real J, Real K, MPO& H, ITensor& Ledge, ITensor& Redge)
    {
        H = MPO(Iss);

        if(H.si(1).m() != SpinOne::Dim) Error("SpinOne::BBChain is only defined for S=1");

        const int k = 11;

        std::vector<Index> links(Ns+1);
        for(int l = 0; l <= Ns; ++l) links[l] = Index(nameint("hl",l),k);

        ITensor W;
        for(int n = 1; n <= Ns; ++n)
        {
            W = ITensor(H.si(n),primed(H.si(n)),links[n-1],links[n]);

            W(1,1,1,1) = 1.0;     W(2,2,1,1) = 1.0;     W(3,3,1,1) = 1.0;      //Ident in ul corner
            W(1,1,2,1) = 1.0;     W(2,2,2,1) = 0.0;     W(3,3,2,1) = -1.0;     //Sz
            W(2,1,3,1) = Sqrt2;   W(3,2,3,1) = Sqrt2;                          //S-
            W(1,2,4,1) = Sqrt2;   W(2,3,4,1) = Sqrt2;                          //S+
            W(1,1,5,1) = 1.0;     W(2,2,5,1) = 0.0;     W(3,3,5,1) = 1.0;      //Sz^2
            W(1,2,6,1) = Sqrt2;                                                //Sz S+
            W(3,2,7,1) = -Sqrt2;                                               //Sz S-
            W(1,3,8,1) = 2;                                                    //S+^2
            W(2,2,9,1) = 2;       W(3,3,9,1) = 2;                              //S- S+
            W(3,1,10,1) = 2;                                                   //S-^2
            W(1,1,11,2) = J;      W(2,2,11,2) = -K;     W(3,3,11,2) = -J-K;    //(J+K) Sz - K/2 S+ S-
            W(1,2,11,3) = J/Sqrt2; W(2,3,11,3) = (J+K)/Sqrt2;                  //(J+K)/2 S+ - K/2 Sz S+
            W(2,1,11,4) = (J+K)/Sqrt2; W(3,2,11,4) = J/Sqrt2;                  //(J+K)/2 S- + K/2 Sz S-
            W(1,1,11,5) = -K;     W(3,3,11,5) = -K;                            //-K Sz^2
            W(2,1,11,6) = -K/Sqrt2; W(3,2,11,6) = K/Sqrt2;                     //-K Sz S- - K/2 S-
            W(1,2,11,7) = -K/Sqrt2; W(2,3,11,7) = K/Sqrt2;                     //-K Sz S+ + K/2 S+
            W(3,1,11,8) = -K/2;                                                 //-K/4 S-^2
            W(1,1,11,9) = -K/2;   W(2,2,11,9) = -K;    W(3,3,11,9) = -K/2;     //-K/2 S+ S- + K/2 Sz
            W(1,3,11,10) = -K/2;                                               //-K/4 S+^2
            W(1,1,11,11) = 1.0;   W(2,2,11,11) = 1.0;   W(3,3,11,11) = 1.0;    //Ident in lr corner

            H.Aref(n) = W;
        }

        Ledge = makeLedge(links[0]);
        Redge = makeRedge(links[Ns]);
    }

    void getMPO(Real J, Real K, MPO& H)
    {
        ITensor Ledge,Redge;
        getIMPO(J,K,H,Ledge,Redge);
        
        H.Aref(1) = Ledge * H.AA(1);
        H.Aref(Ns) = H.AA(Ns) * Redge;
    }

}; //class BBChain
}
*/

#endif
