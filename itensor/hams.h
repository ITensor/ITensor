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
    typedef Model SiteSetT;
private:
    const SiteSetT& iss;
    const int N;
    std::vector<IndexT> currentlinks;
public:

    int NN() const { return N; }
    IndexT si(int i) const { return iss.si(i); }

    HamBuilder(const SiteSetT& iss_) : iss(iss_), N(iss.NN()), currentlinks(N) { }

    void newlinks(std::vector<Index>& currentlinks)
	{
        currentlinks.resize(N);
        static int ver = 0; ++ver;
        for(int i = 1; i < N; i++)
        {
            std::stringstream ss;
            ss << "h" << ver << "-" << i;
            currentlinks[i] = IndexT(ss.str(),1);
        }
	}

    void newlinks(std::vector<IQIndex>& currentlinks)
	{
        currentlinks.resize(N);
        static int ver = 0; ++ver;
        for(int i = 1; i < N; i++)
        {
            std::stringstream ss;
            ss << "h" << ver << "-" << i;
            currentlinks[i] = IndexT(ss.str());
        }
	}

    Tensor unit(int i) const { return Tensor(si(i),primed(si(i)),1); }

    void getidentity(Real factor, MPOt<ITensor>& res)
	{
        newlinks(currentlinks);
        res = MPOt<ITensor>(iss,res.maxm(),res.cutoff());
        res.AAnc(1) = unit(1); res.AAnc(1).addindex1(GET(currentlinks,1));
        res.AAnc(1) *= factor;
        res.AAnc(N) = unit(N); res.AAnc(N).addindex1(GET(currentlinks,N-1));
        for(int i = 2; i < N; ++i)
        {
            res.AAnc(i) = unit(i); 
            res.AAnc(i).addindex1(GET(currentlinks,i-1));
            res.AAnc(i).addindex1(GET(currentlinks,i));
        }
	}

    void getMPO(Real factor, int i, Tensor op, MPOt<ITensor>& res)
	{
        getidentity(1,res);
        res.AAnc(i) = op;
        if(i > 1) res.AAnc(i).addindex1(GET(currentlinks,i-1));
        if(i < N) res.AAnc(i).addindex1(GET(currentlinks,i));
        res *= factor;
	}

    void getMPO(Real factor, int i1, Tensor op1, int i2, Tensor op2, MPOt<ITensor>& res)
	{
        if(i1 == i2) Error("HamBuilder::getMPO: i1 cannot equal i2.");
        getMPO(1,i2,op2,res);
        res.AAnc(i1) = op1;
        if(i1 > 1) res.AAnc(i1).addindex1(GET(currentlinks,i1-1));
        if(i1 < N) res.AAnc(i1).addindex1(GET(currentlinks,i1));
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
            res.AAnc(s) = GET(ops,i);
            assert(GET(ops,i).hasindex(si(sites[i])));
            if(s > 1) res.AAnc(s).addindex1(GET(currentlinks,s-1));
            if(s < N) res.AAnc(s).addindex1(GET(currentlinks,s));
        }
        res *= factor;
	}

};

} //namespace Internal
typedef Internal::HamBuilder<ITensor> HamBuilder;
//typedef Internal::HamBuilder<IQTensor> IQHamBuilder;

class MPOBuilder
    {
    public:

    MPOBuilder(const Model& model_) 
        : 
        model(model_), 
        Ns(model_.NN()) 
        { }

    virtual ~MPOBuilder() { }

    ITensor 
    makeLedge(const Index& L) const
        {
        ITensor res(L); 
        res(L(L.m())) = 1;
        return res;
        }

    ITensor 
    makeRedge(const Index& R) const
        {
        ITensor res(R); 
        res(R(1)) = 1;
        return res;
        }

    IQTensor 
    makeLedge(const IQIndex& L, const std::vector<Index>& start_inds)
        {
        IQTensor res(L);
        Foreach(const Index& ind, start_inds)
            {
            ITensor ledge(ind,0.0);
            ledge(ind(ind.m())) = 1.0; // [0 0 0 ... 1]
            res.insert(ledge);
            }
        return res;
        }

    IQTensor 
    makeRedge(const IQIndex& R, const std::vector<Index>& end_inds)
        {
        IQTensor res(R);
        Foreach(const Index& ind, end_inds)
            {
            ITensor redge(ind,0.0);
            redge(ind(1)) = 1.0; // [1 0 0 ... 0]
            res.insert(redge);
            }
        return res;
        }

    protected:

    const Model& model;
    const int Ns;

    };

/*
namespace SpinHalf 
{

namespace SquareLattice
{

//SpinHalf::SquareLattice::Heisenberg
class Heisenberg : public MPOBuilder
{
public:
    const Model& mod;
    const int ny,nx; 

    Heisenberg(const Model& mod_, int ny_) : MPOBuilder(mod_), mod(mod_), ny(ny_), nx(this->Ns/ny)
    { }

    void getIMPO(Real J, MPO& H, ITensor& Ledge, ITensor& Redge, Real boundary_h = 0)
    {
        H = MPO(mod);

        if(H.si(1).m() != SpinHalf::Dim) { Error("SpinHalf::SquareLattice::Heisenberg is only defined for S=1/2"); }

        const int nop = 3;
        const int k = 2+nop*ny;

        std::vector<Index> links(Ns+1);
        for(int l = 0; l <= Ns; ++l) links[l] = Index(nameint("hl",l),k);

        for(int n = 1; n <= Ns; ++n)
        {
            const int x = (n-1)/ny+1, y = (n-1)%ny+1;
            //cerr << format("n (x,y) = %d (%d,%d)\n")%n%x%y;
            ITensor& W = H.AAnc(n);
            Index &row = links[n-1], &col = links[n];

            W = ITensor(mod.si(n),mod.siP(n),row,col);

            W += mod.id(n) * row(1) * col(1);
            W += mod.id(n) * row(k) * col(k);

            W += mod.sz(n) * row(2) * col(1);
            W += mod.sp(n) * row(3) * col(1);
            W += mod.sm(n) * row(4) * col(1);

            for(int q = 1; q <= (ny-1)*nop; ++q)
            { W += mod.id(n) * row(1+nop+q) * col(1+q); }

            W += mod.sz(n) * row(k) * col(2+nop*(ny-1));
            W += mod.sm(n) * row(k) * col(3+nop*(ny-1)) * 0.5;
            W += mod.sp(n) * row(k) * col(4+nop*(ny-1)) * 0.5;

            if(y == 1 && ny > 1)
            {
            W += mod.sz(n) * row(k) * col(2+nop*(ny-2));
            W += mod.sm(n) * row(k) * col(3+nop*(ny-2)) * 0.5;
            W += mod.sp(n) * row(k) * col(4+nop*(ny-2)) * 0.5;
            }

            if(y != ny)
            {
            W += mod.sz(n) * row(k) * col(2);
            W += mod.sm(n) * row(k) * col(3) * 0.5;
            W += mod.sp(n) * row(k) * col(4) * 0.5;
            }

            if(boundary_h!=0 && (x == 1 || x == nx) )
            {
                Real eff_h = boundary_h;
                if(J > 0) eff_h *= (x%2==1 ? -1 : 1)*(y%2==1 ? -1 : 1);
                //cerr << format("Doing a staggered bf of %.2f at site %d (%d,%d)\n")%eff_h%n%x%y;
                W += mod.sz(n) * ITensor(row(k),col(1)) * eff_h;
            }
            //W.print("W",ShowData);
        }
        Ledge = makeLedge(links[0]);
        Redge = makeRedge(links[Ns]);
        //Ledge.print("Ledge",ShowData);
        //Redge.print("Redge",ShowData);
    }

    void getMPO(Real J, MPO& H, Real boundary_h=0)
    {
        ITensor Ledge,Redge;
        getIMPO(J,H,Ledge,Redge,boundary_h);
        
        H.AAnc(1) = Ledge * H.AA(1);
        H.AAnc(Ns) = H.AA(Ns) * Redge;
    }

    MPO operator()(Real J = 1, Real boundary_h=0)
    {
        MPO H; getMPO(J,H,boundary_h);
        return H;
    }

}; //class SquareLattice::Heisenberg

} //namespace SquareLattice

} //end namespace SpinHalf
*/

/*
namespace SpinOne 
{

class Heisenberg : public MPOBuilder
{
public:
    const Model& model;
    Real J;

    Heisenberg(const Model& model_) : MPOBuilder(model_), model(model_), J(1)
    { }

    void getIMPO(MPO& H, ITensor& Ledge, ITensor& Redge)
    {
        H = MPO(model);

        if(H.si(1).m() != SpinOne::Dim) Error("SpinOne::Heisenberg is only defined for S=1");

        const int k = 5;

        std::vector<Index> links(Ns+1);
        for(int l = 0; l <= Ns; ++l) links[l] = Index(nameint("hl",l),k);

        ITensor W;
        for(int n = 1; n <= Ns; ++n)
        {
            ITensor& W = H.AAnc(n);
            Index &row = links[n-1], &col = links[n];

            W = ITensor(model.si(n),model.siP(n),row,col);

            W += model.id(n) * row(1) * col(1);
            W += model.id(n) * row(k) * col(k);

            W += model.sz(n) * row(2) * col(1);
            W += model.sp(n) * row(3) * col(1);
            W += model.sm(n) * row(4) * col(1);

            W +=  J   *model.sz(n) * row(k) * col(2);
            W += (J/2)*model.sm(n) * row(k) * col(3);
            W += (J/2)*model.sp(n) * row(k) * col(4);
        }

        Ledge = makeLedge(links[0]);
        Redge = makeRedge(links[Ns]);
    }

    void getMPO(MPO& H)
    {
        ITensor Ledge,Redge;
        getIMPO(H,Ledge,Redge);
        
        H.AAnc(1) = Ledge * H.AA(1);
        H.AAnc(Ns) = H.AA(Ns) * Redge;
    }

    MPO operator()(Real J_ = 1)
    {
        const Real Jorig = J;
        J = J_;
        MPO H; getMPO(H);
        J = Jorig;
        return H;
    }

}; //class Heisenberg


} //end namespace SpinOne
*/

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

            H.AAnc(n) = W;
        }

        Ledge = makeLedge(links[0]);
        Redge = makeRedge(links[Ns]);
    }

    void getMPO(Real J, Real K, MPO& H)
    {
        ITensor Ledge,Redge;
        getIMPO(J,K,H,Ledge,Redge);
        
        H.AAnc(1) = Ledge * H.AA(1);
        H.AAnc(Ns) = H.AA(Ns) * Redge;
    }

}; //class BBChain
}
*/

/*
namespace Hubbard
{

class HubbardChain : public MPOBuilder
{
public:
    const Model& mod;

    HubbardChain(const Model& mod_) : MPOBuilder(mod_), mod(mod_) {}

    void getTrotterGates(Real U, Real tau, Real Etot, std::vector<IQTensor>& gates)
	{
	int N = Ns;
	for(int i = 1; i < N; i++)
	    {
	    IQTensor ui(mod.id(i)), uj(mod.id(i+1));
	    IQTensor udi(mod.NupNdn(i)), udj(mod.NupNdn(i+1));
	    IQTensor cui(mod.Cup(i)), cuj(mod.Cup(i+1));
	    IQTensor cdi(mod.Cdn(i)), cdj(mod.Cdn(i+1));
	    IQTensor dui(mod.Cdagup(i)), duj(mod.Cdagup(i+1));
	    IQTensor ddi(mod.Cdagdn(i)), ddj(mod.Cdagdn(i+1));
	    IQTensor fi(mod.FermiPhase(i));
	    IQTensor one = ui * uj;
	    udi *= (i == 1 ? U : U*0.5);
	    udj *= (i == N-1 ? U : U*0.5);
	    IQTensor udbond = udi * uj + ui * udj; 
	    IQTensor ke = multSiteOps(dui,fi) * cuj + duj * multSiteOps(fi,cui) + 
			multSiteOps(ddi,fi) * cdj + ddj * multSiteOps(fi,cdi);
	    ke *= -1.0;
	    IQTensor con = one * (-Etot/(N-1.0));
	    IQTensor ham = (udbond + ke + con);
	    ham *= -tau;
	    //cout << "ham is " << ham;
	    //cout << "one is " << one;
	    IQTensor term = ham, kk;
// exp(x) = 1 + x +  x^2/2! + x^3/3! ..
// = 1 + x * (1 + x/2 *(1 + x/3 * (...
// ~ ((x/3 + 1) * x/2 + 1) * x + 1
	    for(int o = 100; o >= 1; o--)
		{
		if(o != 1) term *= 1.0 / o;
		kk = one + term;
		term = multSiteOps(kk,ham);
		}
	    //kk.reverse_dirs();
	    gates[i] = kk;
	    //if(i == 2) cout << "expbondham 2 is " << expbondham[2];
	    }
	}
    void getIMPO(Real U, MPO& H, ITensor& Ledge, ITensor& Redge)
	{
        H = MPO(mod);

        const int k = 6;

	std::vector<Index> links(Ns+1);
        for(int l = 0; l <= Ns; ++l) links.at(l) = Index(nameint("hl",l),k);

        ITensor W;
        for(int n = 1; n <= Ns; ++n)
	    {
            ITensor& W = H.AAnc(n);
            Index &row = links.at(n-1), &col = links.at(n);

            // MPO is (save this):
            //
            // 1 / I                                  \
            // 2 | -c+up  0                           |
            // 3 | -c+dn  0       0                   |
            // 4 | -cup   0       0    0              |
            // 5 | -cdn   0       0    0       0      | 
            // 6 | Unud   F cup F cdn c+up F c+dn F  I |
            //     1       2      3    4       5     6 
            // F = FermiPhase

            W = ITensor(mod.si(n),mod.siP(n),row,col);
            W += mod.NupNdn(n) * row(6) * col(1) * U;
            W += multSiteOps(mod.FermiPhase(n),mod.Cup(n)) * row(6) * col(2);
            W += multSiteOps(mod.FermiPhase(n),mod.Cdn(n)) * row(6) * col(3);
            W += multSiteOps(mod.Cdagup(n),mod.FermiPhase(n)) * row(6) * col(4);
            W += multSiteOps(mod.Cdagdn(n),mod.FermiPhase(n)) * row(6) * col(5);
            W += mod.id(n) * row(6) * col(6);

            W += mod.id(n) * row(1) * col(1);
            W += mod.Cdagup(n) * row(2) * col(1) * (-1.0);
            W += mod.Cdagdn(n) * row(3) * col(1) * (-1.0);
            W += mod.Cup(n) * row(4) * col(1) * (-1.0);
            W += mod.Cdn(n) * row(5) * col(1) * (-1.0);
	    }

        Ledge = makeLedge(links.at(0));
        Redge = makeRedge(links.at(Ns));
	}

    void getMPO(Real U, MPO& H)
    {
        ITensor Ledge,Redge;
        getIMPO(U,H,Ledge,Redge);
        
        H.AAnc(1) = Ledge * H.AA(1);
        H.AAnc(Ns) = H.AA(Ns) * Redge;
    }

}; //class HubbardChain

} //end namespace Hubbard
*/

/*
namespace Hubbard
{

class HubbardChain : public MPOBuilder
{
public:
    const int Nsp; //number of 'split sites' (== 2 * number of hubbard sites)

    HubbardChain(IndexSiteSet& iss) : MPOBuilder(iss), Nsp(iss.NN()) { }

    void getIMPO(Real U, MPO& H, ITensor& Ledge, ITensor& Redge)
    {
        const int maxm = 5000;
        const Real cutoff = 1E-15;

        H = MPO(Iss,maxm,cutoff);

        const int k = 7;

        std::vector<Index> links(Nsp+1);
        for(int l = 0; l <= Nsp; ++l) links[l] = Index(nameint("hl",l),k);

        const int Emp = 1, Occ = 2;

        ITensor W;
        for(int n = 1; n <= Nsp; ++n)
        {
            W = ITensor(links[n-1],links[n],H.si(n),primed(H.si(n)));

            W(1,1,Emp,Emp) = 1; W(1,1,Occ,Occ) = 1; //Ident at 1,1
            W(k,k,Emp,Emp) = 1; W(k,k,Occ,Occ) = 1; //Ident at k,k (k==7)

            // Rest of MPO is (save this):
            //
            // 1 / I                   \
            // 2 | a  0                |
            // 3 | 0  F  0             |
            // 4 | ad 0  0  0          |
            // 5 | 0  0  0  F  0       | 
            // 6 | n  0  0  0  0  0    |
            // 7 \ 0  0  ad 0  a  Un I /
            //     1  2  3  4  5  6  7
            //
            // Note: the hubbard U should only be put in for up (odd) sites
            //
            //

            if(n%2==1) //Up site
            {

            }
            else       //Dn site
            {

            }

            H.AAnc(n) = W;
        }

        Ledge = makeLedge(links[0]);
        Redge = makeRedge(links[Nsp]);
    }

    void getMPO(Real U, MPO& H)
    {
        ITensor Ledge,Redge;
        getIMPO(U,H,Ledge,Redge);
        
        H.AAnc(1) = Ledge * H.AA(1);
        H.AAnc(Nsp) = H.AA(Nsp) * Redge;
    }

}; //class HubbardChain

} //end namespace Hubbard
*/


#endif
