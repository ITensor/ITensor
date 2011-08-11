#ifndef __HAMS_H
#define __HAMS_H
#include "mps.h"

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

        vector<Index> links(Ns+1);
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
        H = MPO(sst);

        if(H.si(1).m() != SpinOne::Dim) Error("SpinOne::Heisenberg is only defined for S=1");

        const int k = 5;

        vector<Index> links(Ns+1);
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

/*
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

        vector<Index> links(Ns+1);
        for(int l = 0; l <= Ns; ++l) links[l] = Index(nameint("hl",l),k);

        ITensor W;
        for(int n = 1; n <= Ns; ++n)
        {
            W = ITensor(H.si(n),H.si(n).primed(),links[n-1],links[n]);

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
*/

} //end namespace SpinOne

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

        vector<Index> links(Nsp+1);
        for(int l = 0; l <= Nsp; ++l) links[l] = Index(nameint("hl",l),k);

        const int Emp = 1, Occ = 2;

        ITensor W;
        for(int n = 1; n <= Nsp; ++n)
        {
            W = ITensor(links[n-1],links[n],H.si(n),H.si(n).primed());

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
