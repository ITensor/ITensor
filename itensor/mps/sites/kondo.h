// Created 2017/09/16
// Modified from hubbard.h in the ITensor Library.

#ifndef __ITENSOR_KONDO_H
#define __ITENSOR_KONDO_H
#include "itensor/mps/siteset.h"

namespace itensor {

class KondoSite;

using Kondo = BasicSiteSet<KondoSite>;

class KondoSite
    {
    IQIndex s;
    public:

    KondoSite() { }

    KondoSite(IQIndex I) : s(I) { }

    KondoSite(int n, Args const& args = Args::global())
        {
        auto conserveNf = args.getBool("ConserveNf",true);
        auto conserveSz = args.getBool("ConserveSz",true); //total Sz
        int Up = (conserveSz ? +1 : 0),
            Dn = -Up;
        if(conserveNf)
            {
                s = IQIndex{nameint("site=",n),
                    Index(nameint("Emp.Up ",n),1,Site), QN("Sz=", Up,"Nf=",0),    //Sz is the total Sz,
                    Index(nameint("Up.Up ",n),1,Site),  QN("Sz=",2*Up,"Nf=",1),   //including that of 
                    Index(nameint("Dn.Up ",n),1,Site),  QN("Sz=",0,"Nf=",1),      //local spin
                    Index(nameint("UpDn.Up ",n),1,Site),QN("Sz=", Up,"Nf=",2),
                    Index(nameint("Emp.Dn ",n),1,Site), QN("Sz=", Dn,"Nf=",0),
                    Index(nameint("Up.Dn ",n),1,Site),  QN("Sz=",0,"Nf=",1),
                    Index(nameint("Dn.Dn ",n),1,Site),  QN("Sz=",2*Dn,"Nf=",1),
                    Index(nameint("UpDn.Dn ",n),1,Site),QN("Sz=", Dn,"Nf=",2)};
            

            }
        else //don't conserve Nf, only fermion parity
            {
            if(!conserveSz) Error("One of ConserveSz or ConserveNf must be true for Kondo sites");

             s = IQIndex{nameint("site=",n),
                    Index(nameint("Emp.Up ",n),1,Site), QN("Sz=", 1,"Pf=",0),
                    Index(nameint("Up.Up ",n),1,Site),  QN("Sz=",2,"Pf=",1),
                    Index(nameint("Dn.Up ",n),1,Site),  QN("Sz=",0,"Pf=",1),
                    Index(nameint("UpDn.Up ",n),1,Site),QN("Sz=", 1,"Pf=",0),
                    Index(nameint("Emp.Dn ",n),1,Site), QN("Sz=", -1,"Pf=",0),
                    Index(nameint("Up.Dn ",n),1,Site),  QN("Sz=",0,"Pf=",1),
                    Index(nameint("Dn.Dn ",n),1,Site),  QN("Sz=",-2,"Pf=",1),
                    Index(nameint("UpDn.Dn ",n),1,Site),QN("Sz=", -1,"Pf=",0)};
            }
        }

    IQIndex
    index() const { return s; }

    IQIndexVal
    state(std::string const& state)
        {
        if(state == "0+" || state == "Emp.Up") 
            {
            return s(1);
            }
        else 
        if(state == "++" || state == "Up.Up") 
            {
            return s(2);
            }
        else 
        if(state == "-+" || state == "Dn.Up") 
            {
            return s(3);
            }
        else 
        if(state == "S+" || state == "UpDn.Up") 
            {
            return s(4);
            }
        else
        if(state == "0-" || state == "Emp.Dn") 
            {
            return s(5);
            }
        else 
        if(state == "+-" || state == "Up.Dn") 
            {
            return s(6);
            }
        else 
        if(state == "--" || state == "Dn.Dn") 
            {
            return s(7);
            }
        else 
        if(state == "S-" || state == "UpDn.Dn") 
            {
            return s(8);
            }
        else
            {
            Error("State " + state + " not recognized");
            }
        return IQIndexVal{};
        }

	IQTensor
	op(std::string const& opname,
	   Args const& args) const
        {
        auto sP = prime(s);

        IQIndexVal Em_Up(s(1)),    // Up local spin
                   Em_Up_P(sP(1)), // P ~ prime
                   Up_Up(s(2)),
                   Up_Up_P(sP(2)),
                   Dn_Up(s(3)),
                   Dn_Up_P(sP(3)),
                   UD_Up(s(4)),
                   UD_Up_P(sP(4)),
                   Em_Dn(s(5)),    
                   Em_Dn_P(sP(5)), 
                   Up_Dn(s(6)),
                   Up_Dn_P(sP(6)),
                   Dn_Dn(s(7)),
                   Dn_Dn_P(sP(7)),
                   UD_Dn(s(8)),
                   UD_Dn_P(sP(8));

        IQTensor Op(dag(s),sP);

        if(opname == "Nup")
            {
            Op.set(Up_Up,Up_Up_P,1);
            Op.set(UD_Up,UD_Up_P,1);
            Op.set(Up_Dn,Up_Dn_P,1);
            Op.set(UD_Dn,UD_Dn_P,1);
            }
        else
        if(opname == "Ndn")
            {
            Op.set(Dn_Up,Dn_Up_P,1);
            Op.set(UD_Up,UD_Up_P,1);
            Op.set(Dn_Dn,Dn_Dn_P,1);
            Op.set(UD_Dn,UD_Dn_P,1);
            }
        else
        if(opname == "Nupdn")
            {
            Op.set(UD_Up,UD_Up_P,1);
            Op.set(UD_Dn,UD_Dn_P,1);
            }
        else
        if(opname == "Ntot")
            {
            Op.set(Up_Up,Up_Up_P,1);
            Op.set(Dn_Up,Dn_Up_P,1);
            Op.set(UD_Up,UD_Up_P,2);
            Op.set(Up_Dn,Up_Dn_P,1);
            Op.set(Dn_Dn,Dn_Dn_P,1);
            Op.set(UD_Dn,UD_Dn_P,2);
            }
        else
        if(opname == "Cup")
            {
            Op.set(Up_Up,Em_Up_P,1); 
            Op.set(UD_Up,Dn_Up_P,1); 
            Op.set(Up_Dn,Em_Dn_P,1); 
            Op.set(UD_Dn,Dn_Dn_P,1); 
            }
        else
        if(opname == "Cdagup")
            {
            Op.set(Em_Up,Up_Up_P,1); 
            Op.set(Dn_Up,UD_Up_P,1);
            Op.set(Em_Dn,Up_Dn_P,1); 
            Op.set(Dn_Dn,UD_Dn_P,1);
            }
        else
        if(opname == "Cdn")
            {
            Op.set(Dn_Up,Em_Up_P,1); 
            Op.set(UD_Up,Up_Up_P,-1); 
            Op.set(Dn_Dn,Em_Dn_P,1); 
            Op.set(UD_Dn,Up_Dn_P,-1); 
            }
        else
        if(opname == "Cdagdn")
            {
            Op.set(Em_Up,Dn_Up_P,1); 
            Op.set(Up_Up,UD_Up_P,-1);
            Op.set(Em_Dn,Dn_Dn_P,1); 
            Op.set(Up_Dn,UD_Dn_P,-1);
            }
        else
        if(opname == "Aup")
            {
            Op.set(Up_Up,Em_Up_P,1); 
            Op.set(UD_Up,Dn_Up_P,1); 
            Op.set(Up_Dn,Em_Dn_P,1); 
            Op.set(UD_Dn,Dn_Dn_P,1); 
            }
        else
        if(opname == "Adagup")
            {
            Op.set(Em_Up,Up_Up_P,1); 
            Op.set(Dn_Up,UD_Up_P,1);
            Op.set(Em_Dn,Up_Dn_P,1); 
            Op.set(Dn_Dn,UD_Dn_P,1);
            }
        else
        if(opname == "Adn")
            {
            Op.set(Dn_Up,Em_Up_P,1); 
            Op.set(UD_Up,Up_Up_P,1); 
            Op.set(Dn_Dn,Em_Dn_P,1); 
            Op.set(UD_Dn,Up_Dn_P,1); 
            }
        else
        if(opname == "Adagdn")
            {
            Op.set(Em_Up,Dn_Up_P,1); 
            Op.set(Up_Up,UD_Up_P,1);
            Op.set(Em_Dn,Dn_Dn_P,1); 
            Op.set(Up_Dn,UD_Dn_P,1);
            }
        else
        if(opname == "FermiPhase" || opname == "F")
            {
            Op.set(Em_Up,Em_Up_P,+1); 
            Op.set(Up_Up,Up_Up_P,-1);
            Op.set(Dn_Up,Dn_Up_P,-1);
            Op.set(UD_Up,UD_Up_P,+1);

            Op.set(Em_Dn,Em_Dn_P,+1); 
            Op.set(Up_Dn,Up_Dn_P,-1);
            Op.set(Dn_Dn,Dn_Dn_P,-1);
            Op.set(UD_Dn,UD_Dn_P,+1);
            }
        else
        if(opname == "Fup")
            {
            Op.set(Em_Up,Em_Up_P,+1); 
            Op.set(Up_Up,Up_Up_P,-1);
            Op.set(Dn_Up,Dn_Up_P,+1);
            Op.set(UD_Up,UD_Up_P,-1);

            Op.set(Em_Dn,Em_Dn_P,+1); 
            Op.set(Up_Dn,Up_Dn_P,-1);
            Op.set(Dn_Dn,Dn_Dn_P,+1);
            Op.set(UD_Dn,UD_Dn_P,-1);
            }
        else
        if(opname == "Fdn")
            {
            Op.set(Em_Up,Em_Up_P,+1); 
            Op.set(Up_Up,Up_Up_P,+1);
            Op.set(Dn_Up,Dn_Up_P,-1);
            Op.set(UD_Up,UD_Up_P,-1);

            Op.set(Em_Dn,Em_Dn_P,+1); 
            Op.set(Up_Dn,Up_Dn_P,+1);
            Op.set(Dn_Dn,Dn_Dn_P,-1);
            Op.set(UD_Dn,UD_Dn_P,-1);
            }
        else
        if(opname == "Sz")
            {
            Op.set(Up_Up,Up_Up_P,+0.5); 
            Op.set(Dn_Up,Dn_Up_P,-0.5);
            Op.set(Up_Dn,Up_Dn_P,+0.5); 
            Op.set(Dn_Dn,Dn_Dn_P,-0.5);
            }
        else
        if(opname == "S+")
            {
            Op.set(Dn_Up,Up_Up_P,1); 
            Op.set(Dn_Dn,Up_Dn_P,1); 
            }
        else
        if(opname == "S-")
            {
            Op.set(Up_Up,Dn_Up_P,1); 
            Op.set(Up_Dn,Dn_Dn_P,1); 
            }
        else
        if(opname == "S2")
            {
            //S dot S on-site
            Op.set(Up_Up,Up_Up_P,0.75); 
            Op.set(Dn_Up,Dn_Up_P,0.75);
            Op.set(Up_Dn,Up_Dn_P,0.75); 
            Op.set(Dn_Dn,Dn_Dn_P,0.75);
            }
        else
        if(opname == "LSz")    // L ~ local spin
            {
            Op.set(Em_Up,Em_Up_P,+0.5); 
            Op.set(Em_Dn,Em_Dn_P,-0.5);
            Op.set(Up_Up,Up_Up_P,+0.5); 
            Op.set(Up_Dn,Up_Dn_P,-0.5);
            Op.set(Dn_Up,Dn_Up_P,+0.5); 
            Op.set(Dn_Dn,Dn_Dn_P,-0.5);
            Op.set(UD_Up,UD_Up_P,+0.5); 
            Op.set(UD_Dn,UD_Dn_P,-0.5);
            }
        else
        if(opname == "LSx")    
            {
            Op.set(Em_Up,Em_Dn_P,+0.5); 
            Op.set(Em_Dn,Em_Up_P,+0.5);
            Op.set(Up_Up,Up_Dn_P,+0.5); 
            Op.set(Up_Dn,Up_Up_P,+0.5);
            Op.set(Dn_Up,Dn_Dn_P,+0.5); 
            Op.set(Dn_Dn,Dn_Up_P,+0.5);
            Op.set(UD_Up,UD_Dn_P,+0.5); 
            Op.set(UD_Dn,UD_Up_P,+0.5);
            }
        else
        if(opname == "LSy")    
            {
            Op.set(Em_Up,Em_Dn_P,-0.5*Cplx_i); 
            Op.set(Em_Dn,Em_Up_P,+0.5*Cplx_i);
            Op.set(Up_Up,Up_Dn_P,-0.5*Cplx_i); 
            Op.set(Up_Dn,Up_Up_P,+0.5*Cplx_i);
            Op.set(Dn_Up,Dn_Dn_P,-0.5*Cplx_i); 
            Op.set(Dn_Dn,Dn_Up_P,+0.5*Cplx_i);
            Op.set(UD_Up,UD_Dn_P,-0.5*Cplx_i); 
            Op.set(UD_Dn,UD_Up_P,+0.5*Cplx_i);
            }
        else
        if(opname == "LS+")
            {
            Op.set(Em_Dn,Em_Up_P,1); 
            Op.set(Up_Dn,Up_Up_P,1); 
            Op.set(Dn_Dn,Dn_Up_P,1); 
            Op.set(UD_Dn,UD_Up_P,1); 
            }
        else
        if(opname == "LS-")
            {
            Op.set(Em_Up,Em_Dn_P,1); 
            Op.set(Up_Up,Up_Dn_P,1); 
            Op.set(Dn_Up,Dn_Dn_P,1); 
            Op.set(UD_Up,UD_Dn_P,1); 
            }
        else
        if(opname == "LS2")
            {
            //S dot S on-site
            Op.set(Em_Up,Em_Up_P,0.75); 
            Op.set(Em_Dn,Em_Dn_P,0.75);
            Op.set(Up_Up,Up_Up_P,0.75); 
            Op.set(Up_Dn,Up_Dn_P,0.75);
            Op.set(Dn_Up,Dn_Up_P,0.75); 
            Op.set(Dn_Dn,Dn_Dn_P,0.75);
            Op.set(UD_Up,UD_Up_P,0.75); 
            Op.set(UD_Dn,UD_Dn_P,0.75);
            }
        else
            {
            Error("Operator \"" + opname + "\" name not recognized");
            }

        return Op;
        }
    };


} //namespace itensor

#endif
