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
#ifndef __ITENSOR_ELECTRONK_H
#define __ITENSOR_ELECTRONK_H
#include "itensor/mps/siteset.h"
#include "itensor/util/str.h"

namespace itensor {

class ElectronKSite;

using ElectronK = BasicSiteSet<ElectronKSite>;

class ElectronKSite
  {
  Index s;
  public:

  ElectronKSite(Index I) : s(I) { }

  ElectronKSite(Args const& args = Args::global())
    {
    auto ts = TagSet("Site,ElecK");
    if(args.defined("SiteNumber"))
      ts.addTags("n="+str(args.getInt("SiteNumber")));
    auto conserveQNs = args.getBool("ConserveQNs",true);
    auto conserveNf = args.getBool("ConserveNf",conserveQNs);
    auto conserveSz = args.getBool("ConserveSz",conserveQNs);
    auto conserveK = args.getBool("ConserveK",conserveQNs);

    int site = args.getInt("SiteNumber");
    int modulus = args.getInt("Kmod", 0);
    int mod = (site-1) % modulus;

    int Up = (conserveSz ? +1 : 0),
        Dn = -Up;
    if(conserveQNs || conserveNf || conserveSz || conserveK)
      {
      if(conserveK)
        {
        s = Index(QN({"Sz", 0},{"Nf",0,-1},{"K",0,modulus}),  1,
                  QN({"Sz",Up},{"Nf",1,-1},{"K",mod,modulus}),1,
                  QN({"Sz",Dn},{"Nf",1,-1},{"K",mod,modulus}),1,
                  QN({"Sz", 0},{"Nf",2,-1},{"K",(2*mod)%modulus,modulus}),1,Out,ts);
        }
      else if(conserveNf)
        {
        s = Index(QN({"Sz", 0},{"Nf",0,-1}),1,
                  QN({"Sz",Up},{"Nf",1,-1}),1,
                  QN({"Sz",Dn},{"Nf",1,-1}),1,
                  QN({"Sz", 0},{"Nf",2,-1}),1,Out,ts);
        }
      else if(conserveSz) //don't conserve Nf, only fermion parity
        {
        s = Index(QN({"Sz", 0},{"Pf",0,-2}),1,
                  QN({"Sz",+1},{"Pf",1,-2}),1,
                  QN({"Sz",-1},{"Pf",1,-2}),1,
                  QN({"Sz", 0},{"Pf",0,-2}),1,Out,ts);
        }
      else
        {
        s = Index(QN({"Pf",0,-2}),1,
                  QN({"Pf",1,-2}),1,
                  QN({"Pf",1,-2}),1,
                  QN({"Pf",0,-2}),1,Out,ts);
        }
      }
    else
      {
      s = Index(4,ts);
      }
    }

  Index
  index() const { return s; }

  IndexVal
  state(std::string const& state)
    {
    if(state == "0" || state == "Emp") 
      {
      return s(1);
      }
    else 
    if(state == "+" || state == "Up") 
      {
      return s(2);
      }
    else 
    if(state == "-" || state == "Dn") 
      {
      return s(3);
      }
    else 
    if(state == "S" || state == "UpDn") 
      {
      return s(4);
      }
    else
      {
      Error("State " + state + " not recognized");
      }
    return IndexVal{};
    }

  ITensor
  op(std::string const& opname,
     Args const& args) const
    {
    auto sP = prime(s);

    IndexVal Em(s(1)),
             EmP(sP(1)),
             Up(s(2)),
             UpP(sP(2)),
             Dn(s(3)),
             DnP(sP(3)),
             UD(s(4)),
             UDP(sP(4));

    ITensor Op(dag(s),sP);

    if(opname == "Nup")
      {
      Op.set(Up,UpP,1);
      Op.set(UD,UDP,1);
      }
    else
    if(opname == "Ndn")
      {
      Op.set(Dn,DnP,1);
      Op.set(UD,UDP,1);
      }
    else
    if(opname == "Nupdn")
      {
      Op.set(UD,UDP,1);
      }
    else
    if(opname == "Ntot")
      {
      Op.set(Up,UpP,1);
      Op.set(Dn,DnP,1);
      Op.set(UD,UDP,2);
      }
    else
    if(opname == "Cup")
      {
      Op.set(Up,EmP,1); 
      Op.set(UD,DnP,1); 
      }
    else
    if(opname == "Cdagup")
      {
      Op.set(Em,UpP,1); 
      Op.set(Dn,UDP,1);
      }
    else
    if(opname == "Cdn")
      {
      Op.set(Dn,EmP,1); 
      Op.set(UD,UpP,-1); 
      }
    else
    if(opname == "Cdagdn")
      {
      Op.set(Em,DnP,1); 
      Op.set(Up,UDP,-1);
      }
    else
    if(opname == "Aup")
      {
      Op.set(Up,EmP,1); 
      Op.set(UD,DnP,1); 
      }
    else
    if(opname == "Adagup")
      {
      Op.set(Em,UpP,1); 
      Op.set(Dn,UDP,1);
      }
    else
    if(opname == "Adn")
      {
      Op.set(Dn,EmP,1); 
      Op.set(UD,UpP,1); 
      }
    else
    if(opname == "Adagdn")
      {
      Op.set(Em,DnP,1); 
      Op.set(Up,UDP,1);
      }
    else
    if(opname == "FermiPhase" || opname == "F")
      {
      Op.set(Em,EmP,+1); 
      Op.set(Up,UpP,-1);
      Op.set(Dn,DnP,-1);
      Op.set(UD,UDP,+1);
      }
    else
    if(opname == "Fup")
      {
      Op.set(Em,EmP,+1); 
      Op.set(Up,UpP,-1);
      Op.set(Dn,DnP,+1);
      Op.set(UD,UDP,-1);
      }
    else
    if(opname == "Fdn")
      {
      Op.set(Em,EmP,+1); 
      Op.set(Up,UpP,+1);
      Op.set(Dn,DnP,-1);
      Op.set(UD,UDP,-1);
      }
    else
    if(opname == "Sz")
      {
      Op.set(Up,UpP,+0.5); 
      Op.set(Dn,DnP,-0.5);
      }
    else
    if(opname == "S+")
      {
      Op.set(Dn,UpP,1); 
      }
    else
    if(opname == "S-")
      {
      Op.set(Up,DnP,1); 
      }
    else
    if(opname == "S2")
      {
      //S dot S on-site
      Op.set(Up,UpP,0.75); 
      Op.set(Dn,DnP,0.75);
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
