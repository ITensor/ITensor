#ifndef __ITENSOR_MPS_MPO_TEST_HELPER_H
#define __ITENSOR_MPS_MPO_TEST_HELPER_H
#include "itensor/mps/autompo.h"

namespace itensor {

bool inline
checkTags(MPO const& K,
          std::string const& siteTagsBra = "Site,0",
          std::string const& siteTagsKet = "Site,1",
          std::string const& linkTags = "Link,0")
  {
  auto N = length(K);
  for(auto n : range1(N))
      {
      if( n < N )
        {
        auto l = linkIndex(K,n);
        if( !hasTags(l,linkTags+",l="+str(n)) ) return false;
        }
      auto s = siteInds(K,n);
      if( !findIndex(s,siteTagsBra+",n="+str(n)) ) return false;
      if( !findIndex(s,siteTagsKet+",n="+str(n)) ) return false;
      }
  return true;
  }

bool inline
checkTags(MPS const& psi,
          std::string const& siteTags = "Site,0",
          std::string const& linkTags = "Link,0")
  {
  auto N = length(psi);
  for(auto n : range1(N))
      {
      if( n < N )
        {
        auto l = linkIndex(psi,n);
        if( !hasTags(l,linkTags+",l="+str(n)) ) return false;
        }
      auto s = siteIndex(psi,n);
      if( !hasTags(s,siteTags+",n="+str(n)) ) return false;
      }
  return true;
  }

Real inline
diff(MPS const& psi, MPS const& phi)
  {
  auto norm_psi = real(innerC(psi,psi));
  auto norm_phi = real(innerC(phi,phi));
  auto psi_phi = innerC(psi,phi);
  return std::sqrt(std::abs(norm_psi/norm_phi+1.-2.*real(psi_phi)/norm_phi));
  }

MPO inline
randomUnitaryMPO(SiteSet const& sites, Args const& args = Args::global())
  {
  auto tau = args.getReal("Timestep",0.1);
  auto N = length(sites);
  //Use AutoMPO as a trick to get
  //an MPO with bond dimension > 1
  auto ampo = AutoMPO(sites);
  auto J1 = detail::quickran();
  auto J2 = detail::quickran();
  auto J3 = detail::quickran();
  for(auto j : range1(N-1))
      {
      ampo += J1,"Sz",j,"Sz",j+1;
      ampo += 0.5*J2,"S+",j,"S-",j+1;
      ampo += 0.5*J3,"S-",j,"S+",j+1;
      }
  auto U = toExpH(ampo,tau*Cplx_i);
  return U;
  }

}

#endif
