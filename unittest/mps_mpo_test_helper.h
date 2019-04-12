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
  auto norm_psi = inner(psi,psi);
  auto norm_phi = inner(phi,phi);
  auto psi_phi = innerC(psi,phi);
  return std::sqrt(std::abs(norm_psi/norm_phi+1.-2.*real(psi_phi)/norm_phi));
  }

MPO inline
randomMPO(SiteSet const& sites, Args const& args = Args::global())
  {
  auto N = length(sites);
  //Use AutoMPO as a trick to get
  //an MPO with bond dimension > 1
  auto ampo = AutoMPO(sites);
  for(auto j : range1(N-1))
      {
      ampo += "Sz",j,"Sz",j+1;
      ampo += 0.5,"S+",j,"S-",j+1;
      ampo += 0.5,"S-",j,"S+",j+1;
      }
  auto H = toMPO(ampo);
  //Randomize the MPO
  for(auto j : range1(N))
      {
      H.ref(j).randomize(args);
      H.ref(j) /= norm(H(j));
      }
  return H;
  }

}

#endif
