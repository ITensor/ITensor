#include "itensor/all.h"

using namespace itensor;

ITensor
ising(Index const& sh,
      Index const& sv,
      Real beta,
      bool sz = false)
  {
  auto d = dim(sh);
  auto shp = prime(sh);
  auto svp = prime(sv);
  auto T = ITensor(sh, shp, sv, svp);
  for(auto i : range1(d)){ T.set(i, i, i, i, 1.0); }
  if(sz) { T.set(1, 1, 1, 1, -1.0); }
  auto th = sim(sh);
  auto thp = sim(shp);
  auto tv = sim(sv);
  auto tvp = sim(svp);
  auto Tp = T * delta(sh, th) *
                delta(shp, thp) *
                delta(sv, tv) *
                delta(svp, tvp);
  // Analytic square root of the bond matrix:
  // [exp( beta) exp(-beta)
  //  exp(-beta) exp( beta)]
  auto lambda_p = sqrt(exp(beta) + exp(-beta));
  auto lambda_m = sqrt(exp(beta) - exp(-beta));
  auto x_p = 0.5 * (lambda_p + lambda_m);
  auto x_m = 0.5 * (lambda_p - lambda_m);
  auto Xh = ITensor(th, sh);
  Xh.set(1, 1, x_p);
  Xh.set(2, 1, x_m);
  Xh.set(1, 2, x_m);
  Xh.set(2, 2, x_p);
  auto Xhp = replaceInds(Xh, {th, sh}, {thp, shp});
  auto Xv  = replaceInds(Xh, {th, sh}, {tv,  sv });
  auto Xvp = replaceInds(Xh, {th, sh}, {tvp, svp});
  return Tp * Xhp * Xvp * Xh * Xv;
  }

//
// Dual partition function
//

//ITensor
//ising(Index const& sh, Index const& sv,
//      Real beta)
//  {
//  int dim0 = 2;
//  auto T = ITensor(sh, prime(sh), sv, prime(sv));
//
//  // Fill the T tensor with correct Boltzmann weights:
//  auto Sig = [](int s) { return 1. - 2. * (s - 1); };
//  for(auto ssh  : range1(dim0))
//  for(auto ssvp : range1(dim0))
//  for(auto sshp : range1(dim0))
//  for(auto ssv  : range1(dim0))
//    {
//    auto E = Sig(ssh)  * Sig(ssvp) +
//             Sig(ssvp) * Sig(sshp) +
//             Sig(sshp) * Sig(ssv)  +
//             Sig(ssv)  * Sig(ssh);
//    auto P = exp(-beta * E);
//    T.set(sh = ssh, prime(sh) = sshp, sv = ssv, prime(sv) = ssvp, P);
//    }
//  return T;
//  }

