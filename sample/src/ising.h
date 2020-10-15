#include "itensor/all.h"

using namespace itensor;

ITensor
ising(Index const& sh, Index const& sv,
      Real beta)
  {
  int dim0 = 2;
  auto A = ITensor(sh, prime(sh), sv, prime(sv));

  // Fill the A tensor with correct Boltzmann weights:
  auto Sig = [](int s) { return 1. - 2. * (s - 1); };
  for(auto ssh  : range1(dim0))
  for(auto ssvp : range1(dim0))
  for(auto sshp : range1(dim0))
  for(auto ssv  : range1(dim0))
    {
    auto E = Sig(ssh)  * Sig(ssvp) +
             Sig(ssvp) * Sig(sshp) +
             Sig(sshp) * Sig(ssv)  +
             Sig(ssv)  * Sig(ssh);
    auto P = exp(-beta * E);
    A.set(sh = ssh, prime(sh) = sshp, sv = ssv, prime(sv) = ssvp, P);
    }
  return A;
  }

