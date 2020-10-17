#include "itensor/all.h"

using namespace itensor;

ITensor
ising(Index const& l, Index const& r,
      Index const& u, Index const& d,
      Real beta)
  {
  int dim0 = 2;
  auto A = ITensor(l, r, u, d);
  auto T = 1/beta;

  // Fill the A tensor with correct Boltzmann weights:
  auto Sig = [](int s) { return 1.-2.*(s-1); };
  for(auto sl : range1(dim0))
  for(auto sd : range1(dim0))
  for(auto sr : range1(dim0))
  for(auto su : range1(dim0))
    {
    auto E = Sig(sl)*Sig(sd)+Sig(sd)*Sig(sr)
            +Sig(sr)*Sig(su)+Sig(su)*Sig(sl);
    auto P = exp(-E/T);
    A.set(l(sl),r(sr),u(su),d(sd),P);
    }
  return A;
  }

