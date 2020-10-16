#include "itensor/all.h"

using namespace itensor;

std::tuple<ITensor, Real>
trg(ITensor const& A0,
    int maxdim, int topscale,
    Real cutoff = 0.0)
  {
  auto A = A0;
  auto is = findInds(A, "0");
  auto sh = is(1);
  auto sv = is(2);

  // Keep track of partition function per site, z = Z^(1/N)
  Real z = 1.0;

  for(auto scale : range1(topscale))
    {
    //printfln("\n---------- Scale %d -> %d  ----------",scale-1,scale);

    // Get the upper-left and lower-right tensors
    auto [Fh, Fhp] = factor(A, {prime(sh), prime(sv)}, {sh, sv},
                            {"MaxDim = ", maxdim,
                             "Tags = ", "horiz",
                             "SVDMethod = ", "gesdd",
                             "Cutoff = ", cutoff,
                             "ShowEigs = ", false});

    // Grab the new left Index
    auto sh_new = commonIndex(Fh, Fhp);
    Fhp *= delta(sh_new, prime(sh_new));

    // Get the upper-right and lower-left tensors
    auto [Fv, Fvp] = factor(A, {sh, prime(sv)}, {prime(sh), sv},
                            {"MaxDim = ", maxdim,
                             "Tags = ", "vert",
                             "SVDMethod = ", "gesdd",
                             "Cutoff = ", cutoff,
                             "ShowEigs = ", false});

    // Grab the new up Index
    auto sv_new = commonIndex(Fv, Fvp);
    Fvp *= delta(sv_new, prime(sv_new));

    A = (Fh  * delta(prime(sh), sh)) *
        (Fv  * delta(prime(sv), sv)) *
        (Fhp * delta(sh, prime(sh))) *
        (Fvp * delta(sv, prime(sv)));
    
    // Update the indices
    sh = sh_new;
    sv = sv_new;

    // Normalize the current tensor and keep track of
    // the total normalization
    Real TrA = elt(A * delta(sh, prime(sh)) * delta(sv, prime(sv)));
    A /= TrA;
    z *= pow(TrA, 1.0 / pow(2, scale));

    // If using the dual Ising partition function
    //z *= pow(TrA, 1.0 / pow(2, 1 + scale));
    }

  //printfln("log(Z)/N = %.12f",log(z));

  return std::tuple<ITensor, Real>({A, z});
  }

