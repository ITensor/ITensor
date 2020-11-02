#include "itensor/all.h"

using namespace itensor;

std::tuple<ITensor, ITensor>
ctmrg(ITensor T,
      ITensor Clu,
      ITensor Al,
      int maxdim,
      int nsteps,
      Real cutoff = 0.0)
  {
  auto sh = commonIndex(T, Al);
  auto sv = uniqueIndex(T, {Al, prime(Al)}, "0");
  auto lv = commonIndex(Clu, Al);
  auto lh = uniqueIndex(Clu, Al);
  auto Au = replaceInds(Al, {lv, prime(lv), sh},
                            {lh, prime(lh), sv});
  ITensor Uv;
  for(auto i : range1(nsteps))
    {
    // Get the grown corner transfer matrix (CTM)
    auto Clu_new = Al * Clu * Au * T;

    Clu_new.noPrime();
    Clu_new.replaceInds({lh, sh}, {prime(lv), prime(sv)});

    // Diagonalize the grown CTM
    std::tie(Uv, Clu) = diagPosSemiDef(Clu_new,
                                       {"MaxDim = ", maxdim,
                                        "Tags = ", tags(lv),
                                        "Cutoff = ", cutoff});

    lv = commonIndex(Clu, Uv);
    lh = setTags(lv, tags(lh));

    Clu = toDense(Clu);
    Clu.replaceInds({prime(lv)}, {lh});

    // The renormalized CTM is the diagonal matrix of eigenvalues
    // Normalize the CTM
    auto Cl = Clu * prime(dag(Clu), lh);
    auto normC = pow(elt(Cl * dag(Cl)), 0.25);
    Clu /= normC;

    // Calculate the renormalized half row transfer matrix (HRTM)
    Uv.noPrime();
    Al = Al * Uv * T * dag(prime(Uv));
    Al = replaceInds(Al, {prime(sh)}, {sh});

    // Normalize the HRTM
    auto ACl = Al * Clu * prime(dag(Clu));
    auto normA = sqrt(elt(ACl * dag(ACl)));
    Al /= normA;
    Au = replaceInds(Al, {lv, prime(lv), sh},
                         {lh, prime(lh), sv});
    }
  return std::tuple<ITensor, ITensor>({Clu, Al});
  }

