#include "src/ctmrg.h"
#include "src/ising.h"
#include "itensor/util/print_macro.h"

int
main()
  {
  Real betac = 0.5 * log(sqrt(2) + 1.0);
  Real beta = 1.1 * betac;
  int maxdim = 20;
  int nsteps = 20;
  
  auto dim0 = 2;
  
  // Define an initial Index making up
  // the Ising partition function
  auto s = Index(dim0, "Site");
  
  // Define the indices of the scale-0
  // Boltzmann weight tensor "A"
  auto sh = addTags(s, "horiz");
  auto sv = addTags(s, "vert");
  
  auto T = ising(sh, sv, beta);

  auto l = Index(1, "Link");
  auto lh = addTags(l, "horiz");
  auto lv = addTags(l, "vert");
  auto Clu0 = ITensor(lv, lh);
  Clu0.set(1, 1, 1.0);
  auto Al0 = ITensor(lv, prime(lv), sh);
  Al0.set(lv = 1, prime(lv) = 1, sh = 1, 1.0);

  auto [Clu, Al] = ctmrg(T, Clu0, Al0, maxdim, nsteps);

  lv = commonIndex(Clu, Al);
  lh = uniqueIndex(Clu, Al);

  auto Au = replaceInds(Al, {lv, prime(lv), sh},
                            {lh, prime(lh), sv});

  auto ACl = Al * Clu * dag(prime(Clu));

  auto ACTl = prime(ACl * dag(prime(Au)) * T * Au, -1);
  auto kappa = elt(ACTl * dag(ACl));

  auto Tsz = ising(sh, sv, beta, true);
  auto ACTszl = prime(ACl * dag(prime(Au)) * Tsz * Au, -1);
  auto m = elt(ACTszl * dag(ACl)) / kappa;

  printfln("beta = %.12f", beta);
  printfln("beta/betac = %.12f", beta / betac);
  println("maxdim = ", maxdim);
  println("niters = ", nsteps);
  printfln("kappa = %.12f", kappa);
  printfln("m = %.12f", m);

  return 0;
  }

