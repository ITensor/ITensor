#include "src/trg.h"
#include "src/ising.h"
#include "itensor/util/print_macro.h"

int
main()
  {
  Real betac = 0.5 * log(sqrt(2) + 1.0);
  Real beta = 1.1 * betac;
  int maxdim = 20;
  int topscale = 20;
  
  auto dim0 = 2;
  
  // Define an initial Index making up
  // the Ising partition function
  auto s = Index(dim0);
  
  // Define the indices of the scale-0
  // Boltzmann weight tensor "A"
  auto l = addTags(s, "left");
  auto r = addTags(s, "right");
  auto u = addTags(s, "up");
  auto d = addTags(s, "down");
  
  auto A0 = ising(l, r, u, d, beta);
  auto [A, z] = trg(A0, l, r, u, d, maxdim, topscale);
  (void)A;

  printfln("beta = %.12f", beta);
  printfln("beta/betac = %.12f", beta / betac);
  println("maxdim = ", maxdim);
  println("niters = ", topscale);
  printfln("kappa = %.12f", z);

  return 0;
  }

