#include <itensor/all.h>
#include "sample/src/electronk.h"
#include "sample/src/hubbard.h"

using namespace itensor;

int
main(int argc, char *argv[])
  {
  int Nx = 6;
  int Ny = 3;
  double U = 4.0;
  if(argc > 3)
    U = std::stof(argv[3]);
  if(argc > 2)
    Ny = std::stoi(argv[2]);
  if(argc > 1)
    Nx = std::stoi(argv[1]);

  double t = 1.0;

  auto args = Args("Kmod", Ny);
  args.add("ConserveQNs", true);
  args.add("ConserveK", true);
  int N = Nx * Ny;
  
  auto sweeps = Sweeps(15);
  sweeps.maxdim() = 100, 200, 400, 800, 2000, 3000;
  sweeps.cutoff() = 1e-6;
  sweeps.noise() = 1e-6, 1e-7, 1e-8, 0.0;

  SiteSet sites = ElectronK(N, args);
  auto ampo = hubbard_2d_ky(sites, {"Nx = ", Nx, "Ny = ", Ny, "t = ", t, "U = ", U});
  auto H = toMPO(ampo);

  // Create start state
  auto state = InitState(sites);
  for (auto i : range1(N))
    {
    int x = (i-1)/Ny;
    int y = (i-1)%Ny;
    if(x%2==0)
      {
      if(y%2==0) state.set(i,"Up");
      else        state.set(i,"Dn");
      }
    else
      {
      if(y%2==0) state.set(i,"Dn");
      else        state.set(i,"Up");
      }
    }

  PrintData(sweeps);

  auto psi0 = randomMPS(state);
  auto [energy, psi] = dmrg(H, psi0, sweeps, {"Quiet = ", true});

  PrintData(Nx);
  PrintData(Ny);
  PrintData(U);
  PrintData(t);
  PrintData(totalQN(psi));
  PrintData(maxLinkDim(psi));
  PrintData(energy);

  return 0;
  }
