#include <itensor/all.h>
#include "sample/src/electronk.h"

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
  sweeps.maxdim() = 20, 60, 100, 100, 200, 400, 800, 2000, 3000;
  sweeps.cutoff() = 1e-6;
  sweeps.noise() = 1e-6, 1e-7, 1e-8, 0.0, 1e-7, 0.0, 1e-6, 0.0, 1e-6, 0.0;

  /////////////////////////////////
  // K-space Hamiltonian

  SiteSet sitesk = ElectronK(N, args);

  auto ampo = AutoMPO(sitesk);
  // hopping in y-direction
  for(auto x : range(Nx))
    for(auto ky : range(Ny))
      {
      int s = x * Ny + ky + 1;
      double disp = -2 * t * cos((2 * M_PI / Ny) * ky);
      if(std::abs(disp)>1e-12)
        {
        ampo += disp,"Nup",s; 
        ampo += disp,"Ndn",s;
        }
      }

  // hopping in x-direction
  for(auto x : range(Nx-1))
    for(auto ky : range(Ny))
      {
      int s1 = x*Ny + ky + 1;
      int s2 = (x+1)*Ny + ky + 1;
      ampo += -t, "Cdagup", s1, "Cup", s2; 
      ampo += -t, "Cdagup", s2, "Cup", s1;
      ampo += -t, "Cdagdn", s1, "Cdn", s2;
      ampo += -t, "Cdagdn", s2, "Cdn", s1;
      }

  // Hubbard interaction
  for(auto x : range(Nx))
    for(auto ky : range(Ny))
      for(auto py : range(Ny))
        for(auto qy : range(Ny))
          {
          int s1 = x*Ny + (ky+qy+Ny)%Ny + 1;
          int s2 = x*Ny + (py-qy+Ny)%Ny + 1;
          int s3 = x*Ny + py + 1;
          int s4 = x*Ny + ky + 1;
          //ampo += (U/Ny), "Cdagdn", s1, "Cdagup", s2, "Cup", s3, "Cdn", s4;
          if(s1 == s4 && s2 == s3)
            ampo += (U/Ny), "Ndn", s1, "Nup", s2;
          else if(s1 == s4)
            ampo += (U/Ny), "Ndn", s1, "Cdagup", s2, "Cup", s3;
          else if(s2 == s3)
            ampo += (U/Ny), "Cdagdn", s1, "Cdn", s4, "Nup", s2;
          else
            ampo += (U/Ny), "Cdagdn", s1, "Cdagup", s2, "Cup", s3, "Cdn", s4;
          }
  auto H = toMPO(ampo);

  // Create start state
  auto state = InitState(sitesk);
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
