#include "itensor/all.h"

using namespace itensor;

int main(int argc, char *argv[])
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

  auto N = Nx * Ny;
  auto sites = Electron(N);

  auto t = 1.0;

  auto ampo = AutoMPO(sites);
  auto lattice = squareLattice(Nx, Ny, {"YPeriodic = ", true});
  for(auto j : lattice)
      {
      ampo += -t, "Cdagup", j.s1, "Cup", j.s2;
      ampo += -t, "Cdagup", j.s2, "Cup", j.s1;
      ampo += -t, "Cdagdn", j.s1, "Cdn", j.s2;
      ampo += -t, "Cdagdn", j.s2, "Cdn", j.s1;
      }
  for(auto j : range1(N))
      {
      ampo += U, "Nupdn", j;
      }
  auto H = toMPO(ampo);

  auto state = InitState(sites);
  for(auto j : range1(N))
      {
      state.set(j, (j % 2 == 1 ? "Up" : "Dn"));
      }

  auto sweeps = Sweeps(15);
  sweeps.maxdim() = 20, 60, 100, 100, 200, 400, 800, 2000, 3000;
  sweeps.noise() = 1E-7, 1E-8, 1E-10, 0;
  sweeps.cutoff() = 1E-6;

  PrintData(sweeps);

  auto psi0 = randomMPS(state);
  auto [energy,psi] = dmrg(H,psi0,sweeps,{"Quiet=",true});

  PrintData(Nx);
  PrintData(Ny);
  PrintData(U);
  PrintData(t);
  PrintData(totalQN(psi));
  PrintData(maxLinkDim(psi));
  PrintData(energy);


  return 0;
  }
