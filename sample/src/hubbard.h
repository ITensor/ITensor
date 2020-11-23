#include <itensor/all.h>

using namespace itensor;

AutoMPO
hubbard_2d_ky(SiteSet const& sites, Args const& args)
  {
  auto Nx = args.getInt("Nx");
  auto Ny = args.getInt("Ny");
  auto U = args.getReal("U", 0.0);
  auto t = args.getReal("t", 1.0);
  auto ampo = AutoMPO(sites);
  for(auto x : range(Nx))
    for(auto ky : range(Ny))
      {
      int s = x * Ny + ky + 1;
      double disp = -2 * t * cos((2 * M_PI / Ny) * ky);
      if(std::abs(disp) > 1e-12)
        {
        ampo += disp, "Nup", s;
        ampo += disp, "Ndn", s;
        }
      }
  for(auto x : range(Nx-1))
    for(auto ky : range(Ny))
      {
      int s1 = x * Ny + ky + 1;
      int s2 = (x + 1) * Ny + ky + 1;
      ampo += -t, "Cdagup", s1, "Cup", s2;
      ampo += -t, "Cdagup", s2, "Cup", s1;
      ampo += -t, "Cdagdn", s1, "Cdn", s2;
      ampo += -t, "Cdagdn", s2, "Cdn", s1;
      }
  if(U != 0)
    {
    for(auto x : range(Nx))
      for(auto ky : range(Ny))
        for(auto py : range(Ny))
          for(auto qy : range(Ny))
            {
            int s1 = x * Ny + (ky + qy + Ny) % Ny + 1;
            int s2 = x * Ny + (py - qy + Ny) % Ny + 1;
            int s3 = x * Ny + py + 1;
            int s4 = x * Ny + ky + 1;
            ampo += (U / Ny), "Cdagdn", s1, "Cdagup", s2, "Cup", s3, "Cdn", s4;
            }
    }
  return ampo;
  }

