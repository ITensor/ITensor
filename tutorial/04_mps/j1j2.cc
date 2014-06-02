#include "j1j2.h"

using std::vector;

ITensor
B(const SiteSet& sites, int b)
    {
    ITensor B_ = sites.op("Sz",b)*sites.op("Sz",b+1)
              + 0.5*sites.op("Sp",b)*sites.op("Sm",b+1)
              + 0.5*sites.op("Sm",b)*sites.op("Sp",b+1);
    return B_;
    }

int
main(int argc, char* argv[])
    {
    const int N = 20;

    //Model objects represent a collection of 
    //lattice degrees of freedom of a certain type
    SpinHalf sites(N);

    vector<Real> J2s(1,NAN),
                 dimer(1,NAN);

    for(Real J2 = 0.0; J2 <= 0.5; J2 += 0.05)
        {
        //
        // Compute ground state using black box
        //
        MPS psi = computeGroundState(sites,J2);

        Real val = 0;
        //
        // Write code here that
        // measures the dimer order parameter:
        //
        // B(N/2) - 1/2*B(N/2-1) - 1/2*B(N/2+1)
        //

        J2s.push_back(J2);
        dimer.push_back(val);
        }

    for(int j = 1; j < int(dimer.size()); ++j)
        {
        printfln("%.5f %.10f",J2s.at(j),dimer.at(j));
        }

    return 0;
    }
