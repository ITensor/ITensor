#include "j1j2.h"

using std::vector;
using namespace itensor;

ITensor
B(SiteSet const& sites, int b)
    {
    ITensor B_ = sites.op("Sz",b)*sites.op("Sz",b+1)
              + 0.5*sites.op("Sp",b)*sites.op("Sm",b+1)
              + 0.5*sites.op("Sm",b)*sites.op("Sp",b+1);
    return B_;
    }

int main()
    {
    int N = 20;

    auto sites = SpinHalf(N);

    auto J2s = vector<Real>();
    auto dimer = vector<Real>();

    for(Real J2 = 0.0; J2 <= 0.5; J2 += 0.05)
        {
        //
        // Compute ground state using
        // "black box" routine (or have a look at j1j2.h)
        //
        MPS psi = computeGroundState(sites,J2);

        Real val = 0;
        //
        // TODO: ADD CODE here that
        // measures the dimer order parameter:
        //
        // B(N/2) - 1/2*B(N/2-1) - 1/2*B(N/2+1)
        //

        J2s.push_back(J2);
        dimer.push_back(val);
        }

    for(auto j : range(J2s))
        {
        printfln("%.5f %.10f",J2s.at(j),dimer.at(j));
        }

    return 0;
    }
