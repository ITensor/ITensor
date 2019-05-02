#include "j1j2.h"

using std::vector;
using namespace itensor;

ITensor
makeB(SiteSet const& sites, int b)
    {
    auto B_ =     op(sites,"Sz",b)*op(sites,"Sz",b+1) +
              0.5*op(sites,"Sp",b)*op(sites,"Sm",b+1) +
              0.5*op(sites,"Sm",b)*op(sites,"Sp",b+1);
    return B_;
    }

int main()
    {
    int N = 20;

    auto sites = SpinHalf(N,{"ConserveQNs=",false});

    auto J2s = vector<Real>();
    auto dimer = vector<Real>();

    for(Real J2 = 0.0; J2 <= 0.5; J2 += 0.05)
        {
        //
        // Compute ground state using
        // "black box" routine (or have a look at j1j2.h)
        //
        auto psi = computeGroundState(sites,J2);

        Real D = 0.;
        //
        // We will add code below to
        // measure the dimer order parameter:
        //
        // D = B_(N/2) - 1/2*B_(N/2-1) - 1/2*B_(N/2+1)
        //
        // Tip: to prime only physical indices
        // do prime(...,Site)
        //

        //
        // Compute:  <B_(N/2)>
        //
        auto i1 = N/2;
        psi.position(i1);
        auto B1 = makeB(sites,i1);
        auto wf1 = psi(i1)*psi(i1+1);
        // compute <wf1|B1|wf1>
        D +=  elt(dag(prime(wf1,"Site")) * B1 * wf1);

        //
        // Compute: -1/2 <B_(N/2-1)>
        //
        auto i2 = N/2-1;
        psi.position(i2);
        auto B2 = makeB(sites,i2);
        auto wf2 = psi(i2)*psi(i2+1);
        //TODO: ADD CODE to compute <wf2|B2|wf2>
        //      replacing the ... below
        //D += (-0.5) * ...

        //
        // Compute: -1/2 <B_(N/2+1)>
        //
        auto i3 = N/2+1;
        psi.position(i3);
        auto B3 = makeB(sites,i3);
        auto wf3 = psi(i3)*psi(i3+1);
        //TODO: ADD CODE to compute <wf3|B3|wf3>
        //      replacing the ... below
        //D +=  (-0.5) * ...


        J2s.push_back(J2);
        dimer.push_back(D);
        }

    for(auto j : range(J2s))
        {
        printfln("%.5f %.10f",J2s.at(j),dimer.at(j));
        }

    return 0;
    }
