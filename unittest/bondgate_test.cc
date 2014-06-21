#include "test.h"
#include "bondgate.h"
#include "sites/spinhalf.h"

using std::vector;
using namespace itensor;

TEST_CASE("BondGateTest")
{

const int N = 10;
SpinHalf sites(N);

SECTION("Site Accessors")
    {
    Gate g1(sites,1,2);
    CHECK(g1.i1() < g1.i2());

    Gate g2(sites,2,1);
    CHECK(g2.i1() < g2.i2());

    const int s1 = 3, s2 = 4;
    ITensor hh =  sites.op("Sz",s1)*sites.op("Sz",s2);
    hh += sites.op("Sm",s1)*sites.op("Sp",s2) * 0.5;
    hh += sites.op("Sp",s1)*sites.op("Sm",s2) * 0.5;

    Gate g3(sites,s1,s2,Gate::tImag,0.1,hh);
    CHECK(g3.i1() < g3.i2());

    Gate g4(sites,s2,s1,Gate::tImag,0.1,hh);
    CHECK(g4.i1() < g4.i2());
    }

SECTION("SwapGate")
    {
    Gate sw12(sites,1,2);
    ITensor id = multSiteOps(sw12.gate(),sw12.gate());
    CHECK((id-sites.op("Id",1)*sites.op("Id",2)).norm() < 1E-12);

    IQGate qsw12(sites,1,2);
    IQTensor qid = multSiteOps(qsw12.gate(),qsw12.gate());
    CHECK((qid-sites.op("Id",1)*sites.op("Id",2)).norm() < 1E-12);
    }

SECTION("ImagTimeGates")
    {
    vector<IQTensor> H(N+1);
    for(int b = 1; b < N; ++b)
        {
        const int s1 = b;
        const int s2 = s1+1;
        IQTensor& hh = H.at(b);

        //Heisenberg model
        hh =  sites.op("Sz",s1)*sites.op("Sz",s2) * 1;
        hh += sites.op("Sm",s1)*sites.op("Sp",s2) * 0.5;
        hh += sites.op("Sp",s1)*sites.op("Sm",s2) * 0.5;
        }

    vector<IQGate> gates;
    const IQGate::Type type = IQGate::tImag;

    const Real tau = 0.01;
    for(int b = 1; b < N; ++b)
        {
        const int s1 = b;
        const int s2 = s1+1;
        const IQTensor& hh = H.at(b);
        gates.push_back(IQGate(sites,s1,s2,type,tau/2.,hh));
        }
    }

SECTION("RealTimeGate")
    {
    vector<IQTensor> H(N+1);
    for(int b = 1; b < N; ++b)
        {
        const int s1 = b;
        const int s2 = s1+1;
        IQTensor& hh = H.at(b);

        //Heisenberg model
        hh =  sites.op("Sz",s1)*sites.op("Sz",s2) * 1;
        hh += sites.op("Sm",s1)*sites.op("Sp",s2) * 0.5;
        hh += sites.op("Sp",s1)*sites.op("Sm",s2) * 0.5;
        }

    vector<IQGate> gates;
    const IQGate::Type type = IQGate::tReal;

    const Real tau = 0.01;
    for(int b = 1; b < N; ++b)
        {
        const int s1 = b;
        const int s2 = s1+1;
        const IQTensor& hh = H.at(b);
        gates.push_back(IQGate(sites,s1,s2,type,tau/2.,hh));
        }
    }

}
