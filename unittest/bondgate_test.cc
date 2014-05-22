#include "test.h"
#include "bondgate.h"
#include "model/spinhalf.h"

using std::vector;

TEST_CASE("BondGateTest")
{

const int N = 10;
SpinHalf sites(N);

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
