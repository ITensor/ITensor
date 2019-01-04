#include "test.h"
#include "itensor/iterativesolvers.h"
#include "sample/Heisenberg.h"
#include "itensor/mps/sites/spinhalf.h"
#include "itensor/mps/localmpo.h"

using namespace itensor;
using namespace std;

class ITensorMap
    {
    ITensor const* A_;
    mutable long size_;
    public:

    ITensorMap(ITensor const& A)
      : A_(nullptr),
        size_(-1)
        {
        A_ = &A;
        }

    void
    product(ITensor const& x, ITensor& b) const
        {
        b = *A_*x;
        b.mapPrime(1,0);
        }

    long
    size() const
        {
        if(size_ == -1)
            {
            size_ = 1;
            for(auto& I : A_->inds())
                {
                if(I.primeLevel() > 0)
                    size_ *= I.m();
                }
            }
        return size_;
        }

    };

TEST_CASE("EigenSolverTest")
{

SECTION("FourSite")
    {
    //Exact 4 site energy is -1.6160254038 from DMRG

    const int N = 4;
    SpinHalf sites(N);
    MPO H = Heisenberg(sites);

    InitState initState(sites);
    for(int i = 1; i <= N; ++i)
        initState.set(i,i%2==1 ? "Up" : "Dn");

    MPS psi(initState);

    LocalMPO PH(H);

    ITensor phip;
    psi.position(2);
    PH.position(2,psi);

    ITensor phi1 = psi.A(2) * psi.A(3);

    Real En1 = davidson(PH,phi1,"MaxIter=9");
    CHECK_CLOSE(En1,-0.95710678118);

    cout << endl << endl;

    //TODO: should this be put back in?
    //ITensor mpoh = H.A(2)*H.A(3);
    //ITensor phi2 = psi.A(2)*psi.A(3);
    //Real En2 = doDavidson(phi2,mpoh,PH.L(),PH.R(),9,2,1E-4);
    //cout << format("Energy from matrix Davidson (b=2) = %.20f")%En2 << endl;

    //cout << endl << endl;

    //cout << endl << "Overlap = " << Dot(phi1,phi2) << endl;

    //cout << "---------------------------------------" << endl << endl;

    //psi.doSVD(2,phi1,Fromleft);
    //psi.position(3);
    //PH.position(3,psi);
    //ITensor phi3 = psi.A(3) * psi.A(4);
    //Real En3 = d.davidson(PH,phi3);
    //cout << format("Energy from tensor Davidson (b=3) = %.20f")%En3 << endl;

    //cout << endl << endl;

    //mpoh = H.A(3)*H.A(4);
    //ITensor phi3m = psi.A(3)*psi.A(4);
    //Real En3m = doDavidson(phi3m,mpoh,PH.L(),PH.R(),9,2,1E-4);
    //cout << format("Energy from matrix Davidson (b=3) = %.20f")%En3m << endl;

    //cout << "---------------------------------------" << endl << endl;

    //psi.doSVD(3,phi3,Fromright);
    //psi.position(3);
    //PH.position(2,psi);
    //ITensor phi4 = psi.A(2) * psi.A(3);
    //Real En4 = d.davidson(PH,phi4);
    //cout << format("Energy from tensor Davidson (b=2) = %.20f")%En4 << endl;

    //cout << endl << endl;

    //mpoh = H.A(2)*H.A(3);
    //ITensor phi4m = psi.A(2)*psi.A(3);
    //Real En4m = doDavidson(phi4m,mpoh,PH.L(),PH.R(),9,2,1E-4);
    //cout << format("Energy from matrix Davidson (b=2) = %.20f")%En4m << endl;

    //cout << "---------------------------------------" << endl << endl;

    //psi.doSVD(2,phi4,Fromright);
    //psi.position(1);
    //PH.position(1,psi);

    ////With doDavidson
    //mpoh = H.A(1)*H.A(2);
    //ITensor phi5 = psi.A(1) * psi.A(2);
    //Real En5 = doDavidson(phi5,mpoh,PH.L(),PH.R(),9,2,1E-4);
    //cout << format("Energy from matrix Davidson (b=1) = %.20f")%En5 << endl;

    //cout << endl << endl;

    //ITensor phi6 = psi.A(1) * psi.A(2);
    ////Print(phi6);
    ////Print(PH.L());
    ////Print(PH.R());
    ////ITensor AB = phi6 * PH.R();
    ////AB *= H.A(2);
    ////AB *= H.A(1);
    ////AB.noprime();
    ////ITensor AB; PH.product(phi6,AB);
    ////Print(Dot(phi6,AB));
    //Real En6 = d.davidson(PH,phi6);
    //cout << format("Energy from tensor Davidson (b=1) = %.20f")%En6 << endl;

    }

SECTION("IQFourSite")
    {
    //Exact 4 site energy is -1.6160254038 from DMRG

    const int N = 4;
    SpinHalf sites(N);
    MPO H = Heisenberg(sites);

    InitState initState(sites);
    for(int i = 1; i <= N; ++i)
        initState.set(i,i%2==1 ? "Up" : "Dn");

    MPS psi(initState);

    LocalMPO PH(H);

    ITensor phip;
    psi.position(2);
    PH.position(2,psi);

    ITensor phi1 = psi.A(2) * psi.A(3);

    Real En1 = davidson(PH,phi1,"MaxIter=9");
    //cout << format("Energy from tensor Davidson (b=2) = %.20f")%En1 << endl;
    CHECK_CLOSE(En1,-0.95710678118);


    }

SECTION("GMRES (ITensor, Real)")
    {
    auto a1 = Index(3,"Site,a1");
    auto a2 = Index(4,"Site,a2");
    auto a3 = Index(3,"Site,a3");

    auto A = randomITensor(prime(a1),prime(a2),prime(a3),a1,a2,a3);
    auto x = randomITensor(a1,a2,a3);
    auto b = randomITensor(a1,a2,a3);

    // ITensorMap is defined above, it simply wraps an ITensor that is of the
    // form of a matrix (i.e. has indices of the form {i,j,k,...,i',j',k',...})
    gmres(ITensorMap(A),b,x,{"MaxIter",36,"ErrGoal",1e-10});

    CHECK_CLOSE(norm((A*x).mapPrime(1,0)-b)/norm(b),0.0);

    }

SECTION("GMRES (ITensor, Complex)")
    {
    auto a1 = Index(3,"Site,a1");
    auto a2 = Index(4,"Site,a2");
    auto a3 = Index(3,"Site,a3");

    auto A = randomITensorC(prime(a1),prime(a2),prime(a3),a1,a2,a3);
    auto x = randomITensorC(a1,a2,a3);
    auto b = randomITensorC(a1,a2,a3);

    gmres(ITensorMap(A),b,x,{"MaxIter",100,"ErrGoal",1e-10});

    CHECK_CLOSE(norm((A*x).mapPrime(1,0)-b)/norm(b),0.0);

    }

SECTION("GMRES (ITensor, QN)")
    {
    auto i = Index(QN(+1),5,
                   QN(0),5,
                   QN(-1),5);
    auto j = Index(QN(+1),5,
                   QN(0),5,
                   QN(-1),5,
                   In);
    
    auto A = randomITensor(QN(0),prime(dag(i)),prime(dag(j)),i,j);

    auto x = randomITensor(QN(0),dag(i),dag(j));
    auto b = randomITensor(QN(0),dag(i),dag(j));

    gmres(ITensorMap(A),b,x,{"MaxIter",100,"DebugLevel",0,"ErrGoal",1e-10});

    CHECK_CLOSE(norm((A*x).mapPrime(1,0)-b)/norm(b),0.0);

    }

}
