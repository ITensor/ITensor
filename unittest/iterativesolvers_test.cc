#include "test.h"
#include "itensor/iterativesolvers.h"
#include "sample/Heisenberg.h"
#include "itensor/mps/sites/spinhalf.h"
#include "itensor/mps/localmpo.h"
#include "itensor/mps/autompo.h"

using namespace itensor;
using namespace std;

template <class Tensor>
class TensorMap
    {
    Tensor const* A_;
    mutable long size_;
    public:

    TensorMap(Tensor const& A)
      : A_(nullptr),
        size_(-1)
        {
        A_ = &A;
        }

    void
    product(Tensor const& x, Tensor& b) const
        {
        b = *A_*x;
        b.mapprime(1,0);
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

using ITensorMap = TensorMap<ITensor>;
using IQTensorMap = TensorMap<IQTensor>;

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

    LocalMPO<ITensor> PH(H);

    ITensor phip;
    psi.position(2);
    PH.position(2,psi);

    ITensor phi1 = psi.A(2) * psi.A(3);

    Real En1 = davidson(PH,phi1,"MaxIter=9");
    CHECK_CLOSE(En1,-0.95710678118);

    cout << endl << endl;
    /*

    ITensor mpoh = H.A(2)*H.A(3);
    ITensor phi2 = psi.A(2)*psi.A(3);
    Real En2 = doDavidson(phi2,mpoh,PH.L(),PH.R(),9,2,1E-4);
    cout << format("Energy from matrix Davidson (b=2) = %.20f")%En2 << endl;

    cout << endl << endl;

    cout << endl << "Overlap = " << Dot(phi1,phi2) << endl;

    cout << "---------------------------------------" << endl << endl;

    psi.doSVD(2,phi1,Fromleft);
    psi.position(3);
    PH.position(3,psi);
    ITensor phi3 = psi.A(3) * psi.A(4);
    Real En3 = d.davidson(PH,phi3);
    cout << format("Energy from tensor Davidson (b=3) = %.20f")%En3 << endl;

    cout << endl << endl;

    mpoh = H.A(3)*H.A(4);
    ITensor phi3m = psi.A(3)*psi.A(4);
    Real En3m = doDavidson(phi3m,mpoh,PH.L(),PH.R(),9,2,1E-4);
    cout << format("Energy from matrix Davidson (b=3) = %.20f")%En3m << endl;

    cout << "---------------------------------------" << endl << endl;

    psi.doSVD(3,phi3,Fromright);
    psi.position(3);
    PH.position(2,psi);
    ITensor phi4 = psi.A(2) * psi.A(3);
    Real En4 = d.davidson(PH,phi4);
    cout << format("Energy from tensor Davidson (b=2) = %.20f")%En4 << endl;

    cout << endl << endl;

    mpoh = H.A(2)*H.A(3);
    ITensor phi4m = psi.A(2)*psi.A(3);
    Real En4m = doDavidson(phi4m,mpoh,PH.L(),PH.R(),9,2,1E-4);
    cout << format("Energy from matrix Davidson (b=2) = %.20f")%En4m << endl;

    cout << "---------------------------------------" << endl << endl;

    psi.doSVD(2,phi4,Fromright);
    psi.position(1);
    PH.position(1,psi);

    //With doDavidson
    mpoh = H.A(1)*H.A(2);
    ITensor phi5 = psi.A(1) * psi.A(2);
    Real En5 = doDavidson(phi5,mpoh,PH.L(),PH.R(),9,2,1E-4);
    cout << format("Energy from matrix Davidson (b=1) = %.20f")%En5 << endl;

    cout << endl << endl;

    ITensor phi6 = psi.A(1) * psi.A(2);
    //Print(phi6);
    //Print(PH.L());
    //Print(PH.R());
    //ITensor AB = phi6 * PH.R();
    //AB *= H.A(2);
    //AB *= H.A(1);
    //AB.noprime();
    //ITensor AB; PH.product(phi6,AB);
    //Print(Dot(phi6,AB));
    Real En6 = d.davidson(PH,phi6);
    cout << format("Energy from tensor Davidson (b=1) = %.20f")%En6 << endl;
*/

    }

SECTION("IQFourSite")
    {
    //Exact 4 site energy is -1.6160254038 from DMRG

    const int N = 4;
    SpinHalf sites(N);
    IQMPO H = Heisenberg(sites);

    InitState initState(sites);
    for(int i = 1; i <= N; ++i)
        initState.set(i,i%2==1 ? "Up" : "Dn");

    IQMPS psi(initState);

    LocalMPO<IQTensor> PH(H);

    IQTensor phip;
    psi.position(2);
    PH.position(2,psi);

    IQTensor phi1 = psi.A(2) * psi.A(3);

    Real En1 = davidson(PH,phi1,"MaxIter=9");
    //cout << format("Energy from tensor Davidson (b=2) = %.20f")%En1 << endl;
    CHECK_CLOSE(En1,-0.95710678118);


    }

SECTION("Arnoldi (No QN)")
    {
    const int N = 50;
    const int Nc = N/2;
    SpinHalf sites(N);

    // Create random MPO
    auto ampo = AutoMPO(sites);
    for(int j = 1; j < N; ++j)
        {
        ampo += 0.5,"S+",j,"S-",j+1;
        ampo += 0.5,"S-",j,"S+",j+1;
        ampo +=     "Sz",j,"Sz",j+1;
        }
    auto H = MPO(ampo);
    for(int j = 1; j <= N; ++j)
      {
      randomize(H.Aref(j));
      H.Aref(j) /= norm(H.A(j));
      }
    auto normH = sqrt(overlap(H,H));
    for(int j = 1; j <= N; ++j)
      H.Aref(j) /= pow(normH,1.0/N);

    // Create random starting MPS
    auto initState = InitState(sites);
    for(int i = 1; i <= N; ++i)
        initState.set(i,i%2==1 ? "Up" : "Dn");
    auto psi = MPS(initState);
    for(int j = 1; j < 2; ++j)
      {
      psi = applyMPO(H,psi);
      normalize(psi);
      }
    for(int j = 1; j <= N; ++j)
      {
      randomize(psi.Aref(j));
      psi.Aref(j) /= norm(psi.A(j));
      }
    psi.position(Nc);
    normalize(psi);

    LocalMPO<ITensor> PH(H);

    psi.position(Nc);
    PH.position(Nc,psi);

    auto x = psi.A(Nc) * psi.A(Nc+1);
    randomize(x);

    auto lambda = arnoldi(PH,x,{"MaxIter",20,"ErrGoal",1e-14,"DebugLevel",0,"WhichEig","LargestMagnitude"});
    auto PHx = x;
    PH.product(x,PHx);

    CHECK_CLOSE(norm(PHx-lambda*x)/norm(PHx),0.0);
    }

SECTION("Arnoldi (QN)")
    {
    const int N = 50;
    const int Nc = N/2;
    SpinHalf sites(N);

    // Create random MPO
    auto ampo = AutoMPO(sites);
    for(int j = 1; j < N; ++j)
        {
        ampo += 0.5,"S+",j,"S-",j+1;
        ampo += 0.5,"S-",j,"S+",j+1;
        ampo +=     "Sz",j,"Sz",j+1;
        }
    auto H = IQMPO(ampo);
    for(int j = 1; j <= N; ++j)
      {
      randomize(H.Aref(j));
      H.Aref(j) /= norm(H.A(j));
      }
    auto normH = sqrt(overlap(H,H));
    for(int j = 1; j <= N; ++j)
      H.Aref(j) /= pow(normH,1.0/N);

    // Create random starting MPS
    auto initState = InitState(sites);
    for(int i = 1; i <= N; ++i)
        initState.set(i,i%2==1 ? "Up" : "Dn");
    auto psi = IQMPS(initState);
    for(int j = 1; j < 2; ++j)
      {
      psi = applyMPO(H,psi);
      normalize(psi);
      }
    for(int j = 1; j <= N; ++j)
      {
      randomize(psi.Aref(j));
      psi.Aref(j) /= norm(psi.A(j));
      }
    psi.position(Nc);
    normalize(psi);

    LocalMPO<IQTensor> PH(H);

    psi.position(Nc);
    PH.position(Nc,psi);

    auto x = psi.A(Nc) * psi.A(Nc+1);
    randomize(x);

    auto lambda = arnoldi(PH,x,{"MaxIter",20,"ErrGoal",1e-14,"DebugLevel",0,"WhichEig","LargestMagnitude"});
    auto PHx = x;
    PH.product(x,PHx);

    CHECK_CLOSE(norm(PHx-lambda*x)/norm(PHx),0.0);

    }
SECTION("GMRES (ITensor, Real)")
    {
    auto a1 = Index("a1",3,Site);
    auto a2 = Index("a2",4,Site);
    auto a3 = Index("a3",3,Site);

    auto A = randomTensor(prime(a1),prime(a2),prime(a3),a1,a2,a3);
    auto x = randomTensor(a1,a2,a3);
    auto b = randomTensor(a1,a2,a3);

    // ITensorMap is defined above, it simply wraps an ITensor that is of the
    // form of a matrix (i.e. has indices of the form {i,j,k,...,i',j',k',...})
    gmres(ITensorMap(A),b,x,{"MaxIter",36,"ErrGoal",1e-10});

    CHECK_CLOSE(norm((A*x).mapprime(1,0)-b)/norm(b),0.0);

    }

SECTION("GMRES (ITensor, Complex)")
    {
    auto a1 = Index("a1",3,Site);
    auto a2 = Index("a2",4,Site);
    auto a3 = Index("a3",3,Site);

    auto A = randomTensorC(prime(a1),prime(a2),prime(a3),a1,a2,a3);
    auto x = randomTensorC(a1,a2,a3);
    auto b = randomTensorC(a1,a2,a3);

    gmres(ITensorMap(A),b,x,{"MaxIter",100,"ErrGoal",1e-10});

    CHECK_CLOSE(norm((A*x).mapprime(1,0)-b)/norm(b),0.0);

    }

SECTION("GMRES (IQTensor)")
    {
    auto i = IQIndex("i",Index("+1",5),QN(+1),
                          Index("0",5),QN(0),
                          Index("-1",5),QN(-1));
    auto j = IQIndex("j",Index("+1",5),QN(+1),
                         Index("0",5),QN(0),
                         Index("-1",5),QN(-1),
                         In);
    
    auto A = randomTensor(QN(0),prime(dag(i)),prime(dag(j)),i,j);

    auto x = randomTensor(QN(0),dag(i),dag(j));
    auto b = randomTensor(QN(0),dag(i),dag(j));

    gmres(IQTensorMap(A),b,x,{"MaxIter",100,"DebugLevel",0,"ErrGoal",1e-10});

    CHECK_CLOSE(norm((A*x).mapprime(1,0)-b)/norm(b),0.0);

    }

}
