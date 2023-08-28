#include "test.h"
#include "itensor/iterativesolvers.h"
#include "sample/Heisenberg.h"
#include "itensor/mps/sites/spinhalf.h"
#include "itensor/mps/localmpo.h"
#include "itensor/mps/autompo.h"

using namespace itensor;
using namespace std;

class ITensorMap
    {
    ITensor const& A_;
    mutable long size_;
    public:

    ITensorMap(ITensor const& A)
      : A_(A)
        {
        size_ = 1;
        for(auto& I : A_.inds())
            {
            if(I.primeLevel() > 0)
                size_ *= dim(I);
            }
        }

    void
    product(ITensor const& x, ITensor& b) const
        {
        b = A_*x;
        b.noPrime();
        }

    long
    size() const { return size_; }

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

    psi.position(2);
    PH.position(2,psi);

    ITensor phi1 = psi(2) * psi(3);

    Real En1 = davidson(PH,phi1,"MaxIter=9");
    CHECK_CLOSE(En1,-0.95710678118);

    cout << endl << endl;

    //TODO: should this be put back in?
    //ITensor mpoh = H(2)*H(3);
    //ITensor phi2 = psi(2)*psi(3);
    //Real En2 = doDavidson(phi2,mpoh,PH.L(),PH.R(),9,2,1E-4);
    //cout << format("Energy from matrix Davidson (b=2) = %.20f")%En2 << endl;

    //cout << endl << endl;

    //cout << endl << "Overlap = " << Dot(phi1,phi2) << endl;

    //cout << "---------------------------------------" << endl << endl;

    //psi.doSVD(2,phi1,Fromleft);
    //psi.position(3);
    //PH.position(3,psi);
    //ITensor phi3 = psi(3) * psi(4);
    //Real En3 = d.davidson(PH,phi3);
    //cout << format("Energy from tensor Davidson (b=3) = %.20f")%En3 << endl;

    //cout << endl << endl;

    //mpoh = H(3)*H(4);
    //ITensor phi3m = psi(3)*psi(4);
    //Real En3m = doDavidson(phi3m,mpoh,PH.L(),PH.R(),9,2,1E-4);
    //cout << format("Energy from matrix Davidson (b=3) = %.20f")%En3m << endl;

    //cout << "---------------------------------------" << endl << endl;

    //psi.doSVD(3,phi3,Fromright);
    //psi.position(3);
    //PH.position(2,psi);
    //ITensor phi4 = psi(2) * psi(3);
    //Real En4 = d.davidson(PH,phi4);
    //cout << format("Energy from tensor Davidson (b=2) = %.20f")%En4 << endl;

    //cout << endl << endl;

    //mpoh = H(2)*H(3);
    //ITensor phi4m = psi(2)*psi(3);
    //Real En4m = doDavidson(phi4m,mpoh,PH.L(),PH.R(),9,2,1E-4);
    //cout << format("Energy from matrix Davidson (b=2) = %.20f")%En4m << endl;

    //cout << "---------------------------------------" << endl << endl;

    //psi.doSVD(2,phi4,Fromright);
    //psi.position(1);
    //PH.position(1,psi);

    ////With doDavidson
    //mpoh = H(1)*H(2);
    //ITensor phi5 = psi(1) * psi(2);
    //Real En5 = doDavidson(phi5,mpoh,PH.L(),PH.R(),9,2,1E-4);
    //cout << format("Energy from matrix Davidson (b=1) = %.20f")%En5 << endl;

    //cout << endl << endl;

    //ITensor phi6 = psi(1) * psi(2);
    ////Print(phi6);
    ////Print(PH.L());
    ////Print(PH.R());
    ////ITensor AB = phi6 * PH.R();
    ////AB *= H(2);
    ////AB *= H(1);
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

    psi.position(2);
    PH.position(2,psi);

    ITensor phi1 = psi(2) * psi(3);

    Real En1 = davidson(PH,phi1,"MaxIter=9");
    //cout << format("Energy from tensor Davidson (b=2) = %.20f")%En1 << endl;
    CHECK_CLOSE(En1,-0.95710678118);


    }

SECTION("Davidson (Custom Linear Map)")
    {
    auto a1 = Index(3,"Site,a1");
    auto a2 = Index(4,"Site,a2");
    auto a3 = Index(5,"Site,a3");

    auto A = randomITensor(prime(a1),prime(a2),prime(a3),a1,a2,a3);
    A = 0.5*(A + swapPrime(dag(A),0,1));
    auto x = randomITensor(a1,a2,a3);

    // ITensorMap is defined above, it simply wraps an ITensor that is of the
    // form of a matrix (i.e. has indices of the form {i,j,k,...,i',j',k',...})
    auto lambda = davidson(ITensorMap(A),x,{"MaxIter",40,"ErrGoal",1e-14});

    CHECK_CLOSE(norm(noPrime(A*x)-lambda*x)/norm(x),0.0);

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

    CHECK_CLOSE(norm(noPrime(A*x)-b)/norm(b),0.0);

    }

SECTION("GMRES (ITensor, Real), size 1")
    {
    auto a1 = Index(1,"Site,a1");
    auto a2 = Index(1,"Site,a2");
    auto a3 = Index(1,"Site,a3");

    auto A = randomITensor(prime(a1),prime(a2),prime(a3),a1,a2,a3);
    auto x = randomITensor(a1,a2,a3);
    auto b = randomITensor(a1,a2,a3);

    // ITensorMap is defined above, it simply wraps an ITensor that is of the
    // form of a matrix (i.e. has indices of the form {i,j,k,...,i',j',k',...})
    gmres(ITensorMap(A),b,x,{"MaxIter",36,"ErrGoal",1e-10});

    CHECK_CLOSE(norm(noPrime(A*x)-b)/norm(b),0.0);

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

    CHECK_CLOSE(norm(noPrime(A*x)-b)/norm(b),0.0);

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

    CHECK_CLOSE(norm(noPrime(A*x)-b)/norm(b),0.0);

    }

SECTION("Arnoldi (No QN)")
    {
    const int N = 50;
    const int Nc = N/2;
    SpinHalf sites(N,{"ConserveQNs=",false});

    // Create random MPO
    auto ampo = AutoMPO(sites);
    for(int j = 1; j < N; ++j)
        {
        ampo += 0.5,"S+",j,"S-",j+1;
        ampo += 0.5,"S-",j,"S+",j+1;
        ampo +=     "Sz",j,"Sz",j+1;
        }
    auto H = toMPO(ampo);
    for(int j = 1; j <= N; ++j)
      {
      H.ref(j).randomize();
      H.ref(j) /= norm(H(j));
      }
    auto normH = sqrt(trace(H,prime(H)));
    for(int j = 1; j <= N; ++j)
      H.ref(j) /= pow(normH,1.0/N);

    // Create random starting MPS
    auto initState = InitState(sites);
    for(int i = 1; i <= N; ++i)
        initState.set(i,i%2==1 ? "Up" : "Dn");
    auto psi = MPS(initState);
    for(int j = 1; j < 2; ++j)
      {
      psi = applyMPO(H,psi);
      psi.noPrime();
      psi.normalize();
      }
    for(int j = 1; j <= N; ++j)
      {
      psi.ref(j).randomize();
      psi.ref(j) /= norm(psi(j));
      }

    psi.position(Nc);
    psi.normalize();

    LocalMPO PH(H);

    PH.position(Nc,psi);

    auto x = psi(Nc) * psi(Nc+1);
    x.randomize();

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
    auto H = toMPO(ampo);
    for(int j = 1; j <= N; ++j)
      {
      H.ref(j).randomize();
      H.ref(j) /= norm(H(j));
      }
    auto normH = sqrt(trace(H,prime(H)));
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
      psi.noPrime();
      psi.normalize();
      }
    for(int j = 1; j <= N; ++j)
      {
      psi.ref(j).randomize();
      psi.ref(j) /= norm(psi(j));
      }
    psi.position(Nc);
    psi.normalize();

    LocalMPO PH(H);

    psi.position(Nc);
    PH.position(Nc,psi);

    auto x = psi.A(Nc) * psi.A(Nc+1);
    x.randomize();

    auto lambda = arnoldi(PH,x,{"MaxIter",20,"ErrGoal",1e-14,"DebugLevel",0,"WhichEig","LargestMagnitude"});
    auto PHx = x;
    PH.product(x,PHx);

    CHECK_CLOSE(norm(PHx-lambda*x)/norm(PHx),0.0);

    }

SECTION("arnoldi (multiple eigenvectors)")
    {
    auto i = Index(10,"i");
    auto A = randomITensor(prime(i), i);

    auto v1r = randomITensor(i);
    auto lambda1r = arnoldi(ITensorMap(A),v1r,{"ErrGoal=",1E-14,"MaxIter=",20,"MaxRestart=",5});

    auto v1l = randomITensor(i);
    auto lambda1l = arnoldi(ITensorMap(swapPrime(A,0,1)),v1l,{"ErrGoal=",1E-14,"MaxIter=",20,"MaxRestart=",5});

    auto s1 = eltC(v1r*v1l);
    auto v2r = randomITensor(i);
    auto lambda2r = arnoldi(ITensorMap(A-(lambda1r/s1)*prime(v1r)*v1l),v2r,{"ErrGoal=",1E-14,"MaxIter=",20,"MaxRestart=",5});
    auto v2l = randomITensor(i);
    auto lambda2l = arnoldi(ITensorMap(swapPrime(A,0,1)-(lambda1l/s1)*prime(v1l)*v1r),v2l,{"ErrGoal=",1E-14,"MaxIter=",20,"MaxRestart=",5});

    CHECK_CLOSE(real(lambda1l), real(lambda1r));
    CHECK_CLOSE(imag(lambda1l), imag(lambda1r));
    CHECK_CLOSE(real(lambda2l), real(lambda2r));
    CHECK_CLOSE(imag(lambda2l), imag(lambda2r));
    CHECK_CLOSE(norm(noPrime(A*v1r) - lambda1r*v1r), 0.0);
    CHECK_CLOSE(norm(noPrime(A*v2r) - lambda2r*v2r), 0.0);

/*
    auto x = randomITensor(i);
    auto y = randomITensor(i);
    auto vec = std::vector<ITensor>({x,y});

    auto lambda = arnoldi(ITensorMap(A),vec,{"ErrGoal=",1E-14,"MaxIter=",20,"MaxRestart=",5});

    CHECK_CLOSE(real(lambda[0]), real(lambda1r));
    CHECK_CLOSE(imag(lambda[0]), imag(lambda1r));

    // XXX: BROKEN
    //CHECK_CLOSE(real(lambda[1]), real(lambda2r));
    //CHECK_CLOSE(imag(lambda[1]), imag(lambda2r));

    CHECK_CLOSE(norm(noPrime(A*vec[0]) - lambda[0]*vec[0]), 0.0);

    // XXX: BROKEN
    //CHECK_CLOSE(norm(noPrime(A*vec[1]) - lambda[1]*vec[1]), 0.0);
*/

    // Compare to ED
    auto [U,D] = eigen(A);

    auto l = commonIndex(U, D);

    int d1 = (abs(arg(lambda1r) - arg(D.eltC(1,1))) < 1E-12)? 1 : 2;
    CHECK_CLOSE(real(D.eltC(d1,d1)), real(lambda1r));
    CHECK_CLOSE(imag(D.eltC(d1,d1)), imag(lambda1r));
    // These should compare up to a phase
    PrintData(v1r);
    PrintData(U*setElt(l=d1));

    int d2 = (abs(arg(lambda2r) - arg(D.eltC(2,2))) < 1E-12)? 2 : 3;
    CHECK_CLOSE(real(D.eltC(d2,d2)), real(lambda2r));
    CHECK_CLOSE(imag(D.eltC(d2,d2)), imag(lambda2r));
    // These should compare up to a phase
    PrintData(v2r);
    PrintData(U*setElt(l=d2));
    }

SECTION("applyExp (QNs)")
    {
    auto i = Index(QN(-1),10,QN(1),10,"i");

    auto A = randomITensor(QN(0),dag(i),prime(i));
    A += swapPrime(dag(A),0,1);
    A *= 0.5;
    auto x0 = randomITensor(QN(-1),i);

    auto Ac = randomITensorC(QN(0),dag(i),prime(i));
    Ac += swapPrime(dag(Ac),0,1);
    Ac *= 0.5;
    auto x0c = randomITensorC(QN(-1),i);

    SECTION("Real timestep")
        {
        auto t = 0.1;

        auto x = x0;
        applyExp(ITensorMap(A),x,-t,{"ErrGoal=",1E-14,
                                     "MaxIter=",10});

        auto exptA = expHermitian(A,-t);
        auto exptAx = noPrime(exptA*x0);

        CHECK_CLOSE(norm(exptAx - x), 0.);
        }
    
    SECTION("Complex timestep")
        {
        auto t = 0.1*1_i;
        
        auto x = x0;
        applyExp(ITensorMap(A),x,-t,{"ErrGoal=",1E-14,
                                     "MaxIter=",10});
        
        auto exptA = expHermitian(A,-t);
        auto exptAx = noPrime(exptA*x0);
        
        CHECK_CLOSE(norm(exptAx - x), 0.);
        }

    SECTION("Complex tensors, Real timestep")
        {
        auto t = 0.1;

        auto x = x0;
        applyExp(ITensorMap(Ac),x,-t,{"ErrGoal=",1E-14,
                                      "MaxIter=",10});

        auto exptA = expHermitian(Ac,-t);
        auto exptAx = noPrime(exptA*x0);

        CHECK_CLOSE(norm(exptAx - x), 0.);
        }

    SECTION("Complex tensors Complex timestep")
        {
        auto t = 0.1*1_i;

        auto x = x0c;
        applyExp(ITensorMap(Ac),x,-t,{"ErrGoal=",1E-14,
                                     "MaxIter=",10});

        auto exptA = expHermitian(Ac,-t);
        auto exptAx = noPrime(exptA*x0c);

        CHECK_CLOSE(norm(exptAx - x), 0.);
        }

    }

}
