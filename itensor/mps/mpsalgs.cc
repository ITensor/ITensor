//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/mps/mps.h"
//#include "itensor/mps/mpo.h"
#include "itensor/mps/localop.h"
#include "itensor/util/print_macro.h"
#include "itensor/tensor/slicemat.h"

namespace itensor {

using std::istream;
using std::ostream;
using std::cout;
using std::endl;
using std::vector;
using std::find;
using std::pair;
using std::make_pair;
using std::string;
using std::move;

void 
plussers(Index const& l1, 
         Index const& l2, 
         Index      & sumind, 
         ITensor    & first, 
         ITensor    & second)
    {
    if(not hasQNs(l1) && not hasQNs(l2))
        {
        auto m = l1.m()+l2.m();
        if(m <= 0) m = 1;
        sumind = Index(sumind.rawname(),m);

        first = delta(l1,sumind);
        auto S = Matrix(l2.m(),sumind.m());
        for(auto i : range(l2.m()))
            {
            S(i,l1.m()+i) = 1;
            }
        second = matrixTensor(std::move(S),l2,sumind);
        }
    else
        {
        Error("plussers not yet implemented for QN case");
        
        auto siq = stdx::reserve_vector<QNInt>(l1.nblock()+l2.nblock());
        for(auto n : range1(l1.nblock()))
            {
            siq.emplace_back(l1.qn(n),l1.blocksize(n));
            }
        for(auto n : range1(l2.nblock()))
            {
            siq.emplace_back(l2.qn(n),l2.blocksize(n));
            }
#ifdef DEBUG
        if(siq.empty()) Error("siq is empty in plussers");
#endif
        sumind = Index(sumind.rawname(),
                       std::move(siq),
                       sumind.dir(),
                       sumind.type(),
                       sumind.primeLevel());
        //first = ITensor(dag(l1),sumind);
        int n = 1;
        for(auto j : range1(l1.nblock()))
            {
            auto D = Matrix(l1.blocksize(j),sumind.blocksize(n));
            auto minsize = std::min(ncols(D),nrows(D));
            for(auto i : range(minsize)) D(i,i) = 1.0;
            //first += matrixTensor(move(D),iq1.index,s1);
            //TODO: need the ability to add to a certain block of a QN ITensor
            //      form may be that of QDiag...
            ++n;
            }
        //second = ITensor(dag(l2),sumind);
        for(auto j : range1(l2.nblock()))
            {
            auto D = Matrix(l2.blocksize(j),sumind.blocksize(n));
            auto minsize = std::min(ncols(D),nrows(D));
            for(auto i : range(minsize)) D(i,i) = 1.0;
            //second += matrixTensor(move(D),iq2.index,s2);
            ++n;
            }
        }
    }

//
// Adds two MPSs but doesn't attempt to
// orthogonalize them first
//
MPS&
addAssumeOrth(MPS      & L,
              MPS const& R, 
              Args const& args)
    {
    auto N = L.N();
    if(R.N() != N) Error("Mismatched MPS sizes");

    L.primelinks(0,4);

    auto first = vector<ITensor>(N);
    auto second = vector<ITensor>(N);

    for(auto i : range1(N-1))
        {
        auto l1 = rightLinkInd(L,i);
        auto l2 = rightLinkInd(R,i);
        auto r = l1;
        plussers(l1,l2,r,first[i],second[i]);
        }

    L.Aref(1) = L.A(1) * first.at(1) + R.A(1) * second.at(1);
    for(auto i : range1(2,N-1))
        {
        L.Aref(i) = dag(first.at(i-1)) * L.A(i) * first.at(i) 
                     + dag(second.at(i-1)) * R.A(i) * second.at(i);
        }
    L.Aref(N) = dag(first.at(N-1)) * L.A(N) + dag(second.at(N-1)) * R.A(N);

    L.noprimelink();

    L.orthogonalize(args);

    return L;
    }

void 
fitWF(MPS const& psi_basis, MPS & psi_to_fit)
    {
    if(!itensor::isOrtho(psi_basis)) 
        Error("psi_basis must be orthogonolized.");
    if(orthoCenter(psi_basis) != 1) 
        Error("psi_basis must be orthogonolized to site 1.");

    auto N = psi_basis.N();
    if(psi_to_fit.N() != N) 
        Error("Wavefunctions must have same number of sites.");

    auto A = psi_to_fit.A(N) * dag(prime(psi_basis.A(N),Link));
    for(int n = N-1; n > 1; --n)
        {
        A *= dag(prime(psi_basis.A(n),Link));
        A *= psi_to_fit.A(n);
        }
    A = psi_to_fit.A(1) * A;
    A.noprime();

    auto nrm = norm(A);
    if(nrm == 0) Error("Zero overlap of psi_to_fit and psi_basis");
    A /= nrm;

    psi_to_fit = psi_basis;
    psi_to_fit.Anc(1) = A;
    }

bool 
checkQNs(MPS const& psi)
    {
    const int N = psi.N();

    QN Zero;

    int center = findCenter(psi);
    if(center == -1)
        {
        cout << "Did not find an ortho. center\n";
        return false;
        }

    //Check that all IQTensors have zero div
    //except possibly the ortho. center
    for(int i = 1; i <= N; ++i) 
        {
        if(i == center) continue;
        if(!psi.A(i))
            {
            println("A(",i,") null, QNs not well defined");
            return false;
            }
        if(div(psi.A(i)) != Zero)
            {
            cout << "At i = " << i << "\n";
            Print(psi.A(i));
            cout << "IQTensor other than the ortho center had non-zero divergence\n";
            return false;
            }
        }

    //Check arrows from left edge
    for(int i = 1; i < center; ++i)
        {
        if(rightLinkInd(psi,i).dir() != In) 
            {
            println("checkQNs: At site ",i," to the left of the OC, Right side Link not pointing In");
            return false;
            }
        if(i > 1)
            {
            if(leftLinkInd(psi,i).dir() != Out) 
                {
                println("checkQNs: At site ",i," to the left of the OC, Left side Link not pointing Out");
                return false;
                }
            }
        }

    //Check arrows from right edge
    for(int i = N; i > center; --i)
        {
        if(i < N)
        if(rightLinkInd(psi,i).dir() != Out) 
            {
            println("checkQNs: At site ",i," to the right of the OC, Right side Link not pointing Out");
            return false;
            }
        if(leftLinkInd(psi,i).dir() != In) 
            {
            println("checkQNs: At site ",i," to the right of the OC, Left side Link not pointing In");
            return false;
            }
        }

    //Done checking arrows
    return true;
    }

QN
totalQN(MPS const& psi)
    {
    const int center = findCenter(psi);
    if(center == -1)
        Error("Could not find ortho. center");
    return div(psi.A(center));
    }

} //namespace itensor
