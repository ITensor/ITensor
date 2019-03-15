//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/mps/mps.h"
#include "itensor/mps/mpo.h"
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
        auto m = dim(l1)+dim(l2);
        if(m <= 0) m = 1;
        sumind = Index(m,"Link");

        first = delta(l1,sumind);
        auto S = Matrix(dim(l2),dim(sumind));
        for(auto i : range(dim(l2)))
            {
            S(i,dim(l1)+i) = 1;
            }
        second = matrixITensor(std::move(S),l2,sumind);
        }
    else
        {
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
        sumind = Index(std::move(siq),
                       sumind.dir(),
                       sumind.tags()).prime(sumind.primeLevel());
        first = ITensor(dag(l1),sumind);
        int n = 1;
        for(auto j : range1(l1.nblock()))
            {
            auto D = Tensor(l1.blocksize(j),sumind.blocksize(n));
            auto minsize = std::min(D.extent(0),D.extent(1));
            for(auto i : range(minsize)) D(i,i) = 1.0;
            getBlock<Real>(first,{j,n}) &= D;
            ++n;
            }
        second = ITensor(dag(l2),sumind);
        for(auto j : range1(l2.nblock()))
            {
            auto D = Tensor(l2.blocksize(j),sumind.blocksize(n));
            auto minsize = std::min(D.extent(0),D.extent(1));
            for(auto i : range(minsize)) D(i,i) = 1.0;
            getBlock<Real>(second,{j,n}) &= D;
            ++n;
            }
        }
    }

//
// Adds two MPSs but doesn't attempt to
// orthogonalize them first
//
template <class MPSType>
MPSType&
addAssumeOrth(MPSType      & L,
              MPSType const& R, 
              Args const& args)
    {
    auto N = length(L);
    if(length(R) != N) Error("Mismatched MPS sizes");

    // TODO: implement this for MPS, MPO
    L.replaceTags("0","4","Link");

    auto first = vector<ITensor>(N);
    auto second = vector<ITensor>(N);

    for(auto i : range1(N-1))
        {
        auto l1 = rightLinkIndex(L,i);
        auto l2 = rightLinkIndex(R,i);
        auto r = l1;
        plussers(l1,l2,r,first[i],second[i]);
        }

    L.ref(1) = L(1) * first.at(1) + R(1) * second.at(1);
    for(auto i : range1(2,N-1))
        {
        L.ref(i) = dag(first.at(i-1)) * L(i) * first.at(i) 
                     + dag(second.at(i-1)) * R(i) * second.at(i);
        }
    L.ref(N) = dag(first.at(N-1)) * L(N) + dag(second.at(N-1)) * R(N);

    L.noPrime("Link");

    L.orthogonalize(args);

    return L;
    }
template MPS& addAssumeOrth<MPS>(MPS & L,MPS const& R, Args const& args);
template MPO& addAssumeOrth<MPO>(MPO & L,MPO const& R, Args const& args);

void 
fitWF(MPS const& psi_basis, MPS & psi_to_fit)
    {
    if(!itensor::isOrtho(psi_basis)) 
        Error("psi_basis must be orthogonolized.");
    if(orthoCenter(psi_basis) != 1) 
        Error("psi_basis must be orthogonolized to site 1.");

    auto N = length(psi_basis);
    if(length(psi_to_fit) != N) 
        Error("Wavefunctions must have same number of sites.");

    auto A = psi_to_fit(N) * dag(prime(psi_basis(N),"Link"));
    for(int n = N-1; n > 1; --n)
        {
        A *= dag(prime(psi_basis(n),"Link"));
        A *= psi_to_fit(n);
        }
    A = psi_to_fit(1) * A;
    A.noPrime();

    auto nrm = norm(A);
    if(nrm == 0) Error("Zero overlap of psi_to_fit and psi_basis");
    A /= nrm;

    psi_to_fit = psi_basis;
    psi_to_fit.ref(1) = A;
    }

bool 
checkQNs(MPS const& psi)
    {
    const int N = length(psi);

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
        if(!psi(i))
            {
            println("A(",i,") null, QNs not well defined");
            return false;
            }
        if(div(psi(i)) != Zero)
            {
            cout << "At i = " << i << "\n";
            Print(psi(i));
            cout << "IQTensor other than the ortho center had non-zero divergence\n";
            return false;
            }
        }

    //Check arrows from left edge
    for(int i = 1; i < center; ++i)
        {
        if(rightLinkIndex(psi,i).dir() != In) 
            {
            println("checkQNs: At site ",i," to the left of the OC, Right side Link not pointing In");
            return false;
            }
        if(i > 1)
            {
            if(leftLinkIndex(psi,i).dir() != Out) 
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
        if(rightLinkIndex(psi,i).dir() != Out) 
            {
            println("checkQNs: At site ",i," to the right of the OC, Right side Link not pointing Out");
            return false;
            }
        if(leftLinkIndex(psi,i).dir() != In) 
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
    return div(psi(center));
    }

} //namespace itensor
