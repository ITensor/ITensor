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

void 
plussers(IQIndex const& l1, 
         IQIndex const& l2, 
         IQIndex& sumind, 
         IQTensor& first, IQTensor& second)
    {
    auto siq = stdx::reserve_vector<IndexQN>(l1.nindex()+l2.nindex());
    for(auto iq1 : l1)
        {
        auto s1 = Index(iq1.index.rawname(),iq1.m(),iq1.type());
        siq.emplace_back(s1,iq1.qn);
        }
    for(auto iq2 : l2)
        {
        auto s2 = Index(iq2.index.rawname(),iq2.m(),iq2.type());
        siq.emplace_back(s2,iq2.qn);
        }
#ifdef DEBUG
    if(siq.empty()) Error("siq is empty in plussers");
#endif
    sumind = IQIndex(sumind.rawname(),std::move(siq),sumind.dir(),sumind.primeLevel());
    first = IQTensor(dag(l1),sumind);
    int n = 1;
    for(auto iq1 : l1)
        {
        auto s1 = sumind.index(n);
        auto D = Matrix(iq1.index.m(),s1.m());
        auto minsize = std::min(iq1.index.m(),s1.m());
        for(auto i : range(minsize)) D(i,i) = 1.0;
        first += matrixTensor(move(D),iq1.index,s1);
        ++n;
        }
    second = IQTensor(dag(l2),sumind);
    for(auto iq2 : l2)
        {
        auto s2 = sumind.index(n);
        auto D = Matrix(iq2.index.m(),s2.m());
        auto minsize = std::min(iq2.index.m(),s2.m());
        for(auto i : range(minsize)) D(i,i) = 1.0;
        second += matrixTensor(move(D),iq2.index,s2);
        ++n;
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
    using Tensor = typename MPSType::TensorT;

    auto N = L.N();
    if(R.N() != N) Error("Mismatched MPS sizes");

    L.primelinks(0,4);

    auto first = vector<Tensor>(N);
    auto second = vector<Tensor>(N);

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
template MPS& addAssumeOrth(MPS & L,MPS const& R, Args const& args);
template IQMPS& addAssumeOrth(IQMPS & L,IQMPS const& R, Args const& args);
template MPO& addAssumeOrth(MPO & L,MPO const& R, Args const& args);
template IQMPO& addAssumeOrth(IQMPO & L,IQMPO const& R, Args const& args);

template <class Tensor>
void 
fitWF(const MPSt<Tensor>& psi_basis, MPSt<Tensor>& psi_to_fit)
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
template void fitWF(const MPSt<ITensor>& psi_basis, MPSt<ITensor>& psi_to_fit);
template void fitWF(const MPSt<IQTensor>& psi_basis, MPSt<IQTensor>& psi_to_fit);

bool 
checkQNs(const IQMPS& psi)
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
totalQN(const IQMPS& psi)
    {
    const int center = findCenter(psi);
    if(center == -1)
        Error("Could not find ortho. center");
    return div(psi.A(center));
    }

} //namespace itensor
