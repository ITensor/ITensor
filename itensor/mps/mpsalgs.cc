//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
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
        sumind = Index(m,tags(sumind));

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
        auto siq = stdx::reserve_vector<QNInt>(nblock(l1)+nblock(l2));
        for(auto n : range1(nblock(l1)))
            {
            siq.emplace_back(qn(l1,n),blocksize(l1,n));
            }
        for(auto n : range1(nblock(l2)))
            {
            siq.emplace_back(qn(l2,n),blocksize(l2,n));
            }
#ifdef DEBUG
        if(siq.empty()) Error("siq is empty in plussers");
#endif
        sumind = Index(std::move(siq),
                       dir(sumind),
                       tags(sumind));
        first = ITensor(dag(l1),sumind);
        int n = 1;
        for(auto j : range1(nblock(l1)))
            {
            auto D = Tensor(blocksize(l1,j),blocksize(sumind,n));
            auto minsize = std::min(D.extent(0),D.extent(1));
            for(auto i : range(minsize)) D(i,i) = 1.0;
            getBlock<Real>(first,{j,n}) &= D;
            ++n;
            }
        second = ITensor(dag(l2),sumind);
        for(auto j : range1(nblock(l2)))
            {
            auto D = Tensor(blocksize(l2,j),blocksize(sumind,n));
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

    // Make sure there aren't link index clashes between L and R
    // by priming by a random amount
    // TODO: use L.simLinkInds() to avoid clashing insteda of priming
    auto rand_plev = 1254313;
    auto l = linkInds(L);
    L.replaceLinkInds(prime(linkInds(L),rand_plev));

    auto first = vector<ITensor>(N);
    auto second = vector<ITensor>(N);

    for(auto i : range1(N-1))
        {
        auto l1 = linkIndex(L,i);
        auto l2 = linkIndex(R,i);
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

    L.replaceLinkInds(prime(linkInds(L),-rand_plev));
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
    if(nrm == 0) Error("Zero inner of psi_to_fit and psi_basis");
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
            cout << "ITensor other than the ortho center had non-zero divergence\n";
            return false;
            }
        }

    //Check arrows from left edge
    for(int i = 1; i < center; ++i)
        {
        if(dir(rightLinkIndex(psi,i)) != In) 
            {
            println("checkQNs: At site ",i," to the left of the OC, Right side Link not pointing In");
            return false;
            }
        if(i > 1)
            {
            if(dir(leftLinkIndex(psi,i)) != Out) 
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
        if(dir(rightLinkIndex(psi,i)) != Out) 
            {
            println("checkQNs: At site ",i," to the right of the OC, Right side Link not pointing Out");
            return false;
            }
        if(dir(leftLinkIndex(psi,i)) != In) 
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
    auto tq = QN();
    auto sj = psi.leftLim()+1;
    auto ej = psi.rightLim()-1;
    for(int j = sj; j <= ej; ++j)
        {
        tq += flux(psi(j));
        }
    return tq;
    }

} //namespace itensor
