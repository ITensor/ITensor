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
#ifndef __ITENSOR_LOCALMPOSET
#define __ITENSOR_LOCALMPOSET
#include "itensor/mps/localmpo.h"

namespace itensor {

class LocalMPOSet
    {
    std::vector<MPO> const* Op_ = nullptr;
    std::vector<LocalMPO> lmpo_;
    public:

    LocalMPOSet(Args const& args = Args::global()) { }

    LocalMPOSet(std::vector<MPO> const& Op,
                Args const& args = Args::global());

    LocalMPOSet(std::vector<MPO> const& H, 
                std::vector<ITensor> const& LH, 
                int LHlim,
                std::vector<ITensor> const& RH,
                int RHlim,
                Args const& args = Args::global());

    void
    product(ITensor const& phi, 
            ITensor & phip) const;

    Real
    expect(ITensor const& phi) const;

    ITensor
    deltaRho(ITensor const& AA, 
             ITensor const& comb, 
             Direction dir) const;

    ITensor
    diag() const;

    void
    position(int b, 
             MPS const& psi);

    std::vector<ITensor>
    L() const 
        { 
        auto L = std::vector<ITensor>(lmpo_.size());
        for(auto n : range(lmpo_)) L.at(n) = lmpo_.at(n).L();
        return L;
        }
    std::vector<ITensor>
    R() const 
        { 
        auto R = std::vector<ITensor>(lmpo_.size());
        for(auto n : range(lmpo_)) R.at(n) = lmpo_.at(n).R();
        return R;
        }

    void
    L(std::vector<ITensor> const& nL)
        { 
        for(auto n : range(lmpo_)) lmpo_.at(n).L(nL.at(n));
        }
    void
    R(std::vector<ITensor> const& nR)
        { 
        for(auto n : range(lmpo_)) lmpo_.at(n).R(nR.at(n));
        }

    void
    shift(int j, Direction dir, ITensor const& A)
        {
        for(auto n : range(lmpo_)) lmpo_[n].shift(j,dir,A);
        }

    int
    numCenter() const { return lmpo_.front().numCenter(); }
    void
    numCenter(int val);

    size_t
    size() const { return lmpo_.front().size(); }

    explicit
    operator bool() const { return bool(Op_); }

    bool
    doWrite() const { return lmpo_.front().doWrite(); }
    void
    doWrite(bool val, Args const& args = Args::global()) 
        { 
        for(auto& lm : lmpo_) lm.doWrite(val,args);
        }

    };

inline LocalMPOSet::
LocalMPOSet(std::vector<MPO> const& Op,
            Args const& args)
  : Op_(&Op),
    lmpo_(Op.size())
    { 
    for(auto n : range(lmpo_.size()))
        {
        lmpo_[n] = LocalMPO(Op.at(n),args);
        }
    }

inline LocalMPOSet::
LocalMPOSet(std::vector<MPO> const& H, 
            std::vector<ITensor> const& LH, 
            int LHlim,
            std::vector<ITensor> const& RH,
            int RHlim,
            Args const& args)
  : Op_(&H),
    lmpo_(H.size())
    { 
    for(auto n : range(lmpo_.size()))
        {
        lmpo_[n] = LocalMPO(H.at(n),LH.at(n),LHlim,RH.at(n),RHlim,args);
        }
    }

void inline LocalMPOSet::
product(ITensor const& phi, 
        ITensor & phip) const
    {
    lmpo_.front().product(phi,phip);

    ITensor phi_n;
    for(auto n : range(1,lmpo_.size()))
        {
        lmpo_[n].product(phi,phi_n);
        phip += phi_n;
        }
    }

Real inline LocalMPOSet::
expect(ITensor const& phi) const
    {
    Real ex_ = 0;
    for(size_t n = 0; n < lmpo_.size(); ++n)
    for(auto n : range(lmpo_.size()))
        {
        ex_ += lmpo_[n].expect(phi);
        }
    return ex_;
    }

ITensor inline LocalMPOSet::
deltaRho(ITensor const& AA,
         ITensor const& comb, 
         Direction dir) const
    {
    ITensor delta = lmpo_.front().deltaRho(AA,comb,dir);
    for(auto n : range(1,lmpo_.size()))
        {
        delta += lmpo_[n].deltaRho(AA,comb,dir);
        }
    return delta;
    }

ITensor inline LocalMPOSet::
diag() const
    {
    ITensor D = lmpo_.front().diag();
    for(auto n : range(1,lmpo_.size()))
        {
        D += lmpo_[n].diag();
        }
    return D;
    }

void inline LocalMPOSet::
position(int b, 
         MPS const& psi)
    {
    for(auto n : range(lmpo_.size()))
        {
        lmpo_[n].position(b,psi);
        }
    }

void inline LocalMPOSet::
numCenter(int val)
    {
    for(auto n : range(lmpo_.size()))
        {
        lmpo_[n].numCenter(val);
        }
    }

} //namespace itensor

#endif
