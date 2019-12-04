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
#ifndef __ITENSOR_LOCALMPO_MPS
#define __ITENSOR_LOCALMPO_MPS
#include "itensor/mps/localmpo.h"

namespace itensor {

class LocalMPO_MPS
    {
    private:
    MPO const* Op_ = nullptr;
    //std::vector<MPS> const* psis_ = nullptr;
    //LocalMPO object representing projected version
    //of the MPO Op_
    LocalMPO lmpo_; 
    //LocalMPO objects representing projected version
    //of each MPS in psis_
    std::vector<LocalMPO> lmps_;
    Real weight_ = 1;
    public:

    LocalMPO_MPS(Args const& args = Args::global()) { }

    LocalMPO_MPS(MPO const& Op, 
                 std::vector<MPS > const& psis,
                 Args const& args = Args::global());

    LocalMPO_MPS(MPO const& Op, 
                 ITensor const& LOp,
                 ITensor const& ROp,
                 std::vector<MPS> const& psis,
                 std::vector<ITensor> const& Lpsi,
                 std::vector<ITensor> const& Rpsi,
                 Args const& args = Args::global());


    //
    // Sparse matrix methods
    //

    void
    product(ITensor const& phi, 
            ITensor& phip) const;

    Real
    expect(ITensor const& phi) const { return lmpo_.expect(phi); }

    ITensor
    deltaRho(ITensor const& AA, 
             ITensor const& comb, 
             Direction dir) const
        { return lmpo_.deltaRho(AA,comb,dir); }

    ITensor
    diag() const { return lmpo_.diag(); }

    void
    position(int b, MPS const& psi);

    size_t
    size() const { return lmpo_.size(); }

    explicit
    operator bool() const { return bool(Op_); }

    Real
    weight() const { return weight_; }
    void
    weight(Real val) { weight_ = val; }

    bool
    doWrite() const { return lmpo_.doWrite(); }
    void
    doWrite(bool val, Args const& args = Args::global()) { lmpo_.doWrite(val,args); }

    };

inline LocalMPO_MPS::
LocalMPO_MPS(MPO const& Op,
             std::vector<MPS> const& psis,
             Args const& args)
  : Op_(&Op),
    lmps_(psis.size()),
    weight_(args.getReal("Weight",1))
    { 
    lmpo_ = LocalMPO(Op);

    for(auto j : range(lmps_.size()))
        {
        lmps_[j] = LocalMPO(psis[j],args);
        }
    }

inline LocalMPO_MPS::
LocalMPO_MPS(MPO const& Op, 
             ITensor const& LOp,
             ITensor const& ROp,
             std::vector<MPS> const& psis,
             std::vector<ITensor> const& Lpsi,
             std::vector<ITensor> const& Rpsi,
             Args const& args)
  : Op_(&Op),
    lmps_(psis.size()),
    weight_(args.getReal("Weight",1))
    { 
    lmpo_ = LocalMPO(Op,LOp,ROp,args);
#ifdef DEBUG
    if(Lpsi.size() != psis.size()) Error("Lpsi must have same number of elements as psis");
    if(Rpsi.size() != psis.size()) Error("Rpsi must have same number of elements as psis");
#endif

    for(auto j : range(lmps_.size()))
        {
        lmps_[j] = LocalMPO(psis[j],Lpsi[j],Rpsi[j],args);
        }
    }

void inline LocalMPO_MPS::
product(ITensor const& phi, 
        ITensor & phip) const
    {
    lmpo_.product(phi,phip);

    ITensor outer;
    for(auto& M : lmps_)
        {
        M.product(phi,outer);
        outer *= weight_;
        phip += outer;
        }
    }

void inline LocalMPO_MPS::
position(int b, const MPS& psi)
    {
    lmpo_.position(b,psi);
    for(auto& M : lmps_)
        {
        M.position(b,psi);
        }
    }

} //namespace itensor

#endif
