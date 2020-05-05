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
#ifndef __ITENSOR_LOCALMPO
#define __ITENSOR_LOCALMPO
#include "itensor/mps/mpo.h"
#include "itensor/mps/localop.h"
//#include "itensor/util/print_macro.h"

namespace itensor {

//
// The LocalMPO class projects an MPO 
// into the reduced Hilbert space of
// some number of sites of an MPS.
// (The default is 2 sites.)
//
//   .----...---                ----...--.
//   |  |     |      |      |     |      | 
//   W1-W2-..Wj-1 - Wj - Wj+1 -- Wj+2..-WN
//   |  |     |      |      |     |      | 
//   '----...---                ----...--'
//
// 
//  Here the W's are the site tensors
//  of the MPO "Op" and the method position(j,psi)
//  has been called using the MPS 'psi' as a basis 
//  for the projection.
//
//  This results in an unprojected region of
//  num_center sites starting at site j.
//

class LocalMPO
    {
    public:

    //
    // Constructors
    //

    LocalMPO();

    //
    //Regular case where H is an MPO for a finite system
    //
    LocalMPO(MPO const& H, 
             Args const& args = Args::global());

    //
    //Use an MPS instead of an MPO. Equivalent to using an MPO
    //of the outer product |Psi><Psi| but much more efficient.
    //
    LocalMPO(MPS const& Psi, 
             Args const& args = Args::global());
    
    //
    //Use left and right boundary tensors LH and RH. Ok if both boundary 
    //tensors are not default-constructed.
    //
    LocalMPO(const ITensor& LH, 
             const ITensor& RH,
             const Args& args = Args::global());

    //
    //Use an MPO having boundary indices capped off by left and
    //right boundary tensors LH and RH. Ok if one or both boundary 
    //tensors are default-constructed.
    //
    LocalMPO(const MPO& H, 
             const ITensor& LH, 
             const ITensor& RH,
             const Args& args = Args::global());

    //
    //Use an MPS with boundary indices capped off by left and right
    //boundary tensors LP and RP. Ok if one or both boundary tensors 
    //are default-constructed.
    //
    LocalMPO(const MPS& Psi, 
             const ITensor& LP,
             const ITensor& RP,
             const Args& args = Args::global());

    //
    //Use an MPO having boundary indices capped off by left and
    //right boundary tensors LH and RH at positions LHlim and RHlim
    //These positions indicate the site number of the right-most MPO
    //tensor included in LH and the left-most MPO tensor included in RH.
    //
    LocalMPO(MPO const& H, 
             ITensor const& LH, 
             int LHlim,
             ITensor const& RH,
             int RHlim,
             Args const& args = Args::global());

    //
    // Sparse Matrix Methods
    //

    void
    product(const ITensor& phi, ITensor& phip) const;

    Real
    expect(const ITensor& phi) const { return lop_.expect(phi); }

    ITensor
    deltaRho(const ITensor& AA, 
             const ITensor& comb, Direction dir) const
        { return lop_.deltaRho(AA,comb,dir); }

    ITensor
    diag() const { return lop_.diag(); }

    //
    // position(b,psi) uses the MPS psi
    // to adjust the edge tensors such
    // that the MPO tensors at positions
    // b and b+1 are exposed
    //
    void
    position(int b, MPS const& psi);

    int
    position() const;

    void
    shift(int j, Direction dir, ITensor const& A);

    //
    // Accessor Methods
    //

    void
    reset()
        {
        LHlim_ = 0;
        RHlim_ = Op_->length()+1;
        }

    ITensor const&
    L() const { return PH_[LHlim_]; }
    // Replace left edge tensor at current bond
    void
    L(ITensor const& nL) { PH_[LHlim_] = nL; }
    // Replace left edge tensor bordering site j
    // (so that nL includes sites < j)
    void
    L(int j, ITensor const& nL);


    ITensor const&
    R() const { return PH_[RHlim_]; }
    // Replace right edge tensor at current bond
    void
    R(ITensor const& nR) { PH_[RHlim_] = nR; }
    // Replace right edge tensor bordering site j
    // (so that nR includes sites > j)
    void
    R(int j, ITensor const& nR);

    const MPO&
    H() const 
        { 
        if(Op_ == 0)
            Error("LocalMPO is null or contains an MPS");
        return *Op_;
        }

    int
    numCenter() const { return nc_; }
    void
    numCenter(int val) 
        { 
        if(val < 0 || val > 2) Error("numCenter must be set 0 or 1 or 2");
        nc_ = val; 
		lop_.numCenter(val);
        }

    size_t
    size() const { return lop_.size(); }

    explicit operator bool() const { return Op_ != 0 || Psi_ != 0; }

    bool
    doWrite() const { return do_write_; }
    void
    doWrite(bool val,
            Args const& args = Args::global()) 
        { 
        if(Psi_ != 0) Error("Write to disk not yet supported for LocalMPO initialized with an MPS");
        if(!do_write_ && (val == true))
            {
            initWrite(args); 
            }
        do_write_ = val; 
        }

    std::string const&
    writeDir() const { return writedir_; }

    int
    leftLim() const { return LHlim_; }

    int
    rightLim() const { return RHlim_; }

    private:

    /////////////////
    //
    // Data Members
    //

    const MPO* Op_;
    std::vector<ITensor> PH_;
    int LHlim_,RHlim_;
    int nc_;

    LocalOp lop_;

    bool do_write_ = false;
    std::string writedir_ = "./";

    const MPS* Psi_;

    //
    /////////////////

    void
    makeL(const MPS& psi, int k);

    void
    makeR(const MPS& psi, int k);

    void
    setLHlim(int val);

    void
    setRHlim(int val);

    void
    initWrite(Args const& args);

    std::string
    PHFName(int j) const
        {
        return format("%s/PH_%03d",writedir_,j);
        }

    };

inline LocalMPO::
LocalMPO()
    : Op_(0),
      LHlim_(-1),
      RHlim_(-1),
      nc_(2),
      Psi_(0)
    { }

inline LocalMPO::
LocalMPO(const MPO& H, 
         const Args& args)
    : Op_(&H),
      PH_(H.length()+2),
      LHlim_(0),
      RHlim_(H.length()+1),
      nc_(2),
      Psi_(0)
    { 
    if(args.defined("NumCenter"))
        numCenter(args.getInt("NumCenter"));
    }

inline LocalMPO::
LocalMPO(const MPS& Psi, 
         const Args& args)
    : Op_(0),
      PH_(Psi.length()+2),
      LHlim_(0),
      RHlim_(Psi.length()+1),
      nc_(2),
      Psi_(&Psi)
    { 
    if(args.defined("NumCenter"))
        numCenter(args.getInt("NumCenter"));
    }

inline LocalMPO::
LocalMPO(const ITensor& LH, const ITensor& RH,
         const Args& args)
    : Op_(0),
      PH_(2),
      LHlim_(0),
      RHlim_(1),
      nc_(0),
      Psi_(0)
    { 
    PH_[0] = LH;
    PH_[1] = RH;
    lop_.update(L(), R());//nc_ must be set to 0 in this case
    }

inline LocalMPO::
LocalMPO(const MPO& H, 
         const ITensor& LH, const ITensor& RH,
         const Args& args)
    : Op_(&H),
      PH_(H.length()+2),
      LHlim_(0),
      RHlim_(H.length()+1),
      nc_(2),
      Psi_(0)
    { 
    PH_[0] = LH;
    PH_[H.length()+1] = RH;
    if(H.length() == 1)
        lop_.update(Op_->A(1), L(), R());
	else if(H.length() == 2)
        lop_.update(Op_->A(1), Op_->A(2), L(), R());
    if(args.defined("NumCenter"))
        numCenter(args.getInt("NumCenter"));
    }

inline LocalMPO::
LocalMPO(MPS const& Psi, 
         ITensor const& LP,
         ITensor const& RP,
         Args const& args)
    : Op_(0),
      PH_(Psi.length()+2),
      LHlim_(0),
      RHlim_(Psi.length()+1),
      nc_(2),
      Psi_(&Psi)
    { 
    PH_[0] = LP;
    PH_[Psi.length()+1] = RP;
    if(args.defined("NumCenter"))
        numCenter(args.getInt("NumCenter"));
    }

inline LocalMPO::
LocalMPO(MPO const& H, 
         ITensor const& LH, 
         int LHlim,
         ITensor const& RH,
         int RHlim,
         Args const& args)
    : Op_(&H),
      PH_(H.length()+2),
      LHlim_(LHlim),
      RHlim_(RHlim),
      nc_(2),
      Psi_(0)
    { 
    PH_.at(LHlim) = LH;
    PH_.at(RHlim) = RH;
    if(H.length() == 1)
        lop_.update(Op_->A(1), L(), R());
    if(H.length() == 2) 
        lop_.update(Op_->A(1), Op_->A(2), L(), R());
    if(args.defined("NumCenter")) numCenter(args.getInt("NumCenter"));
    }

void inline LocalMPO::
product(ITensor const& phi, 
        ITensor& phip) const
    {
    if(Op_ != 0)
        {
        lop_.product(phi,phip);
        }
    else 
    if(Psi_ != 0)
        {
        int b = position();
 
        ITensor othr;
        if(nc_ == 2)
            {
            othr = (!L() ? dag(prime(Psi_->A(b),"Link")) : L()*dag(prime(Psi_->A(b),"Link")));
            othr *= (!R() ? dag(prime(Psi_->A(b+1),"Link")) : R()*dag(prime(Psi_->A(b+1),"Link")));
            }
        else if(nc_ == 1)
            {
            othr = (!L() ? dag(prime(Psi_->A(b),"Link")) : L()*dag(prime(Psi_->A(b),"Link")));
            if(R()) othr *= R();	
            }
        else if(nc_ == 0)
            {
            if(!L())
                {
                if(!R()) Error("LocalMPO: Empty L() and R() in function product");
                else othr = R();
                }
            else
                {
                othr = L();
                if(R()) othr *= R();
                }
            }
        
        auto z = (othr*phi).eltC();

        phip = dag(othr);
        phip *= z;
        }
    else
        {
        Error("LocalMPO is null");
        }
    }

void inline LocalMPO::
L(int j, ITensor const& nL)
    {
    if(LHlim_ > j-1) setLHlim(j-1);
    PH_[LHlim_] = nL;
    }

void inline LocalMPO::
R(int j, ITensor const& nR)
    {
    if(RHlim_ < j+1) setRHlim(j+1);
    PH_[RHlim_] = nR;
    }

inline void LocalMPO::
position(int b, MPS const& psi)
    {
    if(!(*this)) Error("LocalMPO is null");

    makeL(psi,b-1);
    makeR(psi,b+nc_);

    setLHlim(b-1); //not redundant since LHlim_ could be > b-1
    setRHlim(b+nc_); //not redundant since RHlim_ could be < b+nc_

#ifdef DEBUG
    if(nc_ != 2 && nc_ != 1 && nc_ != 0)
        {
        Error("LocalOp only supports 0 and 1 and 2 center sites currently");
        }
#endif

    if(Op_ != 0) //normal MPO case
        {
        if(nc_ == 2)
            lop_.update(Op_->A(b), Op_->A(b+1), L(), R());
        else if(nc_ == 1)
            lop_.update(Op_->A(b), L(), R());
        else if(nc_ == 0)
            lop_.update(L(),R());
        }
    }

int inline LocalMPO::
position() const
    {
    if(RHlim_-LHlim_ != (nc_+1))
        {
        throw ITError("LocalMPO position not set");
        }
    return LHlim_+1;
    }

inline void LocalMPO::
shift(int j, 
      Direction dir, 
      ITensor const& A)
    {
    if(!(*this)) Error("LocalMPO is null");

#ifdef DEBUG
    if(nc_ != 2 && nc_ != 1 && nc_ != 0)
        {
        Error("LocalOp only supports 0 and 1 and 2 center sites currently");
        }
#endif

    if(dir == Fromleft)
        {
        if((j-1) != LHlim_)
            {
            std::cout << "j-1 = " << (j-1) << ", LHlim = " << LHlim_ << std::endl;
            Error("Can only shift at LHlim");
            }
        auto& E = PH_.at(LHlim_);
        auto& nE = PH_.at(j);
        nE = E * A;
        nE *= Op_->A(j);
        nE *= dag(prime(A));
        setLHlim(j);
        setRHlim(j+nc_+1);

        if(nc_ == 2)
            lop_.update(Op_->A(j+1), Op_->A(j+2), L(), R());
        else if(nc_ == 1)
            lop_.update(Op_->A(j+1), L(), R());
        else if(nc_ == 0)
            lop_.update(L(), R());
        }
    else //dir == Fromright
        {
        if((j+1) != LHlim_)
            {
            std::cout << "j+1 = " << (j+1) << ", RHlim_ = " << RHlim_ << std::endl;
            Error("Can only shift at RHlim_");
            }
        auto& E = PH_.at(RHlim_);
        auto& nE = PH_.at(j);
        nE = E * A;
        nE *= Op_->A(j);
        nE *= dag(prime(A));
        setLHlim(j-nc_-1);
        setRHlim(j);
	
        if(nc_ == 2)
            lop_.update(Op_->A(j-2), Op_->A(j-1), L(), R());
        else if(nc_ == 1)
            lop_.update(Op_->A(j-1), L(), R());
        else if(nc_ == 0)
            lop_.update(L(), R());
        }
    }

inline void LocalMPO::
makeL(MPS const& psi, int k)
    {
    if(!PH_.empty())
        {
        if(Op_ == 0) //Op is actually an MPS
            {
            while(LHlim_ < k)
                {
                auto ll = LHlim_;
                PH_.at(ll+1) = (!PH_.at(ll) ? psi(ll+1) : PH_[ll]*psi(ll+1));
                PH_[ll+1] *= dag(prime(Psi_->A(ll+1),"Link"));
                setLHlim(ll+1);
                }
            }
        else //normal MPO case
            {
            while(LHlim_ < k)
                {
                auto ll = LHlim_;
                if(PH_.at(ll))
                    {
                    PH_.at(ll+1) = PH_.at(ll)*psi(ll+1);
                    }
                else
                    {
                    PH_.at(ll+1) = psi(ll+1);
                    }
                PH_.at(ll+1) *= Op_->A(ll+1);
                PH_.at(ll+1) *= dag(prime(psi(ll+1)));
                setLHlim(ll+1);
                }
            }
        }
    }

inline void LocalMPO::
makeR(MPS const& psi, int k)
    {
    if(!PH_.empty())
        {
        if(Op_ == 0) //Op is actually an MPS
            {
            while(RHlim_ > k)
                {
                const int rl = RHlim_;
                PH_.at(rl-1) = (!PH_.at(rl) ? psi(rl-1) : PH_[rl]*psi(rl-1));
                PH_[rl-1] *= dag(prime(Psi_->A(rl-1),"Link"));
                setRHlim(rl-1);
                }
            }
        else //normal MPO case
            {
            while(RHlim_ > k)
                {
                auto rl = RHlim_;
                //printfln(" Making environment with rl=%d (using H[%d])",rl,rl-1);
                //Print(PH_.at(rl));
                //Print(Op_->A(rl-1));
                //Print(psi(rl-1));
                if(PH_.at(rl))
                    {
                    PH_.at(rl-1) = PH_.at(rl)*psi(rl-1);
                    }
                else
                    {
                    PH_.at(rl-1) = psi(rl-1);
                    }
                PH_.at(rl-1) *= Op_->A(rl-1);
                PH_.at(rl-1) *= dag(prime(psi(rl-1)));
                //printfln("PH[%d] = \n%s",rl-1,PH_.at(rl-1));
                //PAUSE
                setRHlim(rl-1);
                }
            }
        }
    }

void inline LocalMPO::
setLHlim(int val)
    {
    if(!do_write_)
        {
        LHlim_ = val;
        return;
        }

    if(LHlim_ != val && PH_.at(LHlim_))
        {
        writeToFile(PHFName(LHlim_),PH_.at(LHlim_));
        PH_.at(LHlim_) = ITensor();
        }
    LHlim_ = val;
    if(LHlim_ < 1) 
        {
        //Set to null tensor and return
        PH_.at(LHlim_) = ITensor();
        return;
        }
    if(!PH_.at(LHlim_))
        {
        std::string fname = PHFName(LHlim_);
        readFromFile(fname,PH_.at(LHlim_));
        }
    }

void inline LocalMPO::
setRHlim(int val)
    {
    if(!do_write_)
        {
        RHlim_ = val;
        return;
        }

    if(RHlim_ != val && PH_.at(RHlim_))
        {
        writeToFile(PHFName(RHlim_),PH_.at(RHlim_));
        PH_.at(RHlim_) = ITensor();
        }
    RHlim_ = val;
    if(RHlim_ > Op_->length()) 
        {
        //Set to null tensor and return
        PH_.at(RHlim_) = ITensor();
        return;
        }
    if(!PH_.at(RHlim_))
        {
        std::string fname = PHFName(RHlim_);
        readFromFile(fname,PH_.at(RHlim_));
        }
    }

void inline LocalMPO::
initWrite(Args const& args)
    {
    auto basedir = args.getString("WriteDir","./");
    writedir_ = mkTempDir("PH",basedir);
    }

} //namespace itensor


#endif
