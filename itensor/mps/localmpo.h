//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_LOCALMPO
#define __ITENSOR_LOCALMPO
#include "itensor/mps/mpo.h"
#include "itensor/mps/localop.h"
#include "itensor/util/print_macro.h"

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

template <class Tensor>
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
    LocalMPO(const MPOt<Tensor>& H, 
             const Args& args = Global::args());

    //
    //Use an MPS instead of an MPO. Equivalent to using an MPO
    //of the outer product |Psi><Psi| but much more efficient.
    //
    LocalMPO(const MPSt<Tensor>& Psi, 
             const Args& args = Global::args());

    //
    //Use an MPO having boundary indices capped off by left and
    //right boundary tensors LH and RH. Ok if one or both boundary 
    //tensors are default-constructed.
    //
    LocalMPO(const MPOt<Tensor>& H, 
             const Tensor& LH, 
             const Tensor& RH,
             const Args& args = Global::args());

    //
    //Use an MPS with boundary indices capped off by left and right
    //boundary tensors LP and RP. Ok if one or both boundary tensors 
    //are default-constructed.
    //
    LocalMPO(const MPSt<Tensor>& Psi, 
             const Tensor& LP,
             const Tensor& RP,
             const Args& args = Global::args());

    //
    //Use an MPO having boundary indices capped off by left and
    //right boundary tensors LH and RH at positions LHlim and RHlim
    //These positions indicate the site number of the right-most MPO
    //tensor included in LH and the left-most MPO tensor included in RH.
    //
    LocalMPO(MPOt<Tensor> const& H, 
             Tensor const& LH, 
             int LHlim,
             Tensor const& RH,
             int RHlim,
             Args const& args = Args::global());

    //
    // Sparse Matrix Methods
    //

    void
    product(const Tensor& phi, Tensor& phip) const;

    Real
    expect(const Tensor& phi) const { return lop_.expect(phi); }

    Tensor
    deltaRho(const Tensor& AA, 
             const Tensor& comb, Direction dir) const
        { return lop_.deltaRho(AA,comb,dir); }

    Tensor
    diag() const { return lop_.diag(); }

    //
    // position(b,psi) uses the MPS psi
    // to adjust the edge tensors such
    // that the MPO tensors at positions
    // b and b+1 are exposed
    //
    template <class MPSType>
    void
    position(int b, const MPSType& psi);

    int
    position() const;

    void
    shift(int j, Direction dir, const Tensor& A);

    //
    // Accessor Methods
    //

    void
    reset()
        {
        LHlim_ = 0;
        RHlim_ = Op_->N()+1;
        }

    Tensor const&
    L() const { return PH_[LHlim_]; }
    // Replace left edge tensor at current bond
    void
    L(Tensor const& nL) { PH_[LHlim_] = nL; }
    // Replace left edge tensor bordering site j
    // (so that nL includes sites < j)
    void
    L(int j, Tensor const& nL);


    Tensor const&
    R() const { return PH_[RHlim_]; }
    // Replace right edge tensor at current bond
    void
    R(Tensor const& nR) { PH_[RHlim_] = nR; }
    // Replace right edge tensor bordering site j
    // (so that nR includes sites > j)
    void
    R(int j, Tensor const& nR);

    const MPOt<Tensor>&
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
        if(val < 1) Error("numCenter must be set >= 1");
        nc_ = val; 
        }

    long
    size() const { return lop_.size(); }

    explicit operator bool() const { return Op_ != 0 || Psi_ != 0; }

    bool
    doWrite() const { return do_write_; }
    void
    doWrite(bool val) 
        { 
        if(Psi_ != 0)
            Error("Write to disk not yet supported for LocalMPO initialized with an MPS");
        if(!do_write_ && (val == true))
            initWrite(); 
        do_write_ = val; 
        }

    const std::string&
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

    const MPOt<Tensor>* Op_;
    std::vector<Tensor> PH_;
    int LHlim_,RHlim_;
    int nc_;

    LocalOp<Tensor> lop_;

    bool do_write_;
    std::string writedir_;

    const MPSt<Tensor>* Psi_;

    //
    /////////////////

    template <class MPSType>
    void
    makeL(const MPSType& psi, int k);

    template <class MPSType>
    void
    makeR(const MPSType& psi, int k);

    void
    setLHlim(int val);

    void
    setRHlim(int val);

    void
    initWrite();

    std::string
    PHFName(int j) const
        {
        return format("%s/PH_%03d",writedir_,j);
        }

    };

template <class Tensor>
inline LocalMPO<Tensor>::
LocalMPO()
    : Op_(0),
      LHlim_(-1),
      RHlim_(-1),
      nc_(2),
      do_write_(false),
      writedir_("."),
      Psi_(0)
    { }

template <class Tensor>
inline LocalMPO<Tensor>::
LocalMPO(const MPOt<Tensor>& H, 
         const Args& args)
    : Op_(&H),
      PH_(H.N()+2),
      LHlim_(0),
      RHlim_(H.N()+1),
      nc_(2),
      do_write_(false),
      writedir_("."),
      Psi_(0)
    { 
    if(args.defined("NumCenter"))
        numCenter(args.getInt("NumCenter"));
    }

template <class Tensor>
inline LocalMPO<Tensor>::
LocalMPO(const MPSt<Tensor>& Psi, 
         const Args& args)
    : Op_(0),
      PH_(Psi.N()+2),
      LHlim_(0),
      RHlim_(Psi.N()+1),
      nc_(2),
      do_write_(false),
      writedir_("."),
      Psi_(&Psi)
    { 
    if(args.defined("NumCenter"))
        numCenter(args.getInt("NumCenter"));
    }

template <class Tensor>
inline LocalMPO<Tensor>::
LocalMPO(const MPOt<Tensor>& H, 
         const Tensor& LH, const Tensor& RH,
         const Args& args)
    : Op_(&H),
      PH_(H.N()+2),
      LHlim_(0),
      RHlim_(H.N()+1),
      nc_(2),
      do_write_(false),
      writedir_("."),
      Psi_(0)
    { 
    PH_[0] = LH;
    PH_[H.N()+1] = RH;
    if(H.N()==2)
        lop_.update(Op_->A(1),Op_->A(2),L(),R());
    if(args.defined("NumCenter"))
        numCenter(args.getInt("NumCenter"));
    }

template <class Tensor>
inline LocalMPO<Tensor>::
LocalMPO(const MPSt<Tensor>& Psi, 
         const Tensor& LP,
         const Tensor& RP,
         const Args& args)
    : Op_(0),
      PH_(Psi.N()+2),
      LHlim_(0),
      RHlim_(Psi.N()+1),
      nc_(2),
      do_write_(false),
      writedir_("."),
      Psi_(&Psi)
    { 
    PH_[0] = LP;
    PH_[Psi.N()+1] = RP;
    if(args.defined("NumCenter"))
        numCenter(args.getInt("NumCenter"));
    }

template <class Tensor>
inline LocalMPO<Tensor>::
LocalMPO(MPOt<Tensor> const& H, 
         Tensor const& LH, 
         int LHlim,
         Tensor const& RH,
         int RHlim,
         Args const& args)
    : Op_(&H),
      PH_(H.N()+2),
      LHlim_(LHlim),
      RHlim_(RHlim),
      nc_(2),
      do_write_(false),
      writedir_("."),
      Psi_(0)
    { 
    PH_.at(LHlim) = LH;
    PH_.at(RHlim) = RH;
    if(H.N()==2) lop_.update(Op_->A(1),Op_->A(2),L(),R());
    if(args.defined("NumCenter")) numCenter(args.getInt("NumCenter"));
    }

template <class Tensor> inline
void LocalMPO<Tensor>::
product(const Tensor& phi, Tensor& phip) const
    {
    if(Op_ != 0)
        {
        lop_.product(phi,phip);
        }
    else 
    if(Psi_ != 0)
        {
        int b = position();
        auto othr = (!L() ? dag(prime(Psi_->A(b),Link)) : L()*dag(prime(Psi_->A(b),Link)));
        auto othrR = (!R() ? dag(prime(Psi_->A(b+1),Link)) : R()*dag(prime(Psi_->A(b+1),Link)));
        othr *= othrR;
        auto z = (othr*phi).cplx();

        phip = dag(othr);
        phip *= z;
        }
    else
        {
        Error("LocalMPO is null");
        }
    }

template <class Tensor>
void inline LocalMPO<Tensor>::
L(int j, const Tensor& nL)
    {
    if(LHlim_ > j-1) setLHlim(j-1);
    PH_[LHlim_] = nL;
    }

template <class Tensor>
void inline LocalMPO<Tensor>::
R(int j, const Tensor& nR)
    {
    if(RHlim_ < j+1) setRHlim(j+1);
    PH_[RHlim_] = nR;
    }

template <class Tensor>
template <class MPSType> 
inline void LocalMPO<Tensor>::
position(int b, const MPSType& psi)
    {
    if(!(*this)) Error("LocalMPO is null");

    makeL(psi,b-1);
    makeR(psi,b+nc_);

    setLHlim(b-1); //not redundant since LHlim_ could be > b-1
    setRHlim(b+nc_); //not redundant since RHlim_ could be < b+nc_

#ifdef DEBUG
    if(nc_ != 2)
        {
        Error("LocalOp only supports 2 center sites currently");
        }
#endif

    if(Op_ != 0) //normal MPO case
        {
        lop_.update(Op_->A(b),Op_->A(b+1),L(),R());
        }
    }

template <class Tensor>
int inline LocalMPO<Tensor>::
position() const
    {
    if(RHlim_-LHlim_ != (nc_+1))
        {
        throw ITError("LocalMPO position not set");
        }
    return LHlim_+1;
    }

template <class Tensor>
inline void LocalMPO<Tensor>::
shift(int j, Direction dir, const Tensor& A)
    {
    if(!(*this)) Error("LocalMPO is null");

#ifdef DEBUG
    if(nc_ != 2)
        {
        Error("LocalOp only supports 2 center sites currently");
        }
#endif

    if(dir == Fromleft)
        {
        if((j-1) != LHlim_)
            {
            std::cout << "j-1 = " << (j-1) << ", LHlim = " << LHlim_ << std::endl;
            Error("Can only shift at LHlim");
            }
        Tensor& E = PH_.at(LHlim_);
        Tensor& nE = PH_.at(j);
        nE = E * A;
        nE *= Op_->A(j);
        nE *= dag(prime(A));
        setLHlim(j);
        setRHlim(j+nc_+1);

        lop_.update(Op_->A(j+1),Op_->A(j+2),L(),R());
        }
    else //dir == Fromright
        {
        if((j+1) != LHlim_)
            {
            std::cout << "j+1 = " << (j+1) << ", RHlim_ = " << RHlim_ << std::endl;
            Error("Can only shift at RHlim_");
            }
        Tensor& E = PH_.at(RHlim_);
        Tensor& nE = PH_.at(j);
        nE = E * A;
        nE *= Op_->A(j);
        nE *= dag(prime(A));
        setLHlim(j-nc_-1);
        setRHlim(j);

        lop_.update(Op_->A(j-1),Op_->A(j),L(),R());
        }
    }

template <class Tensor>
template <class MPSType> 
inline void LocalMPO<Tensor>::
makeL(const MPSType& psi, int k)
    {
    if(!PH_.empty())
        {
        if(Op_ == 0) //Op is actually an MPS
            {
            while(LHlim_ < k)
                {
                const int ll = LHlim_;
                PH_.at(ll+1) = (!PH_.at(ll) ? psi.A(ll+1) : PH_[ll]*psi.A(ll+1));
                PH_[ll+1] *= dag(prime(Psi_->A(ll+1),Link));
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
                    PH_.at(ll+1) = PH_.at(ll)*psi.A(ll+1);
                    }
                else
                    {
                    PH_.at(ll+1) = psi.A(ll+1);
                    }
                PH_.at(ll+1) *= Op_->A(ll+1);
                PH_.at(ll+1) *= dag(prime(psi.A(ll+1)));
                setLHlim(ll+1);
                }
            }
        }
    }

template <class Tensor>
template <class MPSType> 
inline void LocalMPO<Tensor>::
makeR(const MPSType& psi, int k)
    {
    if(!PH_.empty())
        {
        if(Op_ == 0) //Op is actually an MPS
            {
            while(RHlim_ > k)
                {
                const int rl = RHlim_;
                PH_.at(rl-1) = (!PH_.at(rl) ? psi.A(rl-1) : PH_[rl]*psi.A(rl-1));
                PH_[rl-1] *= dag(prime(Psi_->A(rl-1),Link));
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
                //Print(psi.A(rl-1));
                if(PH_.at(rl))
                    {
                    PH_.at(rl-1) = PH_.at(rl)*psi.A(rl-1);
                    }
                else
                    {
                    PH_.at(rl-1) = psi.A(rl-1);
                    }
                PH_.at(rl-1) *= Op_->A(rl-1);
                PH_.at(rl-1) *= dag(prime(psi.A(rl-1)));
                //printfln("PH[%d] = \n%s",rl-1,PH_.at(rl-1));
                //PAUSE
                setRHlim(rl-1);
                }
            }
        }
    }

template <class Tensor>
void inline LocalMPO<Tensor>::
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
        PH_.at(LHlim_) = Tensor();
        }
    LHlim_ = val;
    if(LHlim_ < 1) 
        {
        //Set to null tensor and return
        PH_.at(LHlim_) = Tensor();
        return;
        }
    if(!PH_.at(LHlim_))
        {
        std::string fname = PHFName(LHlim_);
        readFromFile(fname,PH_.at(LHlim_));
        }
    }

template <class Tensor>
void inline LocalMPO<Tensor>::
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
        PH_.at(RHlim_) = Tensor();
        }
    RHlim_ = val;
    if(RHlim_ > Op_->N()) 
        {
        //Set to null tensor and return
        PH_.at(RHlim_) = Tensor();
        return;
        }
    if(!PH_.at(RHlim_))
        {
        std::string fname = PHFName(RHlim_);
        readFromFile(fname,PH_.at(RHlim_));
        }
    }

template <class Tensor>
void inline LocalMPO<Tensor>::
initWrite()
    {
    std::string global_write_dir = Global::args().getString("WriteDir","./");
    writedir_ = mkTempDir("PH",global_write_dir);
    }

} //namespace itensor


#endif
