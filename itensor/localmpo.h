//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_LOCALMPO
#define __ITENSOR_LOCALMPO
#include "mpo.h"
#include "localop.h"

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

    typedef typename Tensor::CombinerT
    CombinerT;

    //
    // Constructors
    //

    LocalMPO();

    LocalMPO(const MPOt<Tensor>& Op, 
             const OptSet& opts = Global::opts());

    LocalMPO(const MPSt<Tensor>& Psi, 
             const OptSet& opts = Global::opts());

    //
    // Sparse Matrix Methods
    //

    void
    product(const Tensor& phi, Tensor& phip) const;

    Real
    expect(const Tensor& phi) const { return lop_.expect(phi); }

    Tensor
    deltaRho(const Tensor& AA, 
             const CombinerT& comb, Direction dir) const
        { return lop_.deltaRho(AA,comb,dir); }

    Tensor
    deltaPhi(const Tensor& phi) const { return lop_.deltaPhi(phi); }

    void
    diag(Tensor& D) const { lop_.diag(D); }

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
        RHlim_ = Op_->NN()+1;
        }

    const Tensor&
    L() const;

    //
    // Replace left edge tensor
    // at current bond
    //
    void
    L(const Tensor& nL);

    //
    // Replace left edge tensor 
    // bordering site j
    // (so that nL includes sites < j)
    //
    void
    L(int j, const Tensor& nL);


    const Tensor&
    R() const;

    //
    // Replace right edge tensor
    // at current bond
    //
    void
    R(const Tensor& nR);

    //
    // Replace right edge tensor 
    // bordering site j
    // (so that nR includes sites > j)
    //
    void
    R(int j, const Tensor& nR);


    const Tensor&
    bondTensor() const { return lop_.bondTensor(); }

    bool
    combineMPO() const { return lop_.combineMPO(); }
    void
    combineMPO(bool val) { lop_.combineMPO(val); }

    int
    numCenter() const { return nc_; }
    void
    numCenter(int val) 
        { 
        if(val < 1) Error("numCenter must be set >= 1");
        nc_ = val; 
        }

    int
    size() const { return lop_.size(); }

    bool
    isNull() const { return Op_ == 0 && Psi_ == 0; }
    bool
    isNotNull() const { return Op_ != 0 || Psi_ != 0; }

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

    static LocalMPO& Null()
        {
        static LocalMPO Null_;
        return Null_;
        }

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
        return (boost::format("%s/PH_%03d")%writedir_%j).str();
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
LocalMPO(const MPOt<Tensor>& Op, 
         const OptSet& opts)
    : Op_(&Op),
      PH_(Op.NN()+2),
      LHlim_(0),
      RHlim_(Op.NN()+1),
      nc_(2),
      do_write_(false),
      writedir_("."),
      Psi_(0)
    { 
    if(opts.defined("NumCenter"))
        numCenter(opts.intVal("NumCenter"));
    }

template <class Tensor>
inline LocalMPO<Tensor>::
LocalMPO(const MPSt<Tensor>& Psi, 
         const OptSet& opts)
    : Op_(0),
      PH_(Psi.NN()+2),
      LHlim_(0),
      RHlim_(Psi.NN()+1),
      nc_(2),
      do_write_(false),
      writedir_("."),
      Psi_(&Psi)
    { 
    if(opts.defined("NumCenter"))
        numCenter(opts.intVal("NumCenter"));
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
        Tensor othr = (L().isNull() ? primelink(Psi_->AA(b)) : L()*primelink(Psi_->AA(b)));
        othr *= primelink(Psi_->AA(b+1));
        if(R().isNotNull()) 
            othr *= R();

        Real re = 0, im = 0;
        BraKet(othr,phi,re,im);

        phip = othr;
        phip.mapprime(1,0);
        if(fabs(im) < 1E-10) 
            {
            phip *= re;
            }
        else
            {
            phip *= (re*Tensor::Complex_1() + im*Tensor::Complex_i());
            }
        }
    else
        {
        Error("LocalMPO is null");
        }
    }

template <class Tensor>
inline
const Tensor& LocalMPO<Tensor>::
L() const 
    { 
    return PH_[LHlim_];
    }

template <class Tensor>
void inline LocalMPO<Tensor>::
L(const Tensor& nL)
    {
    PH_[LHlim_] = nL;
    }

template <class Tensor>
void inline LocalMPO<Tensor>::
L(int j, const Tensor& nL)
    {
    if(LHlim_ > j-1) setLHlim(j-1);
    PH_[LHlim_] = nL;
    }

template <class Tensor>
inline
const Tensor& LocalMPO<Tensor>::
R() const 
    { 
    return PH_[RHlim_];
    }

template <class Tensor>
void inline LocalMPO<Tensor>::
R(const Tensor& nR)
    {
    PH_[RHlim_] = nR;
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
    if(this->isNull()) Error("LocalMPO is null");

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
        lop_.update(Op_->AA(b),Op_->AA(b+1),L(),R());
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
    if(this->isNull()) Error("LocalMPO is null");

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
        nE *= Op_->AA(j);
        nE *= conj(primed(A));
        setLHlim(j);
        setRHlim(j+nc_+1);

#ifdef DEBUG
        //std::cerr << boost::format("LocalMPO at (%d,%d) \n") % LHlim_ % RHlim_;
        //PrintIndices(L());
        //PrintIndices(R());
#endif

        lop_.update(Op_->AA(j+1),Op_->AA(j+2),L(),R());
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
        nE *= Op_->AA(j);
        nE *= conj(primed(A));
        setLHlim(j-nc_-1);
        setRHlim(j);

        lop_.update(Op_->AA(j-1),Op_->AA(j),L(),R());
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
                PH_.at(ll+1) = (PH_.at(ll).isNull() ? conj(psi.AA(ll+1)) : PH_[ll]*conj(psi.AA(ll+1)));
                PH_[ll+1] *= primelink(Psi_->AA(ll+1));
                setLHlim(LHlim_+1);
                }
            }
        else //normal MPO case
            {
            while(LHlim_ < k)
                {
                const int ll = LHlim_;
                psi.projectOp(ll+1,Fromleft,PH_.at(ll),Op_->AA(ll+1),PH_.at(ll+1));
                setLHlim(LHlim_+1);
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
                PH_.at(rl-1) = (PH_.at(rl).isNull() ? conj(psi.AA(rl-1)) : PH_[rl]*conj(psi.AA(rl-1)));
                PH_[rl-1] *= primelink(Psi_->AA(rl-1));
                setRHlim(RHlim_-1);
                }
            }
        else //normal MPO case
            {
            while(RHlim_ > k)
                {
                const int rl = RHlim_;
                psi.projectOp(rl-1,Fromright,PH_.at(rl),Op_->AA(rl-1),PH_.at(rl-1));
                setRHlim(RHlim_-1);
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

    if(LHlim_ != val && PH_.at(LHlim_).isNotNull())
        {
        //std::cerr << boost::format("Writing PH(%d) to %s\n")%LHlim_%writedir_;
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
    if(PH_.at(LHlim_).isNull())
        {
        std::string fname = PHFName(LHlim_);
        std::ifstream s(fname.c_str());
        if(s.good())
            {
            PH_.at(LHlim_).read(s);
            s.close();
            }
        else
            {
            std::cerr << boost::format("Tried to read file %s\n")%fname;
            Error("Missing file");
            }
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

    if(RHlim_ != val && PH_.at(RHlim_).isNotNull())
        {
        //std::cerr << boost::format("Writing PH(%d) to %s\n")%RHlim_%writedir_;
        writeToFile(PHFName(RHlim_),PH_.at(RHlim_));
        PH_.at(RHlim_) = Tensor();
        }
    RHlim_ = val;
    if(RHlim_ > Op_->NN()) 
        {
        //Set to null tensor and return
        PH_.at(RHlim_) = Tensor();
        return;
        }
    if(PH_.at(RHlim_).isNull())
        {
        std::string fname = PHFName(RHlim_);
        std::ifstream s(fname.c_str());
        if(s.good())
            {
            PH_.at(RHlim_).read(s);
            s.close();
            }
        else
            {
            std::cerr << boost::format("Tried to read file %s\n")%fname;
            Error("Missing file");
            }
        }
    }

template <class Tensor>
void inline LocalMPO<Tensor>::
initWrite()
    {
    std::string global_write_dir = Global::opts().stringOrDefault("WriteDir","./");
    writedir_ = mkTempDir("PH",global_write_dir);
    //std::cout << "Successfully created directory " + writedir_ << std::endl;
    }

#endif
