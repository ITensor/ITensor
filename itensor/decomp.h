//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_DECOMP_H
#define __ITENSOR_DECOMP_H
#include "itensor/spectrum.h"
#include "itensor/itensor.h"


namespace itensor {

//
// Singular value decomposition (SVD)
//
// Factors a tensor AA such that AA=U*D*V
// with D diagonal, real, and non-negative.
//
Spectrum 
svd(ITensor AA, ITensor& U, ITensor& D, ITensor& V, 
    Args args = Args::global());


//
// The "factor" decomposition is based on the SVD,
// but factorizes a tensor T into only two
// tensors T=A*B where A and B share a single
// common index.
//
// If the SVD of T is T=U*S*V where S is a diagonal
// matrix of singular values, then A and B
// are schematically A=U*sqrt(S) and B=sqrt(S)*V.
//
// In addition to the named Args recognized by the 
// svd routine, factor accepts an Arg "IndexName"
// which will be the name of the common index 
// connecting A and B.
// 
Spectrum
factor(ITensor const& T,
       ITensor      & A,
       ITensor      & B,
       Args const& args = Args::global());

//
// Density Matrix Decomposition
// 
// Factors a tensor AA such that AA = A*B.
// Result is equivalent to SVD such that AA = U*D*V where if
// dir==Fromleft, A=U and B=(D*V) or, if dir==Fromright, A=(U*D) and B=V.
// Implementation is faster than SVD, though, and allows the
// noise term to be used.
//
// To determine which indices end up on which factors (i.e. on A versus B),
// the method examines the initial indices of A and B.
// If a given index is present on, say, A, then it will on A 
// upon return (although the elements of A will be overwritten and other indices
// may be added to it). Any indices not initially present on A or B 
// will end up on B if dir==Fromleft or on A if dir==Fromright.
//

//Class that acts as an empty LocalOp
//to provide a default implementation for denmatDecomp
//case where noise is zero / off
class NoOp
    {
    public:

    NoOp() { }

    explicit operator bool() const { return false; }

    ITensor
    deltaRho(ITensor const& rho,
             ITensor const& combine,
             Direction dir) const 
        { 
        Error("Non-zero noise not supported without providing a Hamiltonian");
        return rho; 
        }
    };

//Density matrix decomp with BigMatrixT object supporting the noise term
//The BigMatrixT argument PH has to provide the deltaRho method
//to enable the noise term feature (see localop.h for example)
template<class BigMatrixT>
Spectrum 
denmatDecomp(ITensor const& AA, 
             ITensor & A, 
             ITensor & B, 
             Direction dir, 
             BigMatrixT const& PH,
             Args args = Args::global());

Spectrum inline
denmatDecomp(ITensor const& AA, 
             ITensor & A, 
             ITensor & B, 
             Direction dir, 
             Args const& args = Args::global())
    {
    return denmatDecomp(AA,A,B,dir,NoOp(),args);
    }



//
// Hermitian eigenvalue decomposition / diagonalization
//
// Assumes input is a Hermitian tensor with indices
// i,j,k,.... and i',j',k',...
// (tensor must be conjugate symmetric under
//  exchange primed and unprimed indices)
// Result is unitary tensor U and diagonal sparse tensor D
// such that M == dag(U)*D*prime(U)
//
Spectrum 
diagHermitian(ITensor const& M, 
              ITensor      & U, 
              ITensor      & D,
              Args args = Args::global());


ITensor
expHermitian(ITensor const& T, Cplx t = 1.);

//
//
// Eigenvalues and eigenvectors
//
// Computes eigenvalues V and eigenvectors D of an arbitrary tensor T.
//
// T must be "square-matrix-like" in the sense that
// T has only indices I,J,K,... and indices I',J',K',...
//
// D is a diagonal order 2 tensor (matrix) containing the eigenvalues.
// On return, V has the "column" indices of T and a new index shared with D
// (the index labeled "C" below).
//
// The result is such that V and D give:
//       _         _                 _ 
// I'-<-| |-<-I-<-| |          I'-<-| |     
//      |T|       |V|-<-C  ==       |V|-<-C'-<-(D)-<-C
// J'-<-|_|-<-J-<-|_|          J'-<-|_|    
//
//
void 
eigen(ITensor const& T, 
      ITensor & V, 
      ITensor & D,
      Args const& args = Args::global());

//template<class I>
//void 
//eigDecomp(ITensorT<I> const& T, 
//          ITensorT<I> & R, 
//          ITensorT<I> & D,
//          ITensorT<I> & Rinv, 
//          Args const& args = Args::global());


///////////////////////////
//
// Implementation (non-template parts in decomp.cc)
//
//////////////////////////


Spectrum 
svdOrd2(ITensor const& A, 
         Index const& ui, 
         Index const& vi,
         ITensor & U, 
         ITensor & D, 
         ITensor & V,
         Args args = Args::global());

Spectrum inline
svd(ITensor AA, 
    ITensor & U, 
    ITensor & D, 
    ITensor & V, 
    Args args)
    {
#ifdef DEBUG
    if(!U && !V) 
        Error("U and V default-initialized in svd, must indicate at least one index on U or V");
#endif

    auto noise = args.getReal("Noise",0);
    auto useOrigM = args.getBool("UseOrigM",false);

    if(noise > 0)
        Error("Noise term not implemented for svd");
    
    //if(isZero(AA,Args("Fast"))) 
    //    throw ResultIsZero("svd: AA is zero");


    //Combiners which transform AA
    //into a order 2 tensor
    std::vector<Index> Uinds, 
                        Vinds;
    Uinds.reserve(AA.order());
    Vinds.reserve(AA.order());
    //Divide up indices based on U
    //If U is null, use V instead
    auto &L = (U ? U : V);
    auto &Linds = (U ? Uinds : Vinds),
         &Rinds = (U ? Vinds : Uinds);
    for(const auto& I : AA.inds())
        { 
        if(hasIndex(L,I)) Linds.push_back(I);
        else              Rinds.push_back(I);
        }
    ITensor Ucomb,
            Vcomb;
    if(!Uinds.empty())
        {
        //TODO: Add some tags here?
        Ucomb = combiner(std::move(Uinds)); //,{"IndexName","uc"});
        AA *= Ucomb;
        }
    if(!Vinds.empty())
        {
        //TODO: Add some tags here?
        Vcomb = combiner(std::move(Vinds)); //,{"IndexName","vc"});
        AA *= Vcomb;
        }

    if(useOrigM)
        {
        //Try to determine current m,
        //then set minm_ and maxm_ to this.
        args.add("Cutoff",-1);
        long minm = 1,
             maxm = MAX_M;
        if(D.order() == 0)
            {
            //auto mid = commonIndex(U,V,Link);
            //TODO: check this does the same thing
            auto mid = commonIndex(U,V,"Link");
            if(mid) minm = maxm = dim(mid);
            else    minm = maxm = 1;
            }
        else
            {
            minm = maxm = dim(D.inds().front());
            }
        args.add("Minm",minm);
        args.add("Maxm",maxm);
        }

    auto ui = commonIndex(AA,Ucomb);
    auto vi = commonIndex(AA,Vcomb);

    auto spec = svdOrd2(AA,ui,vi,U,D,V,args);

    U = dag(Ucomb) * U;
    V = V * dag(Vcomb);

    return spec;
    } //svd

Spectrum
diag_hermitian(ITensor  rho, 
               ITensor  & U, 
               ITensor  & D,
               Args const& args);

template<class BigMatrixT>
Spectrum 
denmatDecomp(ITensor const& AA, 
             ITensor & A, 
             ITensor & B, 
             Direction dir, 
             BigMatrixT const& PH,
             Args args)
    {
    auto noise = args.getReal("Noise",0.);

    //TODO: try to avoid using "Link" here
    auto mid = commonIndex(A,B,"Link");

    //If dir==NoDir, put the O.C. on the side
    //that keeps mid's arrow the same
    if(dir == NoDir)
        {
        dir = (mid.dir() == Out ? Fromright : Fromleft);
        }

    auto& to_orth = (dir==Fromleft ? A : B);
    auto& newoc   = (dir==Fromleft ? B : A);
    
    auto& activeInds = (to_orth ? to_orth : AA).inds();

    auto cinds = stdx::reserve_vector<Index>(activeInds.order());
    for(auto& I : activeInds)
        {
        if(!hasIndex(newoc,I)) cinds.push_back(I);
        }

    //Apply combiner
    START_TIMER(8)
    //TODO: decide on a tag convention for denmatDecomp
    auto itagset = getTagSet(args,"Tags","Link");
    args.add("Tags",toString(itagset));
    auto cmb = combiner(std::move(cinds),{"Tags",toString(itagset)});
    auto ci = cmb.inds().front();

    auto AAc = cmb * AA;

    //Form density matrix
    auto rho = AAc*dag(prime(AAc,ci)); 


    //Add noise term if requested
    if(noise > 0 && PH)
        {
        rho += noise*PH.deltaRho(AA,cmb,dir);
        //println("delta(dag(ci),prime(ci)) = ",delta(dag(ci),prime(ci)));
        //print("realPart(rho) = ",realPart(rho));
        auto tr = (delta(dag(ci),prime(ci))*realPart(rho)).real();
        if(tr > 1E-16) rho *= 1./tr;
        }

    STOP_TIMER(8)

    if(args.getBool("UseOrigM",false))
        {
        args.add("Cutoff",-1);
        args.add("Minm",dim(mid));
        args.add("Maxm",dim(mid));
        }

    if(args.getBool("TraceReIm",false))
        {
        rho = realPart(rho);
        }

    ITensor U,D;
    args.add("Truncate",true);
    auto spec = diag_hermitian(rho,U,D,args);

    cmb.dag();

    to_orth = cmb * dag(U);
    newoc = U * AAc;

    return spec;

    } //denmatDecomp




Spectrum inline
diagHermitian(ITensor const& M, 
              ITensor      & U, 
              ITensor      & D,
              Args args)
    {
    //TODO: create tag convention
    if(!args.defined("Tags")) args.add("Tags","Link");

    //
    // Pick an arbitrary index and do some analysis
    // on its prime level spacing
    //
    auto k = M.inds().front();
    auto kps = stdx::reserve_vector<int>(order(M));
    for(auto& i : M.inds()) if(equalsIgnorePrime(i,k)) kps.push_back(i.primeLevel());
    if(kps.size() <= 1ul || kps.size()%2 != 0ul) 
        {
        Error("Input tensor to diagHermitian should have pairs of indices with equally spaced prime levels");
        }
    auto nk = kps.size();
    std::sort(kps.begin(),kps.end());
    //idiff == "inner" difference between cluster of low-prime-level copies
    //         of k, if more than one
    auto idiff = kps.at(nk/2-1)-kps.front();
    //mdiff == max prime-level difference of copies of k
    auto mdiff = kps.back()-kps.front();
    //pdiff == spacing between lower and higher prime level index pairs
    auto pdiff = mdiff-idiff;

    auto inds = stdx::reserve_vector<Index>(order(M)/2);
    for(auto& i : M.inds())
    for(auto& j : M.inds())
        {
        if(equalsIgnorePrime(i,j) && i.primeLevel()+pdiff == j.primeLevel())
            {
            inds.push_back(i);
            }
        }
    if(inds.empty() || order(M)/2 != (long)inds.size())
        {
        Error("Input tensor to diagHermitian should have pairs of indices with equally spaced prime levels");
        }

    auto comb = combiner(std::move(inds),args);
    auto Mc = M*comb;

    auto combP = dag(prime(comb,pdiff));
    try {
        Mc = combP * Mc;
        }
    catch(ITError const& e)
        {
        println("Diagonalize expects opposite arrow directions for primed and unprimed indices.");
        throw e;
        }

    auto spec = diag_hermitian(Mc,U,D,args);

    U = comb * U;

    return spec;
    } //diagHermitian


//Return value is: (trunc_error,docut)
std::tuple<Real,Real>
truncate(Vector & P,
         long maxm,
         long minm,
         Real cutoff,
         bool absoluteCutoff = false,
         bool doRelCutoff = false,
         Args const& args = Args::global());

template<typename V>
MatRefc<V>
toMatRefc(ITensor const& T, 
          Index const& i1, 
          Index const& i2);

template<typename T>
struct GetBlocks
    {
    using value_type = T;
    IndexSet const& is;
    bool transpose = false;

    GetBlocks(IndexSet const& is_, 
              Index const& i1_, 
              Index const& i2_)
      : is(is_)
        { 
        if(is.order() != 2) Error("GetBlocks only supports order 2 currently");
        transpose = (i2_ == is.front());
        }
    };

template<typename T>
struct Ord2Block
    {
    MatRefc<T> M;
    long i1 = 0,
         i2 = 0;
    };

template<typename T>
std::vector<Ord2Block<T>>
doTask(GetBlocks<T> const& G, 
       QDense<T> const& d);

void
showEigs(Vector const& P,
         Real truncerr,
         LogNum const& scale,
         Args const& args);

struct EigQN
    {
    Real eig = 0.;
    QN qn;

    EigQN() { }

    EigQN(Real eg,QN q) : eig(eg),qn(q) { }

    bool
    operator<(EigQN const& o) const { return eig < o.eig; }

    bool
    operator>(EigQN const& o) const { return eig > o.eig; }
    };

} //namespace itensor


#endif
