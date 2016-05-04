//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_DECOMP_H
#define __ITENSOR_DECOMP_H
#include "itensor/iqtensor.h"
#include "itensor/spectrum.h"
#include "itensor/mps/localop.h"


namespace itensor {


//
// Singular value decomposition (SVD)
//
// Factors a tensor AA such that AA=U*D*V
// with D diagonal, real, and non-negative.
//
template<class Tensor>
Spectrum 
svd(Tensor AA, Tensor& U, Tensor& D, Tensor& V, 
    Args args = Global::args());


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
template<typename Tensor>
void
factor(Tensor const& T,
       Tensor      & A,
       Tensor      & B,
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

template<class Tensor>
Spectrum 
denmatDecomp(Tensor const& AA, 
             Tensor & A, 
             Tensor & B, 
             Direction dir, 
             Args const& args = Global::args())
    {
    return denmatDecomp(AA,A,B,dir,LocalOp<Tensor>{},args);
    }

//Density matrix decomp with BigMatrixT object supporting the noise term
//The BigMatrixT argument PH has to provide the deltaRho method
//to enable the noise term feature (see localop.h for example)
template<class Tensor, class BigMatrixT>
Spectrum 
denmatDecomp(Tensor const& AA, 
             Tensor & A, 
             Tensor & B, 
             Direction dir, 
             BigMatrixT const& PH,
             Args args = Global::args());



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
template<class I>
Spectrum 
diagHermitian(ITensorT<I> const& M, 
              ITensorT<I>      & U, 
              ITensorT<I>      & D,
              Args args = Args::global());


template<typename I>
ITensorT<I>
expHermitian(ITensorT<I> const& T, Cplx t = 1.);


//
// Orthogonal decomposition
//
// Given a tensor T, decomposes it into two tensors A and B
// such that T=A*B. If dir==Fromleft, A is guaranteed to be
// real and orthogonal, similar for B if dir==Fromright.
//
template<class Tensor>
Spectrum 
orthoDecomp(Tensor T, Tensor& A, Tensor& B, 
            Direction dir, 
            const Args& args = Global::args());


//
// Inverse Canonical SVD
//
// Factors a tensor AA such that AA=L*V*R
// where V is the inverse of the diagonal tensor
// appearing in the SVD
//

template<class Tensor>
Spectrum 
csvd(const Tensor& AA, Tensor& L, Tensor& V, Tensor& R, 
     const Args& args = Global::args());


//
//
// Eigenvalues and eigenvectors
//
// Computes eigenvalues V and eigenvectors D of an arbitrary tensor T.
//
// T must be "square-matrix-like" in the sense that
// T has only indices I,J,K,... and indices I',J',K',...
//
// D is a diagonal rank 2 tensor (matrix) containing the eigenvalues.
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
template<class I>
void 
eigen(ITensorT<I> const& T, 
      ITensorT<I> & V, 
      ITensorT<I> & D,
      Args const& args = Args::global());

template<class I>
void 
eigDecomp(ITensorT<I> const& T, 
          ITensorT<I> & R, 
          ITensorT<I> & D,
          ITensorT<I> & Rinv, 
          Args const& args = Args::global());


///////////////////////////
//
// Implementation (non-template parts in decomp.cc)
//
//////////////////////////


template<typename IndexT>
Spectrum 
svdRank2(ITensorT<IndexT> const& A, 
         IndexT const& ui, 
         IndexT const& vi,
         ITensorT<IndexT> & U, 
         ITensorT<IndexT> & D, 
         ITensorT<IndexT> & V,
         Args args = Args::global());

template<class Tensor>
Spectrum 
svd(Tensor AA, 
    Tensor & U, 
    Tensor & D, 
    Tensor & V, 
    Args args)
    {
    using IndexT = typename Tensor::index_type;

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
    //into a rank 2 tensor
    std::vector<IndexT> Uinds, 
                        Vinds;
    Uinds.reserve(AA.r());
    Vinds.reserve(AA.r());
    //Divide up indices based on U
    //If U is null, use V instead
    auto &L = (U ? U : V);
    auto &Linds = (U ? Uinds : Vinds),
         &Rinds = (U ? Vinds : Uinds);
    for(const auto& I : AA.inds())
        { 
        if(hasindex(L,I)) Linds.push_back(I);
        else              Rinds.push_back(I);
        }
    Tensor Ucomb,
           Vcomb;
    if(!Uinds.empty())
        {
        Ucomb = combiner(std::move(Uinds),{"IndexName","uc"});
        AA *= Ucomb;
        }
    if(!Vinds.empty())
        {
        Vcomb = combiner(std::move(Vinds),{"IndexName","vc"});
        AA *= Vcomb;
        }

    if(useOrigM)
        {
        //Try to determine current m,
        //then set minm_ and maxm_ to this.
        args.add("Cutoff",-1);
        long minm = 1,
             maxm = MAX_M;
        if(D.r() == 0)
            {
            auto mid = commonIndex(U,V,Link);
            if(mid) minm = maxm = mid.m();
            else    minm = maxm = 1;
            }
        else
            {
            minm = maxm = D.inds().front().m();
            }
        args.add("Minm",minm);
        args.add("Maxm",maxm);
        }

    auto ui = commonIndex(AA,Ucomb);
    auto vi = commonIndex(AA,Vcomb);

    auto spec = svdRank2(AA,ui,vi,U,D,V,args);

    U = dag(Ucomb) * U;
    V = V * dag(Vcomb);

    return spec;
    } //svd

template<class Tensor>
Spectrum 
csvd(const Tensor& AA, Tensor& L, Tensor& V, Tensor& R, 
     const Args& args)
    {
    /*
    Tensor UU(L),VV(R);
    Tensor D(V);
    Spectrum spec = svd(AA,UU,D,VV,args);

    L = UU*D;
    R = D*VV;

    V = dag(D);
    V.pseudoInvert(0);
    return spec;
    */
    //TODO remove this line:
    return Spectrum();
    }

template<class Tensor, class BigMatrixT>
Spectrum 
denmatDecomp(Tensor const& AA, 
             Tensor & A, 
             Tensor & B, 
             Direction dir, 
             BigMatrixT const& PH,
             Args args)
    {
    using IndexT = typename Tensor::index_type;

    auto noise = args.getReal("Noise",0.);

    auto mid = commonIndex(A,B,Link);

    //If dir==NoDir, put the O.C. on the side
    //that keeps mid's arrow the same
    if(dir == NoDir)
        {
        dir = (mid.dir() == Out ? Fromright : Fromleft);
        }

    auto& to_orth = (dir==Fromleft ? A : B);
    auto& newoc   = (dir==Fromleft ? B : A);
    
    auto& activeInds = (to_orth ? to_orth : AA).inds();

    std::vector<IndexT> cinds;
    cinds.reserve(activeInds.r());
    for(auto& I : activeInds)
        if(!hasindex(newoc,I))
            cinds.push_back(I);

    //Apply combiner
    START_TIMER(8)
    auto iname = args.getString("IndexName",mid ? mid.rawname() : "mid");
    auto cmb = combiner(std::move(cinds),iname);
    auto ci = cmb.inds().front();

    auto AAc = cmb * AA;

    //Form density matrix
    auto rho = AAc*dag(prime(AAc,ci)); 


    //Add noise term if requested
    if(noise > 0 && PH)
        {
        rho += noise*PH.deltaRho(AA,cmb,dir);
        auto tr = (delta(dag(ci),prime(ci))*realPart(rho)).real();
        if(tr > 1E-16) rho *= 1./tr;
        }

    STOP_TIMER(8)

    if(args.getBool("UseOrigM",false))
        {
        args.add("Cutoff",-1);
        args.add("Minm",mid.m());
        args.add("Maxm",mid.m());
        }

    if(args.getBool("TraceReIm",false))
        {
        rho = realPart(rho);
        }

    Tensor U,D;
    args.add("Truncate",true);
    auto spec = diag_hermitian(rho,U,D,args);

    cmb.dag();

    to_orth = cmb * dag(U);
    newoc = U * AAc;

    return spec;

    } //denmatDecomp


template<typename I>
Spectrum
diag_hermitian(ITensorT<I>    rho, 
               ITensorT<I>  & U, 
               ITensorT<I>  & D,
               Args const& args);


template<class I>
Spectrum 
diagHermitian(ITensorT<I> const& M, 
              ITensorT<I>      & U, 
              ITensorT<I>      & D,
              Args args)
    {
    if(!args.defined("IndexName")) args.add("IndexName","d");

    auto inds = stdx::reserve_vector<I>(rank(M)/2);
    for(auto& i : M.inds())
        { 
        if(i.primeLevel() == 0) inds.push_back(i);
        }

    auto comb = combiner(std::move(inds),args);
    auto Mc = M*comb;

    auto combP = dag(prime(comb));
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


template<class Tensor>
Spectrum 
orthoDecomp(Tensor T, Tensor& A, Tensor& B, 
            Direction dir, 
            const Args& args)
    {
    /*
    using IndexT = typename Tensor::IndexT;

    if(isZero(T,Args("Fast"))) 
        throw ResultIsZero("orthoDecomp: T is zero");

    const
    bool usedenmat = false;

    Spectrum spec;

    if(usedenmat)
        {
        spec = denmatDecomp(T,A,B,dir,args + Args("TraceReIm",true,"Noise",0));
        }
    else //use svd
        {
        //Combiners which transform T
        //into a rank 2 tensor
        CombinerT Acomb, Bcomb;

        const
        IndexT reim = IQIndex("ReIm",Index("reim",2),QN());

        //Divide up indices based on U
        //If U is null, use V instead
        const Tensor &L = (A ? A : B);
        CombinerT &Lcomb = (A ? Acomb : Bcomb),
                  &Rcomb = (A ? Bcomb : Acomb);
        for(const IndexT& I : T.indices())
            { 
            if(hasindex(L,I))
                Lcomb.addleft(I);
            else
                Rcomb.addleft(I);
            }

        if(dir == Fromleft)
            Rcomb.addleft(reim);
        else
            Lcomb.addleft(reim);

        T = realPart(T)*reim(1) + imagPart(T)*reim(2);

        T = Acomb * T * Bcomb;


        Tensor D;
        spec = svdRank2(T,Acomb.right(),Bcomb.right(),A,D,B,args);

        A = dag(Acomb) * A;
        B = B * dag(Bcomb);

        if(dir==Fromleft) 
            {
            B *= D;
            B = B*dag(reim)(1) + Complex_i*B*dag(reim)(2);
            }
        else              
            {
            A *= D;
            A = A*dag(reim)(1) + Complex_i*A*dag(reim)(2);
            }
        }

    return spec;
    */

    //TODO remove this line:
    return Spectrum();

    } //orthoDecomp

//Return value is: (trunc_error,docut)
std::tuple<Real,Real>
truncate(Vector & P,
         long maxm,
         long minm,
         Real cutoff,
         bool absoluteCutoff = false,
         bool doRelCutoff = false);

template<typename V>
MatRefc<V>
toMatRefc(ITensor const& T, 
          Index const& i1, 
          Index const& i2);

template<typename T>
struct GetBlocks
    {
    using value_type = T;
    IQIndexSet const& is;
    bool transpose = false;

    GetBlocks(IQIndexSet const& is_, 
              IQIndex const& i1_, 
              IQIndex const& i2_)
      : is(is_)
        { 
        if(is.r() != 2) Error("GetBlocks only supports rank 2 currently");
        transpose = (i2_ == is.front());
        }
    };

template<typename T>
struct Rank2Block
    {
    MatRefc<T> M;
    long i1 = 0,
         i2 = 0;
    };

template<typename T>
std::vector<Rank2Block<T>>
doTask(GetBlocks<T> const& G, 
       QDense<T> const& d);

void
showEigs(Vector const& P,
         Real truncerr,
         LogNum const& scale,
         Args const& args);

} //namespace itensor


#endif
