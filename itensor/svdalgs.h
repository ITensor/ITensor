//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SVDALGS_H
#define __ITENSOR_SVDALGS_H
#include "iqcombiner.h"
#include "iqtsparse.h"
#include "localop.h"
#include "spectrum.h"

#define Cout std::cout
#define Endl std::endl
#define Format boost::format


//
// Singular value decomposition (SVD)
//
// Factors a tensor AA such that AA=U*D*V
// with D diagonal, real, and non-negative.
//
template<class Tensor, class SparseT>
Spectrum
svd(const Tensor& AA, Tensor& U, SparseT& D, Tensor& V,
    const OptSet& opts = Global::opts())
    {
    Spectrum spec;
    svd(AA,U,D,V,spec,opts);
    return spec;
    }

//SVD with Spectrum object controlling truncation
//and storing eigs (squares of singular values) on return
template<class Tensor, class SparseT>
void 
svd(Tensor AA, Tensor& U, SparseT& D, Tensor& V, 
    Spectrum& spec,
    const OptSet& opts = Global::opts());

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
denmatDecomp(const Tensor& AA, Tensor& A, Tensor& B, Direction dir, 
             const OptSet& opts = Global::opts())
    {
    Spectrum spec;
    denmatDecomp(AA,A,B,dir,spec,LocalOp<Tensor>::Null(),opts);
    return spec;
    }


//Density matrix decomp with a Spectrum object controlling
//the truncation parameters and storing the eigs on return.
template<class Tensor>
void 
denmatDecomp(const Tensor& AA, Tensor& A, Tensor& B, Direction dir, 
             Spectrum& spec,
             const OptSet& opts = Global::opts())
    {
    denmatDecomp(AA,A,B,dir,spec,LocalOp<Tensor>::Null(),opts);
    }

//Density matrix decomp with a Spectrum object 
//and LocalOpT object supporting the noise term
//The LocalOpT argument PH has to provide the deltaRho method
//to enable the noise term feature (see localop.h for example)
template<class Tensor, class LocalOpT>
void 
denmatDecomp(const Tensor& AA, Tensor& A, Tensor& B, Direction dir, 
             Spectrum& spec, const LocalOpT& PH,
             const OptSet& opts = Global::opts());



//
// Hermitian eigenvalue decomposition / diagonalization
//
// Assumes input is a Hermitian tensor with indices
// i,j,k,.... and i',j',k',...
// (tensor must be conjugate symmetric under
//  exchange primed and unprimed indices)
// Result is unitary tensor U and diagonal sparse tensor D
// such that M == conj(U)*D*primed(U)
//
template<class Tensor, class SparseT>
Spectrum 
diagHermitian(const Tensor& M, Tensor& U, SparseT& D,
              const OptSet& opts = Global::opts())
    {
    Spectrum spec;
    diagHermitian(M,U,D,spec,opts);
    return spec;
    }

//Version accepting a Spectrum object which
//stores eigs on return
template<class Tensor, class SparseT>
void 
diagHermitian(const Tensor& M, Tensor& U, SparseT& D, 
              Spectrum& spec,
              const OptSet& opts = Global::opts());




//
// Orthogonal decomposition
//
// Given a tensor T, decomposes it into two tensors A and B
// such that T=A*B. If dir==Fromleft, A is guaranteed to be
// real and orthogonal, similar for B if dir==Fromright.
//
template<class Tensor>
Spectrum 
orthoDecomp(const Tensor& T, Tensor& A, Tensor& B, Direction dir, 
            const OptSet& opts = Global::opts())
    {
    Spectrum spec;
    orthoDecomp(T,A,B,dir,spec,opts);
    return spec;
    }

template<class Tensor>
void 
orthoDecomp(Tensor T, Tensor& A, Tensor& B, Direction dir, 
            Spectrum& spec,
            const OptSet& opts = Global::opts());




//
// Inverse Canonical SVD
//
// Factors a tensor AA such that AA=L*V*R
// where V is the inverse of the diagonal tensor
// appearing in the SVD
//
template<class Tensor, class SparseT>
Spectrum 
csvd(const Tensor& AA, Tensor& L, SparseT& V, Tensor& R, 
     const OptSet& opts = Global::opts())
    {
    Spectrum spec;
    csvd(AA,L,V,R,spec,opts);
    return spec;
    }

template<class Tensor, class SparseT>
void 
csvd(const Tensor& AA, Tensor& L, SparseT& V, Tensor& R, 
     Spectrum& spec, 
     const OptSet& opts = Global::opts());


///////////////////////////
//
// Implementation (non-template parts in svdalgs.cc)
//
//////////////////////////


void 
svdRank2(ITensor A, const Index& ui, const Index& vi,
         ITensor& U, ITSparse& D, ITensor& V, Spectrum& spec,
         const OptSet& opts = Global::opts());

void 
svdRank2(IQTensor A, const IQIndex& uI, const IQIndex& vI,
         IQTensor& U, IQTSparse& D, IQTensor& V, Spectrum& spec,
         const OptSet& opts = Global::opts());

template<class Tensor, class SparseT>
void 
svd(Tensor AA, Tensor& U, SparseT& D, Tensor& V, 
    Spectrum& spec, 
    const OptSet& opts)
    {
    typedef typename Tensor::IndexT 
    IndexT;
    typedef typename Tensor::CombinerT 
    CombinerT;
    
    if(isZero(AA,Opt("Fast"))) 
        throw ResultIsZero("svd: AA is zero");

    if(spec.noise() > 0)
        Error("Noise term not implemented for svd");

    //if(noise_ > 0 && !PH.isNull())
    //    {
    //    //Add in noise term
    //    Real orig_norm = AA.norm();
    //    AA += noise_*PH.deltaPhi(AA);
    //    AA *= orig_norm/AA.norm();
    //    }

    //Combiners which transform AA
    //into a rank 2 tensor
    CombinerT Ucomb, Vcomb;
    Ucomb.doCondense(true);
    Vcomb.doCondense(true);


    //Divide up indices based on U
    //If U is null, use V instead
    const Tensor &L = (U.isNull() ? V : U);
    CombinerT &Lcomb = (U.isNull() ? Vcomb : Ucomb),
              &Rcomb = (U.isNull() ? Ucomb : Vcomb);
    Foreach(const IndexT& I, AA.indices())
        { 
        if(hasindex(L,I))
            Lcomb.addleft(I);
        else
            Rcomb.addleft(I);
        }

    AA = Ucomb * AA * Vcomb;

    const Real saved_cutoff = spec.cutoff(); 
    const int saved_minm = spec.minm(),
              saved_maxm = spec.maxm(); 
    if(spec.useOrigM())
        {
        //Try to determine current m,
        //then set minm_ and maxm_ to this.
        spec.cutoff(-1);
        if(D.r() == 0)
            {
            IndexT mid = commonIndex(U,V,Link);
            if(!mid.isNull())
                {
                spec.minm(mid.m());
                spec.maxm(mid.m());
                }
            else
                {
                spec.minm(1);
                spec.maxm(1);
                }
            }
        else
            {
            spec.minm(D.index(1).m());
            spec.maxm(D.index(1).m());
            }
        }

    svdRank2(AA,Ucomb.right(),Vcomb.right(),U,D,V,spec,opts);

    spec.cutoff(saved_cutoff);
    spec.minm(saved_minm);
    spec.maxm(saved_maxm);

    U = conj(Ucomb) * U;
    V = V * conj(Vcomb);

    } //svd

template<class Tensor, class SparseT>
void 
csvd(const Tensor& AA, Tensor& L, SparseT& V, Tensor& R, 
     Spectrum& spec, 
     const OptSet& opts)
    {
    typedef typename Tensor::IndexT 
    IndexT;
    typedef typename Tensor::CombinerT 
    CombinerT;

    Tensor UU(L),VV(R);
    SparseT D(V);
    svd(AA,UU,D,VV,spec,opts);

    L = UU*D;
    R = D*VV;

    V = conj(D);
    V.pseudoInvert(0);
    }

Real 
diag_hermitian(ITensor rho, ITensor& U, ITSparse& D, Spectrum& spec,
               const OptSet& opts = Global::opts());

Real 
diag_hermitian(IQTensor rho, IQTensor& U, IQTSparse& D, Spectrum& spec,
               const OptSet& opts = Global::opts());

template<class Tensor, class LocalOpT>
void 
denmatDecomp(const Tensor& AA, Tensor& A, Tensor& B, Direction dir, 
             Spectrum& spec, const LocalOpT& PH,
             const OptSet& opts)
    {
    typedef typename Tensor::IndexT 
    IndexT;
    typedef typename Tensor::CombinerT 
    CombinerT;
    typedef typename Tensor::SparseT 
    SparseT;

    if(isZero(AA,Opt("Fast"))) 
        {
        throw ResultIsZero("denmatDecomp: AA is zero");
        }

    IndexT mid = commonIndex(A,B,Link);

    //If dir==None, put the O.C. on the side
    //that keeps mid's arrow the same
    if(dir == None)
        {
        //Cout << Format("Arrow before = %s")%(mid.dir() == Out ? "Out" : "In") << Endl;
        dir = (mid.dir() == Out ? Fromright : Fromleft);
        }

    Tensor& to_orth = (dir==Fromleft ? A : B);
    Tensor& newoc   = (dir==Fromleft ? B : A);

    CombinerT comb;

    Foreach(const IndexT& I, to_orth.indices())
        { 
        if(!hasindex(newoc,I))
            comb.addleft(I);
        }

    //Apply combiner
    comb.doCondense(true);
    comb.init(mid.isNull() ? "mid" : mid.rawname());

    Tensor AAc; 
    comb.product(AA,AAc);

    //Form density matrix
    Tensor AAcc = conj(AAc); 
    AAcc.prime(comb.right()); 

    Tensor rho = AAc*AAcc; 

    //Add noise term if requested
    if(spec.noise() > 0 && !PH.isNull())
        {
        rho += spec.noise()*PH.deltaRho(AA,comb,dir);
        rho *= 1./trace(realPart(rho));
        }

    const Real saved_cutoff = spec.cutoff(); 
    const int saved_minm = spec.minm(),
              saved_maxm = spec.maxm(); 
    if(spec.useOrigM())
        {
        spec.cutoff(-1);
        spec.minm(mid.m());
        spec.maxm(mid.m());
        }

    if(opts.getBool("TraceReIm",false))
        {
        rho = realPart(rho);
        }

    Tensor U;
    SparseT D;
    Real truncerr = diag_hermitian(rho,U,D,spec,opts);
    spec.truncerr(truncerr);

    spec.cutoff(saved_cutoff);
    spec.minm(saved_minm);
    spec.maxm(saved_maxm);

    comb.conj();
    comb.product(conj(U),to_orth);
    newoc = U * AAc;

    } //denmatDecomp



template<class Tensor, class SparseT>
void 
diagHermitian(const Tensor& M, Tensor& U, SparseT& D, Spectrum& spec,
              const OptSet& opts)
    {
    typedef typename Tensor::IndexT 
    IndexT;
    typedef typename Tensor::CombinerT 
    CombinerT;

    if(isZero(M,Opt("Fast"))) 
        throw ResultIsZero("denmatDecomp: M is zero");

    CombinerT comb;
    Foreach(const IndexT& I, M.indices())
        { 
        if(I.primeLevel() == 0)
            {
            comb.addleft(I);
            }
        }

    //Apply combiner
    comb.doCondense(true);
    comb.init("d");

    Tensor Mc; 
    comb.product(M,Mc);

    CombinerT combP(comb);
    combP.prime();
    combP.conj();

    try {
        Mc = combP * Mc;
        }
    catch(const ITError& e)
        {
        Cout << "Diagonalize expects opposite arrow directions for primed and unprimed indices." << Endl;
        throw e;
        }

    const bool truncate_setting = spec.truncate();
    spec.truncate(false);
    diag_hermitian(Mc,U,D,spec,opts);
    spec.truncerr(0);
    spec.truncate(truncate_setting);

    U = comb * U;

    } //diagHermitian


template<class Tensor>
void 
orthoDecomp(Tensor T, Tensor& A, Tensor& B, Direction dir, 
            Spectrum& spec,
            const OptSet& opts)
    {
    typedef typename Tensor::IndexT 
    IndexT;
    typedef typename Tensor::CombinerT 
    CombinerT;
    typedef typename Tensor::SparseT
    SparseT;

    if(isZero(T,Opt("Fast"))) 
        throw ResultIsZero("orthoDecomp: T is zero");

    const Real orig_noise = spec.noise();
    spec.noise(0);

    const
    bool usedenmat = false;

    if(usedenmat)
        {
        denmatDecomp(T,A,B,dir,spec,opts & Opt("TraceReIm"));
        }
    else //use svd
        {
        //Combiners which transform T
        //into a rank 2 tensor
        CombinerT Acomb, Bcomb;
        Acomb.doCondense(true);
        Bcomb.doCondense(true);

        const
        IndexT reim = IQIndex("ReIm",Index("reim",2),QN());

        //Divide up indices based on U
        //If U is null, use V instead
        const Tensor &L = (A.isNull() ? B : A);
        CombinerT &Lcomb = (A.isNull() ? Bcomb : Acomb),
                  &Rcomb = (A.isNull() ? Acomb : Bcomb);
        Foreach(const IndexT& I, T.indices())
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


        SparseT D;
        svdRank2(T,Acomb.right(),Bcomb.right(),A,D,B,spec,opts);

        A = conj(Acomb) * A;
        B = B * conj(Bcomb);

        if(dir==Fromleft) 
            {
            B *= D;
            B = B*conj(reim)(1) + Complex_i*B*conj(reim)(2);
            }
        else              
            {
            A *= D;
            A = A*conj(reim)(1) + Complex_i*A*conj(reim)(2);
            }
        }

    spec.noise(orig_noise);
    } //orthoDecomp


#undef Cout
#undef Format
#undef Endl


#endif
